// Workflow:
// 1) Parse accession list and ENST variant file
// 2) Fetch UniProt entries (batched) — provides canonical sequence, isoform VAR_SEQ
//    features, CC ALTERNATIVE PRODUCTS, and ENST→isoform cross-references
// 3) Build global ENST→isoform mapping from cross-references; assert coverage
// 4) Reconstruct isoform sequences by applying VAR_SEQ features to canonical
// 5) Fetch EBI variation data (batched) for all canonical + isoform accessions
// 6) Join EBI variants with ENST-derived variants per isoform; write output
// 7) ENSTs unmapped to any UniProt entry → fetch sequence from Ensembl REST;
//    emit synthetic flat-file entries
//
// Warnings and skipped/malformed records are logged via `tracing`, writing
// to both stderr and the exceptions file given via `--exceptions` (see
// `init_logging` below) — there is no more in-process ExceptionLog type.
//
// Usage: map_variants --accessions <accessions.txt> --variants <variants.txt> [--output <out.txt>]

mod cli;
mod ebi;
mod ensembl;
mod format;
mod isoform;
mod parsing;
mod types;
mod uniprot;
mod util;
mod variant;

use cli::Cli;
use ensembl::fetch_ensembl_sequences;
use format::format_entry;
use isoform::collect_and_reconstruct_isoforms;
use tracing::info;
use tracing_subscriber::Layer;
use variant::*;
use parsing::{parse_accession_list, parse_sample_variants};
use types::{EnsemblId, UniProtId};
use uniprot::{fetch_uniprot_entries, UniProtEntry, get_ensembl_mapping};
use util::RetryConfig;

use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, Write};
use std::sync::Arc;
use std::time::Duration;

use crate::ebi::{feature_to_canon_variant, feature_to_iso_variant, fetch_canon_variations, fetch_iso_variations, EBI_MAX_REQUESTS_PER_SECOND};
use crate::types::UniProtCanonId;
use crate::uniprot::UniProtFeature;
use crate::util::RateLimiter;
use futures::stream::{self, StreamExt};
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

/// Caps how many UniProt entries are processed (and thus how many EBI
/// fetches are in flight) concurrently. EBI request *rate* is separately
/// capped by the shared `RateLimiter` passed into every task, so this only
/// bounds concurrency, not throughput.
const MAX_CONCURRENT_ENTRIES: usize = 16;

/// Sets up `tracing` to write every `warn!`/`error!` (and above) both to
/// stderr (for live progress) and to `exceptions_path` (for a persistent,
/// traceable record) — replacing the old hand-rolled `ExceptionLog`.
///
/// Returns the `tracing-appender` worker guard, which must be kept alive
/// for the duration of `main` — the file writer is non-blocking and
/// buffers writes on a background thread; dropping the guard early can
/// lose buffered log lines.
fn init_logging(exceptions_path: &str) -> Result<tracing_appender::non_blocking::WorkerGuard> {
    let log_file = File::create(exceptions_path)
        .with_context(|| format!("failed to create exceptions file '{}'", exceptions_path))?;
    let (non_blocking, guard) = tracing_appender::non_blocking(log_file);

    let file_layer = tracing_subscriber::fmt::layer()
        .with_writer(non_blocking)
        .with_ansi(false)
        .with_target(false);

    let stderr_layer = tracing_subscriber::fmt::layer()
        .with_writer(std::io::stderr)
        .with_ansi(true)
        .with_target(false)
        .with_filter(tracing_subscriber::filter::LevelFilter::WARN);

    tracing_subscriber::registry()
        .with(file_layer)
        .with(stderr_layer)
        .try_init()
        .context("failed to initialize tracing subscriber")?;

    Ok(guard)
}

#[tokio::main]
async fn main() -> anyhow::Result<()>{
    use clap::Parser;
    let args = Cli::parse();
    let source_type_refs: Vec<&str> = args.source_type.iter().map(String::as_str).collect();
    let retry_config = RetryConfig {
        timeout: Duration::from_secs(args.timeout_secs),
        max_attempts: args.max_attempts,
        backoff_base_ms: args.retry_backoff_ms,
    };

    let _log_guard = init_logging(&args.exceptions)
        .context("unable to set up exception logging")?;


    run(
        &args.accessions,
        &args.variants,
        args.output.as_deref(),
        &source_type_refs,
        &retry_config,
        args.ebi_variants,
        args.uniprot_variants,
        args.ensembl_fallback
    )
        .await
        .context("run failed")?;

    Ok(())
}

/// Write one synthetic flat-file entry, followed by a blank line, to `out`.
fn write_entry(
    out: &mut dyn Write,
    accession: &UniProtCanonId,
    sequence: &str,
    var_seq: Option<&[UniProtFeature]>,
    variants: &[Variant],
    entry: Option<&UniProtEntry>,
) -> Result<()> {
    writeln!(out, "{}", format_entry(accession, sequence, var_seq, variants, entry))
        .with_context(|| format!("failed to write output entry for {}", accession.as_str()))
}

/// Attempts to resolve a premature-stop / out-of-bounds condition for `var`
/// against `seq`. On success, pushes `var` into `collection`; on failure,
/// logs the variant id and error via `tracing::info!` and drops the variant.
fn try_add_variant(collection: &mut Vec<Variant>, mut var: Variant, seq: &str) {
    match var.resolve_stop(seq) {
        Ok(_) => collection.push(var),
        Err(e) => tracing::info!("{},{}", var.id.as_str(), e),
    }
}

/// Processes a single UniProt entry: reconstructs isoforms, joins
/// ENST-derived sample variants with UniProt-annotated and EBI-derived
/// variants, and returns the fully formatted output line plus the set of
/// ENSTs it consumed. Does not write to `out` itself — callers collect
/// these and write sequentially, since entries may be processed
/// concurrently but output order doesn't need to match input order.
async fn process_entry(
    canonical_id: UniProtCanonId,
    entry: UniProtEntry,
    global_sample_variants: Arc<HashMap<EnsemblId, Vec<Variant>>>,
    source_types: Arc<Vec<String>>,
    retry_config: RetryConfig,
    ebi_rate_limiter: RateLimiter,
    ebi_client: reqwest::Client,
    ebi_variants_enabled: bool,
    uniprot_variants_enabled: bool,
) -> Result<(String, HashSet<EnsemblId>)> {
    let source_type_refs: Vec<&str> = source_types.iter().map(String::as_str).collect();
    let mut assigned_enst = HashSet::new();

    let (
        canonical_iso_alias,
        isoforms,
        var_seq_features,
        uniprot_variants,
    ) = collect_and_reconstruct_isoforms(&entry)?;

    // collect ensembl variants
    let mut entry_sample_variants: Vec<Variant> = Vec::new();

    let entry_ensembl_map = get_ensembl_mapping(&entry.cross_references)?;

    for (uniprot_id, ensembl_id) in entry_ensembl_map {
        let matches = match &uniprot_id {
            UniProtId::Id(id) => id == &canonical_id,
            UniProtId::Iso(id) => Some(id) == canonical_iso_alias.as_ref(),
        };

        if matches {
            if let Some(variants) = global_sample_variants.get(&ensembl_id) {
                for var in variants {
                    try_add_variant(&mut entry_sample_variants, var.clone(), &entry.sequence.value);
                }
                assigned_enst.insert(ensembl_id);
            }
        } else {
            let iso_id = match uniprot_id {
                UniProtId::Iso(iso_id) => iso_id,
                UniProtId::Id(_) => continue,
            };

            if let Some(iso_variants) = global_sample_variants.get(&ensembl_id) {
                if let Some(iso_seq) = isoforms.get(&iso_id).map(|s| s.as_str()) {
                    for var in iso_variants {
                        let mut var = var.clone();
                        var.isoform = Some(iso_id.clone());
                        try_add_variant(&mut entry_sample_variants, var, iso_seq);
                    }
                    assigned_enst.insert(ensembl_id);
                }
            }
        }
    }

    let mut combined_variants = entry_sample_variants;

    if uniprot_variants_enabled {
        combined_variants.join(uniprot_variants);
    }

    if ebi_variants_enabled {
        let mut ebi_canon_variants = fetch_canon_variations(
            &vec![canonical_id.clone()],
            &source_type_refs,
            &retry_config,
            &ebi_rate_limiter,
            &ebi_client,
        ).await?;

        let mut iso_keys = Vec::new();
        for isoform in isoforms.clone() {
            iso_keys.push(isoform.0);
        }
        let mut ebi_iso_variants = fetch_iso_variations(
            &iso_keys,
            &source_type_refs,
            &retry_config,
            &ebi_rate_limiter,
            &ebi_client,
        ).await?;

        let mut ebi_variants = Vec::new();

        for (_, feat_info) in ebi_canon_variants.drain() {
            for ft in feat_info.features {
                if let Some(canon_variant) = feature_to_canon_variant(&entry.primary_accession, &ft) {
                    try_add_variant(&mut ebi_variants, canon_variant, &entry.sequence.value);
                }
            }
        }

        for (iso_id, feat_info) in ebi_iso_variants.drain() {
            for ft in feat_info.features {
                if let Some(iso_seq) = isoforms.get(&iso_id) {
                    if let Some(iso_variant) = feature_to_iso_variant(&iso_id, &ft) {
                        try_add_variant(&mut ebi_variants, iso_variant, iso_seq.as_str());
                    }
                }
            }
        }

        combined_variants.join(ebi_variants);
    }

    let line = format_entry(&canonical_id, &entry.sequence.value, Some(&var_seq_features), &combined_variants, Some(&entry));

    Ok((line, assigned_enst))
}

async fn run(
    uniprot_accessions_path: &str,
    variants_path: &Option<String>,
    output_path: Option<&str>,
    source_types: &[&str],
    retry_config: &RetryConfig,
    ebi_variants_enabled: bool,
    uniprot_variants_enabled: bool,
    use_ensembl_fallback: bool,
) -> Result<()> {
    let mut out: Box<dyn Write> = match output_path {
        Some(path) => Box::new(
            File::create(path).with_context(|| format!("failed to create output file '{}'", path))?,
        ),
        None => Box::new(io::stdout()),
    };

    // 1. Parse inputs
    let accessions = parse_accession_list(uniprot_accessions_path)?;
    let global_sample_variants = match variants_path {
        Some(path) => parse_sample_variants(path)?,
        None => HashMap::new(),
    };
    let global_sample_variants = Arc::new(global_sample_variants);

    let source_types_owned = Arc::new(
        source_types.iter().map(|s| s.to_string()).collect::<Vec<String>>()
    );

    // 2. Fetch all UniProt entries
    let uniprot_entries = fetch_uniprot_entries(&accessions, retry_config).await?;

    // 3. Build global ENST→isoform mapping from all entry cross-references,
    //    plus a canonical-id -> isoform-ids index.
    let mut global_assinged_enst = HashSet::new();

    let ebi_rate_limiter = RateLimiter::per_second(EBI_MAX_REQUESTS_PER_SECOND);
    let ebi_client = reqwest::Client::new();

    let total_entries = uniprot_entries.len();
    let mut completed = 0usize;

    let mut results = stream::iter(uniprot_entries.into_iter())
        .map(|(canonical_id, entry)| {
            let global_sample_variants = Arc::clone(&global_sample_variants);
            let source_types_owned = Arc::clone(&source_types_owned);
            let retry_config = *retry_config;
            let ebi_rate_limiter = ebi_rate_limiter.clone();
            let ebi_client = ebi_client.clone();
            async move {
                let result = process_entry(
                    canonical_id.clone(),
                    entry,
                    global_sample_variants,
                    source_types_owned,
                    retry_config,
                    ebi_rate_limiter,
                    ebi_client,
                    ebi_variants_enabled,
                    uniprot_variants_enabled,
                ).await;
                (canonical_id, result)
            }
        })
        .buffer_unordered(MAX_CONCURRENT_ENTRIES);

    while let Some((canonical_id, result)) = results.next().await {
        match result {
            Ok((line, assigned_enst)) => {
                writeln!(out, "{}", line)
                    .with_context(|| format!("failed to write output entry for {}", canonical_id.as_str()))?;
                global_assinged_enst.extend(assigned_enst);
            }
            Err(e) => {
                tracing::warn!(
                    entry = canonical_id.as_str(),
                    "skipping entry entirely due to error: {:?}", e
                );
            }
        }

        completed += 1;
        eprint!("\rProcessed entry {}/{}: {}", completed, total_entries, canonical_id.as_str());
        io::stderr().flush().ok();
    }
    eprintln!();

    // 8. ENSTs with no UniProt cross-reference mapping → synthetic entries
    let unmapped_enst_ids: Vec<EnsemblId> = global_sample_variants
        .keys()
        .filter(|ensembl_id| !global_assinged_enst.contains(*ensembl_id))
        .cloned()
        .collect();

    if !unmapped_enst_ids.is_empty() {
        if use_ensembl_fallback {
            let enst_sequences = fetch_ensembl_sequences(&unmapped_enst_ids, retry_config).await?;
            for (enst_id, seq) in &enst_sequences {
                let variants = global_sample_variants.get(enst_id).map(Vec::as_slice).unwrap_or(&[]);
                let synthetic_id = UniProtCanonId(enst_id.as_str().to_string());
                write_entry(&mut *out, &synthetic_id, seq.as_str(), None, &variants, None)?;
            }
        }
        else {
            for enst in unmapped_enst_ids {
                let variants = global_sample_variants.get(&enst).map(Vec::as_slice).unwrap_or(&[]);
                for var in variants {
                    info!("{},Variant {} could not be mapped", &enst.as_str(), &var.id);
                }
            }
        }
    }
    
    Ok(())
}