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
// Usage: map_variants --accessions <accessions.txt> --variants <variants.txt> [--output <out.txt>]

mod cli;
mod ebi;
mod ensembl;
mod exceptions;
mod format;
mod isoform;
mod parsing;
mod types;
mod uniprot;
mod util;
mod variant;

use cli::Cli;
use ensembl::fetch_ensembl_sequences;
use exceptions::ExceptionLog;
use format::format_entry;
use isoform::{
    collect_isoform_ids, 
    collect_and_reconstruct_isoforms,
};
use variant::*;
use parsing::{parse_accession_list, parse_sample_variants};
use types::{EnsemblId, Isoform, Sequence, UniProtId, UniProtIsoId};
use uniprot::{fetch_uniprot_entries, UniProtEntry, get_ensembl_mapping};
use util::RetryConfig;

use anyhow::{anyhow, Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, Write};
use std::process;
use std::time::Duration;

use crate::ebi::{ProteinFeatureInfo, feature_to_canon_variant, feature_to_iso_variant, fetch_canon_variations, fetch_iso_variations};
use crate::types::UniProtCanonId;
use crate::uniprot::UniProtFeature;

fn main() {
    use clap::Parser;
    let args = Cli::parse();
    let source_type_refs: Vec<&str> = args.source_type.iter().map(String::as_str).collect();
    let retry_config = RetryConfig {
        timeout: Duration::from_secs(args.timeout_secs),
        max_attempts: args.max_attempts,
        backoff_base_ms: args.retry_backoff_ms,
    };

    if let Err(e) = run(
        &args.accessions,
        &args.variants,
        args.output.as_deref(),
        &source_type_refs,
        &args.exceptions,
        &retry_config,
    ) {
        eprintln!("Error: {:?}", e);
        process::exit(1);
    }
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

fn run(
    uniprot_accessions_path: &str,
    variants_path: &str,
    output_path: Option<&str>,
    source_types: &[&str],
    exceptions_path: &str,
    retry_config: &RetryConfig,
) -> Result<()> {
    let mut out: Box<dyn Write> = match output_path {
        Some(path) => Box::new(
            File::create(path).with_context(|| format!("failed to create output file '{}'", path))?,
        ),
        None => Box::new(io::stdout()),
    };
    let mut exceptions = ExceptionLog::to_file(exceptions_path)?;

    // 1. Parse inputs
    let accessions = parse_accession_list(uniprot_accessions_path)?;
    let global_sample_variants = parse_sample_variants(variants_path, &mut exceptions)?;

    // 2. Fetch all UniProt entries
    let uniprot_entries = fetch_uniprot_entries(&accessions, retry_config)?;

    // 3. Build global ENST→isoform mapping from all entry cross-references,
    //    plus a canonical-id -> isoform-ids index.
    let mut global_assinged_enst = HashSet::new();

    let total = uniprot_entries.len();
    for (i, (canonical_id, entry)) in uniprot_entries.iter().enumerate() {
        eprint!("\rProcessing entry {}/{}: {}", i + 1, total, canonical_id.as_str());
        io::stderr().flush().ok();
        let (
            canonical_iso_alias, 
            isoforms, 
            var_seq_features, 
            uniprot_variants
        ) = collect_and_reconstruct_isoforms(entry)?;

        // collect ensembl variants
        let mut entry_sample_variants: Vec<Variant> = Vec::new();
        
        let entry_ensembl_map = get_ensembl_mapping(&entry.cross_references)?;     

        for (uniprot_id, ensembl_id) in entry_ensembl_map { 
            let matches = match &uniprot_id {
                UniProtId::Id(id) => id == canonical_id,
                UniProtId::Iso(id) => Some(id) == canonical_iso_alias.as_ref(),
            };

            //variants for canonical sequence
            if matches {
                if let Some(variants) = global_sample_variants.get(&ensembl_id) {
                    for var in variants {
                        let mut var = var.clone();
                        var.resolve_stop(&entry.sequence.value)?;
                        entry_sample_variants.push(var);
                    }
                    global_assinged_enst.insert(ensembl_id);
                }
            }
            
            // variants for isoforms
            else {
                let iso_id = match uniprot_id {
                    UniProtId::Iso(iso_id) => iso_id,
                    UniProtId::Id(_) => continue,
                };

                if let Some(iso_variants) = global_sample_variants.get(&ensembl_id) {
                    if let Some(iso_seq) = isoforms.get(&iso_id).map(|s| s.as_str()) {
                        for var in iso_variants {
                            let mut var = var.clone();
                            var.isoform = Some(iso_id.clone());
                            var.resolve_stop(iso_seq)?;
                            entry_sample_variants.push(var);
                        }
                        global_assinged_enst.insert(ensembl_id);
                    }
                }
            }
        }

        let mut combined_variants = entry_sample_variants;
        combined_variants.join(uniprot_variants);

        // collect ebi variants
        // canonical
        let mut ebi_canon_variants= fetch_canon_variations(&vec![canonical_id.clone()], source_types, retry_config)?;

        // isoform
        let mut iso_keys = Vec::new();
        for isoform in isoforms.clone() {
            iso_keys.push(isoform.0);
        }
        let mut ebi_iso_variants = fetch_iso_variations(&iso_keys, source_types, retry_config)?;

        let mut ebi_variants = Vec::new();

        for (_, feat_info) in ebi_canon_variants.drain() {
            for ft in feat_info.features {
                if let Some(mut canon_variant) = feature_to_canon_variant(&ft) {
                    canon_variant.resolve_stop(&entry.sequence.value)?;
                    ebi_variants.push(canon_variant)
                } 
            }
        }

        for (iso_id, feat_info) in ebi_iso_variants.drain() {
            for ft in feat_info.features {
                if let Some(iso_seq) = isoforms.get(&iso_id) {
                    if let Some(mut iso_variant) = feature_to_iso_variant(&iso_id, &ft) {
                        iso_variant.resolve_stop(iso_seq.as_str())?;
                        ebi_variants.push(iso_variant);   
                    }
                }
            }
        }

        combined_variants.join(ebi_variants);

        write_entry(&mut *out, canonical_id, &entry.sequence.value, Some(&var_seq_features), &combined_variants, Some(entry))?;
    }


    // 8. ENSTs with no UniProt cross-reference mapping → synthetic entries
    let unmapped_enst_ids: Vec<EnsemblId> = global_sample_variants
        .keys()
        .filter(|ensembl_id| !global_assinged_enst.contains(*ensembl_id))
        .cloned()
        .collect();

    if !unmapped_enst_ids.is_empty() {
        let enst_sequences = fetch_ensembl_sequences(&unmapped_enst_ids, &mut exceptions, retry_config)?;
        for (enst_id, seq) in &enst_sequences {
            let variants = global_sample_variants.get(enst_id).map(Vec::as_slice).unwrap_or(&[]);
            let synthetic_id = UniProtCanonId(enst_id.as_str().to_string());
            write_entry(&mut *out, &synthetic_id, seq.as_str(), None, &variants, None)?;
        }
    }

    Ok(())
}
