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
mod join;
mod parsing;
mod types;
mod uniprot;
mod util;

use cli::Cli;
use ebi::{fetch_variations, variants_from_protein_feature_info};
use ensembl::fetch_ensembl_sequences;
use exceptions::ExceptionLog;
use format::format_entry;
use isoform::{
    canonicalize_displayed_isoform, check_enst_isoform_coverage, extract_enst_to_isoform,
    group_enst_variants_by_isoform, reconstruct_isoforms,
};
use join::{append_id_tag, join_variants};
use parsing::{parse_accession_list, parse_variants, resolve_stop_variants};
use types::{EnsemblId, Isoform, Sequence, UniprotId, Variant};
use uniprot::{fetch_uniprot_entries, UniProtEntry};

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use std::process;

fn main() {
    use clap::Parser;
    let args = Cli::parse();
    let source_type_refs: Vec<&str> = args.source_type.iter().map(String::as_str).collect();

    if let Err(e) = run(
        &args.accessions,
        &args.variants,
        args.output.as_deref(),
        &source_type_refs,
        &args.exceptions,
    ) {
        eprintln!("Error: {:?}", e);
        process::exit(1);
    }
}

/// Write one synthetic flat-file entry, followed by a blank line, to `out`.
fn write_entry(
    out: &mut dyn Write,
    accession: &UniprotId,
    sequence: &Sequence,
    variants: &[Variant],
    entry: Option<&UniProtEntry>,
) -> Result<()> {
    writeln!(out, "{}", format_entry(accession, sequence, variants, entry))
        .with_context(|| format!("failed to write output entry for {}", accession.as_str()))
}

fn run(
    accessions_path: &str,
    variants_path: &str,
    output_path: Option<&str>,
    source_types: &[&str],
    exceptions_path: &str,
) -> Result<()> {
    let mut out: Box<dyn Write> = match output_path {
        Some(path) => Box::new(
            File::create(path).with_context(|| format!("failed to create output file '{}'", path))?,
        ),
        None => Box::new(io::stdout()),
    };
    let mut exceptions = ExceptionLog::to_file(exceptions_path)?;

    // 1. Parse inputs
    let accessions = parse_accession_list(accessions_path)?;
    let enst_to_variants = parse_variants(variants_path, &mut exceptions)?;

    // 2. Fetch all UniProt entries
    let uniprot_entries = fetch_uniprot_entries(&accessions);

    // 3. Build global ENST→isoform mapping from all entry cross-references
    let mut global_enst_to_isoform: HashMap<EnsemblId, UniprotId> = HashMap::new();
    for (canonical_id, entry) in &uniprot_entries {
        let mapping = extract_enst_to_isoform(entry, canonical_id);

        // Assertion: ENSTs in variant list that map to this entry should
        // have explicit isoform assignments when the entry has isoforms
        let relevant_enst_variants: HashMap<EnsemblId, Vec<Variant>> = enst_to_variants
            .iter()
            .filter(|(enst, _)| mapping.contains_key(*enst))
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        check_enst_isoform_coverage(&mapping, &relevant_enst_variants, entry, canonical_id);

        let mapping = canonicalize_displayed_isoform(mapping, entry, canonical_id);
        global_enst_to_isoform.extend(mapping);
    }

    // 4. Group ENST variants by their target isoform UniprotId
    let grouped_enst_variants =
        group_enst_variants_by_isoform(&enst_to_variants, &global_enst_to_isoform);

    // 5. Reconstruct isoforms for every entry; collect all UniprotIds for EBI fetch
    let mut entry_isoforms: HashMap<UniprotId, Vec<Isoform>> = HashMap::new();
    let mut all_uniprot_ids: Vec<UniprotId> = Vec::new();
    for (canonical_id, entry) in &uniprot_entries {
        all_uniprot_ids.push(canonical_id.clone());
        let isoforms = reconstruct_isoforms(entry);
        for Isoform(iso_id, _) in &isoforms {
            all_uniprot_ids.push(iso_id.clone());
        }
        entry_isoforms.insert(canonical_id.clone(), isoforms);
    }

    // 6. Fetch EBI variation data for all UniprotIds in one batched pass
    let rest_variants_by_uniprot = fetch_variations(&all_uniprot_ids, source_types);

    // 7. Process each entry: one flat-file entry per canonical accession — the
    // canonical entry itself, matching a real UniProt flat file (ProtGraph builds
    // one graph per canonical accession, folding in isoforms via VAR_SEQ notes and
    // the ALTERNATIVE PRODUCTS comment, not via separate files). Variants specific
    // to an isoform (from EBI or an isoform-tagged ENST cross-reference) are
    // included here too, tagged with `isoform_ref` so their FT line reads
    // "ACCESSION:position" instead of a bare position — ProtGraph's biopython-based
    // FT parser reads that as a remote cross-reference and routes the feature onto
    // that isoform's own vertex chain instead of canonical's.
    for (canonical_id, entry) in &uniprot_entries {
        let canonical_seq = Sequence(entry.sequence.value.clone());
        let rest_variants = rest_variants_by_uniprot
            .get(canonical_id)
            .map(variants_from_protein_feature_info)
            .unwrap_or_default();
        let enst_resolved = resolve_stop_variants(grouped_enst_variants.get(canonical_id), canonical_seq.as_str());
        let mut joined = join_variants(rest_variants, enst_resolved.as_ref());

        if let Some(isoforms) = entry_isoforms.get(canonical_id) {
            for Isoform(iso_id, iso_seq) in isoforms {
                let rest_variants = rest_variants_by_uniprot
                    .get(iso_id)
                    .map(variants_from_protein_feature_info)
                    .unwrap_or_default();
                let enst_resolved = resolve_stop_variants(grouped_enst_variants.get(iso_id), iso_seq.as_str());
                let mut iso_joined = join_variants(rest_variants, enst_resolved.as_ref());
                for v in &mut iso_joined {
                    v.isoform_ref = Some(iso_id.clone());
                }
                joined.extend(iso_joined);
            }
        }

        write_entry(&mut *out, canonical_id, &canonical_seq, &joined, Some(entry))?;
    }

    // 8. ENSTs with no UniProt cross-reference mapping → synthetic entries
    let unmapped_enst_ids: Vec<EnsemblId> = enst_to_variants
        .keys()
        .filter(|enst| !global_enst_to_isoform.contains_key(*enst))
        .cloned()
        .collect();

    if !unmapped_enst_ids.is_empty() {
        let enst_sequences = fetch_ensembl_sequences(&unmapped_enst_ids, &mut exceptions);
        for (enst_id, seq) in &enst_sequences {
            let variants = enst_to_variants.get(enst_id).map(Vec::as_slice).unwrap_or(&[]);
            let tagged_variants: Vec<Variant> = variants
                .iter()
                .enumerate()
                .map(|(idx, v)| {
                    let mut v = v.clone();
                    append_id_tag(&mut v.id, &format!("{}:{}", enst_id.as_str(), idx));
                    v
                })
                .collect();
            let synthetic_id = UniprotId(enst_id.as_str().to_string());
            write_entry(&mut *out, &synthetic_id, seq, &tagged_variants, None)?;
        }
    }

    Ok(())
}
