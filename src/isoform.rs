// ============================================================================
// ENST → isoform mapping, coverage assertion, and isoform sequence
// reconstruction (used to compute isoform-local variant positions; the
// isoforms themselves are no longer written as separate synthetic entries —
// see main.rs)
// ============================================================================

use crate::types::{EnsemblId, Isoform, Sequence, UniprotId, Variant};
use crate::uniprot::{UniProtEntry, UniProtFeature};
use crate::util::strip_version;
use anyhow::{anyhow, Result};
use std::collections::HashMap;

/// The isoform id (e.g. "P04637-1") that UniProt marks "Displayed" in the
/// ALTERNATIVE PRODUCTS comment, if any. Its sequence is by definition
/// exactly the entry's canonical sequence — it is not a distinct protein.
fn displayed_isoform_id(entry: &UniProtEntry) -> Option<UniprotId> {
    entry
        .comments
        .iter()
        .find(|c| c.comment_type == "ALTERNATIVE PRODUCTS")
        .and_then(|c| {
            c.isoforms
                .iter()
                .find(|iso| iso.isoform_sequence_status.as_deref() == Some("Displayed"))
        })
        .and_then(|iso| iso.isoform_ids.first())
        .map(|id| UniprotId(id.clone()))
}

/// Extract the ENST → specific UniprotId (isoform or canonical) mapping from
/// an entry's Ensembl cross-references. ENSTs without an explicit `isoformId`
/// are mapped to the canonical accession.
pub(crate) fn extract_enst_to_isoform(
    entry: &UniProtEntry,
    canonical_id: &UniprotId,
) -> HashMap<EnsemblId, UniprotId> {
    entry
        .cross_references
        .iter()
        .filter(|xr| xr.database == "Ensembl")
        .map(|xr| {
            let enst = EnsemblId(strip_version(&xr.id).to_string());
            let target = xr
                .isoform_id
                .as_ref()
                .map(|iso| UniprotId(iso.clone()))
                .unwrap_or_else(|| canonical_id.clone());
            (enst, target)
        })
        .collect()
}

/// Redirect any isoform id in `mapping` that is merely an alias for the
/// canonical sequence (the "Displayed" isoform) to `canonical_id`. Must run
/// after [`check_enst_isoform_coverage`], which relies on seeing the raw,
/// pre-redirect mapping to tell an explicit "-1" assignment (fine) apart
/// from a defaulted one (ambiguous, when the entry has other isoforms).
pub(crate) fn canonicalize_displayed_isoform(
    mapping: HashMap<EnsemblId, UniprotId>,
    entry: &UniProtEntry,
    canonical_id: &UniprotId,
) -> HashMap<EnsemblId, UniprotId> {
    let displayed_id = displayed_isoform_id(entry);
    mapping
        .into_iter()
        .map(|(enst, target)| {
            let target = if Some(&target) == displayed_id.as_ref() {
                canonical_id.clone()
            } else {
                target
            };
            (enst, target)
        })
        .collect()
}

/// Assert that every ENST in the variant list that maps to this entry has an
/// explicit isoform mapping. Logs a warning (but continues) when an ENST
/// resolves to the canonical accession of an entry that has isoforms — we
/// can't determine which isoform it belongs to in that case.
pub(crate) fn check_enst_isoform_coverage(
    enst_to_isoform: &HashMap<EnsemblId, UniprotId>,
    enst_variants: &HashMap<EnsemblId, Vec<Variant>>,
    entry: &UniProtEntry,
    canonical_id: &UniprotId,
) {
    let has_isoforms = entry
        .comments
        .iter()
        .any(|c| c.comment_type == "ALTERNATIVE PRODUCTS" && !c.isoforms.is_empty());

    if !has_isoforms {
        return;
    }

    for enst in enst_variants.keys() {
        match enst_to_isoform.get(enst) {
            Some(mapped_id) if mapped_id == canonical_id => {
                eprintln!(
                    "Warning: ENST {} maps to canonical {} but the entry has isoforms; \
                     isoform-specific variant assignment may be incorrect",
                    enst.as_str(),
                    canonical_id.as_str()
                );
            }
            None => {
                // Not mapped to this entry at all — handled globally
            }
            _ => {}
        }
    }
}

/// Apply a set of VAR_SEQ features to the canonical sequence to produce an
/// isoform sequence. Features must be applied in descending position order
/// (handled here) to avoid index drift from earlier substitutions/deletions.
fn apply_var_seq_features(canonical: &str, features: &[&UniProtFeature]) -> Result<String> {
    let mut sorted: Vec<&UniProtFeature> = features.to_vec();
    sorted.sort_by(|a, b| b.location.start.value.cmp(&a.location.start.value));

    let mut seq = canonical.to_string();

    for feature in sorted {
        let start = feature.location.start.value; // 1-based
        let end = feature.location.end.value; // 1-based inclusive

        if start == 0 || start > seq.len() || end > seq.len() {
            return Err(anyhow!(
                "VAR_SEQ position {}..{} out of range for sequence length {}",
                start,
                end,
                seq.len()
            ));
        }

        let range_start = start - 1;
        let range_end = end;

        if let Some(orig) = feature
            .alternative_sequence
            .as_ref()
            .and_then(|a| a.original_sequence.as_deref())
        {
            let actual = &seq[range_start..range_end];
            if actual != orig {
                return Err(anyhow!(
                    "VAR_SEQ original sequence mismatch at {}..{}: expected '{}', found '{}'",
                    start,
                    end,
                    orig,
                    actual
                ));
            }
        }

        let replacement = feature
            .alternative_sequence
            .as_ref()
            .and_then(|a| a.alternative_sequences.first())
            .map(|s| s.as_str())
            .unwrap_or(""); // empty = "Missing" → deletion

        seq = format!("{}{}{}", &seq[..range_start], replacement, &seq[range_end..]);
    }

    Ok(seq)
}

/// Reconstruct all isoform sequences for an entry from its VAR_SEQ features
/// and ALTERNATIVE PRODUCTS comment. Returns only isoforms whose sequence
/// could be successfully reconstructed.
pub(crate) fn reconstruct_isoforms(entry: &UniProtEntry) -> Vec<Isoform> {
    let vsp_features: HashMap<&str, &UniProtFeature> = entry
        .features
        .iter()
        .filter(|f| f.feature_type == "Alternative sequence")
        .filter_map(|f| f.feature_id.as_deref().map(|id| (id, f)))
        .collect();

    let alt_products = match entry
        .comments
        .iter()
        .find(|c| c.comment_type == "ALTERNATIVE PRODUCTS")
    {
        Some(c) => c,
        None => return Vec::new(),
    };

    let mut isoforms = Vec::new();

    for iso_info in &alt_products.isoforms {
        let isoform_id = match iso_info.isoform_ids.first() {
            Some(id) => UniprotId(id.clone()),
            None => continue,
        };

        match iso_info.isoform_sequence_status.as_deref().unwrap_or("") {
            // "Displayed" means this isoform's sequence *is* the canonical
            // sequence — it's not a distinct protein, so it gets no synthetic
            // entry of its own; the canonical entry already covers it.
            "Displayed" => {}
            "Described" => {
                let features: Vec<&UniProtFeature> = iso_info
                    .sequence_ids
                    .iter()
                    .filter_map(|id| vsp_features.get(id.as_str()).copied())
                    .collect();
                if features.len() != iso_info.sequence_ids.len() {
                    eprintln!(
                        "Warning: isoform {} references {} VSP IDs but only {} found",
                        isoform_id.as_str(),
                        iso_info.sequence_ids.len(),
                        features.len()
                    );
                }
                match apply_var_seq_features(&entry.sequence.value, &features) {
                    Ok(seq) => isoforms.push(Isoform(isoform_id, Sequence(seq))),
                    Err(e) => eprintln!(
                        "Warning: failed to reconstruct isoform {}: {:?}",
                        isoform_id.as_str(),
                        e
                    ),
                }
            }
            "Not described" => {
                eprintln!(
                    "Warning: isoform {} is 'Not described'; skipping",
                    isoform_id.as_str()
                );
            }
            other => {
                eprintln!(
                    "Warning: unknown isoform status '{}' for {}; skipping",
                    other,
                    isoform_id.as_str()
                );
            }
        }
    }

    isoforms
}

// ============================================================================
// Variant grouping
// ============================================================================

pub(crate) fn group_enst_variants_by_isoform(
    enst_to_variants: &HashMap<EnsemblId, Vec<Variant>>,
    enst_to_isoform: &HashMap<EnsemblId, UniprotId>,
) -> HashMap<UniprotId, Vec<(EnsemblId, usize, Variant)>> {
    let mut grouped: HashMap<UniprotId, Vec<(EnsemblId, usize, Variant)>> = HashMap::new();
    for (enst_id, variants) in enst_to_variants {
        if let Some(uniprot_id) = enst_to_isoform.get(enst_id) {
            let entry = grouped.entry(uniprot_id.clone()).or_default();
            for (idx, variant) in variants.iter().enumerate() {
                entry.push((enst_id.clone(), idx, variant.clone()));
            }
        }
    }
    grouped
}

#[cfg(test)]
#[path = "isoform_tests.rs"]
mod tests;
