// ============================================================================
// ENST → isoform mapping, coverage assertion, and isoform sequence
// reconstruction (used to compute isoform-local variant positions; the
// isoforms themselves are no longer written as separate synthetic entries —
// see main.rs)
// ============================================================================

use crate::types::{EnsemblId, Isoform, Sequence, UniProtIsoId, UniProtCanonId};
use crate::variant::Variant;
use crate::uniprot::{UniProtEntry, UniProtFeature};
use anyhow::{anyhow, Result, Context};
use std::collections::HashMap;

pub(crate) fn collect_isoform_ids(
    entry: &UniProtEntry,
    canonical_id: &UniProtCanonId,
) -> (Option<UniProtIsoId>, HashMap<UniProtCanonId, Vec<UniProtIsoId>>) {
    let mut displayed_id: Option<UniProtIsoId> = None;
    let mut isoform_ids: Vec<UniProtIsoId> = Vec::new();

    for c in entry.comments.iter().filter(|c| c.comment_type == "ALTERNATIVE PRODUCTS") {
        for iso in &c.isoforms {
            if iso.isoform_sequence_status.as_deref() == Some("Displayed") {
                if let Some(first) = iso.isoform_ids.first() {
                    displayed_id = Some(first.parse().unwrap());
                }
                // alias for canonical id is not pushed to isoforms
                continue;
            }
            isoform_ids.extend(iso.isoform_ids.iter().map(|id| id.parse().unwrap()));
        }
    }

    let mut isoform_map = HashMap::new();
    isoform_map.insert(canonical_id.clone(), isoform_ids);

    (displayed_id, isoform_map)
}

/* 
    // Redirect Displayed-isoform aliases to the canonical id.
    let mapping = mapping
        .into_iter()
        .map(|(enst, target)| {
            let target = if Some(&target) == displayed_id.as_ref() {
                canonical_id.clone()
            } else {
                target
            };
            (enst, target)
        })
        .collect();
*/

/// Apply a set of VAR_SEQ features to the canonical sequence to produce an
/// isoform sequence. Features must be applied in descending position order
/// to avoid index drift from earlier substitutions/deletions.
fn apply_var_seq_features(canonical: &str, features: &[&UniProtFeature]) -> Result<String> {
    let mut sorted: Vec<&UniProtFeature> = features.to_vec();
    sorted.sort_by(|a, b| b.location.start.value.cmp(&a.location.start.value));

    let mut seq = canonical.to_string();

    for feature in sorted {
        let start = feature.location.start.value
                .ok_or_else(|| anyhow!("feature location start.value is null"))?;
            let end = feature.location.end.value
                .ok_or_else(|| anyhow!("feature location end.value is null"))?;

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
            .and_then(|a| Some(a.original_sequence.clone()))
        {
            let actual = &seq[range_start..range_end];
            //TODO
            if let Some(original) = orig {
                if actual != original{
                return Err(anyhow!(
                    "VAR_SEQ original sequence mismatch at {}..{}: expected '{}', found '{}'",
                    start,
                    end,
                    original,
                    actual
                ));
            }
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

/// Reconstruct a single "Described" isoform's sequence from its VAR_SEQ
/// features. Returns an error for any other status ("Not described",
/// unknown, etc.).
fn reconstruct_isoform(
    entry: &UniProtEntry,
    isoform_id: UniProtIsoId,
    status: &str,
    sequence_ids: &[String],
    vsp_features: &[UniProtFeature],
) -> Result<(UniProtIsoId, Sequence)> {
    match status {
        "Described" => {
            let features: Vec<&UniProtFeature> = sequence_ids
                .iter()
                .filter_map(|id| {
                    vsp_features
                        .iter()
                        .find(|f| f.feature_id.as_deref() == Some(id.as_str()))
                })
                .collect();

            if features.len() != sequence_ids.len() {
                return Err(anyhow!(
                    "isoform {} references {} VSP IDs but only {} found",
                    isoform_id.as_str(),
                    sequence_ids.len(),
                    features.len()
                ));
            }

            let seq = apply_var_seq_features(&entry.sequence.value, &features)
                .with_context(|| format!("failed to reconstruct isoform {}", isoform_id.as_str()))?;

            Ok((isoform_id, Sequence(seq)))
        }
        "Not described" => Err(anyhow!("isoform {} is 'Not described'", isoform_id.as_str())),
        other => Err(anyhow!(
            "unknown isoform status '{}' for {}",
            other,
            isoform_id.as_str()
        )),
    }
}

/// Collect isoform IDs, reconstruct their sequences, and collect VAR_SEQ
/// features, in a single pass over the ALTERNATIVE PRODUCTS comment(s).
///
/// Returns `(displayed_id, isoform_map, vsp_features)` where `isoform_map`
/// maps the canonical id to the list of reconstructed isoforms (each of
/// which already carries its own `UniProtIsoId`).
pub(crate) fn collect_and_reconstruct_isoforms(
    entry: &UniProtEntry
) -> Result<(
    Option<UniProtIsoId>,
    HashMap<UniProtIsoId, Sequence>,
    Vec<UniProtFeature>,
    Vec<Variant>
)> {
    let mut var_seq_features: Vec<UniProtFeature> = Vec::new();
    let mut uniprot_variants: Vec<Variant> = Vec::new();

    for f in &entry.features {
        if f.feature_type == "Alternative sequence" {
            var_seq_features.push(f.clone());
        } 
        else if f.feature_type == "Natural variant" {
            if let Some(begin) = f.location.start.value
                && let Some(alt_seq) = f.alternative_sequence.as_ref(){
                let mut iso_id: Option<UniProtIsoId> = None;
                if let Some(iso_ref) = &f.location.sequence {
                    iso_id =Some(iso_ref.to_string().parse().unwrap());
                }
                
                let end = f.location.end.value.unwrap_or(begin);
                let aa_ref = 
                    alt_seq.original_sequence.clone().unwrap_or(String::new());
                let aa_new = alt_seq.alternative_sequences.first().cloned().unwrap_or(String::new());
                uniprot_variants.push(
                    Variant {
                        isoform: iso_id,
                        id: f.feature_id.clone().unwrap_or(String::new()),
                        begin: begin,
                        end: end,
                        aa_ref: aa_ref,
                        aa_new: aa_new,
                    }
                );
            }
        }
    }

    let mut displayed_id: Option<UniProtIsoId> = None;
    let mut isoforms  = HashMap::new();

    for c in entry
        .comments
        .iter()
        .filter(|c| c.comment_type == "ALTERNATIVE PRODUCTS")
    {
        for iso in &c.isoforms {
            let isoform_id: Option<UniProtIsoId> =
                iso.isoform_ids.first().map(|id| id.parse().unwrap());

            if iso.isoform_sequence_status.as_deref() == Some("Displayed") {
                if let Some(id) = isoform_id {
                    displayed_id = Some(id);
                }
                // alias for canonical id is not pushed to isoforms
                continue;
            }

            let isoform_id = isoform_id
                .ok_or_else(|| anyhow!("isoform entry has no isoform id"))?;

            let isoform = reconstruct_isoform(
                entry,
                isoform_id,
                iso.isoform_sequence_status.as_deref().unwrap_or(""),
                &iso.sequence_ids,
                &var_seq_features,
            )?;

            isoforms.insert(isoform.0, isoform.1);
        }
    }

    Ok((displayed_id, isoforms, var_seq_features, uniprot_variants))
}

#[cfg(test)]
#[path = "isoform_tests.rs"]
mod tests;
