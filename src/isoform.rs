// ============================================================================
// ENST → isoform mapping, coverage assertion, and isoform sequence
// reconstruction (used to compute isoform-local variant positions; the
// isoforms themselves are no longer written as separate synthetic entries —
// see main.rs)
// ============================================================================

use crate::types::{Sequence, UniProtIsoId};
use crate::variant::Variant;
use crate::uniprot::{UniProtEntry, UniProtFeature};
use anyhow::{anyhow, Result};
use std::collections::HashMap;
/* 
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

/// Reconstruct a single isoform's sequence from its VAR_SEQ features.
fn reconstruct_isoform(
    entry: &UniProtEntry,
    isoform_id: UniProtIsoId,
    status: &str,
    sequence_ids: &[String],
    vsp_features: &[UniProtFeature],
) -> Option<(UniProtIsoId, Sequence)> {
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
                tracing::info!(
                    "{},{},isoform references {} VSP id(s) but only {} were found among this entry's 'Alternative sequence' features; skipping isoform",
                    entry.primary_accession,
                    isoform_id.as_str(),
                    sequence_ids.len(),
                    features.len()
                );
                return None;
            }

            match apply_var_seq_features(&entry.sequence.value, &features) {
                Ok(seq) => Some((isoform_id, Sequence(seq))),
                Err(e) => {
                    tracing::info!(
                        "{},{},failed to reconstruct isoform: {:?}; skipping isoform",
                        entry.primary_accession,
                        isoform_id.as_str(),
                        e
                    );
                    None
                }
            }
        }
        "Not described" => {
        tracing::info!(
            "{},{},isoform sequence status is 'Not described' (no VAR_SEQ features to apply); skipping isoform",
            entry.primary_accession,
            isoform_id.as_str()
        );
        None
        }
        other => {
            tracing::info!(
                "{},{},unknown/unsupported isoform sequence status '{}'; skipping isoform",
                entry.primary_accession,
                isoform_id.as_str(),
                other
            );
            None
        }
    }
}

fn resolve_isoform_reference(
    iso_ref: &Option<String>
) -> Option<UniProtIsoId>{
    let mut iso_id: Option<UniProtIsoId> = None;
        
    if let Some(iso) = iso_ref{
        //infallible unwrap for conversion
        iso_id =Some(iso.to_string().parse().unwrap());
    }
    iso_id
}

fn add_canonical_variant(
    accession: &String,
    feature: &UniProtFeature,
    uniprot_variants: &mut Vec<Variant>
) {
    let isoform_id = feature.feature_id.clone().unwrap_or(
        feature.description.clone().unwrap_or("No Identifier or Description".to_string())
    );

    // position is required
    if let Some(begin) = feature.location.start.value {
        
        
        let iso_id = resolve_isoform_reference(&feature.location.sequence);

        // end defaults to begin if not set
        let end = feature.location.end.value.unwrap_or(begin);

        // ref == None and new == None: Missing sequence
        let mut aa_ref = None;
        let mut aa_new = None;

        if let Some(aa_seq) = feature.alternative_sequence.as_ref() {
            let aa_new_vec = &aa_seq.alternative_sequences;

            // if AlternativeSequence is set, original sequence and alternative sequence must have distinct values
            if let Some(original) = &aa_seq.original_sequence && aa_new_vec.len() > 0 {
                if aa_new_vec.len() > 1 {
                    tracing::info!(
                        "{},Variant [{}] cannot be mapped: Alternative sequence is ambiguous",
                        accession,
                        isoform_id
                    );     
                    return;
                }

                //unwrap safe because of len() check
                aa_new = Some(aa_new_vec.first().cloned().unwrap());

                aa_ref = Some(original);
            }
        } else {
            tracing::info!(
                "{},Variant [{}] cannot be mapped: Malformed format",
                accession,
                isoform_id
            );     
            return;
        }
        uniprot_variants.push(
            Variant {
                isoform: iso_id,
                id: feature.feature_id.clone().unwrap_or(String::new()),
                begin: begin,
                end: end,
                aa_ref: aa_ref.cloned(),
                aa_new: aa_new,
            }
        );
    }
    else {
        tracing::info!(
            "{},Variant [{}] cannot be mapped: No position",
                accession,
                isoform_id
            )
    }
}
/// Collect isoform IDs, reconstruct their sequences, and collect VAR_SEQ
/// features, in a single pass over the ALTERNATIVE PRODUCTS comment(s).
///
/// Returns `(displayed_id, isoform_map, vsp_features, uniprot_variants)`
/// where `isoform_map` maps the canonical id to the list of reconstructed
/// isoforms (each of which already carries its own `UniProtIsoId`). Any
/// isoform that can't be reconstructed is skipped and logged via
/// `tracing::info!` rather than failing the whole entry.
pub(crate) fn collect_and_reconstruct_isoforms(
    entry: &UniProtEntry
) -> Result<(
    Option<UniProtIsoId>,
    HashMap<UniProtIsoId, Sequence>,
    Vec<UniProtFeature>,
    Vec<Variant>,
)> {
    let mut var_seq_features: Vec<UniProtFeature> = Vec::new();
    let mut uniprot_variants: Vec<Variant> = Vec::new();

    for f in &entry.features {
        if f.feature_type == "Alternative sequence" {
            var_seq_features.push(f.clone());
        } 
        else if f.feature_type == "Natural variant" {
            add_canonical_variant(&entry.primary_accession, f, &mut uniprot_variants);
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

            let isoform_id = match isoform_id {
                Some(id) => id,
                None => {
                    tracing::info!("isoform entry has no isoform id; skipping isoform");
                    continue;
                }
            };

            if let Some((id, seq)) = reconstruct_isoform(
                entry,
                isoform_id,
                iso.isoform_sequence_status.as_deref().unwrap_or(""),
                &iso.sequence_ids,
                &var_seq_features,
            ) {
                isoforms.insert(id, seq);
            }
        }
    }

    Ok((displayed_id, isoforms, var_seq_features, uniprot_variants))
}

#[cfg(test)]
#[path = "isoform_tests.rs"]
mod tests;
