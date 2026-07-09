// ============================================================================
// Join logic: merge EBI-sourced and ENST-derived variants for one accession
// ============================================================================

use crate::types::{EnsemblId, Variant};
use std::collections::HashMap;

/// Append a reference tag to a Variant's id, using `|` as separator.
pub(crate) fn append_id_tag(id: &mut String, tag: &str) {
    if id.is_empty() {
        *id = tag.to_string();
    } else {
        id.push('|');
        id.push_str(tag);
    }
}

pub(crate) fn join_variants(
    rest_variants: Vec<Variant>,
    enst_entries: Option<&Vec<(EnsemblId, usize, Variant)>>,
) -> Vec<Variant> {
    let mut joined = rest_variants;
    let mut index_by_variant: HashMap<Variant, usize> = joined
        .iter()
        .enumerate()
        .map(|(i, v)| (v.clone(), i))
        .collect();

    if let Some(entries) = enst_entries {
        for (enst_id, idx, variant) in entries {
            let tag = format!("{}:{}", enst_id.as_str(), idx);
            if let Some(&existing_idx) = index_by_variant.get(variant) {
                append_id_tag(&mut joined[existing_idx].id, &tag);
            } else {
                joined.push(variant.clone());
                let new_idx = joined.len() - 1;
                index_by_variant.insert(variant.clone(), new_idx);
            }
        }
    }

    joined
}
