use super::*;
use crate::uniprot::{IsoformInfo, NameValue, UniProtComment, UniProtCrossReference, UniProtSequence};

// ── "Displayed" isoform aliases the canonical sequence ────────────────────

fn entry_with_displayed_isoform() -> UniProtEntry {
    UniProtEntry {
        primary_accession: "P04637".to_string(),
        protein_description: None,
        genes: None,
        sequence: UniProtSequence { value: "ACDEFGHIKL".to_string() },
        features: Vec::new(),
        comments: vec![UniProtComment {
            comment_type: "ALTERNATIVE PRODUCTS".to_string(),
            isoforms: vec![
                IsoformInfo {
                    name: NameValue { value: "1".to_string() },
                    synonyms: Vec::new(),
                    isoform_ids: vec!["P04637-1".to_string()],
                    sequence_ids: Vec::new(),
                    isoform_sequence_status: Some("Displayed".to_string()),
                },
                IsoformInfo {
                    name: NameValue { value: "2".to_string() },
                    synonyms: Vec::new(),
                    isoform_ids: vec!["P04637-2".to_string()],
                    sequence_ids: vec!["VSP_0001".to_string()],
                    isoform_sequence_status: Some("Described".to_string()),
                },
            ],
            events: vec!["Alternative splicing".to_string()],
            note: None,
        }],
        cross_references: vec![
            // Explicitly cross-referenced to the "Displayed" isoform — should
            // end up routed to canonical, not to a dead "P04637-1" bucket.
            UniProtCrossReference {
                database: "Ensembl".to_string(),
                id: "ENST00000000001".to_string(),
                isoform_id: Some("P04637-1".to_string()),
            },
            UniProtCrossReference {
                database: "Ensembl".to_string(),
                id: "ENST00000000002".to_string(),
                isoform_id: Some("P04637-2".to_string()),
            },
        ],
    }
}

#[test]
fn displayed_isoform_gets_no_synthetic_entry() {
    let entry = entry_with_displayed_isoform();
    let isoforms = reconstruct_isoforms(&entry);
    assert_eq!(isoforms.len(), 1);
    assert_eq!(isoforms[0].0.as_str(), "P04637-2");
}

#[test]
fn enst_cross_referenced_to_displayed_isoform_routes_to_canonical() {
    let entry = entry_with_displayed_isoform();
    let canonical = UniprotId("P04637".to_string());

    let raw = extract_enst_to_isoform(&entry, &canonical);
    // Pre-redirect: still shows the literal "-1" id (this is what
    // check_enst_isoform_coverage must see to avoid a false-positive).
    assert_eq!(
        raw.get(&EnsemblId("ENST00000000001".to_string())),
        Some(&UniprotId("P04637-1".to_string()))
    );

    let canonicalized = canonicalize_displayed_isoform(raw, &entry, &canonical);
    assert_eq!(
        canonicalized.get(&EnsemblId("ENST00000000001".to_string())),
        Some(&canonical)
    );
    assert_eq!(
        canonicalized.get(&EnsemblId("ENST00000000002".to_string())),
        Some(&UniprotId("P04637-2".to_string()))
    );
}
