use super::*;
use crate::isoform::reconstruct_isoforms;
use crate::uniprot::{
    AlternativeSequence, FeatureLocation, IsoformInfo, NameValue, PositionValue, UniProtSequence,
};

fn natural_variant(
    pos: usize,
    orig: &str,
    alt: &str,
    description: &str,
    feature_id: &str,
) -> UniProtFeature {
    UniProtFeature {
        feature_type: "Natural variant".to_string(),
        location: FeatureLocation { start: PositionValue { value: pos }, end: PositionValue { value: pos } },
        description: Some(description.to_string()),
        feature_id: Some(feature_id.to_string()),
        alternative_sequence: Some(AlternativeSequence {
            original_sequence: Some(orig.to_string()),
            alternative_sequences: vec![alt.to_string()],
        }),
        evidences: Vec::new(),
        feature_cross_references: Vec::new(),
    }
}

fn var_seq(start: usize, end: usize, original: &str, feature_id: &str) -> UniProtFeature {
    UniProtFeature {
        feature_type: "Alternative sequence".to_string(),
        location: FeatureLocation { start: PositionValue { value: start }, end: PositionValue { value: end } },
        description: None,
        feature_id: Some(feature_id.to_string()),
        alternative_sequence: Some(AlternativeSequence {
            original_sequence: Some(original.to_string()),
            alternative_sequences: Vec::new(), // empty = deletion ("Missing")
        }),
        evidences: Vec::new(),
        feature_cross_references: Vec::new(),
    }
}

fn isoform_info(name: &str, iso_id: &str, vsp_id: &str) -> IsoformInfo {
    IsoformInfo {
        name: NameValue { value: name.to_string() },
        synonyms: Vec::new(),
        isoform_ids: vec![iso_id.to_string()],
        sequence_ids: vec![vsp_id.to_string()],
        isoform_sequence_status: Some("Described".to_string()),
    }
}

/// A truncated UniProt entry (no comments/features beyond what's under test) with
/// two isoforms and four Natural variant notes. All four now land on the single
/// canonical entry unconditionally, whatever their note text says — ProtGraph
/// itself resolves isoform applicability downstream (VAR_SEQ notes for splicing,
/// position-matching across vertex chains for Natural variants), not us.
fn build_test_entry() -> UniProtEntry {
    let canonical_seq = "ACDEFGHIKLMNPQRSTVWY".to_string(); // 20 residues, position == index

    UniProtEntry {
        primary_accession: "TEST1".to_string(),
        protein_description: None,
        genes: None,
        sequence: UniProtSequence { value: canonical_seq },
        features: vec![
            var_seq(5, 7, "FGH", "VSP_0001"),  // defines isoform 2
            var_seq(15, 17, "RST", "VSP_0002"), // defines isoform 3
            natural_variant(2, "C", "X", "in a disease", "VAR_1001"), // canonical only
            natural_variant(3, "D", "Z", "in isoform 2", "VAR_1002"), // isoform 2 only
            natural_variant(1, "A", "Q", "in isoform 2 and isoform 3", "VAR_1003"), // both
            natural_variant(
                4,
                "E",
                "M",
                "affects interaction between isoform 2 and isoform 3 complexes",
                "VAR_1004",
            ), // prose mention only — no "(in isoform X)" clause → canonical
        ],
        comments: vec![UniProtComment {
            comment_type: "ALTERNATIVE PRODUCTS".to_string(),
            isoforms: vec![
                isoform_info("2", "TEST1-2", "VSP_0001"),
                isoform_info("3", "TEST1-3", "VSP_0002"),
            ],
            events: vec!["Alternative splicing".to_string()],
            note: None,
        }],
        cross_references: Vec::new(),
    }
}

#[test]
fn all_native_features_land_unfiltered_on_the_single_canonical_entry() {
    let entry = build_test_entry();

    // Isoform sequences are still reconstructed — needed elsewhere in the
    // pipeline to compute isoform-local positions for EBI/ENST variants — but
    // no longer get a synthetic entry of their own.
    let isoforms = reconstruct_isoforms(&entry);
    assert_eq!(isoforms.len(), 2);
    let iso2 = isoforms.iter().find(|iso| iso.0.as_str() == "TEST1-2").unwrap();
    let iso3 = isoforms.iter().find(|iso| iso.0.as_str() == "TEST1-3").unwrap();
    assert_eq!(iso2.1.as_str(), "ACDEIKLMNPQRSTVWY");
    assert_eq!(iso3.1.as_str(), "ACDEFGHIKLMNPQVWY");

    let canonical_id = UniprotId("TEST1".to_string());
    let canonical_seq = Sequence(entry.sequence.value.clone());

    // A variant sourced from EBI/ENST specifically for isoform 2, tagged with
    // isoform_ref the way main.rs's pipeline does when folding isoform-specific
    // variants onto the single canonical entry.
    let iso_specific = Variant {
        id: String::new(),
        begin: 10,
        end: 10,
        replaced: "L".to_string(),
        replacement: "V".to_string(),
        isoform_ref: Some(UniprotId("TEST1-2".to_string())),
    };

    let out = format_entry(&canonical_id, &canonical_seq, &[iso_specific], Some(&entry));

    // Every Natural variant note lands here unconditionally, at its native
    // canonical position, regardless of what its note text says — exactly like
    // a real UniProt flat file. ProtGraph resolves isoform applicability itself:
    // a bare canonical position matches every vertex chain (canonical's and any
    // isoform's) that retained that residue.
    assert!(out.contains("VAR_1001"));
    assert!(out.contains("VAR_1002"));
    assert!(out.contains("VAR_1003"));
    assert!(out.contains("VAR_1004"));

    // VAR_SEQ features are present with their native note/id — ProtGraph's
    // execute_var_seq parses the "(in isoform X)" clause straight out of that
    // note text to build each isoform's own vertex chain.
    assert!(out.contains("FT   VAR_SEQ"));
    assert!(out.contains("/id=\"VSP_0001\""));
    assert!(out.contains("/id=\"VSP_0002\""));

    // The ALTERNATIVE PRODUCTS comment survives — ProtGraph's
    // _get_isoforms_of_entry needs it to know each isoform's own accession.
    assert!(out.contains("CC   -!- ALTERNATIVE PRODUCTS:"));
    assert!(out.contains("IsoId=TEST1-2;"));
    assert!(out.contains("IsoId=TEST1-3;"));

    // The isoform-specific joined variant is written with an isoform-qualified
    // position ("TEST1-2:10"), not a bare "10" — biopython's FT location parser
    // reads "ACCESSION:position" as a remote cross-reference, which is what
    // routes this feature onto isoform 2's own vertex chain instead of canonical's.
    assert!(out.contains("TEST1-2:10"));
    assert!(out.contains("/note=\"L -> V\""));
}

#[test]
fn multi_residue_joined_variant_is_keyed_variant_not_var_seq() {
    // A stop-codon truncation (or any multi-AA substitution/deletion) spans
    // begin != end, but it must still be keyed "VARIANT". "VAR_SEQ" is reserved
    // for UniProt's native "Alternative sequence" feature type — a different FT
    // key would route this to ProtGraph's execute_var_seq, which parses isoform
    // parentage from an "(in isoform X)" note clause we don't have (a joined
    // variant's note is plain "X -> Y"/"Missing"), instead of
    // _execute_generic_feature, which reads the isoform_ref cross-reference.
    let entry = build_test_entry();
    let canonical_id = UniprotId("TEST1".to_string());
    let canonical_seq = Sequence(entry.sequence.value.clone());

    let truncation = Variant {
        id: String::new(),
        begin: 10,
        end: 20,
        replaced: "LMNPQRSTVWY".to_string(),
        replacement: String::new(),
        isoform_ref: Some(UniprotId("TEST1-2".to_string())),
    };

    let out = format_entry(&canonical_id, &canonical_seq, &[truncation], Some(&entry));

    assert!(out.contains("FT   VARIANT         TEST1-2:10..20"));
    assert!(!out.contains("FT   VAR_SEQ         TEST1-2:10..20"));
    assert!(out.contains("/note=\"Missing\""));
}
