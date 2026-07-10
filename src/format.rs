// ============================================================================
// FT / CC formatting and synthetic flat-file entry generation
// ============================================================================

use crate::types::{Sequence, UniprotId, Variant};
use crate::uniprot::{EvidenceSource, UniProtComment, UniProtEntry, UniProtEvidence, UniProtFeature};
use std::collections::HashSet;

/// Column prefix for all FT qualifier continuation lines (21 chars).
const FT_INDENT: &str = "FT                   ";
/// Max characters of content per FT line after the 21-char prefix (79 - 21 = 58).
const FT_CONTENT_WIDTH: usize = 58;

/// 1-based inclusive `begin..end` slice of `seq`, or `None` if out of range.
fn residues_at(seq: &str, begin: usize, end: usize) -> Option<String> {
    if begin == 0 || begin > end || end > seq.len() {
        return None;
    }
    seq.get(begin - 1..end).map(str::to_string)
}

/// `isoform_ref`, when set, is written as an "ACCESSION:position" cross-reference
/// instead of a bare position — biopython's feature-location parser (which
/// ProtGraph is built on) reads that form as a remote reference, routing the
/// feature onto that isoform's own vertex chain. See the `Variant` doc comment.
fn format_ft_position_line(key: &str, isoform_ref: Option<&str>, begin: usize, end: usize) -> String {
    let pos = if begin == end {
        format!("{}", begin)
    } else {
        format!("{}..{}", begin, end)
    };
    let pos = match isoform_ref {
        Some(accession) => format!("{}:{}", accession, pos),
        None => pos,
    };
    format!("FT   {:<16}{}", key, pos)
}

/// Wrap a `/qualifier="value"` across multiple FT continuation lines, breaking
/// at word boundaries where possible.
fn wrap_ft_qualifier(qualifier: &str, value: &str) -> Vec<String> {
    let full = format!("/{}=\"{}\"", qualifier, value);
    let mut lines: Vec<String> = Vec::new();
    let chars: Vec<char> = full.chars().collect();
    let mut pos = 0;

    while pos < chars.len() {
        let chunk_end = (pos + FT_CONTENT_WIDTH).min(chars.len());
        let break_at = if chunk_end < chars.len() {
            chars[pos..chunk_end]
                .iter()
                .rposition(|&c| c == ' ')
                .map(|i| pos + i)
                .unwrap_or(chunk_end)
        } else {
            chunk_end
        };
        let chunk: String = chars[pos..break_at].iter().collect();
        lines.push(format!("{}{}", FT_INDENT, chunk));
        pos = break_at;
        while pos < chars.len() && chars[pos] == ' ' {
            pos += 1;
        }
    }

    lines
}

fn format_uniprot_evidence_string(evidences: &[UniProtEvidence]) -> String {
    evidences
        .iter()
        .map(|e| match &e.source {
            Some(EvidenceSource::Full { name, id }) => format!("{}|{}:{}", e.evidence_code, name, id),
            Some(EvidenceSource::NameOnly(name)) => format!("{}|{}", e.evidence_code, name),
            None => e.evidence_code.clone(),
        })
        .collect::<Vec<_>>()
        .join(", ")
}

/// Format a UniProt entry feature (Alternative sequence or Natural variant)
/// into FT lines.
fn format_ft_from_uniprot_feature(feature: &UniProtFeature) -> Vec<String> {
    let key = match feature.feature_type.as_str() {
        "Alternative sequence" => "VAR_SEQ",
        "Natural variant" => "VARIANT",
        _ => return Vec::new(),
    };

    let begin = feature.location.start.value;
    let end = feature.location.end.value;

    let mut lines = Vec::new();
    lines.push(format_ft_position_line(key, None, begin, end));

    let note = match &feature.alternative_sequence {
        Some(alt) => {
            let desc = feature.description.as_deref().unwrap_or("");
            match (
                alt.original_sequence.as_deref(),
                alt.alternative_sequences.first(),
            ) {
                (Some(orig), Some(alt_seq)) => {
                    if desc.is_empty() {
                        format!("{} -> {}", orig, alt_seq)
                    } else {
                        format!("{} -> {} ({})", orig, alt_seq, desc)
                    }
                }
                _ => {
                    if desc.is_empty() {
                        "Missing".to_string()
                    } else {
                        format!("Missing ({})", desc)
                    }
                }
            }
        }
        None => feature.description.clone().unwrap_or_default(),
    };

    if !note.is_empty() {
        lines.extend(wrap_ft_qualifier("note", &note));
    }

    if !feature.evidences.is_empty() {
        lines.extend(wrap_ft_qualifier(
            "evidence",
            &format_uniprot_evidence_string(&feature.evidences),
        ));
    }

    // /id: prefer featureId, then cross-references, then evidence string
    let id_str = feature
        .feature_id
        .as_deref()
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .or_else(|| {
            if !feature.feature_cross_references.is_empty() {
                Some(
                    feature
                        .feature_cross_references
                        .iter()
                        .map(|xr| format!("{}:{}", xr.database, xr.id))
                        .collect::<Vec<_>>()
                        .join("|"),
                )
            } else if !feature.evidences.is_empty() {
                Some(format_uniprot_evidence_string(&feature.evidences))
            } else {
                None
            }
        });

    if let Some(id) = id_str {
        lines.extend(wrap_ft_qualifier("id", &id));
    }

    lines
}

/// Format a joined Variant (EBI + ENST) into FT lines.
///
/// Always keyed "VARIANT", regardless of span — never "VAR_SEQ". ProtGraph
/// dispatches by FT type to two unrelated handlers: "VAR_SEQ" goes to
/// execute_var_seq, which finds isoform parentage by parsing "(in isoform X)"
/// out of the note text (and never looks at the location's cross-reference);
/// "VARIANT" (also MUTAGEN/CONFLICT) goes to _execute_generic_feature, which
/// uses the location's cross-reference (our `isoform_ref`) and ignores the note
/// text for that purpose. A joined variant's note is plain "X -> Y"/"Missing"
/// with no isoform clause, so keying it "VAR_SEQ" would send it to the wrong
/// handler — it isn't a splice-defining feature (that's exclusively UniProt's
/// native "Alternative sequence" type, handled separately by
/// format_ft_from_uniprot_feature), so "VARIANT" is correct at any span length.
fn format_ft_from_variant(variant: &Variant) -> Vec<String> {
    let mut lines = Vec::new();
    let isoform_ref = variant.isoform_ref.as_ref().map(UniprotId::as_str);
    lines.push(format_ft_position_line("VARIANT", isoform_ref, variant.begin, variant.end));
    let note = if variant.replacement.is_empty() {
        "Missing".to_string()
    } else {
        format!("{} -> {}", variant.replaced, variant.replacement)
    };
    lines.extend(wrap_ft_qualifier("note", &note));
    if !variant.id.is_empty() {
        lines.extend(wrap_ft_qualifier("id", &variant.id));
    }
    lines
}

fn format_cc_alternative_products(comment: &UniProtComment) -> String {
    let mut lines = Vec::new();
    lines.push("CC   -!- ALTERNATIVE PRODUCTS:".to_string());
    lines.push(format!(
        "CC       Event={}; Named isoforms={};",
        comment.events.join(", "),
        comment.isoforms.len()
    ));

    if let Some(note) = &comment.note {
        for text in note.texts() {
            lines.push(format!("CC         Comment={}", text));
        }
    }

    for iso in &comment.isoforms {
        let synonyms: Vec<&str> = iso.synonyms.iter().map(|s| s.value.as_str()).collect();
        if synonyms.is_empty() {
            lines.push(format!("CC       Name={};", iso.name.value));
        } else {
            lines.push(format!(
                "CC       Name={}; Synonyms={};",
                iso.name.value,
                synonyms.join(", ")
            ));
        }
        let iso_id = iso.isoform_ids.first().map(String::as_str).unwrap_or("");
        let seq_str = match iso.isoform_sequence_status.as_deref() {
            Some("Displayed") => "Displayed".to_string(),
            Some("Not described") => "Not described".to_string(),
            _ => {
                if iso.sequence_ids.is_empty() {
                    "Not described".to_string()
                } else {
                    iso.sequence_ids.join(", ")
                }
            }
        };
        lines.push(format!(
            "CC         IsoId={}; Sequence={};",
            iso_id, seq_str
        ));
    }

    lines.join("\n")
}

fn format_seq_body(seq: &str) -> String {
    seq.as_bytes()
        .chunks(60)
        .map(|chunk| {
            let line = std::str::from_utf8(chunk).unwrap_or("");
            let grouped = line
                .as_bytes()
                .chunks(10)
                .map(|g| std::str::from_utf8(g).unwrap_or(""))
                .collect::<Vec<_>>()
                .join(" ");
            format!("     {}", grouped)
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// Build the synthetic Swiss-Prot-style flat file entry text.
/// `entry` is `Some` for UniProt-backed accessions, `None` for synthetic ENST entries.
pub(crate) fn format_entry(
    accession: &UniprotId,
    sequence: &Sequence,
    variants: &[Variant],
    entry: Option<&UniProtEntry>,
) -> String {
    let date = "01-JAN-1970";

    let (protein_desc, gene) = match entry {
        Some(e) => (e.protein_name().to_string(), e.gene_name().to_string()),
        None => (accession.as_str().to_string(), accession.as_str().to_string()),
    };

    let id_ac = format!(
        "ID   {:<24}Reviewed;         {} AA.\nAC   {};",
        accession.as_str(),
        sequence.as_str().len(),
        accession.as_str()
    );

    let mut ft_lines: Vec<String> = Vec::new();

    // VAR_SEQ and Natural variant features from the UniProt entry. `accession` is
    // always the canonical accession now — there's exactly one entry per UniProt
    // record — so every feature is emitted unconditionally, at its native
    // (canonical) position, exactly as UniProt's own flat file does. Isoform
    // applicability is resolved downstream by ProtGraph itself: VAR_SEQ notes
    // already carry "(in isoform X)" verbatim, and a bare canonical position on a
    // Natural variant matches that same position on every isoform chain that
    // retained the residue, so canonical-only variants reach isoforms for free
    // without any remapping here.
    let mut uniprot_variants: HashSet<Variant> = HashSet::new();
    let mut seen_native_feature_ids: HashSet<&str> = HashSet::new();
    if let Some(e) = entry {
        for feature in e
            .features
            .iter()
            .filter(|f| matches!(f.feature_type.as_str(), "Alternative sequence" | "Natural variant"))
        {
            // A feature id already seen means UniProt's own JSON repeated the
            // same feature (has happened for some entries) — render it once.
            if let Some(id) = feature.feature_id.as_deref().filter(|s| !s.is_empty()) {
                if !seen_native_feature_ids.insert(id) {
                    continue;
                }
            }
            let begin = feature.location.start.value;
            let end = feature.location.end.value;
            uniprot_variants.insert(Variant {
                id: String::new(),
                begin,
                end,
                // UniProt sometimes omits `alternativeSequence` entirely for
                // older annotations (an empty `{}`, no original residues
                // recorded at all) even though the position is present. Fall
                // back to reading the actual residues off the canonical
                // sequence so this still matches an EBI-sourced record of the
                // same amino-acid change instead of comparing "" to "V".
                replaced: feature
                    .alternative_sequence
                    .as_ref()
                    .and_then(|a| a.original_sequence.clone())
                    .or_else(|| residues_at(sequence.as_str(), begin, end))
                    .unwrap_or_default(),
                replacement: feature
                    .alternative_sequence
                    .as_ref()
                    .and_then(|a| a.alternative_sequences.first().cloned())
                    .unwrap_or_default(),
                isoform_ref: None,
            });
            ft_lines.extend(format_ft_from_uniprot_feature(feature));
        }
    }

    // Joined variants (EBI + ENST) — skip those already covered by a UniProt
    // feature, and skip exact repeats within `variants` itself (EBI's
    // variation feed can list the same amino-acid change twice under
    // different sourceTypes). Both checks rely purely on position + residue
    // change (`Variant`'s `Eq`/`Hash`, which normalize "del"/"*" to "" so a
    // deletion compares equal across sources) — reconciling variants from
    // different sources is the whole point here, so provenance (ids) plays
    // no part in the match.
    let mut seen_joined: HashSet<&Variant> = HashSet::new();
    for variant in variants {
        if !uniprot_variants.contains(variant) && seen_joined.insert(variant) {
            ft_lines.extend(format_ft_from_variant(variant));
        }
    }

    let ft_block = if ft_lines.is_empty() {
        String::new()
    } else {
        format!("\n{}", ft_lines.join("\n"))
    };

    let cc_alt_products = entry
        .and_then(|e| {
            e.comments
                .iter()
                .find(|c| c.comment_type == "ALTERNATIVE PRODUCTS")
        })
        .map(|c| format!("{}\n", format_cc_alternative_products(c)))
        .unwrap_or_default();

    let sq_line = format!("SQ   SEQUENCE   {} AA;  0 MW;  0000000000000000 CRC64;", sequence.as_str().len());
    let seq_body = format_seq_body(sequence.as_str());
    let acc = accession.as_str();

    format!(
        "{id_ac}
DT   {date}, integrated into UTA.
DT   {date}, sequence version 1.
DT   {date}, entry version 1.
DE   RecName: Full={protein_desc};
GN   Name={gene}; Synonyms={gene}; ORFNames={gene};
OS   Autogenerated.
OC   Autogenerated.
OX   NCBI_TaxID=0000;
RN   [1]
RP   AUTOGENERATED
RA   This file was autogenerated;
RT   \"This is an autogenerated file.\";
RL   (er) Autogenerated.
CC   -!- CAUTION: This file was autogenerated and sequences
CC       can differ from UniProt!
{cc_alt_products}DR   UTA; {acc}; Universal Transcript Archive
PE   4: Predicted;
KW   Generated for ProtGraph{ft_block}
{sq_line}
{seq_body}
//"
    )
}

#[cfg(test)]
#[path = "format_tests.rs"]
mod tests;
