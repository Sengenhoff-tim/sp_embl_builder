// ============================================================================
// FT / CC formatting and synthetic flat-file entry generation
// ============================================================================

use crate::types::{Sequence, UniProtIsoId, UniProtCanonId};
use crate::variant::Variant;
use crate::uniprot::{EvidenceSource, UniProtComment, UniProtEntry, UniProtEvidence, UniProtFeature};
use std::collections::HashSet;
use std::ffi::OsString;

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
 
fn format_ft_position_line(
    key: &str,
    isoform: &Option<UniProtIsoId>,
    begin: usize,
    end: usize,
) -> String {
    let range = if begin == end {
        begin.to_string()
    } else {
        format!("{}..{}", begin, end)
    };
    
    let pos = isoform
        .as_ref()
        .map(|iso| format!("{}:{}", iso.as_str(), range))
        .unwrap_or(range);
    
    format!("FT   {:<16}{}", key, pos)
}


fn wrap_ft_qualifier(qualifier: &str, value: &str) -> Vec<String> {
    let full = format!("/{}=\"{}\"", qualifier, value);
    let mut lines: Vec<String> = Vec::new();
    let chars: Vec<char> = full.chars().collect();
    let mut pos = 0;

    while pos < chars.len() {
        let chunk_end = (pos + FT_CONTENT_WIDTH).min(chars.len());
        let chunk: String = chars[pos..chunk_end].iter().collect();
        lines.push(format!("{}{}", FT_INDENT, chunk));
        pos = chunk_end;
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

/// Format a UniProt Alternative sequence feature into VAR_SEQ FT lines.
fn format_ft_from_uniprot_var_seq(feature: &UniProtFeature) -> Option<Vec<String>> {
    let Some(begin) = feature.location.start.value else {
        return None;
    };
    let end = feature.location.end.value.unwrap_or(begin);
    
    let mut lines = Vec::new();
    lines.push(format_ft_position_line("VAR_SEQ", &None, begin, end));

    let note = match &feature.alternative_sequence {
        Some(alt) => {
            let desc = feature.description.as_deref().unwrap_or("");
            let orig = &alt.original_sequence;
            let alt_seq = alt.alternative_sequences.first();

            match alt_seq {
                Some(alt_seq) => {
                    if desc.is_empty() {
                        format!("{} -> {}", orig.clone().unwrap_or("".to_string()), alt_seq)
                    } else {
                        format!("{} -> {} ({})", orig.clone().unwrap_or("".to_string()), alt_seq, desc)
                    }
                }
                None => {
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

    Some(lines)
}


/// Format a Variant into FT lines.
fn format_variant(variant: &Variant) -> Vec<String> {
    let mut lines = Vec::new();
    lines.push(format_ft_position_line("VARIANT", &variant.isoform, variant.begin, variant.end));
    let note = if variant.aa_new.is_empty() {
        "Missing".to_string()
    } else {
        format!("{} -> {}", variant.aa_ref, variant.aa_new)
    };
    lines.extend(wrap_ft_qualifier("note", &note));
    lines.extend(wrap_ft_qualifier("id", &variant.id));
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
        if let Some(iname) = iso.name.as_ref() {
            if synonyms.is_empty() {
                {
                    lines.push(format!("CC       Name={};", iname.value));
                }
            } else {
                lines.push(format!(
                    "CC       Name={}; Synonyms={};",
                    iname.value,
                    synonyms.join(", ")
                ));
            }
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
    accession: &UniProtCanonId,
    sequence: &str,
    var_seq: Option<&[UniProtFeature]>,
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
        sequence.len(),
        accession.as_str()
    );

    let mut ft_lines: Vec<String> = Vec::new();

    // VAR_SEQ features
    if let Some(feats) = var_seq {
        for feature in feats {
            if let Some(var_seq) = format_ft_from_uniprot_var_seq(feature) {
                ft_lines.extend(var_seq);
            }
        }
    }
    

    // VARIANT features 
    for variant in variants {
        ft_lines.extend(format_variant(variant));
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

    let sq_line = format!("SQ   SEQUENCE   {} AA;  0 MW;  0000000000000000 CRC64;", sequence.len());
    let seq_body = format_seq_body(sequence);
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
