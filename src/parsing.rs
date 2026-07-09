// ============================================================================
// Parsing of the accession list and ENST variant input files
// ============================================================================

use crate::exceptions::ExceptionLog;
use crate::types::{EnsemblId, UniprotId, Variant};
use crate::util::strip_version;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs;

pub(crate) fn parse_accession_list(path: &str) -> Result<Vec<UniprotId>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("failed to read accession list '{}'", path))?;
    Ok(content
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty())
        .map(|l| UniprotId(l.to_string()))
        .collect())
}

fn split_position_and_residues(s: &str) -> Option<(usize, &str)> {
    let digit_count = s.chars().take_while(|c| c.is_ascii_digit()).count();
    if digit_count == 0 {
        return None;
    }
    let (pos_str, residues) = s.split_at(digit_count);
    let pos: usize = pos_str.parse().ok()?;
    Some((pos, residues))
}

fn parse_variant(variant_str: &str) -> Option<Variant> {
    let (left, right) = variant_str.split_once('>')?;
    let (begin, raw_wt) = split_position_and_residues(left)?;
    let (_rpos, raw_mt) = split_position_and_residues(right)?;
    // Stop variants (* in wt or mt) are resolved later with the sequence;
    // end is naive here and will be corrected by resolve_stop_variants.
    let end = begin + raw_wt.len().checked_sub(1)?;
    Some(Variant {
        id: String::new(),
        begin,
        end,
        replaced: raw_wt.to_string(),
        replacement: raw_mt.to_string(),
        isoform_ref: None,
    })
}

/// Atomically strip stop markers and correct begin/end/replaced/replacement using
/// the target sequence. Must be called before join_variants so matching sees clean
/// variants.
///   replaced == "*"          → 1*>1BB*: begin=len, end=len, replaced=last, replacement=last+BB
///   replacement ends with '*' → 1A>1* or 1A>1BB*: end=len, replaced expands to suffix, replacement stripped
pub(crate) fn resolve_stop_variants(
    entries: Option<&Vec<(EnsemblId, usize, Variant)>>,
    seq: &str,
) -> Option<Vec<(EnsemblId, usize, Variant)>> {
    entries.map(|e| {
        e.iter()
            .map(|(enst, idx, v)| {
                let v = if v.replaced == "*" {
                    let last = seq.chars().last().unwrap_or_default().to_string();
                    let mt = v.replacement.trim_end_matches('*');
                    Variant {
                        id: v.id.clone(),
                        begin: seq.len(),
                        end: seq.len(),
                        replaced: last.clone(),
                        replacement: format!("{}{}", last, mt),
                        isoform_ref: v.isoform_ref.clone(),
                    }
                } else if v.replacement.ends_with('*') {
                    let suffix_start = v.begin - 1 + v.replaced.len();
                    let suffix = seq.get(suffix_start..).unwrap_or("");
                    Variant {
                        id: v.id.clone(),
                        begin: v.begin,
                        end: seq.len(),
                        replaced: format!("{}{}", v.replaced, suffix),
                        replacement: v.replacement.trim_end_matches('*').to_string(),
                        isoform_ref: v.isoform_ref.clone(),
                    }
                } else {
                    v.clone()
                };
                (enst.clone(), *idx, v)
            })
            .collect()
    })
}

pub(crate) fn parse_variants(
    path: &str,
    exceptions: &mut ExceptionLog,
) -> Result<HashMap<EnsemblId, Vec<Variant>>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("failed to read variants file '{}'", path))?;
    let mut variants: HashMap<EnsemblId, Vec<Variant>> = HashMap::new();
    for (line_no, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let identifier = format!("{}:{}", path, line_no + 1);

        let mut parts = line.split_whitespace();
        let enst = match parts.next() {
            Some(e) => e,
            None => continue,
        };
        if !enst.starts_with("ENST") {
            exceptions.log(&identifier, &format!("'{}' does not look like an ENST id; skipping line", enst));
            continue;
        }
        let variant_string = match parts.next() {
            Some(v) => v,
            None => {
                exceptions.log(&identifier, &format!("missing variant field for {}", enst));
                continue;
            }
        };
        if variant_string == "." {
            // Explicit "no variant" marker — not malformed input.
            continue;
        }
        let variant = match parse_variant(variant_string) {
            Some(v) => v,
            None => {
                exceptions.log(
                    &identifier,
                    &format!("could not parse variant '{}' for {}", variant_string, enst),
                );
                continue;
            }
        };
        variants
            .entry(EnsemblId(strip_version(enst).to_string()))
            .or_default()
            .push(variant);
    }
    Ok(variants)
}

#[cfg(test)]
#[path = "parsing_tests.rs"]
mod tests;
