// ============================================================================
// Parsing of the accession list and ENST variant input files
// ============================================================================

use crate::types::{EnsemblId, UniProtCanonId};
use crate::variant::Variant;
use crate::util::strip_version;
use anyhow::{Context, Result};
use tracing::info;
use std::collections::HashMap;
use std::fs;

pub(crate) fn parse_accession_list(path: &str) -> Result<Vec<UniProtCanonId>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("failed to read accession list '{}'", path))?;
    Ok(content
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty())
        .map(|l| UniProtCanonId(l.to_string()))
        .collect())
}

/// Parses `<pos><ref_aa>` style position/residue notation.
/// e.g. "726G" -> (726, "G")
fn split_position_and_residues(s: &str) -> Option<(usize, &str)> {
    let digit_count = s.chars().take_while(|c| c.is_ascii_digit()).count();
    if digit_count == 0 {
        return None;
    }
    let (pos_str, residues) = s.split_at(digit_count);
    let pos: usize = pos_str.parse().ok()?;
    Some((pos, residues))
}

/// Checks that the position on the left of `>` matches the position on the right.
/// e.g. "726G>726S" -> Some(true) (positions match: 726 == 726)
/// e.g. "726G>731S" -> Some(false) (positions differ: 726 != 731)
/// e.g. "garbage!!" -> None (couldn't parse either side)
fn positions_match(variant_str: &str) -> Option<bool> {
    let (left, right) = variant_str.split_once('>')?;
    let (begin, _) = split_position_and_residues(left)?;
    let (rpos, _) = split_position_and_residues(right)?;
    Some(begin == rpos)
}

/// Parses a protein variant string of the form `<pos><ref_aa>><pos><alt_aa>`.
/// e.g. "726G>726S" -> Some(Variant { begin: 726, end: 726, aa_ref: "G", aa_new: "S" })
fn to_variant(id: &str, variant_str: &str) -> Option<Variant> {
    let (left, right) = variant_str.split_once('>')?;
    let (begin, raw_wt) = split_position_and_residues(left)?;
    let (_rpos, raw_mt) = split_position_and_residues(right)?;
    let end = begin + raw_wt.len().checked_sub(1)?;
    Some(Variant {
        isoform: None,
        id: id.to_string(),
        begin,
        end,
        aa_ref: Some(raw_wt.to_string()),
        aa_new: Some(raw_mt.to_string()),
    })
}

/// parse a ensembl transcript mapped variant input file, returning a map of Ensembl-ID → list of variants.
pub(crate) fn parse_sample_variants(
    path: &str,
) -> Result<HashMap<EnsemblId, Vec<Variant>>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("failed to read variants file '{}'", path))?;
    let mut variants: HashMap<EnsemblId, Vec<Variant>> = HashMap::new();
    for (line_no, line) in content.lines().enumerate() {
        // e.g. "  ENST00000367770.5    726G>726S  " -> "ENST00000367770.5    726G>726S"
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let identifier = format!("{}:{}", path, line_no + 1);

        let mut parts = line.split_whitespace();

        // e.g. "ENST00000367770.5    726G>726S" -> "ENST00000367770.5"
        let enst = match parts.next() {
            Some(e) => e,
            None => continue,
        };
        // e.g. "ENST00000367770.5" starts with "ENST" -> ok; "NM_001256789.1" -> not ok, logged and skipped
        if !enst.starts_with("ENST") {
            info!("{},'{}' does not look like an ENST id; skipping line", &identifier, enst);
            continue;
        }

        // e.g. "ENST00000367770.5    726G>726S" -> "726G>726S"
        let variant_string = match parts.next() {
            Some(v) => v,
            None => {
                // e.g. "ENST00000367770.5" alone on the line, no second field
                info!("{},missing variant field for {}", &identifier, enst);
                continue;
            }
        };

        // e.g. "726G>726S" -> Some(true); "726G>731S" -> Some(false), logged and skipped
        match positions_match(variant_string) {
            Some(true) => {}
            Some(false) => {
                info!("{},position mismatch in variant '{}' for {}", &identifier, variant_string, enst);
                continue;
            }
            None => {
            }
        }

        //let ensembl_id_versionless = strip_version(enst);

        // get entry for variant index
        let entry = variants
            .entry(EnsemblId(strip_version(enst).to_string()))
            .or_default();
        let idx = entry.len();

        // e.g. "726G>726S" -> Some(Variant { .. }); "garbage!!" -> None, logged and skipped
        // variant id is set to "{ensembl_id_versionless}:{idx}" so that each variant has a unique id
        let variant = match to_variant(&format!("{}_{}", enst, idx), variant_string) {
            Some(v) => v,
            None => {
                info!("{},could not parse variant '{}' for {}", &identifier, variant_string, enst);
                continue;
            }
        };

        entry.push(variant);
    }
    Ok(variants)
}

#[cfg(test)]
#[path = "parsing_tests.rs"]
mod tests;
