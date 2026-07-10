// ============================================================================
// EBI Proteins API — variation model, fetch, and conversion to Variant
// ============================================================================

use crate::types::{UniprotId, Variant};
use crate::util::{with_retries, RateLimiter, RetryConfig};
use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;

// Only fields actually read elsewhere in this module are kept; the EBI response
// carries plenty more, but serde ignores unmapped JSON keys by default (no
// `deny_unknown_fields`), so trimming these doesn't affect parsing.
#[derive(Debug, Clone, serde::Deserialize)]
pub(crate) struct ProteinFeatureInfo {
    accession: String,
    features: Vec<EbiFeature>,
}

#[derive(Debug, Clone, serde::Deserialize)]
struct EbiFeature {
    #[serde(rename = "ftId")]
    ft_id: Option<String>,
    begin: String,
    xrefs: Option<Vec<EbiDbRef>>,
    #[serde(default)]
    evidences: Vec<EbiEvidence>,
    #[serde(rename = "wildType")]
    wild_type: Option<String>,
    #[serde(rename = "mutatedType")]
    mutated_type: Option<String>,
}

#[derive(Debug, Clone, serde::Deserialize)]
struct EbiDbRef {
    name: String,
    id: String,
}

#[derive(Debug, Clone, serde::Deserialize)]
struct EbiEvidence {
    code: String,
    source: Option<EbiDbRef>,
}

const MAX_EBI_BATCH_SIZE: usize = 100;
const EBI_MAX_REQUESTS_PER_SECOND: u32 = 200;

const ALLOWED_SOURCE_TYPES: [&str; 10] = [
    "uniprot",
    "large scale study",
    "mixed",
    "clinvar",
    "nci-tcga",
    "cosmic curated",
    "ensembl",
    "gnomad",
    "topmed",
    "exac",
];

fn validate_source_types(source_types: &[&str]) -> Result<()> {
    if source_types.len() > 2 {
        return Err(anyhow!(
            "at most 2 sourcetype values accepted, got {}",
            source_types.len()
        ));
    }
    for st in source_types {
        if !ALLOWED_SOURCE_TYPES
            .iter()
            .any(|a| a.eq_ignore_ascii_case(st))
        {
            return Err(anyhow!(
                "invalid sourcetype '{}', must be one of {:?}",
                st,
                ALLOWED_SOURCE_TYPES
            ));
        }
    }
    Ok(())
}

fn build_variation_batch_url(uniprot_ids: &[String], source_types: &[&str]) -> Result<String> {
    validate_source_types(source_types)?;
    if uniprot_ids.is_empty() {
        return Err(anyhow!("build_variation_batch_url called with no accessions"));
    }
    if uniprot_ids.len() > MAX_EBI_BATCH_SIZE {
        return Err(anyhow!(
            "at most {} accessions per batch, got {}",
            MAX_EBI_BATCH_SIZE,
            uniprot_ids.len()
        ));
    }
    let mut url = format!(
        "https://www.ebi.ac.uk/proteins/api/variation?accession={}",
        uniprot_ids.join(",")
    );
    if !source_types.is_empty() {
        url.push_str("&sourcetype=");
        url.push_str(&source_types.join(","));
    }
    Ok(url)
}

fn fetch_variation_batch(
    uniprot_ids: &[String],
    source_types: &[&str],
    retry_config: &RetryConfig,
) -> Result<Vec<ProteinFeatureInfo>> {
    with_retries(
        &format!("fetching EBI variation for {} accessions", uniprot_ids.len()),
        retry_config,
        || {
            let url = build_variation_batch_url(uniprot_ids, source_types)?;
            let body = ureq::get(&url)
                .timeout(retry_config.timeout)
                .call()
                .with_context(|| format!("EBI variation request to {} failed", url))?
                .into_string()
                .context("failed to read EBI variation response body")?;
            let info: Vec<ProteinFeatureInfo> = serde_json::from_str(&body)
                .context("failed to parse EBI variation response as JSON")?;
            Ok(info)
        },
    )
}

/// Fetches EBI variation data for all `uniprot_ids`. A batch that still
/// fails after retries is treated as fatal rather than being dropped: a
/// silently missing batch of variants would otherwise pass through as
/// entries with no variation data, indistinguishable from accessions that
/// genuinely have none.
pub(crate) fn fetch_variations(
    uniprot_ids: &[UniprotId],
    source_types: &[&str],
    retry_config: &RetryConfig,
) -> Result<HashMap<UniprotId, ProteinFeatureInfo>> {
    let id_strings: Vec<String> = uniprot_ids
        .iter()
        .map(|id| id.as_str().to_string())
        .collect();
    let mut by_accession: HashMap<String, ProteinFeatureInfo> = HashMap::new();
    let mut rate_limiter = RateLimiter::per_second(EBI_MAX_REQUESTS_PER_SECOND);

    for chunk in id_strings.chunks(MAX_EBI_BATCH_SIZE) {
        rate_limiter.throttle();
        let infos = fetch_variation_batch(chunk, source_types, retry_config)?;
        for info in infos {
            by_accession.insert(info.accession.clone(), info);
        }
    }

    Ok(uniprot_ids
        .iter()
        .filter_map(|id| by_accession.remove(id.as_str()).map(|info| (id.clone(), info)))
        .collect())
}

// ============================================================================
// Variant building from EBI features
// ============================================================================

fn format_ebi_evidence_id(evidences: &[EbiEvidence]) -> String {
    evidences
        .iter()
        .map(|e| match &e.source {
            Some(s) => format!("{}|{}:{}", e.code, s.name, s.id),
            None => e.code.clone(),
        })
        .collect::<Vec<_>>()
        .join("|")
}

/// Build the `id` field for a Variant from an EBI feature.
/// Format: `Uniprot:{ftId}|{xref.name}:{xref.id}|...`
/// Falls back to the evidence string if neither ftId nor xrefs are present.
fn build_variant_id(feature: &EbiFeature) -> String {
    let mut parts: Vec<String> = Vec::new();

    if let Some(ft_id) = &feature.ft_id {
        if !ft_id.is_empty() {
            parts.push(format!("Uniprot:{}", ft_id));
        }
    }

    if let Some(xrefs) = &feature.xrefs {
        for xref in xrefs {
            parts.push(format!("{}:{}", xref.name, xref.id));
        }
    }

    if parts.is_empty() && !feature.evidences.is_empty() {
        parts.push(format_ebi_evidence_id(&feature.evidences));
    }

    parts.join("|")
}

fn feature_to_variant(feature: &EbiFeature) -> Option<Variant> {
    let replaced = feature.wild_type.clone()?;
    let replacement = feature.mutated_type.clone()?;
    let begin: usize = feature.begin.parse().ok()?;
    let end = begin + replaced.len().checked_sub(1)?;
    let id = build_variant_id(feature);
    Some(Variant { id, begin, end, replaced, replacement, isoform_ref: None })
}

pub(crate) fn variants_from_protein_feature_info(info: &ProteinFeatureInfo) -> Vec<Variant> {
    info.features.iter().filter_map(feature_to_variant).collect()
}
