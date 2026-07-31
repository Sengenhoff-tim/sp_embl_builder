// ============================================================================
// EBI Proteins API — variation model, fetch, and conversion to Variant
// ============================================================================

use crate::types::{UniProtCanonId, UniProtId, UniProtIsoId};
use crate::variant::Variant;
use crate::util::{with_retries, RateLimiter, RetryConfig};
use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;

// Only fields read elsewhere in this module are kept
#[derive(Debug, Clone, serde::Deserialize)]
pub(crate) struct ProteinFeatureInfo {
    accession: String,
    pub(crate) features: Vec<EbiFeature>,
}

#[derive(Debug, Clone, serde::Deserialize)]
pub(crate) struct EbiFeature {
    #[serde(rename = "ftId")]
    ft_id: Option<String>,
    begin: String,
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

pub(crate) fn fetch_iso_variations(
    uniprot_ids: &[UniProtIsoId],
    source_types: &[&str],
    retry_config: &RetryConfig,
) -> Result<HashMap<UniProtIsoId, ProteinFeatureInfo>> {
    let mut result: HashMap<UniProtIsoId, ProteinFeatureInfo> = HashMap::new();
    let mut rate_limiter = RateLimiter::per_second(EBI_MAX_REQUESTS_PER_SECOND);

    for chunk in uniprot_ids.chunks(MAX_EBI_BATCH_SIZE) {
        rate_limiter.throttle();
        let id_strings: Vec<String> = chunk.iter().map(|id| id.as_str().to_string()).collect();
        let infos = fetch_variation_batch(&id_strings, source_types, retry_config)?;
        for (i, info) in infos.into_iter().enumerate() {
            result.insert(chunk[i].clone(), info);
        }
    }

    Ok(result)
}

pub(crate) fn fetch_canon_variations(
    uniprot_ids: &[UniProtCanonId],
    source_types: &[&str],
    retry_config: &RetryConfig,
) -> Result<HashMap<UniProtCanonId, ProteinFeatureInfo>> {
    let mut result: HashMap<UniProtCanonId, ProteinFeatureInfo> = HashMap::new();
    let mut rate_limiter = RateLimiter::per_second(EBI_MAX_REQUESTS_PER_SECOND);

    for chunk in uniprot_ids.chunks(MAX_EBI_BATCH_SIZE) {
        rate_limiter.throttle();
        let id_strings: Vec<String> = chunk.iter().map(|id| id.as_str().to_string()).collect();
        let infos = fetch_variation_batch(&id_strings, source_types, retry_config)?;
        for (i, info) in infos.into_iter().enumerate() {
            result.insert(chunk[i].clone(), info);
        }
    }

    Ok(result)
}


// ============================================================================
// Variant building from EBI features
// ============================================================================


/// Build the `id` field for a Variant from an EBI feature.
/// Format: `EBI:{UniProtId}:{position}:{mutatedType}`
/// Falls back to the evidence string if neither ftId nor xrefs are present.
fn build_variant_id(feature: &EbiFeature) -> String {
    format!(
        "EBI:{}:{}:{}",
        feature
            .ft_id
            .as_ref()
            .map(|s| s.as_str())
            .unwrap_or("unknown"),
        feature.begin,
        feature
            .mutated_type
            .as_ref()
            .map(|s| s.as_str())
            .unwrap_or("unknown")
    )
}

pub(crate) fn feature_to_iso_variant(accession: &UniProtIsoId, feature: &EbiFeature) -> Option<Variant> {
    let replaced = feature.wild_type.clone()?;
    let replacement = feature.mutated_type.clone()?;
    let begin: usize = feature.begin.parse().ok()?;
    let end = begin + replaced.len().checked_sub(1)?;
    let id = build_variant_id(feature);
    let isoform = Some(accession.clone());
    Some(Variant {isoform, id, begin, end, aa_ref: replaced, aa_new: replacement })
}

pub(crate) fn feature_to_canon_variant(feature: &EbiFeature) -> Option<Variant> {
    let replaced = feature.wild_type.clone()?;
    let replacement = feature.mutated_type.clone()?;
    let begin: usize = feature.begin.parse().ok()?;
    let end = begin + replaced.len().checked_sub(1)?;
    let id = build_variant_id(feature);
    let isoform = None;
    Some(Variant {isoform, id, begin, end, aa_ref: replaced, aa_new: replacement })
}