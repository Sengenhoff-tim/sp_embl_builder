// ============================================================================
// EBI Proteins API — variation model, fetch, and conversion to Variant
// ============================================================================

use crate::types::{UniProtCanonId, UniProtIsoId};
use crate::variant::Variant;
use crate::util::{with_retries, RateLimiter, RetryConfig};
use anyhow::{anyhow, Context, Result};
use futures::stream::{self, StreamExt};
use tracing::info;

const MAX_CONCURRENT_EBI_REQUESTS: usize = 10;
use std::collections::HashMap;

// Only fields read elsewhere in this module are kept
#[derive(Debug, Clone, serde::Deserialize)]
pub(crate) struct ProteinFeatureInfo {
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

pub(crate) const EBI_MAX_REQUESTS_PER_SECOND: u32 = 200;

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

/// Builds the URL for a single accession using EBI's singular path-style
/// endpoint (`/variation/{accession}`), rather than the query-param batch
/// endpoint (`/variation?accession=A,B,C`). We always fetch one accession
/// at a time now, and the batch endpoint is known to reject/400 on certain
/// single-accession requests, so the singular endpoint is used uniformly.
fn build_variation_url(accession: &str, source_types: &[&str]) -> Result<String> {
    validate_source_types(source_types)?;
    let mut url = format!("https://www.ebi.ac.uk/proteins/api/variation/{}", accession);
    if !source_types.is_empty() {
        url.push_str("?sourcetype=");
        url.push_str(&source_types.join(","));
    }
    Ok(url)
}

async fn fetch_variation_single(
    accession: &str,
    source_types: &[&str],
    retry_config: &RetryConfig,
    client: &reqwest::Client,
) -> Result<Option<ProteinFeatureInfo>> {
    with_retries(
        &format!("fetching EBI variation for {}", accession),
        retry_config,
        || async {
            let url = build_variation_url(accession, source_types)?;
            let response = client
                .get(&url)
                .timeout(retry_config.timeout)
                .send()
                .await
                .with_context(|| format!("EBI variation request to {} failed", url))?;

            /* 
            if response.status() == reqwest::StatusCode::NOT_FOUND {
                // No variation data for this accession - not an error.
                return Ok(None);
            }
            */

            if !response.status().is_success() {
                let status = response.status();
                let body = response.text().await.unwrap_or_default();
                return Err(anyhow!("EBI variation request failed: URL: [{}]; HTTP {}",url, status));
            }

            let body = response
                .text()
                .await
                .context("failed to read EBI variation response body")?;
            let info: ProteinFeatureInfo = serde_json::from_str(&body)
                .context("failed to parse EBI variation response as JSON")?;
            Ok(Some(info))
        },
    )
    .await
}

async fn fetch_variations_generic<K>(
    keys: &[K],
    as_str: impl Fn(&K) -> String,
    source_types: &[&str],
    retry_config: &RetryConfig,
    rate_limiter: &RateLimiter,
    client: &reqwest::Client,
) -> Result<HashMap<K, ProteinFeatureInfo>>
where
    K: Clone + Eq + std::hash::Hash,
{
    let mut result: HashMap<K, ProteinFeatureInfo> = HashMap::new();

    let source_types_owned: Vec<String> = source_types.iter().map(|s| s.to_string()).collect();

    tracing::info!(
        accession_count = keys.len(),
        "Starting per-accession EBI variation fetch"
    );

    let mut results = stream::iter(keys.iter().cloned())
        .map(|key| {
            let accession = as_str(&key);
            let rate_limiter = rate_limiter.clone();
            let retry_config = *retry_config;
            let client = client.clone();
            let source_types_owned = source_types_owned.clone();

            async move {
                rate_limiter.throttle().await;
                let source_type_refs: Vec<&str> = source_types_owned.iter().map(String::as_str).collect();

                tracing::info!(accession = %accession, source_types = ?source_type_refs, "Fetching variation");

                let info = fetch_variation_single(&accession, &source_type_refs, &retry_config, &client).await?;

                Ok::<_, anyhow::Error>((key, info))
            }
        })
        .buffer_unordered(MAX_CONCURRENT_EBI_REQUESTS);

    while let Some(res) = results.next().await {
        let (key, info) = res?;
        if let Some(info) = info {
            result.insert(key, info);
        }
    }

    tracing::info!(
        total_results = result.len(),
        "Per-accession EBI variation fetch completed"
    );

    Ok(result)
}


pub(crate) async fn fetch_iso_variations(
    uniprot_ids: &[UniProtIsoId],
    source_types: &[&str],
    retry_config: &RetryConfig,
    rate_limiter: &RateLimiter,
    client: &reqwest::Client,
) -> Result<HashMap<UniProtIsoId, ProteinFeatureInfo>> {
    fetch_variations_generic(uniprot_ids, |id| id.as_str().to_string(), source_types, retry_config, rate_limiter, client).await
}

pub(crate) async fn fetch_canon_variations(
    uniprot_ids: &[UniProtCanonId],
    source_types: &[&str],
    retry_config: &RetryConfig,
    rate_limiter: &RateLimiter,
    client: &reqwest::Client,
) -> Result<HashMap<UniProtCanonId, ProteinFeatureInfo>> {
    fetch_variations_generic(uniprot_ids, |id| id.as_str().to_string(), source_types, retry_config, rate_limiter, client).await
}


// ============================================================================
// Variant building from EBI features
// ============================================================================

pub(crate) fn feature_to_iso_variant(accession: &UniProtIsoId, feature: &EbiFeature) -> Option<Variant> {
    if let Some(replaced) = &feature.wild_type {
        let begin: usize = feature.begin.parse().ok()?;
        let end = begin + replaced.len().checked_sub(1)?;
        let id = "EBI".to_string();
        let isoform = Some(accession.clone());
        return Some(Variant {isoform, id, begin, end, aa_ref: Some(replaced.to_string()), aa_new: feature.mutated_type.clone() })
    }
    info!("{},Variant at [{}] replaces ambigues sequence", accession.as_str(), feature.begin);
    None
}

pub(crate) fn feature_to_canon_variant(accession: &str, feature: &EbiFeature) -> Option<Variant> {
    if let Some(replaced) = &feature.wild_type {
        let begin: usize = feature.begin.parse().ok()?;
        let end = begin + replaced.len().checked_sub(1)?;
        let id = format!("EBI:{}:{}",begin, replaced );
        let isoform = None;
        return Some(Variant {isoform, id, begin, end, aa_ref: Some(replaced.to_string()), aa_new: feature.mutated_type.clone() })
    }
    info!("{},Variant at [{}] replaces ambigues sequence", accession, feature.begin);
    None
}