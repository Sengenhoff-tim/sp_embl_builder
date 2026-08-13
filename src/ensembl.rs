// ============================================================================
// Ensembl REST — batch sequence fetch via POST /sequence/id
// ============================================================================

use crate::types::EnsemblId;
use crate::util::{with_retries, RateLimiter, RetryConfig};
use anyhow::{anyhow, Context, Result};
use futures::stream::{self, StreamExt};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

const ENSEMBL_SEQUENCE_URL: &str = "https://rest.ensembl.org/sequence/id";
// Ensembl's POST endpoint caps batches at 50 ids per request.
const MAX_IDS_PER_BATCH: usize = 50;
const MAX_CONCURRENT_ENSEMBL_REQUESTS: usize = 5;

// Ensembl's default per-client rate limit; adjust to match whatever
// your RateLimiter expects (e.g. requests per second).
pub(crate) const ENSEMBL_MAX_REQUESTS_PER_SECOND: u32 = 15;

#[derive(Debug, Serialize)]
struct SequenceIdPostRequest<'a> {
    ids: &'a [String],
}

#[derive(Debug, Deserialize)]
struct SequenceIdPostResponseItem {
    // Present on success
    id: Option<String>,
    seq: Option<String>,
    // Present on per-id failure (Ensembl returns {"error": "..."} entries
    // inline in the array rather than failing the whole batch)
    error: Option<String>,
}

async fn fetch_ensembl_sequence_batch(
    ids: &[String],
    retry_config: &RetryConfig,
    client: &reqwest::Client,
) -> Result<Vec<SequenceIdPostResponseItem>> {
    with_retries(
        &format!("fetching Ensembl protein sequences for {} ids", ids.len()),
        retry_config,
        || async {
            let body = SequenceIdPostRequest { ids };
            let response = client
                .post(ENSEMBL_SEQUENCE_URL)
                .query(&[("type", "protein")])
                .header("Content-Type", "application/json")
                .header("Accept", "application/json")
                .json(&body)
                .timeout(retry_config.timeout)
                .send()
                .await
                .with_context(|| format!("Ensembl sequence POST request failed ({} ids)", ids.len()))?;

            if !response.status().is_success() {
                let status = response.status();
                let text = response.text().await.unwrap_or_default();
                return Err(anyhow!(
                    "Ensembl sequence POST request failed: HTTP {} body: {}",
                    status,
                    text
                ));
            }

            let items: Vec<SequenceIdPostResponseItem> = response
                .json()
                .await
                .context("failed to parse Ensembl sequence POST response as JSON")?;
            Ok(items)
        },
    )
    .await
}

pub(crate) async fn fetch_ensembl_sequences(
    enst_ids: &[EnsemblId],
    retry_config: &RetryConfig,
    rate_limiter: &RateLimiter,
) -> Result<HashMap<EnsemblId, String>> {
    fetch_ensembl_sequences_with_client(enst_ids, retry_config, rate_limiter, &reqwest::Client::new()).await
}

pub(crate) async fn fetch_ensembl_sequences_with_client(
    enst_ids: &[EnsemblId],
    retry_config: &RetryConfig,
    rate_limiter: &RateLimiter,
    client: &reqwest::Client,
) -> Result<HashMap<EnsemblId, String>> {
    let mut result: HashMap<EnsemblId, String> = HashMap::new();
    if enst_ids.is_empty() {
        return Ok(result);
    }

    // Preserve a lookup from the string form back to the original EnsemblId,
    // since Ensembl echoes ids back as plain strings.
    let id_lookup: HashMap<String, EnsemblId> = enst_ids
        .iter()
        .map(|id| (id.as_str().to_string(), id.clone()))
        .collect();

    let chunks: Vec<Vec<String>> = enst_ids
        .chunks(MAX_IDS_PER_BATCH)
        .map(|c| c.iter().map(|id| id.as_str().to_string()).collect())
        .collect();

    tracing::info!(
        total_ids = enst_ids.len(),
        num_batches = chunks.len(),
        "Starting batched Ensembl sequence fetch"
    );

    let mut results = stream::iter(chunks.into_iter())
        .map(|chunk| {
            let retry_config = *retry_config;
            let client = client.clone();
            let rate_limiter = rate_limiter.clone();
            async move {
                // Throttle per-batch (i.e. per HTTP request), not per id,
                // since each chunk is a single POST request.
                rate_limiter.throttle().await;
                fetch_ensembl_sequence_batch(&chunk, &retry_config, &client).await
            }
        })
        .buffer_unordered(MAX_CONCURRENT_ENSEMBL_REQUESTS);

    while let Some(batch_result) = results.next().await {
        let items = batch_result?;
        for item in items {
            match (item.id, item.seq, item.error) {
                (Some(id), Some(seq), _) => {
                    if let Some(enst_id) = id_lookup.get(&id) {
                        result.insert(enst_id.clone(), seq);
                    } else {
                        tracing::warn!(id = %id, "Ensembl returned sequence for unrequested id");
                    }
                }
                (id, _, Some(err)) => {
                    tracing::info!(
                        id = id.as_deref().unwrap_or("<unknown>"),
                        error = %err,
                        "Ensembl reported error for id, skipping"
                    );
                }
                _ => {
                    tracing::warn!("Ensembl returned malformed entry with no id/seq/error");
                }
            }
        }
    }

    tracing::info!(
        total_results = result.len(),
        requested = enst_ids.len(),
        "Batched Ensembl sequence fetch completed"
    );

    Ok(result)
}