// ============================================================================
// Ensembl REST — batch sequence fetch via POST /sequence/id
// ============================================================================

use crate::types::EnsemblId;
use crate::util::{with_retries, RateLimiter, RetryConfig};
use anyhow::{Context, Result};
use futures::stream::{self, StreamExt};
use serde::{Deserialize, Serialize};
use core::fmt;
use std::collections::{HashMap, HashSet};

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
    // The id as originally requested — this is what we match against,
    // NOT `id`, which can differ (e.g. Ensembl may return a resolved
    // protein id here even though we requested a transcript id).
    query: Option<String>,
    seq: Option<String>,
    error: Option<String>,
}

#[derive(Debug)]
struct EnsemblSequenceError {
    status: reqwest::StatusCode,
    body: String,
    ids: String,
}

impl fmt::Display for EnsemblSequenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Ensembl sequence POST request failed: HTTP {} for id batch [{}]",
            self.status, self.ids
        )
    }
}

impl std::error::Error for EnsemblSequenceError {}


impl EnsemblSequenceError {
    fn status(&self) -> reqwest::StatusCode {
        self.status
    }
}

fn is_server_error(err: &anyhow::Error) -> bool {
    err.downcast_ref::<EnsemblSequenceError>()
        .map(|e| e.status().is_server_error())
        .unwrap_or(false)
}

/// Sends a single POST request for exactly `ids` (no bisection, no
/// batch-size assumptions beyond Ensembl's hard cap). Retries transient
/// failures per `retry_config`.
async fn fetch_ensembl_sequence_batch_raw(
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
                .with_context(|| {
                    format!(
                        "Ensembl sequence POST request failed for id batch [{}]",
                        ids.join(",")
                    )
                })?;

            if !response.status().is_success() {
                let status = response.status();
                let text = response.text().await.unwrap_or_default();
                return Err(EnsemblSequenceError {
                    status,
                    body: text,
                    ids: ids.join(","),
                }
                .into());
            }

            let items: Vec<SequenceIdPostResponseItem> = response
                .json()
                .await
                .context("failed to parse Ensembl sequence POST response as JSON")?;

            // NOTE: Ensembl may return fewer entries than requested ids —
            // it can silently omit ids with no data instead of including
            // an {"error": ...} placeholder. That's expected, not a
            // transport failure, so we don't error/retry on count here.
            Ok(items)
        },
    )
    .await
}

/// Fetches a batch of ids, bisecting on server errors (5xx) to isolate
/// whichever id(s) are actually causing the failure, so one bad id can't
/// poison an entire 50-id batch forever. Non-server-error failures (e.g.
/// timeouts, connection errors) are NOT bisected — they're assumed
/// transient and are left to `with_retries` inside the raw call.
fn fetch_ensembl_sequence_batch<'a>(
    ids: &'a [String],
    retry_config: &'a RetryConfig,
    client: &'a reqwest::Client,
) -> std::pin::Pin<Box<dyn std::future::Future<Output = Result<Vec<SequenceIdPostResponseItem>>> + 'a>> {
    Box::pin(async move {
        // Base case: give up bisecting past a single id — either it
        // succeeds, or we log and drop it so one bad id can't fail the
        // whole run.
        if ids.len() == 1 {
            return Ok(
                fetch_ensembl_sequence_batch_raw(ids, retry_config, client)
                    .await
                    .unwrap_or_else(|err| {
                        tracing::info!(
                            id = %ids[0],
                            error = %err,
                            "Ensembl request failed for single id, dropped"
                        );
                        Vec::new()
                    }),
            );
        }

        match fetch_ensembl_sequence_batch_raw(ids, retry_config, client).await {
            Ok(items) => Ok(items),
            Err(err) if is_server_error(&err) => {
                tracing::warn!(
                    batch_size = ids.len(),
                    error = %err,
                    "Ensembl batch failed with server error, bisecting to isolate bad id(s)"
                );
                let mid = ids.len() / 2;
                let (left, right) = ids.split_at(mid);
                let mut left_items =
                    fetch_ensembl_sequence_batch(left, retry_config, client).await?;
                let right_items =
                    fetch_ensembl_sequence_batch(right, retry_config, client).await?;
                left_items.extend(right_items);
                Ok(left_items)
            }
            // Non-server-error failures (timeouts, connection errors, JSON
            // parse errors, etc.) are not bisected — they've already been
            // retried inside fetch_ensembl_sequence_batch_raw and are
            // propagated as-is.
            Err(err) => Err(err),
        }
    })
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
    // since Ensembl echoes ids back as plain strings (in the `query` field).
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
                // Throttle per-batch (i.e. per top-level HTTP request), not
                // per id. Note: if a batch gets bisected due to a server
                // error, the resulting sub-requests are NOT separately
                // throttled here — bisection is expected to be rare.
                rate_limiter.throttle().await;
                let items = fetch_ensembl_sequence_batch(&chunk, &retry_config, &client).await?;
                Ok::<_, anyhow::Error>((chunk, items))
            }
        })
        .buffer_unordered(MAX_CONCURRENT_ENSEMBL_REQUESTS);

    while let Some(batch_result) = results.next().await {
        let (chunk, items) = batch_result?;

        let mut seen: HashSet<&str> = HashSet::new();

        for item in &items {
            let Some(requested_id) = item.query.as_deref() else {
                tracing::warn!("Ensembl response entry missing 'query' field, skipping");
                continue;
            };
            seen.insert(requested_id);

            match (&item.seq, &item.error) {
                (Some(seq), _) => {
                    if let Some(enst_id) = id_lookup.get(requested_id) {
                        result.insert(enst_id.clone(), seq.clone());
                    } else {
                        tracing::warn!(id = %requested_id, "no lookup entry for requested id (bug)");
                    }
                }
                (_, Some(err)) => {
                    tracing::info!(
                        id = %requested_id,
                        error = %err,
                        "Ensembl reported error for id, skipping"
                    );
                }
                _ => {
                    tracing::warn!(
                        id = %requested_id,
                        "Ensembl returned entry with neither seq nor error"
                    );
                }
            }
        }

        // Ids Ensembl silently dropped from the response entirely (no
        // data available), or that were dropped by bisection after
        // repeated single-id failures.
        for requested_id in &chunk {
            if !seen.contains(requested_id.as_str()) {
                tracing::info!(
                    id = %requested_id,
                    "Ensembl returned no entry for id (no data or dropped after failure), skipping"
                );
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
