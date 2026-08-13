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


#[derive(Debug, Serialize)]
struct SequenceIdPostRequest<'a> {
    ids: &'a [String],
}

#[derive(Debug, Deserialize)]
struct SequenceIdPostResponseItem {
    query: Option<String>,
    seq: Option<String>,
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
                    rate_limiter.throttle().await;
                    let items = fetch_ensembl_sequence_batch(&chunk, &retry_config, &client).await?;
                    Ok::<_, anyhow::Error>((chunk, items))
                }
            })
            .buffer_unordered(MAX_CONCURRENT_ENSEMBL_REQUESTS);

    while let Some(batch_result) = results.next().await {
        let (chunk, items) = batch_result?;

        let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();

        for item in &items {
            let Some(requested_id) = item.query.as_deref() else {
                tracing::warn!("Ensembl response entry missing 'query' field, skipping");
                continue;
            };
            seen.insert(requested_id);

            match (&item.seq, &item.error) {
                (Some(seq), _) => {
                    if let Some(enst_id) = id_lookup.get(requested_id) {
                        let validated_seq = strip_stop_codon(enst_id.as_str(), seq.clone());
                        result.insert(enst_id.clone(), validated_seq);
                    } else {
                        tracing::warn!(id = %requested_id, "no lookup entry for requested id (bug)");
                    }
                }
                (_, Some(err)) => {
                    tracing::info!(id = %requested_id, error = %err, "Ensembl reported error for id, skipping");
                }
                _ => {
                    tracing::warn!(id = %requested_id, "Ensembl returned entry with neither seq nor error");
                }
            }
        }

        // Ids Ensembl silently dropped from the response entirely.
        for requested_id in &chunk {
            if !seen.contains(requested_id.as_str()) {
                tracing::info!(id = %requested_id, "Ensembl returned no entry for id (no data), skipping");
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

    /// Ensembl occasionally returns a translated protein sequence with containing '*' stop-codon marker (e.g. for select readthrough transcripts). 
    /// Sequence is assumed to stop at the first stop codon
    fn strip_stop_codon(enst: &str, seq: String) -> String {
        match seq.find('*') {
            Some(pos) => {
                tracing::info!(
                    enst = enst,
                    "{}, sequence contained '*' stop codon at position {}; amino acids after are discarded",
                    enst,
                    pos
                );
                seq[..pos].to_string()
            }
            None => seq,
        }
    }