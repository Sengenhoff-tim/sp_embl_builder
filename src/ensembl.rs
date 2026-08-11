// ============================================================================
// Ensembl sequence API model and fetch
// ============================================================================

use crate::types::{EnsemblId, Sequence};
use crate::util::{with_retries, RetryConfig};
use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;

const MAX_ENSEMBL_BATCH_SIZE: usize = 50;
const ENSEMBL_SEQUENCE_URL: &str = "https://rest.ensembl.org/sequence/id?type=protein";

#[derive(serde::Serialize)]
struct EnsemblSequenceRequestBody<'a> {
    ids: &'a [String],
}

#[derive(Debug, serde::Deserialize)]
struct EnsemblSequenceResponse {
    molecule: String,
    seq: String,
    query: String,
    #[serde(rename = "id")]
    #[allow(dead_code)]
    response_id: String,
}

async fn fetch_ensembl_sequence_batch(
    ids: &[String],
    retry_config: &RetryConfig,
    client: &reqwest::Client,
) -> Result<Vec<EnsemblSequenceResponse>> {
    with_retries(
        &format!("fetching Ensembl sequences for {} ids", ids.len()),
        retry_config,
        || async {
            let response = client
                .post(ENSEMBL_SEQUENCE_URL)
                .timeout(retry_config.timeout)
                .json(&EnsemblSequenceRequestBody { ids })
                .send()
                .await
                .context("Ensembl sequence request failed")?;

            if !response.status().is_success() {
                let status = response.status();
                return Err(anyhow!("{},Ensembl sequence request failed: HTTP {}", ids.join("|"), status));
            }

            let responses: Vec<EnsemblSequenceResponse> = response
                .json()
                .await
                .context("failed to parse Ensembl sequence response")?;
            Ok(responses)
        },
    )
    .await
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

pub(crate) async fn fetch_ensembl_sequences(
    enst_ids: &[EnsemblId],
    retry_config: &RetryConfig,
) -> Result<Vec<(EnsemblId, Sequence)>> {
    let id_strings: Vec<String> = enst_ids
        .iter()
        .map(|id| id.as_str().to_string())
        .collect();
    let mut sequences: HashMap<String, Sequence> = HashMap::new();
    let client = reqwest::Client::new();

    for chunk in id_strings.chunks(MAX_ENSEMBL_BATCH_SIZE) {
        let responses = fetch_ensembl_sequence_batch(chunk, retry_config, &client).await?;
        for record in responses {
            if record.molecule != "protein" {
                tracing::info!(
                    enst = record.query.as_str(),
                    "expected molecule type 'protein', got '{}'; skipping",
                    record.molecule
                );
                continue;
            }
            let seq = strip_stop_codon(&record.query, record.seq);
            sequences.insert(record.query, Sequence(seq));
        }
    }

    Ok(enst_ids
        .iter()
        .filter_map(|id| {
            sequences
                .remove(id.as_str())
                .map(|seq| (id.clone(), seq))
        })
        .collect())
}
