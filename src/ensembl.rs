// ============================================================================
// Ensembl sequence API model and fetch
// ============================================================================

use crate::exceptions::ExceptionLog;
use crate::types::{EnsemblId, Sequence};
use crate::util::with_retries;
use anyhow::{Context, Result};
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

fn fetch_ensembl_sequence_batch(ids: &[String]) -> Result<Vec<EnsemblSequenceResponse>> {
    with_retries(
        &format!("fetching Ensembl sequences for {} ids", ids.len()),
        || {
            let responses: Vec<EnsemblSequenceResponse> =
                ureq::post(ENSEMBL_SEQUENCE_URL)
                    .send_json(EnsemblSequenceRequestBody { ids })
                    .context("Ensembl sequence request failed")?
                    .into_json()
                    .context("failed to parse Ensembl sequence response")?;
            Ok(responses)
        },
    )
}

/// Ensembl occasionally returns a translated protein sequence with a
/// trailing '*' stop-codon marker (e.g. for select readthrough transcripts).
/// Strip it so downstream position math, which assumes a pure amino-acid
/// sequence, stays correct.
fn strip_stop_codon(enst: &str, seq: String, exceptions: &mut ExceptionLog) -> String {
    match seq.strip_suffix('*') {
        Some(stripped) => {
            exceptions.log(enst, "retrieved sequence had a trailing '*' stop codon; stripped it");
            stripped.to_string()
        }
        None => seq,
    }
}

pub(crate) fn fetch_ensembl_sequences(
    enst_ids: &[EnsemblId],
    exceptions: &mut ExceptionLog,
) -> Vec<(EnsemblId, Sequence)> {
    let id_strings: Vec<String> = enst_ids
        .iter()
        .map(|id| id.as_str().to_string())
        .collect();
    let mut sequences: HashMap<String, Sequence> = HashMap::new();

    for chunk in id_strings.chunks(MAX_ENSEMBL_BATCH_SIZE) {
        match fetch_ensembl_sequence_batch(chunk) {
            Ok(responses) => {
                for record in responses {
                    if record.molecule != "protein" {
                        exceptions.log(
                            &record.query,
                            &format!("expected molecule type 'protein', got '{}'; skipping", record.molecule),
                        );
                        continue;
                    }
                    let seq = strip_stop_codon(&record.query, record.seq, exceptions);
                    sequences.insert(record.query, Sequence(seq));
                }
            }
            Err(e) => {
                for id in chunk {
                    exceptions.log(id, &format!("Ensembl sequence batch fetch failed: {:?}", e));
                }
            }
        }
    }

    enst_ids
        .iter()
        .filter_map(|id| {
            sequences
                .remove(id.as_str())
                .map(|seq| (id.clone(), seq))
        })
        .collect()
}
