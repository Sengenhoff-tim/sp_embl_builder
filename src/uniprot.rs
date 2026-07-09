// ============================================================================
// UniProt REST JSON model and fetch
// ============================================================================

use crate::types::UniprotId;
use crate::util::{with_retries, RateLimiter};
use anyhow::{Context, Result};
use std::collections::HashMap;

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct UniProtEntry {
    pub(crate) primary_accession: String,
    pub(crate) protein_description: Option<ProteinDescription>,
    pub(crate) genes: Option<Vec<GeneEntry>>,
    pub(crate) sequence: UniProtSequence,
    #[serde(default)]
    pub(crate) features: Vec<UniProtFeature>,
    #[serde(default)]
    pub(crate) comments: Vec<UniProtComment>,
    #[serde(rename = "uniProtKBCrossReferences", default)]
    pub(crate) cross_references: Vec<UniProtCrossReference>,
}

impl UniProtEntry {
    pub(crate) fn protein_name(&self) -> &str {
        self.protein_description.as_ref()
            .and_then(|pd| pd.recommended_name.as_ref())
            .and_then(|rn| rn.full_name.as_ref())
            .map(|n| n.value.as_str())
            .unwrap_or("Unknown")
    }

    pub(crate) fn gene_name(&self) -> &str {
        self.genes.as_ref()
            .and_then(|g| g.first())
            .and_then(|g| g.gene_name.as_ref())
            .map(|n| n.value.as_str())
            .unwrap_or("Unknown")
    }
}

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct ProteinDescription {
    pub(crate) recommended_name: Option<RecommendedName>,
}

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct RecommendedName {
    pub(crate) full_name: Option<NameValue>,
}

#[derive(Debug, serde::Deserialize)]
pub(crate) struct NameValue {
    pub(crate) value: String,
}

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct GeneEntry {
    pub(crate) gene_name: Option<NameValue>,
}

#[derive(Debug, serde::Deserialize)]
pub(crate) struct UniProtSequence {
    pub(crate) value: String,
}

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct UniProtFeature {
    #[serde(rename = "type")]
    pub(crate) feature_type: String,
    pub(crate) location: FeatureLocation,
    pub(crate) description: Option<String>,
    #[serde(rename = "featureId")]
    pub(crate) feature_id: Option<String>,
    pub(crate) alternative_sequence: Option<AlternativeSequence>,
    #[serde(default)]
    pub(crate) evidences: Vec<UniProtEvidence>,
    #[serde(rename = "featureCrossReferences", default)]
    pub(crate) feature_cross_references: Vec<FeatureCrossRef>,
}

#[derive(Debug, serde::Deserialize)]
pub(crate) struct FeatureLocation {
    pub(crate) start: PositionValue,
    pub(crate) end: PositionValue,
}

#[derive(Debug, serde::Deserialize)]
pub(crate) struct PositionValue {
    pub(crate) value: usize,
}

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct AlternativeSequence {
    pub(crate) original_sequence: Option<String>,
    #[serde(default)]
    pub(crate) alternative_sequences: Vec<String>,
}

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct UniProtEvidence {
    pub(crate) evidence_code: String,
    pub(crate) source: Option<EvidenceSource>,
}

#[derive(Debug, serde::Deserialize)]
#[serde(untagged)]
pub(crate) enum EvidenceSource {
    Full { name: String, id: String },
    NameOnly(String),
}

#[derive(Debug, serde::Deserialize)]
pub(crate) struct FeatureCrossRef {
    pub(crate) database: String,
    pub(crate) id: String,
}

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct UniProtComment {
    pub(crate) comment_type: String,
    #[serde(default)]
    pub(crate) isoforms: Vec<IsoformInfo>,
    #[serde(default)]
    pub(crate) events: Vec<String>,
    pub(crate) note: Option<CommentNote>,
}

#[derive(Debug, serde::Deserialize)]
#[serde(untagged)]
pub(crate) enum CommentNote {
    Structured { texts: Vec<NoteText> },
    Plain(String),
}

impl CommentNote {
    pub(crate) fn texts(&self) -> Vec<&str> {
        match self {
            CommentNote::Structured { texts } => texts.iter().map(|t| t.value.as_str()).collect(),
            CommentNote::Plain(s) => vec![s.as_str()],
        }
    }
}

#[derive(Debug, serde::Deserialize)]
pub(crate) struct NoteText {
    pub(crate) value: String,
}

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct IsoformInfo {
    pub(crate) name: NameValue,
    #[serde(default)]
    pub(crate) synonyms: Vec<NameValue>,
    #[serde(default)]
    pub(crate) isoform_ids: Vec<String>,
    #[serde(default)]
    pub(crate) sequence_ids: Vec<String>,
    pub(crate) isoform_sequence_status: Option<String>,
}

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct UniProtCrossReference {
    pub(crate) database: String,
    pub(crate) id: String,
    pub(crate) isoform_id: Option<String>,
}

#[derive(Debug, serde::Deserialize)]
struct UniProtSearchResponse {
    results: Vec<UniProtEntry>,
}

/// Conservative rate limit — UniProt's published limit is not confirmed in this
/// environment; lower this if you observe 429 responses.
const UNIPROT_REQUESTS_PER_SECOND: u32 = 50;
/// Batch size for UniProt stream requests. The stream endpoint handles large OR
/// queries but very long URLs may be rejected; lower if you see 414 errors.
const MAX_UNIPROT_BATCH_SIZE: usize = 200;

fn build_uniprot_stream_url(accessions: &[String]) -> String {
    let query = accessions
        .iter()
        .map(|a| format!("accession:{}", a))
        .collect::<Vec<_>>()
        .join("+OR+");
    format!(
        "https://rest.uniprot.org/uniprotkb/stream?query={}&format=json",
        query
    )
}

fn fetch_uniprot_batch(accessions: &[String]) -> Result<Vec<UniProtEntry>> {
    with_retries(
        &format!("fetching {} UniProt entries", accessions.len()),
        || {
            let url = build_uniprot_stream_url(accessions);
            let body = ureq::get(&url)
                .call()
                .with_context(|| {
                    format!("UniProt request failed for {} accessions", accessions.len())
                })?
                .into_string()
                .context("failed to read UniProt response body")?;
            let response: UniProtSearchResponse = serde_json::from_str(&body)
                .context("failed to parse UniProt stream response")?;
            Ok(response.results)
        },
    )
}

pub(crate) fn fetch_uniprot_entries(accessions: &[UniprotId]) -> HashMap<UniprotId, UniProtEntry> {
    let id_strings: Vec<String> = accessions
        .iter()
        .map(|id| id.as_str().to_string())
        .collect();
    let mut by_accession: HashMap<String, UniProtEntry> = HashMap::new();
    let mut rate_limiter = RateLimiter::per_second(UNIPROT_REQUESTS_PER_SECOND);

    for chunk in id_strings.chunks(MAX_UNIPROT_BATCH_SIZE) {
        rate_limiter.throttle();
        match fetch_uniprot_batch(chunk) {
            Ok(entries) => {
                for entry in entries {
                    by_accession.insert(entry.primary_accession.clone(), entry);
                }
            }
            Err(e) => eprintln!(
                "Warning: failed to fetch UniProt batch of {} accessions: {:?}",
                chunk.len(),
                e
            ),
        }
    }

    accessions
        .iter()
        .filter_map(|id| by_accession.remove(id.as_str()).map(|e| (id.clone(), e)))
        .collect()
}
