// ============================================================================
// UniProt REST JSON model and fetch
// ============================================================================

use crate::types::{UniProtCanonId, UniProtIsoId, UniProtId, EnsemblId};
use crate::util::{with_retries, RateLimiter, RetryConfig};
use anyhow::{anyhow, Context, Result};
use futures::stream::{self, StreamExt};

/// Caps the number of UniProt batch requests in flight at once. Admission
/// is already paced by `RateLimiter`, but that only limits how fast new
/// requests are *started* — without this cap, every admitted request stays
/// open concurrently until it completes, which under load can pile up
/// enough simultaneous connections to make individual requests time out.
const MAX_CONCURRENT_UNIPROT_REQUESTS: usize = 8;
use std::collections::HashMap;

#[derive(Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct UniProtEntry {
    // Maps JSON "primaryAccession". This is the stable accession UniProt
    // returns for the entry, used below as the HashMap key when we
    // reassemble batch results back into per-accession lookups.
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

    /// Returns the isoforms declared across this entry's `ALTERNATIVE
    /// PRODUCTS` comments
    pub(crate) fn alternative_product_isoforms(&self) -> impl Iterator<Item = &IsoformInfo> {
        self.comments
            .iter()
            .filter(|c| c.comment_type == "ALTERNATIVE PRODUCTS")
            .flat_map(|c| c.isoforms.iter())
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

#[derive(Debug, Clone, serde::Deserialize)]
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

#[derive(Debug, Clone, serde::Deserialize)]
pub(crate) struct FeatureLocation {
    pub(crate) start: PositionValue,
    pub(crate) end: PositionValue,
    pub(crate) sequence: Option<String>,
}

#[derive(Debug, Clone, serde::Deserialize)]
pub(crate) struct PositionValue {
    // UniProt sometimes reports an unresolved/unknown position as
    // `"value": null` (typically alongside a `"modifier": "UNKNOWN"`
    // field), so this must be optional rather than a bare `usize`.
    pub(crate) value: Option<usize>,
}

#[derive(Debug, Clone, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct AlternativeSequence {
    pub(crate) original_sequence: Option<String>,
    #[serde(default)]
    pub(crate) alternative_sequences: Vec<String>,
}

#[derive(Debug, Clone, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct UniProtEvidence {
    pub(crate) evidence_code: String,
    pub(crate) source: Option<EvidenceSource>,
}

#[derive(Debug, Clone, serde::Deserialize)]
#[serde(untagged)]
pub(crate) enum EvidenceSource {
    Full { name: String, id: String },
    NameOnly(String),
}

#[derive(Debug, Clone, serde::Deserialize)]
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

    // Nested-object-or-plain-string form of the text, used by
    // DiseaseComment / CofactorComment / AlternativeProductsComment /
    // RnaEditingComment (object) and SequenceCautionComment /
    // WebResourceComment (plain string).
    pub(crate) note: Option<CommentNote>,

    #[serde(default)]
    pub(crate) texts: Vec<NoteText>,
}
/* 
impl UniProtComment {
    /// Unified accessor: returns the comment's text
    pub(crate) fn all_texts(&self) -> Vec<&str> {
        let mut out: Vec<&str> = self.texts.iter().map(|t| t.value.as_str()).collect();
        if let Some(note) = &self.note {
            out.extend(note.texts());
        }
        out
    }
}
    */

#[derive(Debug, serde::Deserialize)]
#[serde(untagged)]
pub(crate) enum CommentNote {
    // DiseaseComment / CofactorComment / AlternativeProductsComment /
    // RnaEditingComment: note is `{ valid: bool, texts: [...] }`.
    Structured { texts: Vec<NoteText> },
    // SequenceCautionComment / WebResourceComment: note is a bare string.
    // This variant must come second: a JSON string would fail to match
    // `Structured` (which requires an object), and correctly fall through
    // to here.
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
    pub(crate) name: Option<NameValue>,

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

const UNIPROT_REQUESTS_PER_SECOND: u32 = 50;

// Backstop on the number of accessions per batch. UniProt's query parser
// enforces a hard maximum of 100 `OR` conditions per query (confirmed via
// their error message: "Too many OR conditions in query. Maximum allowed
// is 100."), so this must never exceed 100.
const MAX_UNIPROT_BATCH_SIZE: usize = 100;

// Conservative cap on the *length* of the generated stream URL, kept as
// a secondary safeguard against pathologically long accession lists
// (many proxies/load balancers reject URLs above ~8KB).
const MAX_UNIPROT_URL_LEN: usize = 3000;

const UNIPROT_STREAM_BASE_URL: &str = "https://rest.uniprot.org/uniprotkb/stream?query=";
const UNIPROT_STREAM_SUFFIX: &str = "&format=json";

fn build_uniprot_stream_url(accessions: &[String]) -> String {
    let query = accessions
        .iter()
        .map(|a| format!("accession:{}", a))
        .collect::<Vec<_>>()
        .join("+OR+");
    format!("{}{}{}", UNIPROT_STREAM_BASE_URL, query, UNIPROT_STREAM_SUFFIX)
}

/// Splits `accessions` into batches that respect both `MAX_UNIPROT_BATCH_SIZE`
/// and `MAX_UNIPROT_URL_LEN`, so batches stay safely within UniProt's query
/// complexity limits and produce reasonably sized URLs.
fn batch_accessions(accessions: &[String]) -> Vec<Vec<String>> {
    let fixed_len = UNIPROT_STREAM_BASE_URL.len() + UNIPROT_STREAM_SUFFIX.len();
    let mut batches: Vec<Vec<String>> = Vec::new();
    let mut current: Vec<String> = Vec::new();
    let mut current_len = fixed_len;

    for accession in accessions {
        let term_len = "accession:".len() + accession.len();
        let added_len = term_len + if current.is_empty() { 0 } else { "+OR+".len() };

        let would_exceed_len = current_len + added_len > MAX_UNIPROT_URL_LEN;
        let would_exceed_count = current.len() >= MAX_UNIPROT_BATCH_SIZE;

        if !current.is_empty() && (would_exceed_len || would_exceed_count) {
            batches.push(std::mem::take(&mut current));
            current_len = fixed_len;
        }

        current_len += term_len + if current.is_empty() { 0 } else { "+OR+".len() };
        current.push(accession.clone());
    }

    if !current.is_empty() {
        batches.push(current);
    }

    batches
}

async fn fetch_uniprot_batch(
    accessions: &[String],
    retry_config: &RetryConfig,
    client: &reqwest::Client,
) -> Result<Vec<UniProtEntry>> {
    with_retries(
        &format!("fetching {} UniProt entries", accessions.len()),
        retry_config,
        || async {
            let url = build_uniprot_stream_url(accessions);
            let response = client
                .get(&url)
                .timeout(retry_config.timeout)
                .send()
                .await
                .with_context(|| {
                    format!("UniProt request failed for {} accessions", accessions.len())
                })?;

            if !response.status().is_success() {
                let status = response.status();
                let body = response.text().await.unwrap_or_default();
                return Err(anyhow!("UniProt request failed: HTTP {}: {}", status, body));
            }

            let body = response
                .text()
                .await
                .context("failed to read UniProt response body")?;
            let parsed: UniProtSearchResponse = serde_json::from_str(&body)
                .context("failed to parse UniProt stream response")?;
            Ok(parsed.results)
        },
    )
    .await
}

/// Fetches all UniProt entries for `accessions`, concurrently, while
/// enforcing a single global `UNIPROT_REQUESTS_PER_SECOND` cap shared by
/// every in-flight batch. A batch that still fails after retries is
/// treated as fatal to avoid silently missing entries.
pub(crate) async fn fetch_uniprot_entries(
    accessions: &[UniProtCanonId],
    retry_config: &RetryConfig,
) -> Result<HashMap<UniProtCanonId, UniProtEntry>> {
    let id_strings: Vec<String> = accessions
        .iter()
        .map(|id| id.as_str().to_string())
        .collect();

    let batches = batch_accessions(&id_strings);
    let total = batches.len();

    // One RateLimiter, cloned into every task below, so the per-second cap
    // is enforced globally across all concurrent requests to this API.
    let rate_limiter = RateLimiter::per_second(UNIPROT_REQUESTS_PER_SECOND);
    let client = reqwest::Client::new();

    let mut by_accession: HashMap<String, UniProtEntry> = HashMap::new();
    let mut done = 0usize;

    let mut results = stream::iter(batches)
        .map(|chunk| {
            let rate_limiter = rate_limiter.clone();
            let retry_config = *retry_config;
            let client = client.clone();
            async move {
                rate_limiter.throttle().await;
                fetch_uniprot_batch(&chunk, &retry_config, &client).await
            }
        })
        .buffer_unordered(MAX_CONCURRENT_UNIPROT_REQUESTS);

    while let Some(result) = results.next().await {
        let entries = result?;
        for entry in entries {
            by_accession.insert(entry.primary_accession.clone(), entry);
        }
        done += 1;
        eprint!("\rFetched UniProt batch {}/{}", done, total);
        use std::io::Write;
        std::io::stderr().flush().ok();
    }
    eprintln!();

    Ok(accessions
        .iter()
        .filter_map(|id| by_accession.remove(id.as_str()).map(|e| (id.clone(), e)))
        .collect())
}

pub(crate) fn get_ensembl_mapping(cross_refernces: &[UniProtCrossReference]) -> Result<Vec<(UniProtId, EnsemblId)>>{
    let mut result = Vec::new();
    for cross_ref in cross_refernces {
        if cross_ref.database == "Ensembl" {
            if let Some(iso_id) = &cross_ref.isoform_id {
                let iso_id: UniProtId = iso_id.parse().unwrap();
                let id: EnsemblId = cross_ref.id.parse().unwrap();
                result.push((iso_id, id))
            }
        }
    }
    Ok(result)
}
