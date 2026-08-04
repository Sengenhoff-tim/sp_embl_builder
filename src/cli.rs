// ============================================================================
// CLI argument definition
// ============================================================================

use crate::util::{DEFAULT_MAX_ATTEMPTS, DEFAULT_RETRY_BACKOFF_MS, DEFAULT_TIMEOUT_SECS};

#[derive(clap::Parser)]
#[command(about = "Map variants from UniProt/EBI onto isoform sequences")]
pub(crate) struct Cli {
    /// Path to accession list file
    #[arg(long)]
    pub(crate) accessions: String,

    /// Path to ENST variant file. If omitted, no ENST-derived sample variants
    /// are included (only UniProt/EBI variants, per --uniprot-variants /
    /// --ebi-variants).
        #[arg(long)]
    pub(crate) variants: Option<String>,

    /// Path to write the synthetic flat-file output to. Defaults to stdout.
    #[arg(long)]
    pub(crate) output: Option<String>,

    /// Filter EBI variants by sourceType (comma-separated, up to 2).
    /// Allowed values: uniprot, large scale study, mixed, clinvar, nci-tcga,
    /// cosmic curated, ensembl, gnomad, topmed, exac.
    #[arg(long, value_delimiter = ',')]
    pub(crate) source_type: Vec<String>,

    /// Whether to fetch and include EBI-derived variants. If false, EBI is
    /// never queried and only UniProt/ENST-derived variants are used.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub(crate) ebi_variants: bool,

    /// Whether to include UniProt-annotated variants. If false, UniProt's own
    /// variant annotations are excluded (ENST/EBI-derived variants from other
    /// sources are unaffected by this flag).
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub(crate) uniprot_variants: bool,

    /// Path to write skipped/malformed-input exceptions to, as tab-separated
    /// `<identifier>\t<message>` lines.
    #[arg(long, default_value = "exceptions.log")]
    pub(crate) exceptions: String,

    /// Per-request timeout, in seconds, for outgoing REST calls (UniProt,
    /// EBI, Ensembl). Covers connect + write + read for a single attempt;
    /// each attempt is retried independently on failure.
    #[arg(long, default_value_t = DEFAULT_TIMEOUT_SECS)]
    pub(crate) timeout_secs: u64,

    /// Number of attempts (including the first) before a failed fetch is
    /// treated as fatal and the program exits.
    #[arg(long, default_value_t = DEFAULT_MAX_ATTEMPTS)]
    pub(crate) max_attempts: u32,

    /// Base backoff, in milliseconds, between retry attempts. Doubles after
    /// each failed attempt (e.g. 500ms -> 1s -> 2s -> ...).
    #[arg(long, default_value_t = DEFAULT_RETRY_BACKOFF_MS)]
    pub(crate) retry_backoff_ms: u64,
}
