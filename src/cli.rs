// ============================================================================
// CLI argument definition
// ============================================================================

#[derive(clap::Parser)]
#[command(about = "Map variants from UniProt/EBI onto isoform sequences")]
pub(crate) struct Cli {
    /// Path to accession list file
    #[arg(long)]
    pub(crate) accessions: String,
    
    /// Path to ENST variant file
    #[arg(long)]
    pub(crate) variants: String,

    /// Path to write the synthetic flat-file output to. Defaults to stdout.
    #[arg(long)]
    pub(crate) output: Option<String>,

    /// Filter EBI variants by sourceType (comma-separated, up to 2).
    /// Allowed values: uniprot, large scale study, mixed, clinvar, nci-tcga,
    /// cosmic curated, ensembl, gnomad, topmed, exac
    #[arg(long, value_delimiter = ',')]
    pub(crate) source_type: Vec<String>,

    /// Path to write skipped/malformed-input exceptions to, as tab-separated
    /// `<identifier>\t<message>` lines.
    #[arg(long, default_value = "exceptions.log")]
    pub(crate) exceptions: String,
}
