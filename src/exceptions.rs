// ============================================================================
// Exception log: a persistent, tab-separated record of skipped or malformed
// input (bad variant lines, failed Ensembl fetches, ...), each tagged with an
// identifier (source line, ENST id, ...) so it can be traced back later.
// ============================================================================

use anyhow::{Context, Result};
use std::fs::File;
use std::io::Write;

pub(crate) struct ExceptionLog {
    writer: Box<dyn Write>,
}

impl ExceptionLog {
    pub(crate) fn to_file(path: &str) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("failed to create exceptions file '{}'", path))?;
        Ok(ExceptionLog { writer: Box::new(file) })
    }

    pub(crate) fn log(&mut self, identifier: &str, message: &str) {
        eprintln!("Warning: [{}] {}", identifier, message);
        let _ = writeln!(self.writer, "{}\t{}", identifier, message);
    }
}
