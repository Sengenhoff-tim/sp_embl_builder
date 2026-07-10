// ============================================================================
// Utility: retry/backoff helper and a simple rate limiter
// ============================================================================

use anyhow::{Context, Result};
use std::thread;
use std::time::Duration;

pub(crate) fn strip_version(enst: &str) -> &str {
    enst.split('.').next().unwrap_or(enst)
}

/// Default per-request timeout (connect + write + read) applied to every
/// outgoing REST call, overridable via `--timeout-secs`. Generous on
/// purpose: these APIs are occasionally slow, and a batch that exhausts
/// `with_retries` now aborts the whole run, so it's worth waiting a few
/// extra seconds before writing a request off as dead.
pub(crate) const DEFAULT_TIMEOUT_SECS: u64 = 5;

/// Default number of attempts (including the first) before a fetch is
/// treated as fatal, overridable via `--max-attempts`.
pub(crate) const DEFAULT_MAX_ATTEMPTS: u32 = 5;

/// Default base backoff (doubled after each failed attempt), overridable
/// via `--retry-backoff-ms`.
pub(crate) const DEFAULT_RETRY_BACKOFF_MS: u64 = 500;

/// Bundles the retry/timeout knobs exposed on the CLI so fetch functions
/// take one value instead of a growing list of parameters.
#[derive(Clone, Copy)]
pub(crate) struct RetryConfig {
    pub(crate) timeout: Duration,
    pub(crate) max_attempts: u32,
    pub(crate) backoff_base_ms: u64,
}

impl Default for RetryConfig {
    fn default() -> Self {
        RetryConfig {
            timeout: Duration::from_secs(DEFAULT_TIMEOUT_SECS),
            max_attempts: DEFAULT_MAX_ATTEMPTS,
            backoff_base_ms: DEFAULT_RETRY_BACKOFF_MS,
        }
    }
}

pub(crate) fn with_retries<T>(
    description: &str,
    config: &RetryConfig,
    mut make_request: impl FnMut() -> Result<T>,
) -> Result<T> {
    let mut attempt = 0;
    loop {
        attempt += 1;
        match make_request() {
            Ok(value) => return Ok(value),
            Err(err) if attempt < config.max_attempts => {
                let backoff_ms = config.backoff_base_ms * 2u64.pow(attempt - 1);
                eprintln!(
                    "Warning: {} failed (attempt {}/{}), retrying in {}ms: {:?}",
                    description, attempt, config.max_attempts, backoff_ms, err
                );
                thread::sleep(Duration::from_millis(backoff_ms));
            }
            Err(err) => {
                return Err(err).with_context(|| {
                    format!("{} failed after {} attempts", description, attempt)
                });
            }
        }
    }
}

pub(crate) struct RateLimiter {
    min_interval: Duration,
    last_call: Option<std::time::Instant>,
}

impl RateLimiter {
    pub(crate) fn per_second(max_per_second: u32) -> Self {
        RateLimiter {
            min_interval: Duration::from_secs_f64(1.0 / max_per_second as f64),
            last_call: None,
        }
    }

    pub(crate) fn throttle(&mut self) {
        if let Some(last) = self.last_call {
            let elapsed = last.elapsed();
            if elapsed < self.min_interval {
                thread::sleep(self.min_interval - elapsed);
            }
        }
        self.last_call = Some(std::time::Instant::now());
    }
}
