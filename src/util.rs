// ============================================================================
// Utility: retry/backoff helper and a global async rate limiter
// ============================================================================

use anyhow::{Context, Result};
use std::future::Future;
use std::time::Duration;
use tokio::sync::Mutex;

pub(crate) fn strip_version(enst: &str) -> &str {
    enst.split('.').next().unwrap_or(enst)
}

/// Default per-request timeout (connect + write + read) applied to every
/// outgoing REST call, overridable via `--timeout-secs`.
/// Generous on purpose, as an invalid fetch is irrecoverable.
pub(crate) const DEFAULT_TIMEOUT_SECS: u64 = 30;

/// Default number of attempts (including the first) before a fetch is
/// treated as fatal, overridable via `--max-attempts`.
pub(crate) const DEFAULT_MAX_ATTEMPTS: u32 = 5;

/// Default base backoff (doubled after each failed attempt), overridable
/// via `--retry-backoff-ms`.
pub(crate) const DEFAULT_RETRY_BACKOFF_MS: u64 = 500;

/// Bundles the retry/timeout knobs exposed on the CLI so fetch functions
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

/// Async retry helper. `make_request` is called fresh on every attempt
/// (it must be `FnMut` returning a new future each time, since futures
/// are single-use and can't be re-awaited).
pub(crate) async fn with_retries<T, Fut>(
    description: &str,
    config: &RetryConfig,
    mut make_request: impl FnMut() -> Fut,
) -> Result<T>
where
    Fut: Future<Output = Result<T>>,
{
    let mut attempt = 0;
    loop {
        attempt += 1;
        match make_request().await {
            Ok(value) => return Ok(value),
            Err(err) if attempt < config.max_attempts => {
                let backoff_ms = config.backoff_base_ms * 2u64.pow(attempt - 1);
                eprintln!(
                    "Warning: {} failed (attempt {}/{}), retrying in {}ms: {:?}",
                    description, attempt, config.max_attempts, backoff_ms, err
                );
                tokio::time::sleep(Duration::from_millis(backoff_ms)).await;
            }
            Err(err) => {
                return Err(err).with_context(|| {
                    format!("{} failed after {} attempts", description, attempt)
                });
            }
        }
    }
}

/// Global, shareable async rate limiter: enforces at most `max_per_second`
/// admissions per second, *in total*, across every caller that holds a
/// clone of it. Cloning is cheap (bumps an `Arc` refcount) — construct one
/// `RateLimiter` per API and clone it into every concurrently-spawned task
/// that calls that API, so the limit is enforced globally rather than
/// per-task.
#[derive(Clone)]
pub(crate) struct RateLimiter {
    inner: std::sync::Arc<Mutex<RateLimiterState>>,
}

struct RateLimiterState {
    min_interval: Duration,
    next_allowed: tokio::time::Instant,
}

impl RateLimiter {
    pub(crate) fn per_second(max_per_second: u32) -> Self {
        RateLimiter {
            inner: std::sync::Arc::new(Mutex::new(RateLimiterState {
                min_interval: Duration::from_secs_f64(1.0 / max_per_second as f64),
                next_allowed: tokio::time::Instant::now(),
            })),
        }
    }

    /// Reserves the next available slot and waits for it. The reservation
    /// (bumping `next_allowed`) happens while the lock is held, so
    /// concurrently-waiting callers are queued strictly in order and each
    /// gets a slot exactly `min_interval` after the previous one — no two
    /// callers can slip through the same window even under a burst of
    /// simultaneous calls.
    pub(crate) async fn throttle(&self) {
        let wait_until = {
            let mut state = self.inner.lock().await;
            let now = tokio::time::Instant::now();
            let scheduled = state.next_allowed.max(now);
            state.next_allowed = scheduled + state.min_interval;
            scheduled
        };
        tokio::time::sleep_until(wait_until).await;
    }
}
