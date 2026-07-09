// ============================================================================
// Utility: retry/backoff helper and a simple rate limiter
// ============================================================================

use anyhow::{Context, Result};
use std::thread;
use std::time::Duration;

pub(crate) fn strip_version(enst: &str) -> &str {
    enst.split('.').next().unwrap_or(enst)
}

const MAX_FETCH_ATTEMPTS: u32 = 3;

pub(crate) fn with_retries<T>(description: &str, mut make_request: impl FnMut() -> Result<T>) -> Result<T> {
    let mut attempt = 0;
    loop {
        attempt += 1;
        match make_request() {
            Ok(value) => return Ok(value),
            Err(err) if attempt < MAX_FETCH_ATTEMPTS => {
                let backoff_ms = 100u64 * 2u64.pow(attempt - 1);
                eprintln!(
                    "Warning: {} failed (attempt {}/{}), retrying in {}ms: {:?}",
                    description, attempt, MAX_FETCH_ATTEMPTS, backoff_ms, err
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
