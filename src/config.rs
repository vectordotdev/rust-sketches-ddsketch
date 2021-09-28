const DEFAULT_MAX_BINS: u32 = 2048;
const DEFAULT_ALPHA: f64 = 0.01;
const DEFAULT_MIN_VALUE: f64 = 1.0e-9;

const DDAGENT_DEFAULT_MAX_BINS: u32 = 4096;
const DDAGENT_DEFAULT_ALPHA: f64 = 1.0 / 128.0;
const DDAGENT_DEFAULT_MIN_VALUE: f64 = 1.0e-9;

/// The configuration struct for constructing a `DDSketch`
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Config {
    pub(crate) max_num_bins: u32,
    pub(crate) gamma: f64,
    gamma_ln: f64,
    pub(crate) min_value: f64,
    pub(crate) offset: i32,
}

#[inline]
fn log_gamma(value: f64, gamma_ln: f64) -> f64 {
    value.ln() / gamma_ln
}

impl Config {
    /// Construct a new `Config` struct with specific parameters. If you are unsure of how to
    /// configure this, the `defaults` method constructs a `Config` with built-in defaults.
    ///
    /// `max_num_bins` is the max number of bins the DDSketch will grow to, in steps of 128 bins.
    pub fn new(alpha: f64, max_num_bins: u32, min_value: f64) -> Self {
        let gamma_ln = (2.0 * alpha) / (1.0 - alpha);
        let gamma_ln = gamma_ln.ln_1p();

        Self {
            max_num_bins,
            gamma: 1.0 + (2.0 * alpha) / (1.0 - alpha),
            gamma_ln,
            min_value,
            offset: 1 - (log_gamma(min_value, gamma_ln) as i32),
        }
    }

    /// Return a `Config` using built-in default settings
    pub fn defaults() -> Self {
        Self::new(DEFAULT_ALPHA, DEFAULT_MAX_BINS, DEFAULT_MIN_VALUE)
    }

    /// Return a `Config` using built-in default settings that the open-source Datadog Agent uses.
    /// These are subtly different from the parameters in the open-source implementations of the
    /// white paper.
    ///
    /// Use this to more closely model what Datadog uses themselves.
    pub fn agent_defaults() -> Self {
        Self::new(
            DDAGENT_DEFAULT_ALPHA,
            DDAGENT_DEFAULT_MAX_BINS,
            DDAGENT_DEFAULT_MIN_VALUE,
        )
    }

    #[inline]
    pub fn key(&self, v: f64) -> i32 {
        if v < -self.min_value {
            -(self.log_gamma(-v).ceil() as i32) - self.offset
        } else if v > self.min_value {
            (self.log_gamma(v).ceil() as i32) + self.offset
        } else {
            0
        }
    }

    pub fn min_key(&self) -> i32 {
        self.key(self.min_value)
    }

    #[inline]
    pub fn lower_bound(&self, key: i32) -> f64 {
        if key < 0 {
            return -self.lower_bound(-key);
        }

        if key == i32::max_value() {
            return f64::INFINITY;
        }

        if key == 0 {
            return 0.0;
        }

        self.pow_gamma(key - self.offset)
    }

    #[inline]
    pub fn log_gamma(&self, value: f64) -> f64 {
        log_gamma(value, self.gamma_ln)
    }

    #[inline]
    pub fn pow_gamma(&self, k: i32) -> f64 {
        ((k as f64) * self.gamma_ln).exp()
    }
}
