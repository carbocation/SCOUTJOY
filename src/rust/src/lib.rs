//! Zero-intercept York regression, ported from yorktools::york() as SCOUTJOY
//! calls it in its analytic outlier test:
//!     yorktools::york(data, intercept = 0, gridstartvals = TRUE)
//! (least-squares York, tol = 1e-12, maxit = 5000). Only the slope and chi2
//! are returned, which is all the outlier test uses.
//!
//! This is a like-for-like port: the same starting values, the same iteration
//! and stopping rule, with arithmetic in the same order as the R expressions.
//! The R code remains the reference; tests/testthat/test-rust-york.R compares
//! the two.
//!
//! The outlier test needs one regression per variant, each on the data with
//! that one variant added or removed. All of them are fitted in one call, and
//! since they are independent they can be spread over threads without changing
//! any result.

// Index loops are kept as they are: they mirror the R expressions.
#![allow(clippy::needless_range_loop)]

use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;
use std::thread;

const TOL: f64 = 1e-12;
const MAXIT: usize = 5000;

const CONVERGED: i32 = 0;
const NOT_CONVERGED: i32 = 1;
const CANNOT_FIT: i32 = 2;

struct Data<'a> {
    x: &'a [f64],
    sx: &'a [f64],
    y: &'a [f64],
    sy: &'a [f64],
    re: &'a [f64],
}

struct Fit {
    b: f64,
    chi2: f64,
    status: i32,
}

const FAILED: Fit = Fit {
    b: f64::NAN,
    chi2: f64::NAN,
    status: CANNOT_FIT,
};

/// Slope of a regression of y on x through the origin with the given weights:
/// what stats::lm(Y ~ 0 + X, weights = w) returns.
fn origin_slope(d: &Data, rows: &[usize], weight: impl Fn(usize) -> f64) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for &i in rows {
        let w = weight(i);
        num += w * d.x[i] * d.y[i];
        den += w * d.x[i] * d.x[i];
    }
    num / den
}

/// yorktools::york(data[rows, ], intercept = 0, gridstartvals = TRUE)
fn york_fit(d: &Data, rows: &[usize]) -> Fit {
    if rows.is_empty() {
        return FAILED;
    }

    // starting values: five weighted regressions through the origin
    let starts = [
        origin_slope(d, rows, |_| 1.0),
        origin_slope(d, rows, |i| 1.0 / (d.sy[i] * d.sy[i])),
        origin_slope(d, rows, |i| 1.0 / (d.sx[i] * d.sx[i])),
        origin_slope(d, rows, |i| 1.0 / (d.sx[i] * d.sx[i] + d.sy[i] * d.sy[i])),
        origin_slope(d, rows, |i| {
            1.0 / (d.sx[i] * d.sx[i] + d.sy[i] * d.sy[i] - 2.0 * (d.re[i] * d.sx[i] * d.sy[i]))
        }),
    ];
    if starts.iter().any(|v| !v.is_finite()) {
        return FAILED;
    }
    let mut maxb = 0.0_f64;
    for v in starts {
        maxb = maxb.max(v.abs());
    }

    // qq <- c(seq(-5*maxb, 5*maxb, length = 10), 0, the five starts)
    let from = -5.0 * maxb;
    let to = 5.0 * maxb;
    let by = (to - from) / 9.0;
    let mut qq = [0.0_f64; 16];
    qq[0] = from;
    for (i, q) in qq.iter_mut().enumerate().take(9).skip(1) {
        *q = from + i as f64 * by;
    }
    qq[9] = to;
    qq[10] = 0.0;
    qq[11..16].copy_from_slice(&starts);

    // which.min(init_loss): the first smallest value, ignoring NaN
    let mut best: Option<(usize, f64)> = None;
    for (c, &q) in qq.iter().enumerate() {
        let mut loss = 0.0;
        for &i in rows {
            let wx = 1.0 / (d.sx[i] * d.sx[i]);
            let wy = 1.0 / (d.sy[i] * d.sy[i]);
            let w = wx * wy / (wx + (q * q) * wy - 2.0 * q * d.re[i] * (wx * wy).sqrt());
            let resid = d.y[i] - d.x[i] * q;
            loss += w * (resid * resid);
        }
        if loss.is_nan() {
            continue;
        }
        match best {
            Some((_, value)) if loss >= value => {}
            _ => best = Some((c, loss)),
        }
    }
    let mut b = match best {
        Some((c, _)) => qq[c],
        None => return FAILED,
    };

    let mut chi2 = f64::NAN;
    let mut converged = false;
    for _ in 0..MAXIT {
        let bold = b;
        let mut num = 0.0;
        let mut den = 0.0;
        for &i in rows {
            let sx2 = d.sx[i] * d.sx[i];
            let sy2 = d.sy[i] * d.sy[i];
            let cov = d.re[i] * d.sx[i] * d.sy[i];
            let w = 1.0 / (sy2 + (bold * bold) * sx2 - 2.0 * bold * cov);
            let (x, v) = (d.x[i], d.y[i]);
            num += w * w * v * (x * sy2 + bold * v * sx2 - v * cov);
            den += w * w * x * (x * sy2 + bold * v * sx2 - bold * x * cov);
        }
        b = num / den;
        // loss <- sum(W * (Y - X*b - a)^2), with W from the previous slope
        let mut loss = 0.0;
        for &i in rows {
            let sx2 = d.sx[i] * d.sx[i];
            let sy2 = d.sy[i] * d.sy[i];
            let cov = d.re[i] * d.sx[i] * d.sy[i];
            let w = 1.0 / (sy2 + (bold * bold) * sx2 - 2.0 * bold * cov);
            let resid = d.y[i] - d.x[i] * b - 0.0;
            loss += w * (resid * resid);
        }
        chi2 = loss;
        let change = (bold / b - 1.0) * (bold / b - 1.0);
        if change.is_nan() {
            // R: `if ((bold/b - 1)^2 < tol)` fails on NA
            return FAILED;
        }
        if change < TOL {
            converged = true;
            break;
        }
    }

    if converged {
        Fit {
            b,
            chi2,
            status: CONVERGED,
        }
    } else {
        Fit {
            b,
            chi2: f64::NAN,
            status: NOT_CONVERGED,
        }
    }
}

fn worker_count(n_fits: usize, requested_threads: usize) -> usize {
    let available = thread::available_parallelism()
        .map(|value| value.get())
        .unwrap_or(1);
    requested_threads.max(1).min(available).min(n_fits.max(1))
}

/// Fit the base subset (rows not excluded) and that subset with one row's
/// membership reversed, once per element of `flip` (0 = base subset, i > 0 =
/// reverse row i, 1-based). Returns 0, or 1 if the routine panicked.
///
/// # Safety
/// All pointers must be valid for the lengths implied by n and n_fits.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn scoutjoy_york_flip(
    x: *const f64,
    sx: *const f64,
    y: *const f64,
    sy: *const f64,
    re: *const f64,
    n: usize,
    excluded: *const i32,
    flip: *const i32,
    n_fits: usize,
    requested_threads: usize,
    b_out: *mut f64,
    chi2_out: *mut f64,
    status_out: *mut i32,
) -> i32 {
    let result = catch_unwind(AssertUnwindSafe(|| {
        let data = Data {
            x: slice::from_raw_parts(x, n),
            sx: slice::from_raw_parts(sx, n),
            y: slice::from_raw_parts(y, n),
            sy: slice::from_raw_parts(sy, n),
            re: slice::from_raw_parts(re, n),
        };
        let excluded = slice::from_raw_parts(excluded, n);
        let flip = slice::from_raw_parts(flip, n_fits);
        let data = &data;

        let fit_one = |f: usize| -> Fit {
            let target = flip[f];
            let rows: Vec<usize> = (0..n)
                .filter(|&i| {
                    let out = excluded[i] != 0;
                    let reversed = target > 0 && i + 1 == target as usize;
                    out == reversed
                })
                .collect();
            york_fit(data, &rows)
        };

        let workers = worker_count(n_fits, requested_threads);
        let mut fits: Vec<Fit> = Vec::with_capacity(n_fits);
        if workers <= 1 {
            for f in 0..n_fits {
                fits.push(fit_one(f));
            }
        } else {
            let chunk = n_fits.div_ceil(workers);
            let fit_one = &fit_one;
            let pieces: Vec<Vec<Fit>> = thread::scope(|scope| {
                let handles: Vec<_> = (0..workers)
                    .map(|w| {
                        scope.spawn(move || {
                            let from = (w * chunk).min(n_fits);
                            let to = ((w + 1) * chunk).min(n_fits);
                            (from..to).map(fit_one).collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .map(|handle| handle.join().expect("York worker panicked"))
                    .collect()
            });
            for piece in pieces {
                fits.extend(piece);
            }
        }

        let b_out = slice::from_raw_parts_mut(b_out, n_fits);
        let chi2_out = slice::from_raw_parts_mut(chi2_out, n_fits);
        let status_out = slice::from_raw_parts_mut(status_out, n_fits);
        for (f, fit) in fits.iter().enumerate() {
            b_out[f] = fit.b;
            chi2_out[f] = fit.chi2;
            status_out[f] = fit.status;
        }
    }));
    match result {
        Ok(()) => 0,
        Err(_) => 1,
    }
}
