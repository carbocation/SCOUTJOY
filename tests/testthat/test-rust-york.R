# The rust engine is a like-for-like port of yorktools::york(intercept = 0,
# gridstartvals = TRUE); the R code is the reference.

make_york_inputs <- function(n = 60L, slope = 1.2, re = 0.2, n_outliers = 0L, seed = 1L) {
  set.seed(seed)
  true_x <- stats::runif(n, 3, 9)
  sx <- stats::runif(n, 0.7, 1.3)
  sy <- stats::runif(n, 0.7, 1.3)
  ex <- stats::rnorm(n)
  ey <- re * ex + sqrt(1 - re^2) * stats::rnorm(n)
  d <- data.frame(x = true_x + sx * ex, sx = sx, y = slope * true_x + sy * ey, sy = sy)
  if (n_outliers > 0L) {
    shift <- rep(c(8, -7, 9, -8, 7, -9), length.out = n_outliers)
    d$y[seq_len(n_outliers)] <- d$y[seq_len(n_outliers)] + shift
  }
  rownames(d) <- paste0("v", seq_len(n))
  d
}

run_scoutjoy <- function(d, engine, ...) {
  suppressWarnings(scoutjoy("y", "x", "sy", "sx", d, CovIntercept = 0.2, printProgress = FALSE,
                            engine = engine, ...))
}

test_that("rust York fits agree with yorktools on arbitrary subsets", {
  d <- make_york_inputs(n_outliers = 4L, seed = 21L)
  d$re <- 0.2
  excluded <- rep(FALSE, nrow(d))
  excluded[c(1, 2, 7)] <- TRUE
  flip <- 0:nrow(d)
  reference <- SCOUTJOY:::york_flip_fits(d, excluded, flip, engine = "R")
  rust <- SCOUTJOY:::york_flip_fits(d, excluded, flip, engine = "rust")
  expect_equal(rust$b, reference$b, tolerance = 1e-10)
  expect_equal(rust$chi2, reference$chi2, tolerance = 1e-10)
})

test_that("scoutjoy gives the same results with either engine", {
  for (n_outliers in c(0L, 2L, 6L)) {
    for (null in c("estimate", "fitted")) {
      d <- make_york_inputs(n_outliers = n_outliers, seed = 22L + n_outliers)
      expect_equal(run_scoutjoy(d, "rust", OutlierNull = null),
                   run_scoutjoy(d, "R", OutlierNull = null), tolerance = 1e-8)
    }
  }
})

test_that("rust engine results do not depend on the number of threads", {
  d <- make_york_inputs(n = 120L, n_outliers = 4L, seed = 30L)
  expect_identical(run_scoutjoy(d, "rust", threads = 1), run_scoutjoy(d, "rust", threads = 4))
})

test_that("engine is validated", {
  d <- make_york_inputs(seed = 31L)
  expect_error(run_scoutjoy(d, "fortran"), "engine must be")
})
