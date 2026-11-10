#!/usr/bin/env Rscript
# Time scoutjoy() with the R and rust engines on simulated effect sizes and
# report how closely they agree. Usage: Rscript benchmarks/benchmark-rust-york.R [n_variants]

suppressPackageStartupMessages(library(SCOUTJOY))

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1) as.integer(args[1]) else 200L

for (n_outliers in c(0L, 6L)) {
  set.seed(1L)
  true_x <- stats::runif(n, 3, 9)
  d <- data.frame(x = true_x + stats::rnorm(n), sx = 1, y = 1.2 * true_x + stats::rnorm(n), sy = 1)
  if (n_outliers > 0L) d$y[seq_len(n_outliers)] <- d$y[seq_len(n_outliers)] + c(8, -7, 9, -8, 7, -9)

  run <- function(engine) {
    seconds <- system.time(fit <- suppressWarnings(
      scoutjoy("y", "x", "sy", "sx", d, printProgress = FALSE, engine = engine)
    ))[["elapsed"]]
    list(fit = fit, seconds = seconds)
  }
  r <- run("R")
  rust <- run("rust")
  p_r <- r$fit[["Outlier Tests"]]$Pvalue
  p_rust <- rust$fit[["Outlier Tests"]]$Pvalue
  cat(sprintf(
    "%d variants, %d planted outliers: R %.3fs, rust %.3fs (%.0fx); largest difference in outlier p-values %.1e\n",
    n, n_outliers, r$seconds, rust$seconds, r$seconds / rust$seconds, max(abs(p_r - p_rust), na.rm = TRUE)
  ))
}
