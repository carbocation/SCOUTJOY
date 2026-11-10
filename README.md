# SCOUTJOY

Significant Cross-trait OUtliers and Trends in JOint York regression

# Installation:

This R package can be installed using the `devtools` package as follows:

```
require(devtools)
devtools::install_github("aelliott08/SCOUTJOY", ref="main")
```

## Rust engine branch

The `dev-jp-rust` branch of `carbocation/SCOUTJOY` adds an optional compiled
engine for the analytic outlier test (`scoutjoy(..., engine = "rust")`). The
default engine is still the R code, and the rust engine is a like-for-like port
of the regression it runs, so the two give the same results. Building this
branch from source requires a Rust toolchain with `cargo` and `rustc` on `PATH`
or in `~/.cargo/bin` (see https://rustup.rs/). Rust is needed to build the
package but not to run it.

```
devtools::install_github("carbocation/SCOUTJOY", ref="dev-jp-rust")
```

`benchmarks/benchmark-rust-york.R` times the two engines against each other.

# Usage

See the documentation in R using `?scoutjoy` 

# Citation

Methods details can be found in:

> Elliott, A., Walters, R.K., Pirinen, M. et al. Distinct and shared genetic architectures of gestational diabetes mellitus and type 2 diabetes. Nat Genet 56, 377–382 (2024). https://doi.org/10.1038/s41588-023-01607-4

We ask that you cite this paper when using SCOUTJOY.
