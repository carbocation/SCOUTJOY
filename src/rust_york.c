#include <R.h>
#include <Rinternals.h>
#include <stddef.h>

extern int scoutjoy_york_flip(
    const double *x,
    const double *sx,
    const double *y,
    const double *sy,
    const double *re,
    size_t n,
    const int *excluded,
    const int *flip,
    size_t n_fits,
    size_t requested_threads,
    double *b_out,
    double *chi2_out,
    int *status_out);

static void check_real_vector(SEXP value, R_xlen_t length, const char *name) {
    if (!isReal(value) || XLENGTH(value) != length) {
        error("'%s' must be a double vector of length %lld", name, (long long) length);
    }
}

/* Zero-intercept York regressions on subsets of the data. `excluded` marks the
   rows left out of the base subset. Each element of `flip` asks for one fit:
   0 is the base subset, and i > 0 is the base subset with row i's membership
   reversed (left out if it was in, put back if it was out).
   Returns list(b, chi2, status); status 0 = converged, 1 = did not converge,
   2 = a straight line cannot be fitted. */
SEXP C_scoutjoy_york_flip(
    SEXP x, SEXP sx, SEXP y, SEXP sy, SEXP re, SEXP excluded, SEXP flip, SEXP threads) {

    if (!isReal(x) || XLENGTH(x) < 1) {
        error("'x' must be a non-empty double vector");
    }
    R_xlen_t n = XLENGTH(x);
    check_real_vector(sx, n, "sx");
    check_real_vector(y, n, "y");
    check_real_vector(sy, n, "sy");
    check_real_vector(re, n, "re");
    if (!isLogical(excluded) || XLENGTH(excluded) != n) {
        error("'excluded' must be a logical vector with one element per row");
    }
    if (!isInteger(flip) || XLENGTH(flip) < 1) {
        error("'flip' must be a non-empty integer vector");
    }
    R_xlen_t n_fits = XLENGTH(flip);
    const int *excluded_values = LOGICAL(excluded);
    for (R_xlen_t i = 0; i < n; i++) {
        if (excluded_values[i] == NA_LOGICAL) {
            error("'excluded' must not contain missing values");
        }
    }
    const int *flip_values = INTEGER(flip);
    for (R_xlen_t i = 0; i < n_fits; i++) {
        if (flip_values[i] == NA_INTEGER || flip_values[i] < 0 || flip_values[i] > n) {
            error("'flip' must hold row numbers, or 0 for the base subset");
        }
    }
    int threads_value = asInteger(threads);
    if (threads_value == NA_INTEGER || threads_value < 1) {
        error("'threads' must be a positive whole number");
    }

    SEXP b_out = PROTECT(allocVector(REALSXP, n_fits));
    SEXP chi2_out = PROTECT(allocVector(REALSXP, n_fits));
    SEXP status_out = PROTECT(allocVector(INTSXP, n_fits));

    int failed = scoutjoy_york_flip(
        REAL(x), REAL(sx), REAL(y), REAL(sy), REAL(re), (size_t) n,
        excluded_values, flip_values, (size_t) n_fits, (size_t) threads_value,
        REAL(b_out), REAL(chi2_out), INTEGER(status_out));
    if (failed != 0) {
        UNPROTECT(3);
        error("the Rust York routine failed unexpectedly");
    }

    const char *names[] = {"b", "chi2", "status", ""};
    SEXP output = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(output, 0, b_out);
    SET_VECTOR_ELT(output, 1, chi2_out);
    SET_VECTOR_ELT(output, 2, status_out);
    UNPROTECT(4);
    return output;
}
