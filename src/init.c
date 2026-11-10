#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP C_scoutjoy_york_flip(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"scoutjoy_york_flip", (DL_FUNC) &C_scoutjoy_york_flip, 8},
    {NULL, NULL, 0}
};

void attribute_visible R_init_SCOUTJOY(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
