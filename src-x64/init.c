#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

void ADprobApproxInf(double *, int *, double *);
void ADprobExactInf(double *, int *, int *, double *);
void ADprobN(double *, int *, double *);
void ADtestR(double *, int *, double *, double *);

static const R_CallMethodDef CallEntries[] = {
  {"ADprobApproxInf",      (DL_FUNC) &ADprobApproxInf,      3},
  {"ADprobExactInf",       (DL_FUNC) &ADprobExactInf,       3},
  {"ADprobN",              (DL_FUNC) &ADprobN,              4},
  {"ADtestR",              (DL_FUNC) &ADtestR,              4},
  {NULL, NULL, 0}
};


void R_init_ECToolsad(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
