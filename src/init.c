#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP dantzig_ladm_scr(SEXP, SEXP, SEXP, SEXP, SEXP,
                             SEXP, SEXP, SEXP, SEXP, SEXP,
                             SEXP, SEXP, SEXP, SEXP, SEXP,
                             SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP slim_dantzig_ladm_scr(SEXP, SEXP, SEXP, SEXP, SEXP,
                                  SEXP, SEXP, SEXP, SEXP, SEXP,
                                  SEXP, SEXP, SEXP, SEXP, SEXP,
                                  SEXP, SEXP, SEXP, SEXP, SEXP,
                                  SEXP);


extern SEXP slim_dantzig_ladm_scr2(SEXP, SEXP, SEXP, SEXP, SEXP,
                                  SEXP, SEXP, SEXP, SEXP, SEXP,
                                  SEXP, SEXP, SEXP, SEXP, SEXP,
                                  SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"dantzig_ladm_scr",            (DL_FUNC) &dantzig_ladm_scr,           20},
  {"slim_dantzig_ladm_scr",    (DL_FUNC) &slim_dantzig_ladm_scr,    21},
  {"slim_dantzig_ladm_scr2",    (DL_FUNC) &slim_dantzig_ladm_scr,    18},
  {NULL, NULL, 0}
};

void R_init_sitmo(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

