#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP EBS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP g_ref(SEXP);
extern SEXP isref(SEXP, SEXP);
extern SEXP MCP(SEXP, SEXP, SEXP);
extern SEXP pat(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"EBS",   (DL_FUNC) &EBS,   9},
    {"g_ref", (DL_FUNC) &g_ref, 1},
    {"isref", (DL_FUNC) &isref, 2},
    {"MCP",   (DL_FUNC) &MCP,   3},
    {"pat",   (DL_FUNC) &pat,   1},
    {NULL, NULL, 0}
};

void R_init_scDDboost(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
