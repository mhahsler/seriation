#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ar(SEXP, SEXP, SEXP);
extern SEXP bar(SEXP, SEXP, SEXP);
extern SEXP gradient(SEXP, SEXP, SEXP);
extern SEXP inertia_criterion(SEXP, SEXP);
extern SEXP lazy_path_length(SEXP, SEXP);
extern SEXP least_squares_criterion(SEXP, SEXP);
extern SEXP order_greedy(SEXP);
extern SEXP order_length(SEXP, SEXP);
extern SEXP order_optimal(SEXP, SEXP);
extern SEXP pathdist_floyd(SEXP);
extern SEXP reorder_dist(SEXP, SEXP);
extern SEXP rgar(SEXP, SEXP, SEXP, SEXP);
extern SEXP stress(SEXP, SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(arsa)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(bburcg)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(bbwrcg)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cbea)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(energy)(void *, void *, void *, void *);
extern void F77_NAME(rbea)(void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"ar",                      (DL_FUNC) &ar,                      3},
    {"bar",                     (DL_FUNC) &bar,                     3},
    {"gradient",                (DL_FUNC) &gradient,                3},
    {"inertia_criterion",       (DL_FUNC) &inertia_criterion,       2},
    {"lazy_path_length",        (DL_FUNC) &lazy_path_length,        2},
    {"least_squares_criterion", (DL_FUNC) &least_squares_criterion, 2},
    {"order_greedy",            (DL_FUNC) &order_greedy,            1},
    {"order_length",            (DL_FUNC) &order_length,            2},
    {"order_optimal",           (DL_FUNC) &order_optimal,           2},
    {"pathdist_floyd",          (DL_FUNC) &pathdist_floyd,          1},
    {"reorder_dist",            (DL_FUNC) &reorder_dist,            2},
    {"rgar",                    (DL_FUNC) &rgar,                    4},
    {"stress",                  (DL_FUNC) &stress,                  4},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"arsa",   (DL_FUNC) &F77_NAME(arsa),   15},
    {"bburcg", (DL_FUNC) &F77_NAME(bburcg), 10},
    {"bbwrcg", (DL_FUNC) &F77_NAME(bbwrcg), 10},
    {"cbea",   (DL_FUNC) &F77_NAME(cbea),    7},
    {"energy", (DL_FUNC) &F77_NAME(energy),  4},
    {"rbea",   (DL_FUNC) &F77_NAME(rbea),    7},
    {NULL, NULL, 0}
};

void R_init_seriation(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
