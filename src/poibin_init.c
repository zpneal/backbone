#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* .C calls */
extern void multi_bin_bh(void *, void *, void *, void *, void *, void *);
extern void multi_bin_dft_cf(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"multi_bin_bh",     (DL_FUNC) &multi_bin_bh,      6},
    {"multi_bin_dft_cf", (DL_FUNC) &multi_bin_dft_cf, 11},
    {NULL, NULL, 0}
};

void R_init_backbone(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
