#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* Native routines from oregMclust.cc (vendored from the 'edci' package).
   Registered under the treeX package DLL so that
   .C("c_oregMcirc", ...) resolves on every platform, including macOS. */

/* .C calls */
extern void c_oregMcirc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void c_oregMclust(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"c_oregMcirc",  (DL_FUNC) &c_oregMcirc,  25},
  {"c_oregMclust", (DL_FUNC) &c_oregMclust, 14},
  {NULL, NULL, 0}
};

void R_init_treeX(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
