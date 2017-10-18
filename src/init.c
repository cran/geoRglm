#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "geoRglm.h"


static const R_CMethodDef cMethods[] = {
  {"binitprod", (DL_FUNC) &binitprod, 9},
  {"mcmcrun4binom", (DL_FUNC) &mcmcrun4binom, 9},
  {"mcmcrun4boxcox", (DL_FUNC) &mcmcrun4boxcox, 9},
  {"mcmcrun5", (DL_FUNC) &mcmcrun5, 9},
  {"mcmcrun5binom", (DL_FUNC) &mcmcrun5binom, 9},
  {"mcmcrun5boxcox", (DL_FUNC) &mcmcrun5boxcox, 9},
  {"mcmc1binom", (DL_FUNC) &mcmc1binom, 9},
  {"mcmc1poislog", (DL_FUNC) &mcmc1poislog, 9},
  {"mcmc1poisboxcox", (DL_FUNC) &mcmc1poisboxcox, 9},
  {NULL, NULL, 0}
};


void R_init_mypkg(DllInfo *dll)
{
   R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
   R_useDynamicSymbols(dll, FALSE);
}
