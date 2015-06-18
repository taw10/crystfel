

#ifndef ASDF_H
#define ASDF_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "index.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int run_asdf(struct image *image, IndexingPrivate *ipriv);

extern IndexingPrivate *asdf_prepare(IndexingMethod *indm,
                                      UnitCell *cell, struct detector *det,
                                      float *ltl);

extern void asdf_cleanup(IndexingPrivate *pp);

#ifdef __cplusplus
}
#endif

#endif	/* DIRAX_H */
