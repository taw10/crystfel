/*
 * integer_matrix.h
 *
 * A small integer matrix library
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2012 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef INTEGER_MATRIX_H
#define INTEGER_MATRIX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/**
 * \file integer_matrix.h
 * Matrix type containing only integers
 */

/**
 * The \p IntegerMatrix is an opaque data structure representing an integer matrix.
 **/
typedef struct _integermatrix IntegerMatrix;

#include "rational.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Alloc/dealloc */
extern IntegerMatrix *intmat_new(unsigned int rows, unsigned int cols);
extern IntegerMatrix *intmat_copy(const IntegerMatrix *m);
extern IntegerMatrix *intmat_identity(int size);
extern void intmat_free(IntegerMatrix *m);

/* Get/set */
extern void intmat_size(const IntegerMatrix *m, unsigned int *rows,
                        unsigned int *cols);

extern void intmat_set(IntegerMatrix *m, unsigned int i, unsigned int j,
                       signed int v);
extern signed int intmat_get(const IntegerMatrix *m,
                             unsigned int i, unsigned int j);

extern void intmat_zero(IntegerMatrix *m);

extern IntegerMatrix *intmat_create_3x3(signed int m11, signed int m12, signed int m13,
                                        signed int m21, signed int m22, signed int m23,
                                        signed int m31, signed int m32, signed int m33);

/* Matrix-vector multiplication */
extern signed int *transform_indices(const IntegerMatrix *P, const signed int *hkl);

/* Matrix-matrix multiplication */
extern IntegerMatrix *intmat_times_intmat(const IntegerMatrix *a,
                                          const IntegerMatrix *b);

/* Inverse */
extern IntegerMatrix *intmat_inverse(const IntegerMatrix *m);

/* Determinant */
extern signed int intmat_det(const IntegerMatrix *m);

/* Is identity? */
extern int intmat_is_identity(const IntegerMatrix *m);

/* Is inversion? */
extern int intmat_is_inversion(const IntegerMatrix *m);

/* Comparison */
extern int intmat_equals(const IntegerMatrix *a, const IntegerMatrix *b);

/* Diagnostics */
extern void intmat_print(const IntegerMatrix *m);

#ifdef __cplusplus
}
#endif

#endif	/* INTEGER_MATRIX_H */
