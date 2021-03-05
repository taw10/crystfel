/*
 * rational.h
 *
 * A small rational number library
 *
 * Copyright © 2019-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019 Thomas White <taw@physics.org>
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

#ifndef RATIONAL_H
#define RATIONAL_H

/**
 * \file rational.h
 * %Rational numbers (including rational matrices)
 */

/**
 * The Rational is an opaque-ish data structure representing a rational number.
 *
 * "Opaque-ish" means that the structure isn't technically opaque, allowing you
 * to assign and allocate them easily.  But you shouldn't look at or set its
 * contents, except by using the accessor functions.
 **/
typedef struct {
	/* Private, don't modify */
	signed long long int num;
	signed long long int den;
} Rational;


/**
 * The RationalMatrix is an opaque data structure representing a matrix of
 * rational numbers.
 **/
typedef struct _rationalmatrix RationalMatrix;

#include "integer_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

extern Rational rtnl_zero(void);
extern Rational rtnl(signed long long int num, signed long long int den);
extern double rtnl_as_double(Rational r);

extern Rational rtnl_mul(Rational a, Rational b);
extern Rational rtnl_div(Rational a, Rational b);
extern Rational rtnl_add(Rational a, Rational b);
extern Rational rtnl_sub(Rational a, Rational b);

extern signed int rtnl_cmp(Rational a, Rational b);
extern Rational rtnl_abs(Rational a);

extern char *rtnl_format(Rational rt);

extern Rational *rtnl_list(signed int num_min, signed int num_max,
                           signed int den_min, signed int den_max,
                           int *pn);

extern RationalMatrix *rtnl_mtx_new(unsigned int rows, unsigned int cols);
extern RationalMatrix *rtnl_mtx_copy(const RationalMatrix *m);
extern Rational rtnl_mtx_get(const RationalMatrix *m, int i, int j);
extern void rtnl_mtx_set(const RationalMatrix *m, int i, int j, Rational v);
extern RationalMatrix *rtnl_mtx_from_intmat(const IntegerMatrix *m);
extern RationalMatrix *rtnl_mtx_identity(int rows);
extern IntegerMatrix *intmat_from_rtnl_mtx(const RationalMatrix *m);
extern void rtnl_mtx_free(RationalMatrix *mtx);
extern RationalMatrix *rtnlmtx_times_rtnlmtx(const RationalMatrix *A,
                                             const RationalMatrix *B);
extern RationalMatrix *rtnlmtx_times_intmat(const RationalMatrix *A,
                                            const IntegerMatrix *B);
extern RationalMatrix *intmat_times_rtnlmtx(const IntegerMatrix *a,
                                            const RationalMatrix *b);

extern int transform_fractional_coords_rtnl(const RationalMatrix *P,
                                            const Rational *ivec,
                                            Rational *ans);
extern void transform_fractional_coords_rtnl_inverse(const RationalMatrix *P,
                                                     const Rational *vec,
                                                     Rational *ans);
extern void rtnl_mtx_print(const RationalMatrix *m);
extern Rational rtnl_mtx_det(const RationalMatrix *m);
extern int rtnl_mtx_is_identity(const RationalMatrix *m);
extern int rtnl_mtx_is_perm(const RationalMatrix *m);

#ifdef __cplusplus
}
#endif

#endif	/* RATIONAL_H */
