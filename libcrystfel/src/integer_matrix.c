/*
 * integer_matrix.c
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "rational.h"
#include "integer_matrix.h"
#include "utils.h"


/**
 * SECTION:integer_matrix
 * @short_description: Integer matrices
 * @title: Integer matrices
 * @section_id:
 * @see_also:
 * @include: "integer_matrix.h"
 * @Image:
 *
 * An integer matrix library
 */


struct _integermatrix
{
	unsigned int rows;
	unsigned int cols;

	signed int *v;
};


/**
 * intmat_new:
 * @rows: Number of rows that the new matrix is to have
 * @cols: Number of columns that the new matrix is to have
 *
 * Allocates a new %IntegerMatrix with all elements set to zero.
 *
 * Returns: a new %IntegerMatrix, or NULL on error.
 **/
IntegerMatrix *intmat_new(unsigned int rows, unsigned int cols)
{
	IntegerMatrix *m;

	m = malloc(sizeof(IntegerMatrix));
	if ( m == NULL ) return NULL;

	m->v = calloc(rows*cols, sizeof(signed int));
	if ( m->v == NULL ) {
		free(m);
		return NULL;
	}

	m->rows = rows;
	m->cols = cols;

	return m;
}


/**
 * intmat_copy:
 * @m: An %IntegerMatrix
 *
 * Returns: a newly allocated copy of @m, or NULL on error/
 **/
IntegerMatrix *intmat_copy(const IntegerMatrix *m)
{
	IntegerMatrix *p;
	int i, j;

	p = intmat_new(m->rows, m->cols);
	if ( p == NULL ) return NULL;

	for ( i=0; i<m->rows; i++ ) {
	for ( j=0; j<m->rows; j++ ) {
		intmat_set(p, i, j, intmat_get(m, i, j));
	}
	}

	return p;
}


/**
 * intmat_free:
 * @m: An %IntegerMatrix
 *
 * Frees @m, unless @m is NULL in which case nothing is done.
 **/
void intmat_free(IntegerMatrix *m)
{
	if ( m == NULL ) return;
	free(m->v);
	free(m);
}


void intmat_size(const IntegerMatrix *m, unsigned int *rows, unsigned int *cols)
{
	if ( m == NULL ) {
		*rows = 0;
		*cols = 0;
		return;
	}

	*rows = m->rows;
	*cols = m->cols;
}


/**
 * intmat_set:
 * @m: An %IntegerMatrix
 * @i: row number to set
 * @j: column number to set
 * @v: value to set to
 *
 * Sets the @i,@j element of @m to @v.
 **/
void intmat_set(IntegerMatrix *m, unsigned int i, unsigned int j, signed int v)
{
	assert(i < m->rows);
	assert(j < m->cols);
	m->v[j + m->cols*i] = v;
}


/**
 * intmat_get:
 * @m: An %IntegerMatrix
 * @i: column number to set
 * @j: row number to set
 *
 * Gets the @i,@j element of @m.
 *
 * Returns: the @i,@j element of @m.
 **/
signed int intmat_get(const IntegerMatrix *m, unsigned int i, unsigned int j)
{
	assert(i < m->rows);
	assert(j < m->cols);
	return m->v[j + m->cols*i];
}


/* intmat_intvec_mult:
 * @m: An %IntegerMatrix
 * @vec: An array of floats
 * @ans: An array of floats in which to store the result
 *
 * Multiplies the matrix @m by the vector @vec.  The size of @vec must equal the
 * number of columns in @m, and the size of the result equals the number of rows
 * in @m.
 *
 * Returns: non-zero on error
 **/
int intmat_floatvec_mult(const IntegerMatrix *m, const float *vec, float *ans)
{
	unsigned int i;

	for ( i=0; i<m->rows; i++ ) {

		unsigned int j;

		ans[i] = 0.0;
		for ( j=0; j<m->cols; j++ ) {
			ans[i] += intmat_get(m, i, j) * vec[j];
		}

	}

	return 0;
}


/* intmat_rationalvec_mult:
 * @m: An %IntegerMatrix
 * @vec: An array of %Rational
 * @ans: An array of %Rational in which to store the result
 *
 * Multiplies the matrix @m by the vector @vec.  The size of @vec must equal the
 * number of columns in @m, and the size of the result equals the number of rows
 * in @m.
 *
 * Returns: non-zero on error
 **/
int intmat_rationalvec_mult(const IntegerMatrix *m, const Rational *vec,
                            Rational *ans)
{
	unsigned int i;


	for ( i=0; i<m->rows; i++ ) {

		unsigned int j;

		ans[i] = rtnl_zero();
		for ( j=0; j<m->cols; j++ ) {
			Rational t;
			t = rtnl_mul(vec[j], rtnl(intmat_get(m, i, j), 1));
			ans[i] = rtnl_add(ans[i], t);
		}

	}

	return 0;
}


/* intmat_solve_rational:
 * @m: An %IntegerMatrix
 * @vec: An array of %Rational
 * @ans: An array of %Rational in which to store the result
 *
 * Solves the matrix equation m*ans = vec, where @ans and @vec are
 * vectors of %Rational.
 *
 * This is just a convenience function which creates a %RationalMatrix out of
 * @m and then calls rtnl_mtx_solve().
 *
 * Returns: non-zero on error
 **/
int intmat_solve_rational(const IntegerMatrix *m, const Rational *vec,
                          Rational *ans)
{
	RationalMatrix *rm;
	int r;

	rm = rtnl_mtx_from_intmat(m);
	r = rtnl_mtx_solve(rm, vec, ans);
	rtnl_mtx_free(rm);
	return r;
}


/**
 * intmat_intvec_mult:
 * @m: An %IntegerMatrix
 * @vec: An array of signed integers
 *
 * Multiplies the matrix @m by the vector @vec.  The size of @vec must equal the
 * number of columns in @m, and the size of the result equals the number of rows
 * in @m.
 *
 * Returns: a newly allocated array of signed integers containing the answer,
 * or NULL on error.
 **/
signed int *intmat_intvec_mult(const IntegerMatrix *m, const signed int *vec)
{
	signed int *ans;
	unsigned int i;

	ans = malloc(m->rows * sizeof(signed int));
	if ( ans == NULL ) return NULL;

	for ( i=0; i<m->rows; i++ ) {

		unsigned int j;

		ans[i] = 0;
		for ( j=0; j<m->cols; j++ ) {
			ans[i] += intmat_get(m, i, j) * vec[j];
		}

	}

	return ans;
}


/**
 * intmat_intmat_mult:
 * @a: An %IntegerMatrix
 * @b: An %IntegerMatrix
 *
 * Multiplies the matrix @a by the matrix @b.
 *
 * Returns: a newly allocated %IntegerMatrix containing the answer, or NULL on
 * error.
 **/
IntegerMatrix *intmat_intmat_mult(const IntegerMatrix *a,
                                  const IntegerMatrix *b)
{
	unsigned int i, j;
	IntegerMatrix *ans;

	if ( a->cols != b->rows ) return NULL;

	ans = intmat_new(a->rows, a->cols);
	if ( ans == NULL ) return NULL;

	for ( i=0; i<ans->rows; i++ ) {
	for ( j=0; j<ans->cols; j++ ) {

		unsigned int k;
		signed int r = 0;

		for ( k=0; k<a->cols; k++ ) {  /* a->cols == b->rows */
			r += intmat_get(a, i, k) * intmat_get(b, k, j);
		}
		intmat_set(ans, i, j, r);

	}
	}

	return ans;
}


static IntegerMatrix *delete_row_and_column(const IntegerMatrix *m,
                                            unsigned int di, unsigned int dj)
{
	IntegerMatrix *n;
	unsigned int i, j;

	n = intmat_new(m->rows-1, m->cols-1);
	if ( n == NULL ) return NULL;

	for ( i=0; i<n->rows; i++ ) {
	for ( j=0; j<n->cols; j++ ) {

		signed int val;
		unsigned int gi, gj;

		gi = (i>=di) ? i+1 : i;
		gj = (j>=dj) ? j+1 : j;
		val = intmat_get(m, gi, gj);
		intmat_set(n, i, j, val);

	}
	}

	return n;
}


static signed int cofactor(const IntegerMatrix *m,
                           unsigned int i, unsigned int j)
{
	IntegerMatrix *n;
	signed int t, C;

	n = delete_row_and_column(m, i, j);
	if ( n == NULL ) {
		fprintf(stderr, "Failed to allocate matrix.\n");
		return 0;
	}

	/* -1 if odd, +1 if even */
	t = (i+j) & 0x1 ? -1 : +1;

	C = t * intmat_det(n);
	intmat_free(n);

	return C;
}


/**
 * intmat_det:
 * @m: An %IntegerMatrix
 *
 * Calculates the determinant of @m.  Inefficiently.
 *
 * Returns: the determinant of @m.
 **/
signed int intmat_det(const IntegerMatrix *m)
{
	unsigned int i, j;
	signed int det = 0;

	assert(m->rows == m->cols);  /* Otherwise determinant doesn't exist */

	if ( m->rows == 2 ) {
		return intmat_get(m, 0, 0)*intmat_get(m, 1, 1)
		     - intmat_get(m, 0, 1)*intmat_get(m, 1, 0);
	}

	i = 0;  /* Fixed */
	for ( j=0; j<m->cols; j++ ) {

		det += intmat_get(m, i, j) * cofactor(m, i, j);

	}

	return det;
}


static IntegerMatrix *intmat_cofactors(const IntegerMatrix *m)
{
	IntegerMatrix *n;
	signed int i, j;

	n = intmat_new(m->rows, m->cols);
	if ( n == NULL ) return NULL;

	for ( i=0; i<n->rows; i++ ) {
	for ( j=0; j<n->cols; j++ ) {

		intmat_set(n, i, j, cofactor(m, i, j));

	}
	}

	return n;
}


/**
 * intmat_inverse:
 * @m: An %IntegerMatrix
 *
 * Calculates the inverse of @m.  Inefficiently.
 *
 * Works only if the inverse of the matrix is also an integer matrix,
 * i.e. if the determinant of @m is +/- 1.
 *
 * Returns: the inverse of @m, or NULL on error.
 **/
IntegerMatrix *intmat_inverse(const IntegerMatrix *m)
{
	IntegerMatrix *adjugateT;
	IntegerMatrix *inverse;
	unsigned int i, j;
	signed int det;

	det = intmat_det(m);
	if ( (det != +1) && (det != -1) ) {
		fprintf(stderr,
		        "Inverse matrix not an integer matrix (det = %i).\n",
		        det);
		return NULL;
	}

	adjugateT = intmat_cofactors(m);
	if ( adjugateT == NULL ) return NULL;

	inverse = intmat_new(m->cols, m->rows);  /* The other way round */
	if ( inverse == NULL ) return NULL;

	for ( i=0; i<inverse->rows; i++ ) {
	for ( j=0; j<inverse->cols; j++ ) {

		signed int v;

		v = intmat_get(adjugateT, j, i);

		/* 1/-1 = -1 and 1/+1 = +1, and these are the only two cases */
		intmat_set(inverse, i, j, v*det);

	}
	}

	intmat_free(adjugateT);

	return inverse;
}


/**
 * intmat_print
 * @m: An %IntegerMatrix
 *
 * Prints @m to stderr.
 *
 */
void intmat_print(const IntegerMatrix *m)
{
	unsigned int i, j;

	for ( i=0; i<m->rows; i++ ) {

		fprintf(stderr, "[ ");
		for ( j=0; j<m->cols; j++ ) {
			fprintf(stderr, "%4i ", intmat_get(m, i, j));
		}
		fprintf(stderr, "]\n");
	}
}


/**
 * intmat_is_identity
 * @m: An %IntegerMatrix
 *
 * Returns: true if @m is an identity matrix.
 *
 */
int intmat_is_identity(const IntegerMatrix *m)
{
	int i, j;

	if ( m->rows != m->cols ) return 0;

	for ( i=0; i<m->rows; i++ ) {
	for ( j=0; j<m->cols; j++ ) {

		signed int v;

		v = intmat_get(m, i, j);

		if ( i == j ) {
			if ( v != 1 ) return 0;
		} else {
			if ( v != 0 ) return 0;
		}

	}
	}

	return 1;
}


/**
 * intmat_is_inversion
 * @m: An %IntegerMatrix
 *
 * Returns: true if @m = -I, where I is an identity matrix.
 *
 */
int intmat_is_inversion(const IntegerMatrix *m)
{
	int i, j;

	if ( m->rows != m->cols ) return 0;

	for ( i=0; i<m->rows; i++ ) {
	for ( j=0; j<m->cols; j++ ) {

		signed int v;

		v = intmat_get(m, i, j);

		if ( i == j ) {
			if ( v != -1 ) return 0;
		} else {
			if ( v != 0 ) return 0;
		}

	}
	}

	return 1;
}


/**
 * intmat_equals
 * @a: An %IntegerMatrix
 * @b: An %IntegerMatrix
 *
 * Returns: true if @a = @b.
 *
 */
int intmat_equals(const IntegerMatrix *a, const IntegerMatrix *b)
{
	int i, j;

	if ( a->rows != b->rows ) return 0;
	if ( a->cols != b->cols ) return 0;

	for ( i=0; i<a->rows; i++ ) {
	for ( j=0; j<b->cols; j++ ) {

		signed int v;

		v = intmat_get(a, i, j);

		if ( v != intmat_get(b, i, j) ) return 0;

	}
	}

	return 1;
}


/**
 * intmat_identity
 * @size: The size of the (square) matrix
 *
 * Returns: an identity %IntegerMatrix with side length @size, or NULL on error.
 *
 */
IntegerMatrix *intmat_identity(int size)
{
	IntegerMatrix *m;
	int i, j;

	m = intmat_new(size, size);
	if ( m == NULL ) return NULL;

	for ( i=0; i<size; i++ ) {
	for ( j=0; j<size; j++ ) {

		if ( i == j ) {
			intmat_set(m, i, j, 1);
		} else {
			intmat_set(m, i, j, 0);
		}

	}
	}

	return m;
}


/**
 * intmat_set_all_3x3
 * @m: An %IntegerMatrix
 *
 * Returns: an identity %IntegerMatrix with side length @size, or NULL on error.
 *
 */
void intmat_set_all_3x3(IntegerMatrix *m,
                        signed int m11, signed int m12, signed int m13,
                        signed int m21, signed int m22, signed int m23,
                        signed int m31, signed int m32, signed int m33)
{
	intmat_set(m, 0, 0, m11);
	intmat_set(m, 0, 1, m12);
	intmat_set(m, 0, 2, m13);

	intmat_set(m, 1, 0, m21);
	intmat_set(m, 1, 1, m22);
	intmat_set(m, 1, 2, m23);

	intmat_set(m, 2, 0, m31);
	intmat_set(m, 2, 1, m32);
	intmat_set(m, 2, 2, m33);
}


IntegerMatrix *intmat_create_3x3(signed int m11, signed int m12, signed int m13,
                                 signed int m21, signed int m22, signed int m23,
                                 signed int m31, signed int m32, signed int m33)
{
	IntegerMatrix *t = intmat_new(3, 3);
	intmat_set_all_3x3(t, m11, m12, m13, m21, m22, m23, m31, m32, m33);
	return t;
}
