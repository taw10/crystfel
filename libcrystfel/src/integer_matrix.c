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


/** \file integer_matrix.h */


struct _integermatrix
{
	unsigned int rows;
	unsigned int cols;

	signed int *v;
};


/**
 * \param rows Number of rows that the new matrix is to have
 * \param cols Number of columns that the new matrix is to have
 *
 * Allocates a new \ref IntegerMatrix with all elements set to zero.
 *
 * \returns A new \ref IntegerMatrix, or NULL on error.
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
 * \param m An \ref IntegerMatrix
 *
 * \returns A newly allocated copy of \p m, or NULL on error
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
 * \param m An \ref IntegerMatrix
 *
 * Frees \p m, unless \p m is NULL in which case nothing is done.
 **/
void intmat_free(IntegerMatrix *m)
{
	if ( m == NULL ) return;
	free(m->v);
	free(m);
}


/**
 * \param m An \ref IntegerMatrix
 * \param rows Location to store number of rows
 * \param cols Location to store number of columns
 *
 * Sets \p rows and \p cols to the size of \p m.
 */
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
 * \param m An \ref IntegerMatrix
 * \param i row number to set
 * \param j column number to set
 * \param v value to set to
 *
 * Sets the \p i,\p j element of \p m to \p v.
 **/
void intmat_set(IntegerMatrix *m, unsigned int i, unsigned int j, signed int v)
{
	assert(i < m->rows);
	assert(j < m->cols);
	m->v[j + m->cols*i] = v;
}


/**
 * \param m An \ref IntegerMatrix
 * \param i column number to set
 * \param j row number to set
 *
 * Gets the \p i,\p j element of \p m.
 *
 * \returns The \p i,\p j element of \p m.
 **/
signed int intmat_get(const IntegerMatrix *m, unsigned int i, unsigned int j)
{
	assert(i < m->rows);
	assert(j < m->cols);
	return m->v[j + m->cols*i];
}


/**
 * \param P An \ref IntegerMatrix
 * \param hkl An array of signed integers
 *
 * Apply transformation matrix P to a set of reciprocal space Miller indices.
 *
 * In other words:
 * Multiplies the matrix \p P by the row vector \p hkl.  The size of \p vec must equal
 * the number of columns in \p P, and the size of the result equals the number of
 * rows in \p P.
 *
 * The multiplication looks like this:
 *    (a1, a2, a3) = (hkl1, hkl2, hkl3) P
 * Therefore matching the notation in ITA chapter 5.1.
 *
 * \returns A newly allocated array of signed integers containing the answer,
 * or NULL on error.
 **/
signed int *transform_indices(const IntegerMatrix *P, const signed int *hkl)
{
	signed int *ans;
	unsigned int j;

	ans = malloc(P->rows * sizeof(signed int));
	if ( ans == NULL ) return NULL;

	for ( j=0; j<P->cols; j++ ) {

		unsigned int i;
		ans[j] = 0;
		for ( i=0; i<P->rows; i++ ) {
			ans[j] += intmat_get(P, i, j) * hkl[i];
		}

	}

	return ans;
}


/**
 * \param a An \ref IntegerMatrix
 * \param b An \ref IntegerMatrix
 *
 * Multiplies the matrix \p a by the matrix \p b.
 *
 * \returns A newly allocated \ref IntegerMatrix containing the answer, or NULL on
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
 * \param m An \ref IntegerMatrix
 *
 * Calculates the determinant of \p m.  Inefficiently.
 *
 * \returns The determinant of \p m.
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
 * \param m An \ref IntegerMatrix
 *
 * Calculates the inverse of \p m.  Inefficiently.
 *
 * Works only if the inverse of the matrix is also an integer matrix,
 * i.e. if the determinant of \p m is +/- 1.
 *
 * \returns The inverse of \p m, or NULL on error.
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
 * \param m An \ref IntegerMatrix
 *
 * Prints \param m to stderr.
 *
 */
void intmat_print(const IntegerMatrix *m)
{
	unsigned int i, j;

	if ( m == NULL ) {
		fprintf(stderr, "(NULL matrix)\n");
		return;
	}

	for ( i=0; i<m->rows; i++ ) {

		fprintf(stderr, "[ ");
		for ( j=0; j<m->cols; j++ ) {
			fprintf(stderr, "%4i ", intmat_get(m, i, j));
		}
		fprintf(stderr, "]\n");
	}
}


/**
 * \param m An \ref IntegerMatrix
 *
 * \returns True if \p m is an identity matrix.
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
 * \param m An \ref IntegerMatrix
 *
 * \returns True if \p m = -I, where I is an identity matrix.
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
 * \param a An \ref IntegerMatrix
 * \param b An \ref IntegerMatrix
 *
 * \returns True if \p a = \p b.
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
 * \param size The size of the (square) matrix
 *
 * \returns An identity \ref IntegerMatrix with side length \p size, or NULL on error.
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
 * \param m11 Matrix element
 * \param m12 Matrix element
 * \param m13 Matrix element
 * \param m21 Matrix element
 * \param m22 Matrix element
 * \param m23 Matrix element
 * \param m31 Matrix element
 * \param m32 Matrix element
 * \param m33 Matrix element
 *
 * \returns A newly allocated 3x3 \ref IntegerMatrix with the given values.
 */
IntegerMatrix *intmat_create_3x3(signed int m11, signed int m12, signed int m13,
                                 signed int m21, signed int m22, signed int m23,
                                 signed int m31, signed int m32, signed int m33)
{
	IntegerMatrix *m = intmat_new(3, 3);
	if ( m == NULL ) return NULL;

	intmat_set(m, 0, 0, m11);
	intmat_set(m, 0, 1, m12);
	intmat_set(m, 0, 2, m13);

	intmat_set(m, 1, 0, m21);
	intmat_set(m, 1, 1, m22);
	intmat_set(m, 1, 2, m23);

	intmat_set(m, 2, 0, m31);
	intmat_set(m, 2, 1, m32);
	intmat_set(m, 2, 2, m33);
	return m;
}
