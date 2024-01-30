/*
 * rational.c
 *
 * A small rational number library
 *
 * Copyright © 2019-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019-2020 Thomas White <taw@physics.org>
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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <locale.h>

#include "rational.h"
#include "integer_matrix.h"
#include "utils.h"

/** \file rational.h */

/* Euclidean algorithm for finding greatest common divisor */
static signed int gcd(signed int a, signed int b)
{
	while ( b != 0 ) {
		signed int t = b;
		b = a % b;
		a = t;
	}
	return a;
}


static void squish(Rational *rt)
{
	signed int g;

	if ( rt->num == 0 ) {
		rt->den = 1;
		return;
	}

	g = gcd(rt->num, rt->den);
	assert(g != 0);
	rt->num /= g;
	rt->den /= g;

	if ( rt->den < 0 ) {
		rt->num = -rt->num;
		rt->den = -rt->den;
	}
}


Rational rtnl_zero()
{
	Rational r;
	r.num = 0;
	r.den = 1;
	return r;
}


Rational rtnl(signed long long int num, signed long long int den)
{
	Rational r;
	r.num = num;
	r.den = den;
	squish(&r);
	return r;
}


double rtnl_as_double(Rational r)
{
	return (double)r.num/r.den;
}


static void overflow(long long int c, long long int a, long long int b)
{
	setlocale(LC_ALL, "");
	ERROR("Overflow detected in rational number library.\n");
	ERROR("%'lli < %'lli * %'lli\n", c, a, b);
	abort();
}


static void check_overflow(long long int c, long long int a, long long int b)
{
	if ( (a==0) || (b==0) ) {
		if ( c != 0 ) overflow(c,a,b);
	} else if ( (llabs(c) < llabs(a)) || (llabs(c) < llabs(b)) ) {
		overflow(c,a,b);
	}
}


Rational rtnl_mul(Rational a, Rational b)
{
	Rational r;
	r.num = a.num * b.num;
	r.den = a.den * b.den;
	check_overflow(r.num, a.num, b.num);
	check_overflow(r.den, a.den, b.den);
	squish(&r);
	return r;
}


Rational rtnl_div(Rational a, Rational b)
{
	signed int t = b.num;
	b.num = b.den;
	b.den = t;
	return rtnl_mul(a, b);
}


Rational rtnl_add(Rational a, Rational b)
{
	Rational r;
	long long int common_den;
	long long int a_num, b_num;

	a_num = a.num * b.den;
	b_num = b.num * a.den;
	check_overflow(a_num, a.num, b.den);
	check_overflow(b_num, b.num, a.den);

	common_den = a.den * b.den;
	check_overflow(common_den, a.den, b.den);

	r.num = a_num + b_num;
	r.den = common_den;
	squish(&r);
	return r;
}


Rational rtnl_sub(Rational a, Rational b)
{
	b.num = -b.num;
	return rtnl_add(a, b);
}


/* -1, 0 +1 respectively for a<b, a==b, a>b */
signed int rtnl_cmp(Rational a, Rational b)
{
	long long int a_num, b_num;

	a_num = a.num * b.den;
	b_num = b.num * a.den;

	if ( a_num > b_num ) return +1;
	if ( a_num < b_num ) return -1;
	return 0;
}


Rational rtnl_abs(Rational a)
{
	Rational r = a;
	squish(&r);
	if ( r.num < 0 ) r.num = -r.num;
	return r;
}


/**
 * \param rt A \ref Rational
 *
 * Formats \p rt as a string
 *
 * \returns A string which should be freed by the caller
 */
char *rtnl_format(Rational rt)
{
	char *v = cfmalloc(32);
	if ( v == NULL ) return NULL;
	if ( rt.den == 1 ) {
		snprintf(v, 31, "%lli", rt.num);
	} else {
		snprintf(v, 31, "%lli/%lli", rt.num, rt.den);
	}
	return v;
}


Rational *rtnl_list(signed int num_min, signed int num_max,
                    signed int den_min, signed int den_max,
                    int *pn)
{
	signed int num, den;
	Rational *list;
	int n = 0;

	list = cfmalloc((1+num_max-num_min)*(1+den_max-den_min)*sizeof(Rational));
	if ( list == NULL ) return NULL;

	for ( num=num_min; num<=num_max; num++ ) {
		for ( den=den_min; den<=den_max; den++ ) {

			Rational r = rtnl(num, den);

			/* Denominator zero? */
			if ( den == 0 ) continue;

			/* Same as last entry? */
			if ( (n>0) && (rtnl_cmp(list[n-1], r)==0) ) continue;

			/* Can be reduced? */
			if ( abs(gcd(num, den)) != 1 ) continue;

			list[n++] = r;
		}
	}
	*pn = n;
	return list;
}


struct _rationalmatrix
{
	unsigned int rows;
	unsigned int cols;
	Rational *v;
};


/**
 * \param rows Number of rows that the new matrix is to have
 * \param cols Number of columns that the new matrix is to have
 *
 * Allocates a new \ref RationalMatrix with all elements set to zero.
 *
 * \returns A new \ref RationalMatrix, or NULL on error.
 **/
RationalMatrix *rtnl_mtx_new(unsigned int rows, unsigned int cols)
{
	RationalMatrix *m;
	int i;

	m = cfmalloc(sizeof(RationalMatrix));
	if ( m == NULL ) return NULL;

	m->v = cfcalloc(rows*cols, sizeof(Rational));
	if ( m->v == NULL ) {
		cffree(m);
		return NULL;
	}

	m->rows = rows;
	m->cols = cols;

	for ( i=0; i<m->rows*m->cols; i++ ) {
		m->v[i] = rtnl_zero();
	}

	return m;
}


RationalMatrix *rtnl_mtx_identity(int rows)
{
	int i;
	RationalMatrix *m = rtnl_mtx_new(rows, rows);
	for ( i=0; i<rows; i++ ) {
		rtnl_mtx_set(m, i, i, rtnl(1,1));
	}
	return m;
}


RationalMatrix *rtnl_mtx_copy(const RationalMatrix *m)
{
	RationalMatrix *n;
	int i;

	n = rtnl_mtx_new(m->rows, m->cols);
	if ( n == NULL ) return NULL;

	for ( i=0; i<m->rows*m->cols; i++ ) {
		n->v[i] = m->v[i];
	}

	return n;
}


Rational rtnl_mtx_get(const RationalMatrix *m, int i, int j)
{
	assert(m != NULL);
	return m->v[j+m->cols*i];
}


void rtnl_mtx_set(const RationalMatrix *m, int i, int j, Rational v)
{
	assert(m != NULL);
	m->v[j+m->cols*i] = v;
}


RationalMatrix *rtnl_mtx_from_intmat(const IntegerMatrix *m)
{
	RationalMatrix *n;
	unsigned int rows, cols;
	int i, j;

	intmat_size(m, &rows, &cols);
	n = rtnl_mtx_new(rows, cols);
	if ( n == NULL ) return NULL;

	for ( i=0; i<rows; i++ ) {
		for ( j=0; j<cols; j++ ) {
			rtnl_mtx_set(n, i, j, rtnl(intmat_get(m, i, j), 1));
		}
	}

	return n;
}


IntegerMatrix *intmat_from_rtnl_mtx(const RationalMatrix *m)
{
	IntegerMatrix *n;
	int i, j;

	n = intmat_new(m->rows, m->cols);
	if ( n == NULL ) return NULL;

	for ( i=0; i<m->rows; i++ ) {
		for ( j=0; j<m->cols; j++ ) {
			Rational v = rtnl_mtx_get(m, i, j);
			squish(&v);
			if ( v.den != 1 ) {
				ERROR("Rational matrix can't be converted to integers\n");
				intmat_free(n);
				return NULL;
			}
			intmat_set(n, i, j, v.num);
		}
	}

	return n;
}


void rtnl_mtx_free(RationalMatrix *mtx)
{
	if ( mtx == NULL ) return;
	cffree(mtx->v);
	cffree(mtx);
}


/* rtnl_mtx_solve:
 * @P: A %RationalMatrix
 * @vec: An array of %Rational
 * @ans: An array of %Rational in which to store the result
 *
 * Solves the matrix equation m*ans = vec, where @ans and @vec are
 * vectors of %Rational.
 *
 * In this version, @m must be square.
 *
 * The number of columns in @m must equal the length of @ans, and the number
 * of rows in @m must equal the length of @vec, but note that there is no way
 * for this function to check that this is the case.
 *
 * Given that P is transformation of the unit cell axes (see ITA chapter 5.1),
 * this function calculates the fractional coordinates of a point in the
 * transformed axes, given its fractional coordinates in the original axes.
 *
 * Returns: non-zero on error
 **/
int transform_fractional_coords_rtnl(const RationalMatrix *P,
                                     const Rational *ivec, Rational *ans)
{
	RationalMatrix *cm;
	Rational *vec;
	int h, k;
	int i;

	if ( P->rows != P->cols ) return 1;

	/* Copy the matrix and vector because the calculation will
	 * be done in-place */
	cm = rtnl_mtx_copy(P);
	if ( cm == NULL ) return 1;

	vec = cfmalloc(cm->rows*sizeof(Rational));
	if ( vec == NULL ) return 1;
	for ( h=0; h<cm->rows; h++ ) vec[h] = ivec[h];

	/* Gaussian elimination with partial pivoting */
	h = 0;
	k = 0;
	while ( h<=cm->rows && k<=cm->cols ) {

		int prow = 0;
		Rational pval = rtnl_zero();
		Rational t;

		/* Find the row with the largest value in column k */
		for ( i=h; i<cm->rows; i++ ) {
			if ( rtnl_cmp(rtnl_abs(rtnl_mtx_get(cm, i, k)), pval) > 0 ) {
				pval = rtnl_abs(rtnl_mtx_get(cm, i, k));
				prow = i;
			}
		}

		if ( rtnl_cmp(pval, rtnl_zero()) == 0 ) {
			k++;
			continue;
		}

		/* Swap 'prow' with row h */
		for ( i=0; i<cm->cols; i++ ) {
			t = rtnl_mtx_get(cm, h, i);
			rtnl_mtx_set(cm, h, i, rtnl_mtx_get(cm, prow, i));
			rtnl_mtx_set(cm, prow, i, t);
		}
		t = vec[prow];
		vec[prow] = vec[h];
		vec[h] = t;

		/* Divide and subtract rows below */
		for ( i=h+1; i<cm->rows; i++ ) {

			int j;
			Rational dval;

			dval = rtnl_div(rtnl_mtx_get(cm, i, k),
			                rtnl_mtx_get(cm, h, k));

			for ( j=0; j<cm->cols; j++ ) {
				Rational tcol = rtnl_mtx_get(cm, i, j);
				Rational p = rtnl_mul(dval, rtnl_mtx_get(cm, h, j));
				tcol = rtnl_sub(tcol, p);
				rtnl_mtx_set(cm, i, j, tcol);
			}

			/* Divide the right hand side as well */
			Rational trhs = vec[i];
			Rational p = rtnl_mul(dval, vec[h]);
			vec[i] = rtnl_sub(trhs, p);
		}

		h++;
		k++;


	}

	/* Back-substitution */
	for ( i=cm->rows-1; i>=0; i-- ) {
		int j;
		Rational sum = rtnl_zero();
		for ( j=i+1; j<cm->cols; j++ ) {
			Rational av;
			av = rtnl_mul(rtnl_mtx_get(cm, i, j), ans[j]);
			sum = rtnl_add(sum, av);
		}
		sum = rtnl_sub(vec[i], sum);
		ans[i] = rtnl_div(sum, rtnl_mtx_get(cm, i, i));
	}

	cffree(vec);
	rtnl_mtx_free(cm);

	return 0;
}


/**
 * \param m A \ref RationalMatrix
 *
 * Prints \p m to stderr.
 *
 */
void rtnl_mtx_print(const RationalMatrix *m)
{
	unsigned int i, j;

	if ( m == NULL ) {
		fprintf(stderr, "(NULL matrix)\n");
		return;
	}

	for ( i=0; i<m->rows; i++ ) {

		fprintf(stderr, "[ ");
		for ( j=0; j<m->cols; j++ ) {
			char *v = rtnl_format(rtnl_mtx_get(m, i, j));
			fprintf(stderr, "%4s ", v);
			cffree(v);
		}
		fprintf(stderr, "]\n");
	}
}


RationalMatrix *rtnlmtx_times_intmat(const RationalMatrix *A,
                                     const IntegerMatrix *B)
{
	int i, j;
	RationalMatrix *ans;
	unsigned int B_rows, B_cols;

	intmat_size(B, &B_rows, &B_cols);
	assert(A->cols == B_rows);

	ans = rtnl_mtx_new(A->rows, B_cols);
	if ( ans == NULL ) return NULL;

	for ( i=0; i<ans->rows; i++ ) {
		for ( j=0; j<ans->cols; j++ ) {
			int k;
			Rational sum = rtnl_zero();
			for ( k=0; k<A->rows; k++ ) {
				Rational add;
				add = rtnl_mul(rtnl_mtx_get(A, i, k),
				               rtnl(intmat_get(B, k, j), 1));
				sum = rtnl_add(sum, add);
			}
			rtnl_mtx_set(ans, i, j, sum);
		}
	}

	return ans;
}


RationalMatrix *rtnlmtx_times_rtnlmtx(const RationalMatrix *A,
                                      const RationalMatrix *B)
{
	int i, j;
	RationalMatrix *ans;

	assert(A->cols == B->rows);

	ans = rtnl_mtx_new(A->rows, B->cols);
	if ( ans == NULL ) return NULL;

	for ( i=0; i<ans->rows; i++ ) {
		for ( j=0; j<ans->cols; j++ ) {
			int k;
			Rational sum = rtnl_zero();
			for ( k=0; k<A->rows; k++ ) {
				Rational add;
				add = rtnl_mul(rtnl_mtx_get(A, i, k),
				               rtnl_mtx_get(B, k, j));
				sum = rtnl_add(sum, add);
			}
			rtnl_mtx_set(ans, i, j, sum);
		}
	}

	return ans;
}


RationalMatrix *intmat_times_rtnlmtx(const IntegerMatrix *A, const RationalMatrix *B)
{
	int i, j;
	RationalMatrix *ans;
	unsigned int A_rows, A_cols;

	intmat_size(A, &A_rows, &A_cols);

	ans = rtnl_mtx_new(A_rows, B->cols);
	if ( ans == NULL ) return NULL;

	for ( i=0; i<ans->rows; i++ ) {
		for ( j=0; j<ans->cols; j++ ) {
			int k;
			Rational sum = rtnl_zero();
			for ( k=0; k<A_rows; k++ ) {
				Rational add;
				add = rtnl_mul(rtnl(intmat_get(A, i, k), 1),
				               rtnl_mtx_get(B, k, j));
				sum = rtnl_add(sum, add);
			}
			rtnl_mtx_set(ans, i, j, sum);
		}
	}

	return ans;
}


/* Given a "P-matrix" (see ITA chapter 5.1), calculate the fractional
 * coordinates of point "vec" in the original axes, given its fractional
 * coordinates in the transformed axes.
 * To transform in the opposite direction, we would multiply by Q = P^-1.
 * Therefore, this direction is the "easy way" and needs just a matrix
 * multiplication. */
void transform_fractional_coords_rtnl_inverse(const RationalMatrix *P,
                                              const Rational *vec,
                                              Rational *ans)
{
	int i, j;

	for ( i=0; i<P->rows; i++ ) {
		ans[i] = rtnl_zero();
		for ( j=0; j<P->cols; j++ ) {
			Rational add;
			add = rtnl_mul(rtnl_mtx_get(P, i, j), vec[j]);
			ans[i] = rtnl_add(ans[i], add);
		}
	}
}


static RationalMatrix *rational_delete_row_and_column(const RationalMatrix *m,
                                                      unsigned int di,
                                                      unsigned int dj)
{
	RationalMatrix *n;
	unsigned int i, j;

	n = rtnl_mtx_new(m->rows-1, m->cols-1);
	if ( n == NULL ) return NULL;

	for ( i=0; i<n->rows; i++ ) {
	for ( j=0; j<n->cols; j++ ) {

		Rational val;
		unsigned int gi, gj;

		gi = (i>=di) ? i+1 : i;
		gj = (j>=dj) ? j+1 : j;
		val = rtnl_mtx_get(m, gi, gj);
		rtnl_mtx_set(n, i, j, val);

	}
	}

	return n;
}


static Rational rational_cofactor(const RationalMatrix *m,
                                  unsigned int i, unsigned int j)
{
	RationalMatrix *n;
	Rational t, C;

	n = rational_delete_row_and_column(m, i, j);
	if ( n == NULL ) {
		fprintf(stderr, "Failed to allocate matrix.\n");
		return rtnl_zero();
	}

	/* -1 if odd, +1 if even */
	t = (i+j) & 0x1 ? rtnl(-1, 1) : rtnl(1, 1);

	C = rtnl_mul(t, rtnl_mtx_det(n));
	rtnl_mtx_free(n);

	return C;
}



Rational rtnl_mtx_det(const RationalMatrix *m)
{
	unsigned int i, j;
	Rational det;

	assert(m->rows == m->cols);  /* Otherwise determinant doesn't exist */

	if ( m->rows == 2 ) {
		Rational a, b;
		a = rtnl_mul(rtnl_mtx_get(m, 0, 0), rtnl_mtx_get(m, 1, 1));
		b = rtnl_mul(rtnl_mtx_get(m, 0, 1), rtnl_mtx_get(m, 1, 0));
		return rtnl_sub(a, b);
	}

	i = 0;  /* Fixed */
	det = rtnl_zero();
	for ( j=0; j<m->cols; j++ ) {
		Rational a;
		a = rtnl_mul(rtnl_mtx_get(m, i, j),
		             rational_cofactor(m, i, j));
		det = rtnl_add(det, a);
	}

	return det;
}


int rtnl_mtx_is_identity(const RationalMatrix *m)
{
	int i, j;

	if ( m->rows != m->cols ) return 0;

	for ( i=0; i<m->rows; i++ ) {
	for ( j=0; j<m->cols; j++ ) {

		Rational v;

		v = rtnl_mtx_get(m, i, j);

		if ( i == j ) {
			if ( rtnl_cmp(v, rtnl(1,1)) != 0 ) return 0;
		} else {
			if ( rtnl_cmp(v, rtnl_zero()) != 0 ) return 0;
		}

	}
	}

	return 1;
}


int rtnl_mtx_is_perm(const RationalMatrix *m)
{
	Rational det;
	int i, j;

	/* Must be square */
	if ( m->rows != m->cols ) return 0;

	/* Determinant must be +1 or -1 */
	det = rtnl_mtx_det(m);
	if ( (rtnl_cmp(det, rtnl(1,1)) != 0)
	  && (rtnl_cmp(det, rtnl(-1,1)) != 0) ) return 0;

	/* All components must be +1, -1 or 0 */
	for ( i=0; i<m->rows; i++ ) {
	for ( j=0; j<m->cols; j++ ) {

		Rational v = rtnl_mtx_get(m, i, j);

		if ( (rtnl_cmp(v, rtnl(1,1)) != 0)
		  && (rtnl_cmp(v, rtnl(-1,1)) != 0)
		  && (rtnl_cmp(v, rtnl_zero()) != 0) ) return 0;

	}
	}

	return 1;
}
