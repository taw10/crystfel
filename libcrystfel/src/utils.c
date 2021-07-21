/*
 * utils.c
 *
 * Utility stuff
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2020 Thomas White <taw@physics.org>
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

#include <libgen.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "utils.h"
#include "image.h"

/** \file utils.h */

/**
 * \param M A matrix
 * \param v A vector
 *
 * Displays a matrix equation of the form \p M.a = \p v.
 **/
void show_matrix_eqn(gsl_matrix *M, gsl_vector *v)
{
	int i, j;

	if ( M->size1 != v->size ) {
		ERROR("Matrix and vector sizes don't agree.\n");
		return;
	}

	for ( i=0; i<M->size1; i++ ) {
		STATUS("[ ");
		for ( j=0; j<M->size2; j++ ) {
			STATUS("%+9.3e ", gsl_matrix_get(M, i, j));
		}
		if ( i < M->size2 ) {
			STATUS("][ a%2i ] = [ %+9.3e ]\n", i,
			       gsl_vector_get(v, i));
		} else {
			STATUS("]        = [ +%9.3e ]\n", gsl_vector_get(v, i));
		}
	}
}


/**
 * \param M A matrix
 *
 * Displays a matrix.
 **/
void show_matrix(gsl_matrix *M)
{
	int i, j;

	for ( i=0; i<M->size1; i++ ) {
		STATUS("[ ");
		for ( j=0; j<M->size2; j++ ) {
			STATUS("%+9.3e ", gsl_matrix_get(M, i, j));
		}
		STATUS("]\n");
	}
}


static int check_eigen(gsl_vector *e_val, int verbose)
{
	int i;
	double vmax, vmin;
	const int n = e_val->size;
	const double max_condition = 1e6;
	int n_filt = 0;

	if ( verbose ) STATUS("Eigenvalues:\n");
	vmin = +INFINITY;
	vmax = 0.0;
	for ( i=0; i<n; i++ ) {
		double val = gsl_vector_get(e_val, i);
		if ( verbose ) STATUS("%i: %e\n", i, val);
		if ( val > vmax ) vmax = val;
		if ( val < vmin ) vmin = val;
	}

	for ( i=0; i<n; i++ ) {
		double val = gsl_vector_get(e_val, i);
		if ( val < vmax/max_condition ) {
			gsl_vector_set(e_val, i, 0.0);
			n_filt++;
		}
	}

	vmin = +INFINITY;
	vmax = 0.0;
	for ( i=0; i<n; i++ ) {
		double val = gsl_vector_get(e_val, i);
		if ( val == 0.0 ) continue;
		if ( val > vmax ) vmax = val;
		if ( val < vmin ) vmin = val;
	}
	if ( verbose ) {
		STATUS("Condition number: %e / %e = %5.2f\n",
		       vmax, vmin, vmax/vmin);
		STATUS("%i out of %i eigenvalues filtered.\n", n_filt, n);
	}

	return n_filt;
}


/**
 * \param v a gsl_vector
 * \param M a gsl_matrix
 * \param pn_filt pointer to store the number of filtered eigenvalues
 * \param verbose flag for verbosity on the terminal
 *
 * Solves the matrix equation M.x = v, returning x.
 * Performs rescaling and eigenvalue filtering.
 **/
gsl_vector *solve_svd(gsl_vector *v, gsl_matrix *M, int *pn_filt, int verbose)
{
	gsl_matrix *s_vec;
	gsl_vector *s_val;
	int err, n;
	gsl_vector *shifts;
	gsl_vector *SB;
	gsl_vector *SinvX;
	gsl_matrix *S;  /* rescaling matrix due to Bricogne */
	gsl_matrix *AS;
	gsl_matrix *SAS;
	int i;
	int n_filt;
	gsl_matrix *SAS_copy;

	n = v->size;
	if ( v->size != M->size1 ) return NULL;
	if ( v->size != M->size2 ) return NULL;

	/* Calculate the rescaling matrix S */
	S = gsl_matrix_calloc(n, n);
	for ( i=0; i<n; i++ ) {
		double sii = pow(gsl_matrix_get(M, i, i), -0.5);
		gsl_matrix_set(S, i, i, sii);
	}

	/* Calculate the matrix SAS, which we will be (not) inverting */
	AS = gsl_matrix_calloc(n, n);
	SAS = gsl_matrix_calloc(n, n);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, M, S, 0.0, AS);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, S, AS, 0.0, SAS);
	gsl_matrix_free(AS);

	/* Calculate the vector SB, which is the RHS of the equation */
	SB = gsl_vector_calloc(n);
	gsl_blas_dgemv(CblasNoTrans, 1.0, S, v, 0.0, SB);

	if ( verbose ) {
		STATUS("The equation after rescaling:\n");
		show_matrix_eqn(SAS, SB);
	}

	SAS_copy = gsl_matrix_alloc(SAS->size1, SAS->size2);
	gsl_matrix_memcpy(SAS_copy, SAS);

	for ( i=0; i<n; i++ ) {
		int j;
		if ( isnan(gsl_vector_get(SB, i)) ) gsl_vector_set(SB, i, 0.0);
		for ( j=0; j<n; j++ ) {
			if ( isnan(gsl_matrix_get(SAS, i, j)) ) {
					gsl_matrix_set(SAS, i, j, 0.0);
			}
		}
	}

	/* Do the SVD */
	s_val = gsl_vector_calloc(n);
	s_vec = gsl_matrix_calloc(n, n);
	err = gsl_linalg_SV_decomp_jacobi(SAS, s_vec, s_val);
	if ( err ) {
		if ( verbose ) ERROR("SVD failed: %s\n", gsl_strerror(err));
		gsl_matrix_free(s_vec);
		gsl_vector_free(s_val);
		gsl_matrix_free(SAS);
		gsl_matrix_free(S);
		return NULL;
	}
	/* "SAS" is now "U" */

	/* Filter the eigenvalues */
	n_filt = check_eigen(s_val, verbose);
	if ( pn_filt != NULL ) *pn_filt = n_filt;

	gsl_matrix_free(SAS_copy);

	/* Solve the equation SAS.SinvX = SB */
	SinvX = gsl_vector_calloc(n);
	err = gsl_linalg_SV_solve(SAS, s_vec, s_val, SB, SinvX);
	gsl_vector_free(SB);
	gsl_matrix_free(SAS);
	gsl_matrix_free(s_vec);
	gsl_vector_free(s_val);

	if ( err ) {
		ERROR("Matrix solution failed: %s\n", gsl_strerror(err));
		gsl_matrix_free(S);
		gsl_vector_free(SinvX);
		return NULL;
	}

	/* Calculate S.SinvX to get X, the shifts */
	shifts = gsl_vector_calloc(n);
	gsl_blas_dgemv(CblasNoTrans, 1.0, S, SinvX, 0.0, shifts);

	gsl_matrix_free(S);
	gsl_vector_free(SinvX);

	return shifts;
}


/* ------------------------------ Message logging ---------------------------- */

/* Lock to keep lines serialised on the terminal */
pthread_mutex_t stderr_lock = PTHREAD_MUTEX_INITIALIZER;


static void log_to_stderr(enum log_msg_type type, const char *msg,
                          void *vp)
{
	int error_print_val = get_status_label();
	pthread_mutex_lock(&stderr_lock);
	if ( error_print_val >= 0 ) {
		fprintf(stderr, "%3i: ", error_print_val);
	}
	fprintf(stderr, "%s", msg);
	pthread_mutex_unlock(&stderr_lock);
}


/* Function to call with ERROR/STATUS messages */
LogMsgFunc log_msg_func = log_to_stderr;
void *log_msg_vp = NULL;


void set_log_message_func(LogMsgFunc new_log_msg_func, void *vp)
{
	log_msg_func = new_log_msg_func;
	log_msg_vp = vp;
}


void STATUS(const char *format, ...)
{
	va_list args;
	char tmp[1024];
	va_start(args, format);
	vsnprintf(tmp, 1024, format, args);
	va_end(args);
	log_msg_func(LOG_MSG_STATUS, tmp, log_msg_vp);
}


void ERROR(const char *format, ...)
{
	va_list args;
	char tmp[1024];
	va_start(args, format);
	vsnprintf(tmp, 1024, format, args);
	va_end(args);
	log_msg_func(LOG_MSG_ERROR, tmp, log_msg_vp);
}


/* ------------------------------ Useful functions ---------------------------- */

int convert_int(const char *str, int *pval)
{
	long int val;
	char *rval;

	val = strtol(str, &rval, 10);
	if ( *rval != '\0' ) {
		return 1;
	} else {
		*pval = val;
		return 0;
	}
}


int convert_float(const char *str, double *pval)
{
	double val;
	char *rval;

	val = strtod(str, &rval);
	if ( *rval != '\0' ) {
		return 1;
	} else {
		*pval = val;
		return 0;
	}
}


size_t notrail(char *s)
{
	ssize_t i;
	size_t munched = 0;

	for ( i=strlen(s)-1; i>=0; i-- ) {
		if ( (s[i] == ' ') || (s[i] == '\t') ) {
			s[i] = '\0';
			munched++;
		} else {
			return munched;
		}
	}

	return munched;
}


void chomp(char *s)
{
	size_t i;
	size_t len;

	if ( s == NULL ) return;
	len = strlen(s);

	for ( i=0; i<len; i++ ) {
		if ( (s[i] == '\n') || (s[i] == '\r') ) {
			s[i] = '\0';
			return;
		}
	}
}


void progress_bar(int val, int total, const char *text)
{
	double frac;
	int n, i;
	char s[1024];
	const int width = 50;

	if ( total == 0 ) return;

	if ( !isatty(STDERR_FILENO) ) return;
	if ( tcgetpgrp(STDERR_FILENO) != getpgrp() ) return;

	frac = (double)val/total;
	n = (int)(frac*width);

	for ( i=0; i<n; i++ ) s[i] = '=';
	for ( i=n; i<width; i++ ) s[i] = '.';
	s[width] = '\0';

	pthread_mutex_lock(&stderr_lock);
	fprintf(stderr, "\r%s: |%s|", text, s);
	if ( val == total ) fprintf(stderr, "\n");
	pthread_mutex_unlock(&stderr_lock);

	fflush(stdout);
}


double random_flat(gsl_rng *rng, double max)
{
	return max * gsl_rng_uniform(rng);
}


double flat_noise(gsl_rng *rng, double expected, double width)
{
	double noise = random_flat(rng, 2.0*width);
	return expected+noise-width;
}


double gaussian_noise(gsl_rng *rng, double expected, double stddev)
{
	double x1, x2, noise;

	/* Generate two uniformly distributed random numbers between 0 and 1,
	 * including 1 but not 0. */
	x1 = 1.0 - gsl_rng_uniform(rng);
	x2 = 1.0 - gsl_rng_uniform(rng);

	noise = sqrt(-2.0*log(x1)) * cos(2.0*M_PI*x2);

	return expected + noise*stddev;
}


static int fake_poisson_noise(gsl_rng *rng, double expected)
{
	double rf = gaussian_noise(rng, expected, sqrt(expected));
	return (int)rf;
}


int poisson_noise(gsl_rng *rng, double expected)
{
	double L;
	int k = 0;
	double p = 1.0;

	/* For large values of the mean, we get big problems with arithmetic.
	 * In such cases, fall back on a Gaussian with the right variance. */
	if ( expected > 100.0 ) return fake_poisson_noise(rng, expected);

	L = exp(-expected);

	do {

		double r;

		k++;
		r = gsl_rng_uniform(rng);
		p *= r;

	} while ( p > L );

	return k - 1;
}


/* Return non-zero if c is in delims */
static int assplode_isdelim(const char c, const char *delims)
{
	size_t i;
	for ( i=0; i<strlen(delims); i++ ) {
		if ( c == delims[i] ) return 1;
	}
	return 0;
}


static int assplode_extract(char ***pbits, int n, size_t n_captured,
                            size_t start, const char *a)
{
	char **bits = *pbits;
	bits = realloc(bits, sizeof(char *)*(n+1));
	bits[n] = malloc(n_captured+1);
	memcpy(bits[n], a+start, n_captured);
	bits[n][n_captured] = '\0';
	n++;
	*pbits = bits;
	return n;
}


/**
 * Split the string 'a' using 'delims' as a zero-terminated list of
 *  deliminators.
 * Store each segment in bits[0...n] where n is the number of segments and is
 *  the return value.  pbits = &bits
 * Each segment needs to be freed with free() when finished with.
 * The array of bits also needs to be freed with free() when finished with,
 *  unless n=0 in which case bits==NULL
 */
int assplode(const char *a, const char *delims, char ***pbits,
             AssplodeFlag flags)
{
	size_t i, start, n_captured;
	int n, last_was_delim;
	char **bits;

	n = 0;
	i = 0;
	n_captured = 0;
	start = 0;
	last_was_delim = 0;
	bits = NULL;
	while ( i < strlen(a) ) {

		if ( assplode_isdelim(a[i], delims) ) {

			if ( n_captured > 0 ) {
				/* This is a deliminator after a sequence of
				 * non-deliminator chars */
				n = assplode_extract(&bits, n, n_captured,
				                     start, a);
			}

			n_captured = 0;
			if ( (flags & ASSPLODE_DUPS) && last_was_delim ) {
				n = assplode_extract(&bits, n, 0, start, a);
			}
			last_was_delim = 1;

		} else {

			if ( n_captured == 0 ) {
				/* No characters currently found, so this is the
				 * start */
				start = i;
			}
			n_captured++;
			last_was_delim = 0;

		}

		i++;

	}
	/* Left over characters at the end? */
	if ( n_captured > 0 ) {
		n = assplode_extract(&bits, n, n_captured, start, a);
	}

	*pbits = bits;
	return n;
}


char *check_prefix(char *prefix)
{
	int r;
	struct stat statbuf;
	char *new;
	size_t len;

	if ( prefix[0] == '\0' ) return prefix;

	/* Is "prefix" a directory? */
	r = stat(prefix, &statbuf);
	if ( r != 0 ) {
		/* "prefix" probably doesn't exist.  This is fine - assume
		 * the user knows what they're doing, and that "prefix"
		 * suffixed with the actual filename will produce something
		 * sensible. */
		return prefix;
	}

	if ( !S_ISDIR(statbuf.st_mode) ) {
		/* Also fine, as above. */
		return prefix;
	}

	/* Does the prefix end in a slash? */
	if ( prefix[strlen(prefix)-1] == '/' ) {
		/* This looks sensible. */
		return prefix;
	}

	STATUS("Your prefix ('%s') is a directory, but doesn't end"
	       " with a slash.  I'm going to add it for you.\n", prefix);
	STATUS("If this isn't what you want, run with --no-check-prefix.\n");
	len = strlen(prefix)+2;
	new = malloc(len);
	snprintf(new, len, "%s/", prefix);
	free(prefix);
	return new;
}


char *safe_strdup(const char *in)
{
	if ( in == NULL ) return NULL;
	return strdup(in);
}


char *safe_basename(const char *in)
{
	int i;
	char *cpy;
	char *res;

	cpy = strdup(in);

	/* Get rid of any trailing slashes */
	for ( i=strlen(cpy)-1; i>0; i-- ) {
		if ( cpy[i] == '/' ) {
			cpy[i] = '\0';
		} else {
			break;
		}
	}

	/* Find the base name */
	for ( i=strlen(cpy)-1; i>0; i-- ) {
		if ( cpy[i] == '/' ) {
			i++;
			break;
		}
	}

	res = strdup(cpy+i);
	/* If we didn't find a previous slash, i==0 so res==cpy */

	free(cpy);

	return res;
}


void strip_extension(char *bfn)
{
	size_t r = strlen(bfn)-1;
	while ( r > 0 ) {
		if ( bfn[r] == '.') {
			bfn[r] = '\0';
			return;
		}
		r--;
	}
}


const char *filename_extension(const char *fn, const char **pext2)
{
	const char *ext = NULL;
	const char *ext2 = NULL;
	size_t r = strlen(fn)-1;

	while ( r > 0 ) {
		if ( fn[r] == '.' ) {
			if ( ext != NULL ) {
				ext2 = fn+r;
				break;
			} else {
				ext = fn+r;
			}
		}
		r--;
	}

	if ( pext2 != NULL ) *pext2 = ext2;
	return ext;
}


/* Force the linker to bring in CBLAS to make GSL happy */
void utils_fudge_gslcblas()
{
        STATUS("%p\n", cblas_sgemm);
}


/**
 * \param q A \ref quaternion
 *
 * If a quaternion represents a pure rotation, its modulus should be unity.
 *
 * \returns The modulus of the given quaternion.
 **/
double quaternion_modulus(struct quaternion q)
{
	return sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
}


/**
 * \param q A \ref quaternion
 *
 * Rescales the quaternion such that its modulus is unity.
 *
 * \returns The normalised version of \p q
 **/
struct quaternion normalise_quaternion(struct quaternion q)
{
	double mod;
	struct quaternion r;

	mod = quaternion_modulus(q);

	r.w = q.w / mod;
	r.x = q.x / mod;
	r.y = q.y / mod;
	r.z = q.z / mod;

	return r;
}


/**
 * \param rng A GSL random number generator to use
 *
 * \returns A randomly generated, normalised, quaternion.
 **/
struct quaternion random_quaternion(gsl_rng *rng)
{
	struct quaternion q;

	q.w = 2.0*gsl_rng_uniform(rng) - 1.0;
	q.x = 2.0*gsl_rng_uniform(rng) - 1.0;
	q.y = 2.0*gsl_rng_uniform(rng) - 1.0;
	q.z = 2.0*gsl_rng_uniform(rng) - 1.0;
	q = normalise_quaternion(q);

	return q;
}


/**
 * \param q A \ref quaternion
 *
 * Checks if the given quaternion is normalised.
 *
 * This function performs a nasty floating point comparison of the form
 * <code>(modulus > 0.999) && (modulus < 1.001)</code>, and so should not be
 * relied upon to spot anything other than the most obvious input error.
 *
 * \returns 1 if the quaternion is normalised, 0 if not.
 **/
int quaternion_valid(struct quaternion q)
{
	double qmod;

	qmod = quaternion_modulus(q);

	/* Modulus = 1 to within some tolerance?
	 * Nasty allowance for floating-point accuracy follows... */
	if ( (qmod > 0.999) && (qmod < 1.001) ) return 1;

	return 0;
}


/**
 * \param q A vector (in the form of an \ref rvec struct)
 * \param z A \ref quaternion
 *
 * Rotates a vector according to a quaternion.
 *
 * \returns a rotated version of \p p.
 **/
struct rvec quat_rot(struct rvec q, struct quaternion z)
{
	struct rvec res;
	double t01, t02, t03, t11, t12, t13, t22, t23, t33;

	t01 = z.w*z.x;
	t02 = z.w*z.y;
	t03 = z.w*z.z;
	t11 = z.x*z.x;
	t12 = z.x*z.y;
	t13 = z.x*z.z;
	t22 = z.y*z.y;
	t23 = z.y*z.z;
	t33 = z.z*z.z;

	res.u = (1.0 - 2.0 * (t22 + t33)) * q.u
	            + (2.0 * (t12 + t03)) * q.v
	            + (2.0 * (t13 - t02)) * q.w;

	res.v =       (2.0 * (t12 - t03)) * q.u
	      + (1.0 - 2.0 * (t11 + t33)) * q.v
	            + (2.0 * (t01 + t23)) * q.w;

	res.w =       (2.0 * (t02 + t13)) * q.u
	            + (2.0 * (t23 - t01)) * q.v
	      + (1.0 - 2.0 * (t11 + t22)) * q.w;

	return res;
}


char *load_entire_file(const char *filename)
{
	struct stat statbuf;
	int r;
	char *contents;
	FILE *fh;

	r = stat(filename, &statbuf);
	if ( r != 0 ) {
		ERROR("File '%s' not found\n", filename);
		return NULL;
	}

	contents = malloc(statbuf.st_size+1);
	if ( contents == NULL ) {
		ERROR("Failed to allocate memory for file\n");
		return NULL;
	}

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		ERROR("Failed to open file '%s'\n", filename);
		free(contents);
		return NULL;
	}

	if ( fread(contents, 1, statbuf.st_size, fh) != statbuf.st_size ) {
		ERROR("Failed to read file '%s'\n", filename);
		fclose(fh);
		free(contents);
		return NULL;
	}
	contents[statbuf.st_size] = '\0';

	fclose(fh);

	return contents;
}


int file_exists(const char *filename)
{
	struct stat statbuf;
	int r;

	r = stat(filename, &statbuf);
	if ( r != 0 ) {
		return 0;
	}

	return 1;
}


int compare_double(const void *av, const void *bv)
{
	double a = *(double *)av;
	double b = *(double *)bv;
	if ( a > b ) return 1;
	if ( a < b ) return -1;
	return 0;
}


void *srealloc(void *arr, size_t new_size)
{
	void *new_arr = realloc(arr, new_size);
	if ( new_arr == NULL ) {
		free(arr);
		return NULL;
	} else {
		return new_arr;
	}
}


/* -------------------------- libcrystfel features  ------------------------ */

int crystfel_has_peakfinder9()
{
#ifdef HAVE_FDIP
	return 1;
#else
	return 0;
#endif
}
