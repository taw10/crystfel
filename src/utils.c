/*
 * utils.c
 *
 * Utility stuff
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#include <libgen.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "utils.h"
#include "image.h"


/**
 * SECTION:utils
 * @short_description: Miscellaneous utilities
 * @title: Utilities
 * @section_id:
 * @see_also:
 * @include: "utils.h"
 * @Image:
 *
 * Wibble
 */

/**
 * show_matrix_eqn:
 * @M: A matrix
 * @v: A vector
 * @r: The number of elements in @v and the side length of @M.
 *
 * Displays a matrix equation of the form @M.a = @v.
 *
 * @M must be square.
 **/
void show_matrix_eqn(gsl_matrix *M, gsl_vector *v, int r)
{
	int i, j;

	for ( i=0; i<r; i++ ) {
		STATUS("[ ");
		for ( j=0; j<r; j++ ) {
			STATUS("%+9.3e ", gsl_matrix_get(M, i, j));
		}
		STATUS("][ a%2i ] = [ %+9.3e ]\n", i, gsl_vector_get(v, i));
	}
}


size_t notrail(char *s)
{
	size_t i;
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

	if ( !s ) return;

	for ( i=0; i<strlen(s); i++ ) {
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


static int fake_poisson_noise(double expected)
{
	double x1, x2, w;
	double noise, rf;

	do {

		x1 = 2.0 * ((double)random()/RAND_MAX) - 1.0;
		x2 = 2.0 * ((double)random()/RAND_MAX) - 1.0;
		w = pow(x1, 2.0) + pow(x2, 2.0);

	} while ( w >= 1.0 );

	w = sqrt((-2.0*log(w))/w);
	noise = w * x1;

	rf = expected + noise*sqrt(expected);

	return (int)rf;
}


int poisson_noise(double expected)
{
	double L;
	int k = 0;
	double p = 1.0;

	/* For large values of the mean, we get big problems with arithmetic.
	 * In such cases, fall back on a Gaussian with the right variance. */
	if ( expected > 100.0 ) return fake_poisson_noise(expected);

	L = exp(-expected);

	do {

		double r;

		k++;
		r = (double)random()/(double)RAND_MAX;
		p *= r;

	} while ( p > L );

	return k - 1;
}


/**
 * SECTION:quaternion
 * @short_description: Simple quaternion handling
 * @title: Quaternion
 * @section_id:
 * @see_also:
 * @include: "utils.h"
 * @Image:
 *
 * There is a simple quaternion structure in CrystFEL.  At the moment, it is
 * only used when simulating patterns, as an argument to %cell_rotate to
 * orient the unit cell.
 */

/**
 * quaternion_modulus:
 * @q: A %quaternion
 *
 * If a quaternion represents a pure rotation, its modulus should be unity.
 *
 * Returns: the modulus of the given quaternion.
 **/
double quaternion_modulus(struct quaternion q)
{
	return sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
}


/**
 * normalise_quaternion:
 * @q: A %quaternion
 *
 * Rescales the quaternion such that its modulus is unity.
 *
 * Returns: the normalised version of @q
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
 * random_quaternion:
 *
 * Returns: a randomly generated, normalised, quaternion.
 **/
struct quaternion random_quaternion()
{
	struct quaternion q;

	q.w = 2.0*(double)random()/RAND_MAX - 1.0;
	q.x = 2.0*(double)random()/RAND_MAX - 1.0;
	q.y = 2.0*(double)random()/RAND_MAX - 1.0;
	q.z = 2.0*(double)random()/RAND_MAX - 1.0;
	q = normalise_quaternion(q);

	return q;
}


/**
 * quaternion_valid:
 * @q: A %quaternion
 *
 * Checks if the given quaternion is normalised.
 *
 * This function performs a nasty floating point comparison of the form
 * <code>(modulus > 0.999) && (modulus < 1.001)</code>, and so should not be
 * relied upon to spot anything other than the most obvious input error.
 *
 * Returns: 1 if the quaternion is normalised, 0 if not.
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
 * quat_rot
 * @q: A vector (in the form of an %rvec)
 * @z: A %quaternion
 *
 * Rotates a vector according to a quaternion.
 *
 * Returns: A rotated version of @p.
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


/* Split the string 'a' using 'delims' as a zero-terminated list of
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

/**
 * SECTION:reflitemlist
 * @short_description: The index list and indexed arrays
 * @title: ReflItemList
 * @section_id:
 * @see_also:
 * @include: "utils.h"
 * @Image:
 *
 * Wibble
 */

struct _reflitemlist {
	struct refl_item *items;
	int n_items;
	int max_items;
};


void clear_items(ReflItemList *items)
{
	items->n_items = 0;
}


static void alloc_items(ReflItemList *items)
{
	items->items = realloc(items->items,
	                       items->max_items*sizeof(struct refl_item));
}


/**
 * new_items:
 *
 * Creates a new %ReflItemList.
 *
 * Returns: The new list, or NULL.
 **/
ReflItemList *new_items()
{
	ReflItemList *new;
	new = malloc(sizeof(ReflItemList));
	if ( new == NULL ) return NULL;
	new->max_items = 1024;
	new->n_items = 0;
	new->items = NULL;
	alloc_items(new);
	return new;
}


void delete_items(ReflItemList *items)
{
	if ( items == NULL ) return;
	if ( items->items != NULL ) free(items->items);
	free(items);
}


void add_item_with_op(ReflItemList *items, signed int h, signed int k,
                      signed int l, int op)
{
	if ( items->n_items == items->max_items ) {
		items->max_items += 1024;
		alloc_items(items);
	}

	items->items[items->n_items].h = h;
	items->items[items->n_items].k = k;
	items->items[items->n_items].l = l;
	items->items[items->n_items].op = op;
	items->n_items++;
}


void add_item(ReflItemList *items, signed int h, signed int k, signed int l)
{
	add_item_with_op(items, h, k, l, 0);
}


int find_item(ReflItemList *items,
                     signed int h, signed int k, signed int l)
{
	int i;

	for ( i=0; i<items->n_items; i++ ) {
		if ( items->items[i].h != h ) continue;
		if ( items->items[i].k != k ) continue;
		if ( items->items[i].l != l ) continue;
		return 1;
	}
	return 0;
}


static int find_op(ReflItemList *items, int op)
{
	int i;

	for ( i=0; i<items->n_items; i++ ) {
		if ( items->items[i].op == op ) return 1;
	}
	return 0;
}


struct refl_item *get_item(ReflItemList *items, int i)
{
	if ( i >= items->n_items ) return NULL;
	return &items->items[i];
}


int num_items(const ReflItemList *items)
{
	return items->n_items;
}


void union_op_items(ReflItemList *items, ReflItemList *newi)
{
	int n, i;

	n = num_items(newi);
	for ( i=0; i<n; i++ ) {

		struct refl_item *r = get_item(newi, i);
		if ( find_op(items, r->op) ) continue;

		add_item_with_op(items, r->h, r->k, r->l, r->op);

	}
}


void union_items(ReflItemList *items, ReflItemList *newi)
{
	int n, i;

	n = num_items(newi);
	for ( i=0; i<n; i++ ) {

		struct refl_item *r = get_item(newi, i);
		if ( find_item(items, r->h, r->k, r->l) ) continue;

		add_item_with_op(items, r->h, r->k, r->l, r->op);

	}
}


ReflItemList *intersection_items(ReflItemList *i1, ReflItemList *i2)
{
	int n, i;
	ReflItemList *res = new_items();

	n = num_items(i1);
	for ( i=0; i<n; i++ ) {

		struct refl_item *r = get_item(i1, i);
		if ( find_item(i2, r->h, r->k, r->l) ) {
			add_item_with_op(res, r->h, r->k, r->l, r->op);
		}

	}

	return res;
}


char *check_prefix(char *prefix)
{
	int r;
	struct stat statbuf;
	char *new;
	size_t len;

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
