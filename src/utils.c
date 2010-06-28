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

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

#include "utils.h"
#include "image.h"


size_t skipspace(const char *s)
{
	size_t i;

	for ( i=0; i<strlen(s); i++ ) {
		if ( (s[i] != ' ') && (s[i] != '\t') ) return i;
	}

	return strlen(s);
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

	if ( !isatty(STDERR_FILENO) ) return;
	if ( tcgetpgrp(STDERR_FILENO) != getpgrp() ) return;

	frac = (double)val/total;
	n = (int)(frac*width);

	for ( i=0; i<n; i++ ) s[i] = '=';
	for ( i=n; i<width; i++ ) s[i] = '.';
	s[width] = '\0';
	STATUS("\r%s: |%s|", text, s);

	if ( val == total ) STATUS("\n");

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


double quaternion_modulus(struct quaternion q)
{
	return sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
}


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


int quaternion_valid(struct quaternion q)
{
	double qmod;

	qmod = quaternion_modulus(q);

	/* Modulus = 1 to within some tolerance?
	 * Nasty allowance for floating-point accuracy follows... */
	if ( (qmod > 0.999) && (qmod < 1.001) ) return 1;

	return 0;
}


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

ReflItemList *new_items()
{
	ReflItemList *new;
	new = malloc(sizeof(ReflItemList));
	new->max_items = 1024;
	new->n_items = 0;
	new->items = NULL;
	alloc_items(new);
	return new;
}

void delete_items(ReflItemList *items)
{
	if ( items->items != NULL ) free(items->items);
	free(items);
}

void add_item(ReflItemList *items,
                     signed int h, signed int k, signed int l)
{
	if ( items->n_items == items->max_items ) {
		items->max_items += 1024;
		alloc_items(items);
	}

	items->items[items->n_items].h = h;
	items->items[items->n_items].k = k;
	items->items[items->n_items].l = l;
	items->n_items++;
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

struct refl_item *get_item(ReflItemList *items, int i)
{
	if ( i >= items->n_items ) return NULL;
	return &items->items[i];
}

int num_items(const ReflItemList *items)
{
	return items->n_items;
}

unsigned int *items_to_counts(ReflItemList *items)
{
	unsigned int *c;
	int i;

	c = new_list_count();

	for ( i=0; i<num_items(items); i++ ) {
		struct refl_item *r;
		r = get_item(items, i);
		set_count(c, r->h, r->k, r->l, 1);
	}

	return c;
}
