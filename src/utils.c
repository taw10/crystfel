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
