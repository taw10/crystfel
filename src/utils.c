/*
 * utils.c
 *
 * Utility stuff
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "utils.h"


/* Angle between two vectors.  Answer in radians */
double angle_between(double x1, double y1, double z1,
                     double x2, double y2, double z2)
{
	double mod1 = modulus(x1, y1, z1);
	double mod2 = modulus(x2, y2, z2);
	return acos( (x1*x2 + y1*y2 + z1*z2) / (mod1*mod2) );
}


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

	frac = (double)val/total;
	n = (int)(frac*width);

	for ( i=0; i<n; i++ ) s[i] = '=';
	for ( i=n; i<width; i++ ) s[i] = '.';
	s[width] = '\0';
	printf("\r%s: |%s|", text, s);

	if ( val == total ) printf("\n");

	fflush(stdout);
}


int poisson_noise(double expected)
{
	double L;
	int k = 0;
	double p = 1.0;

	L = exp(-expected);

	do {

		double r;

		k++;
		r = (double)random()/RAND_MAX;
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

	q.w = (double)random()/RAND_MAX;
	q.x = (double)random()/RAND_MAX;
	q.y = (double)random()/RAND_MAX;
	q.z = (double)random()/RAND_MAX;
	normalise_quaternion(q);

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
