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


void progress_bar(int val, int total)
{
	double frac;
	int n, i;
	char s[1024];

	frac = (double)val/total;
	n = (int)(frac*80);

	for ( i=0; i<n; i++ ) s[i] = '=';
	for ( i=n; i<80; i++ ) s[i] = '.';
	s[80] = '\0';
	printf("\r|%s|", s);

	if ( val == total ) printf("\n");

	fflush(stdout);
}
