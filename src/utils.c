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


/* Return the MOST POSITIVE of two numbers */
unsigned int biggest(signed int a, signed int b)
{
	if ( a>b ) {
		return a;
	}
	return b;
}


/* Return the LEAST POSITIVE of two numbers */
unsigned int smallest(signed int a, signed int b)
{
	if ( a<b ) {
		return a;
	}
	return b;
}


double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}


double modulus(double x, double y, double z)
{
	return sqrt(x*x + y*y + z*z);
}


double modulus_squared(double x, double y, double z) {
	return x*x + y*y + z*z;
}


double distance3d(double x1, double y1, double z1,
                  double x2, double y2, double z2)
{
	return modulus(x1-x2, y1-y2, z1-z2);
}


/* Angle between two vectors.  Answer in radians */
double angle_between(double x1, double y1, double z1,
                     double x2, double y2, double z2)
{
	double mod1 = modulus(x1, y1, z1);
	double mod2 = modulus(x2, y2, z2);
	return acos( (x1*x2 + y1*y2 + z1*z2) / (mod1*mod2) );
}


/* As above, answer in degrees */
double angle_between_d(double x1, double y1, double z1,
                       double x2, double y2, double z2)
{
	return rad2deg(angle_between(x1, y1, z1, x2, y2, z2));
}


/* Wavelength of an electron (in m) given accelerating potential (in V) */
double lambda(double V)
{
	double m = 9.110E-31;
	double h = 6.625E-34;
	double e = 1.60E-19;
	double c = 2.998E8;

	return h / sqrt(2*m*e*V*(1+((e*V) / (2*m*c*c))));
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


int sign(double a)
{
	if ( a < 0 ) return -1;
	if ( a > 0 ) return +1;
	return 0;
}


void mapping_rotate(double x, double y, double z,
                    double *ddx, double *ddy, double *ddz,
                    double omega, double tilt)
{
	double nx, ny, nz;
	double x_temp, y_temp, z_temp;

	/* First: rotate image clockwise until tilt axis is aligned
	 * horizontally. */
	nx = x*cos(omega) + y*sin(omega);
	ny = -x*sin(omega) + y*cos(omega);
	nz = z;

	/* Now, tilt about the x-axis ANTICLOCKWISE around +x, i.e. the
	 * "wrong" way. This is because the crystal is rotated in the
	 * experiment, not the Ewald sphere. */
	x_temp = nx; y_temp = ny; z_temp = nz;
	nx = x_temp;
	ny = cos(tilt)*y_temp + sin(tilt)*z_temp;
	nz = -sin(tilt)*y_temp + cos(tilt)*z_temp;

	/* Finally, reverse the omega rotation to restore the location of the
	 * image in 3D space */
	x_temp = nx; y_temp = ny; z_temp = nz;
	nx = x_temp*cos(-omega) + y_temp*sin(-omega);
	ny = -x_temp*sin(-omega) + y_temp*cos(-omega);
	nz = z_temp;

	*ddx = nx;
	*ddy = ny;
	*ddz = nz;
}


void progress_bar(int val, int total)
{
	double frac;
	int n, i;
	char s[1024];

	frac = (double)val/total;
	n = frac*80;

	for ( i=0; i<n; i++ ) s[i] = '=';
	for ( i=n; i<80; i++ ) s[i] = '.';
	s[80] = '\0';
	printf("\r|%s|", s);

	if ( val == total ) printf("\n");
	
	fflush(stdout);
}
