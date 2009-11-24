/*
 * utils.h
 *
 * Utility stuff
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */

#ifndef UTILS_H
#define UTILS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>


/* Electron charge in C */
#define ELECTRON_CHARGE (1.6021773e-19)

/* Planck's constant (Js) */
#define PLANCK (6.62606896e-34)

/* Speed of light in vacuo (m/s) */
#define C_VACUO (299792458)

/* Thomson scattering length (m) */
#define THOMSON_LENGTH (2.81794e-15)

/* Maxmimum index to go up to */
#define INDMAX 20
#define IDIM (INDMAX*2 +1)


extern unsigned int biggest(signed int a, signed int b);
extern unsigned int smallest(signed int a, signed int b);
extern double distance(double x1, double y1, double x2, double y2);
extern double modulus(double x, double y, double z);
extern double modulus_squared(double x, double y, double z);
extern double angle_between(double x1, double y1, double z1,
                            double x2, double y2, double z2);
extern double angle_between_d(double x1, double y1, double z1,
                              double x2, double y2, double z2);
extern double lambda(double voltage);
extern double distance3d(double x1, double y1, double z1,
                         double x2, double y2, double z2);
extern size_t skipspace(const char *s);
extern void chomp(char *s);
extern int sign(double a);
extern void mapping_rotate(double x, double y, double z,
                           double *ddx, double *ddy, double *ddz,
                           double omega, double tilt);
extern void progress_bar(int val, int total);

#define rad2deg(a) ((a)*180/M_PI)
#define deg2rad(a) ((a)*M_PI/180)

#define is_odd(a) ((a)%2==1)

/* Photon energy (J) to wavelength (m) */
#define ph_en_to_lambda(a) ((PLANCK*C_VACUO)/(a))

/* Photon wavelength (m) to energy (J) */
#define ph_lambda_to_en(a) ((PLANCK*C_VACUO)/(a))

/* eV to Joules */
#define eV_to_J(a) ((a)*ELECTRON_CHARGE)

/* Joules to eV */
#define J_to_eV(a) ((a)/ELECTRON_CHARGE)


static inline void integrate_reflection(double complex *ref, signed int h,
                                        signed int k, signed int l,
                                        double complex i)
{
	int idx;

	/* Not interested in central beam */
	if ( (h==0) && (k==0) && (l==0) ) return;

	if ( h < 0 ) h += IDIM;
	if ( k < 0 ) k += IDIM;
	if ( l < 0 ) l += IDIM;

	idx = h + (IDIM*k) + (IDIM*IDIM*l);
	ref[idx] += i;
}


static inline double complex get_integral(double complex *ref, signed int h,
                                          signed int k, signed int l)
{
	int idx;

	if ( h < 0 ) h += IDIM;
	if ( k < 0 ) k += IDIM;
	if ( l < 0 ) l += IDIM;

	idx = h + (IDIM*k) + (IDIM*IDIM*l);
	return ref[idx];
}


static inline double complex *reflist_new(void)
{
	double complex *r;
	size_t r_size;
	r_size = IDIM*IDIM*IDIM*sizeof(double complex);
	r = malloc(r_size);
	memset(r, 0, r_size);
	return r;
}


#endif	/* UTILS_H */
