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


/* -------------------------- Fundamental constants  ------------------------ */

/* Electron charge in C */
#define ELECTRON_CHARGE (1.6021773e-19)

/* Planck's constant (Js) */
#define PLANCK (6.62606896e-34)

/* Speed of light in vacuo (m/s) */
#define C_VACUO (299792458)

/* Thomson scattering length (m) */
#define THOMSON_LENGTH (2.81794e-15)

/* Density of water in kg/m^3 */
#define WATER_DENSITY (1.0e6)

/* Molar mass of water, in kg/mol */
#define WATER_MOLAR_MASS (18.01528e3)

/* Avogadro's number */
#define AVOGADRO (6.022e23)


/* ------------------------------ Quaternions ------------------------------- */

struct quaternion
{
	double w;
	double x;
	double y;
	double z;
};

extern struct quaternion normalise_quaternion(struct quaternion q);
extern double quaternion_modulus(struct quaternion q);
extern struct quaternion random_quaternion(void);
extern int quaternion_valid(struct quaternion q);


/* --------------------------- Useful functions ----------------------------- */

extern double angle_between(double x1, double y1, double z1,
                            double x2, double y2, double z2);
extern size_t skipspace(const char *s);
extern void chomp(char *s);
extern void progress_bar(int val, int total, const char *text);
extern double poisson_noise(double expected);

/* Keep these ones inline, to avoid function call overhead */
static inline struct quaternion invalid_quaternion(void)
{
	struct quaternion quat;
	quat.w = 0.0;
	quat.x = 0.0;
	quat.y = 0.0;
	quat.z = 0.0;
	return quat;
}

static inline double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}

static inline double modulus(double x, double y, double z)
{
	return sqrt(x*x + y*y + z*z);
}

static inline double modulus_squared(double x, double y, double z) {
	return x*x + y*y + z*z;
}

static inline double distance3d(double x1, double y1, double z1,
                  double x2, double y2, double z2)
{
	return modulus(x1-x2, y1-y2, z1-z2);
}


/* ----------------------------- Useful macros ------------------------------ */

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


/* -------------------- Indexed lists for specified types ------------------- */

/* Maxmimum index to hold values up to (can be increased if necessary) */
#define INDMAX 40

/* Array size */
#define IDIM (INDMAX*2 +1)

/* Create functions for storing reflection intensities indexed as h,k,l */
#define LABEL(x) x##_intensity
#define TYPE double
#include "list_tmp.h"

/* As above, but for complex structure factors */
#define LABEL(x) x##_sfac
#define TYPE double complex
#include "list_tmp.h"

/* As above, but for (unsigned) integer counts */
#define LABEL(x) x##_count
#define TYPE unsigned int
#include "list_tmp.h"

#endif	/* UTILS_H */
