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


/* Electron charge in C */
#define ELECTRON_CHARGE (1.6021773e-19)

/* Planck's constant (Js) */
#define PLANCK (6.62606896e-34)

/* Speed of light in vacuo (m/s) */
#define C_VACUO (299792458)

/* Thomson scattering length (m) */
#define THOMSON_LENGTH (2.81794e-15)


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

#endif	/* UTILS_H */
