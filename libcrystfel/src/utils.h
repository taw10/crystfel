/*
 * utils.h
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

#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <complex.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>


#include "thread-pool.h"


/**
 * \file utils.h
 * Miscellaneous utility functions
 */

/* -------------------------- Fundamental constants  ------------------------ */

/* Electron charge (Coulombs) */
#define ELECTRON_CHARGE (1.6021773e-19)

/* Electron rest mass (kg) */
#define ELECTRON_MASS (9.1093837015e-31)

/* Planck's constant (Js) */
#define PLANCK (6.62606896e-34)

/* Speed of light in vacuo (m/s) */
#define C_VACUO (299792458)

/* Thomson scattering length (m) */
#define THOMSON_LENGTH (2.81794e-15)


#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------- Useful functions ----------------------------- */

extern void show_matrix_eqn(gsl_matrix *M, gsl_vector *v);
extern void show_matrix(gsl_matrix *M);
extern void show_vector(gsl_vector *M);
extern gsl_vector *solve_svd(gsl_vector *v, gsl_matrix *M, int *n_filt,
                            int verbose);
extern gsl_matrix *matrix_mult2(gsl_matrix *A, gsl_matrix *B);
extern gsl_matrix *matrix_mult3(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C);
extern gsl_matrix *matrix_invert(gsl_matrix *m);

extern size_t notrail(char *s);
extern int convert_int(const char *str, int *pval);
extern int convert_float(const char *str, double *pval);
extern void chomp(char *s);
extern void *srealloc(void *arr, size_t new_size);

#define CLEAR_BIT(val, bit) (((val) | (bit)) ^ (bit))

/**
 * Controls the behaviour of \ref assplode.
 **/
typedef enum {
	ASSPLODE_NONE	= 0,    /**< Nothing */
	ASSPLODE_DUPS	= 1<<0  /**< Don't merge deliminators */
} AssplodeFlag;
extern int assplode(const char *a, const char *delims, char ***pbits,
                    AssplodeFlag flags);

extern void progress_bar(int val, int total, const char *text);
extern double random_flat(gsl_rng *rng, double max);
extern double flat_noise(gsl_rng *rng, double expected, double width);
extern double gaussian_noise(gsl_rng *rng, double expected, double stddev);
extern int poisson_noise(gsl_rng *rng, double expected);

static inline double sq(double a)
{
	return a*a;
}

static inline double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}

static inline double modulus(double x, double y, double z)
{
	return sqrt(x*x + y*y + z*z);
}

static inline double modulus2d(double x, double y)
{
	return sqrt(x*x + y*y);
}

static inline double modulus_squared(double x, double y, double z) {
	return x*x + y*y + z*z;
}

static inline double distance3d(double x1, double y1, double z1,
                                double x2, double y2, double z2)
{
	double d = modulus(x1-x2, y1-y2, z1-z2);
	return d;
}

/* Answer in radians */
static inline double angle_between(double x1, double y1, double z1,
                                   double x2, double y2, double z2)
{
	double mod1 = modulus(x1, y1, z1);
	double mod2 = modulus(x2, y2, z2);
	double cosine = (x1*x2 + y1*y2 + z1*z2) / (mod1*mod2);

	/* Fix domain if outside due to rounding */
	if ( cosine > 1.0 ) cosine = 1.0;
	if ( cosine < -1.0 ) cosine = -1.0;

	return acos(cosine);
}

/* Answer in radians */
static inline double angle_between_2d(double x1, double y1,
                                      double x2, double y2)
{
	double mod1 = modulus2d(x1, y1);
	double mod2 = modulus2d(x2, y2);
	double cosine = (x1*x2 + y1*y2) / (mod1*mod2);

	/* Fix domain if outside due to rounding */
	if ( cosine > 1.0 ) cosine = 1.0;
	if ( cosine < -1.0 ) cosine = -1.0;

	return acos(cosine);
}

static inline int within_tolerance(double a, double b, double percent)
{
	double tol = fabs(a) * (percent/100.0);
	if ( fabs(b-a) < tol ) return 1;
	return 0;
}

extern void rotate2d(double *x, double *y, double cx, double cy, double ang);


/* ----------------------------- Useful macros ------------------------------ */

#define rad2deg(a) ((a)*180/M_PI)
#define deg2rad(a) ((a)*M_PI/180)

#define is_odd(a) ((a)%2==1)

#define biggest(a,b) ((a>b) ? (a) : (b))
#define smallest(a,b) ((a<b) ? (a) : (b))


/* Photon energy (J) to wavelength (m) */
#define ph_en_to_lambda(a) ((PLANCK*C_VACUO)/(a))

/* Photon wavelength (m) to energy (J) */
#define ph_lambda_to_en(a) ((PLANCK*C_VACUO)/(a))

/* eV to Joules */
#define eV_to_J(a) ((a)*ELECTRON_CHARGE)

/* Joules to eV */
#define J_to_eV(a) ((a)/ELECTRON_CHARGE)

/* Photon wavelength (m) to energy (eV) */
#define ph_lambda_to_eV(a) J_to_eV(ph_lambda_to_en(a))

/* Photon energy (eV) to wavelength (m) */
#define ph_eV_to_lambda(a) ph_en_to_lambda(eV_to_J(a))

/* Photon energy (eV) to k (1/m) */
#define ph_eV_to_k(a) ((a)*ELECTRON_CHARGE/PLANCK/C_VACUO)

/* Electron accelerating voltage (V) to wavelength (m) */
static inline double el_V_to_lambda(double E)
{
	double Estar;

	/* Relativistically corrected accelerating voltage */
	Estar = E * (1.0 + E * ELECTRON_CHARGE/(2.0*ELECTRON_MASS*C_VACUO*C_VACUO));

	return PLANCK / sqrt(2.0*ELECTRON_MASS*ELECTRON_CHARGE*Estar);
}


/* ------------------------------ Message logging ---------------------------- */

extern pthread_mutex_t stderr_lock;

enum log_msg_type {
                   LOG_MSG_STATUS,
                   LOG_MSG_ERROR
};


extern void STATUS(const char *format, ...);
extern void ERROR(const char *format, ...);

typedef void (*LogMsgFunc)(enum log_msg_type type, const char *msg, void *vp);

extern void set_log_message_func(LogMsgFunc new_log_msg_func,
                                 void *vp);


/* ---------------------------- Memory management --------------------------- */

extern void *cfmalloc(size_t size);
extern void cffree(void *ptr);
extern void *cfcalloc(size_t nmemb, size_t size);
extern void *cfrealloc(void *ptr, size_t size);
extern char *cfstrdup(const char *s);
extern char *cfstrndup(const char *s, size_t n);
extern int set_mm_funcs(void *(*cfmalloc)(size_t size),
                        void (*cffree)(void *ptr),
                        void *(*cfcalloc)(size_t nmemb, size_t size),
                        void *(*cfrealloc)(void *ptr, size_t size));


/* ------------------------------ File handling ----------------------------- */

extern char *check_prefix(char *prefix);
extern char *safe_basename(const char *in);
extern char *safe_strdup(const char *in);
extern void strip_extension(char *bfn);
extern char *load_entire_file(const char *filename);
extern int file_exists(const char *filename);
extern int is_dir(const char *filename);
extern const char *filename_extension(const char *fn, const char **ext2);


/* ------------------------------ Useful stuff ------------------------------ */

#if __GNUC__ >= 3
#define UNUSED __attribute__((unused))
#define likely(x) __builtin_expect (!!(x), 1)
#define unlikely(x) __builtin_expect (!!(x), 0)
#else
#define UNUSED
#define likely(x) (x)
#define unlikely(x) (x)
#endif


/* ------------------------------ Quaternions ------------------------------- */

/**
 * A structure representing a quaternion.
 **/
struct quaternion;

struct quaternion {
	double w;
	double x;
	double y;
	double z;
};

extern struct quaternion normalise_quaternion(struct quaternion q);
extern double quaternion_modulus(struct quaternion q);
extern struct quaternion random_quaternion(gsl_rng *rng);
extern int quaternion_valid(struct quaternion q);
extern struct rvec quat_rot(struct rvec q, struct quaternion z);

extern int compare_double(const void *av, const void *bv);

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


/**
 * \param x value
 * \param w weight
 * \param sumw pointer to accumulator variable for the sum of weights
 * \param mean pointer to online mean value
 * \param M2 pointer to online variance times sumw
 *
 * Function to compute mean and variance stably
 */
static inline void mean_variance(const double x, const double w,
                                 double *sumw, double *mean, double *M2)
{
	const double temp  = w + *sumw;
	const double delta = x - *mean;
	const double R = delta * w / temp;

	if ( w < DBL_MIN ) return;

	*mean += R;
	*M2   += *sumw*delta*R;
	*sumw  = temp;
}


/* -------------------------- libcrystfel features  ------------------------ */

extern int crystfel_has_peakfinder9(void);


#ifdef __cplusplus
}
#endif

#endif	/* UTILS_H */
