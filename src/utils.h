/*
 * utils.h
 *
 * Utility stuff
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
extern struct rvec quat_rot(struct rvec q, struct quaternion z);


/* --------------------------- Useful functions ----------------------------- */

extern size_t notrail(char *s);
extern void chomp(char *s);

typedef enum {
	ASSPLODE_NONE	= 0,
	ASSPLODE_DUPS	= 1<<0
} AssplodeFlag;
extern int assplode(const char *a, const char *delims, char ***pbits,
                    AssplodeFlag flags);

extern void progress_bar(int val, int total, const char *text);
extern int poisson_noise(double expected);

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


/* -------------------- Indexed lists for specified types ------------------- */

#include "defs.h"

#define LIST_SIZE (IDIM*IDIM*IDIM)

/* Create functions for storing reflection intensities indexed as h,k,l */
#define LABEL(x) x##_intensity
#define TYPE double
#include "list_tmp.h"

/* CAs above, but for phase values */
#define LABEL(x) x##_phase
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

/* As above, but for sigmas */
#define LABEL(x) x##_sigma
#define TYPE double
#include "list_tmp.h"

/* As above, but for simple flags */
#define LABEL(x) x##_flag
#define TYPE unsigned char
#include "list_tmp.h"


/* ----------- Reflection lists indexed by sequence (not indices) ----------- */

typedef struct _reflitemlist ReflItemList;  /* Opaque */

struct refl_item {
	signed int h;
	signed int k;
	signed int l;
	int op;
};

extern void clear_items(ReflItemList *items);
extern ReflItemList *new_items(void);
extern void delete_items(ReflItemList *items);
extern void add_item(ReflItemList *items,
                     signed int h, signed int k, signed int l);
extern void add_item_with_op(ReflItemList *items,
                             signed int h, signed int k, signed int l, int op);
extern int find_item(ReflItemList *items,
                     signed int h, signed int k, signed int l);
extern struct refl_item *get_item(ReflItemList *items, int i);
extern int num_items(const ReflItemList *items);
extern unsigned int *items_to_counts(ReflItemList *items);
extern void union_op_items(ReflItemList *items, ReflItemList *newi);
extern void union_items(ReflItemList *items, ReflItemList *newi);
extern ReflItemList *intersection_items(ReflItemList *i1, ReflItemList *i2);


/* ------------------------------ Message macros ---------------------------- */

#define ERROR(...) fprintf(stderr, __VA_ARGS__)
#define STATUS(...) fprintf(stderr, __VA_ARGS__)


/* ------------------------------ File handling ----------------------------- */

extern char *check_prefix(char *prefix);
extern char *safe_basename(const char *in);


#endif	/* UTILS_H */
