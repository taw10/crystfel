/*
 * diffraction.c
 *
 * Calculate diffraction patterns by Fourier methods
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <assert.h>

#include "image.h"
#include "utils.h"
#include "cell.h"
#include "diffraction.h"
#include "sfac.h"
#include "parameters-lcls.tmp"


#define SAMPLING (4)
#define BWSAMPLING (1)
#define BANDWIDTH (0.0 / 100.0)


static double lattice_factor(struct rvec q, double ax, double ay, double az,
                                            double bx, double by, double bz,
                                            double cx, double cy, double cz,
                                            int na, int nb, int nc)
{
	struct rvec Udotq;
	double f1, f2, f3;

	Udotq.u = ax*q.u + ay*q.v + az*q.w;
	Udotq.v = bx*q.u + by*q.v + bz*q.w;
	Udotq.w = cx*q.u + cy*q.v + cz*q.w;

	/* At exact Bragg condition, f1 = na */
	if ( na > 1 ) {
		f1 = sin(M_PI*(double)na*Udotq.u) / sin(M_PI*Udotq.u);
	} else {
		f1 = 1.0;
	}

	/* At exact Bragg condition, f2 = nb */
	if ( nb > 1 ) {
		f2 = sin(M_PI*(double)nb*Udotq.v) / sin(M_PI*Udotq.v);
	} else {
		f2 = 1.0;
	}

	/* At exact Bragg condition, f3 = nc */
	if ( nc > 1 ) {
		f3 = sin(M_PI*(double)nc*Udotq.w) / sin(M_PI*Udotq.w);
	} else {
		f3 = 1.0;
	}

	/* At exact Bragg condition, this will multiply the molecular
	 * part of the structure factor by the number of unit cells,
	 * as desired (more scattering from bigger crystal!) */
	return f1 * f2 * f3;
}


static double interpolate_linear(const double *ref,
                                 const unsigned int *counts,
                                 float hd, signed int k, signed int l)
{
	signed int h;
	double val1, val2;
	float f;
	unsigned int c1, c2;

	h = (signed int)hd;
	if ( hd < 0.0 ) h -= 1;
	f = hd - (float)h;
	assert(f >= 0.0);

	val1 = lookup_intensity(ref, h, k, l);
	val2 = lookup_intensity(ref, h+1, k, l);
	c1 = lookup_count(counts, h, k, l);
	c2 = lookup_count(counts, h+1, k, l);

	if ( c1 == 0 ) {
		ERROR("Needed intensity for %i %i %i, but don't have it.\n",
		      h, k, l);
		return 1.0e20;
	}

	if ( c2 == 0 ) {
		ERROR("Needed intensity for %i %i %i, but don't have it.\n",
		      h+1, k, l);
		return 1.0e20;
	}

	val1 = val1 / (double)c1;
	val2 = val2 / (double)c2;

	return (1.0-f)*val1 + f*val2;
}


static double interpolate_bilinear(const double *ref,
                                   const unsigned int *counts,
                                   float hd, float kd, signed int l)
{
	signed int k;
	double val1, val2;
	float f;

	k = (signed int)kd;
	if ( kd < 0.0 ) k -= 1;
	f = kd - (float)k;
	assert(f >= 0.0);

	val1 = interpolate_linear(ref, counts, hd, k, l);
	val2 = interpolate_linear(ref, counts, hd, k+1, l);

	return (1.0-f)*val1 + f*val2;
}


static double interpolate_intensity(const double *ref,
                                    const unsigned int *counts,
                                    float hd, float kd, float ld)
{
	signed int l;
	double val1, val2;
	float f;

	l = (signed int)ld;
	if ( ld < 0.0 ) l -= 1;
	f = ld - (float)l;
	assert(f >= 0.0);

	val1 = interpolate_bilinear(ref, counts, hd, kd, l);
	val2 = interpolate_bilinear(ref, counts, hd, kd, l+1);

	return (1.0-f)*val1 + f*val2;
}


/* Look up the structure factor for the nearest Bragg condition */
static double molecule_factor(const double *intensities,
                              const unsigned int *counts, struct rvec q,
                              double ax, double ay, double az,
                              double bx, double by, double bz,
                              double cx, double cy, double cz,
                              GradientMethod m)
{
	float hd, kd, ld;
	signed int h, k, l;
	double r;

	hd = q.u * ax + q.v * ay + q.w * az;
	kd = q.u * bx + q.v * by + q.w * bz;
	ld = q.u * cx + q.v * cy + q.w * cz;

	switch ( m ) {
	case GRADIENT_MOSAIC :
		h = (signed int)rint(hd);
		k = (signed int)rint(kd);
		l = (signed int)rint(ld);
		if ( lookup_count(counts, h, k, l) == 0 ) {
			ERROR("Needed intensity for %i %i %i, but don't have it.\n",
			      h, k, l);
			return 1.0e20;
		}
		r = lookup_intensity(intensities, h, k, l);
		break;
	case GRADIENT_INTERPOLATE :
		r = interpolate_intensity(intensities, counts, hd, kd, ld);
		break;
	case GRADIENT_PHASED :
	default:
		ERROR("This gradient method not implemented yet.\n");
		exit(1);
	}

	return r;
}


double water_diffraction(struct rvec q, double en,
                                 double beam_r, double water_r)
{
	double complex fH, fO;
	double s, modq;
	double width;
	double complex ifac;

	/* Interatomic distances in water molecule */
	const double rOH = 0.09584e-9;
	const double rHH = 0.1515e-9;

	/* Volume of water column, approximated as:
	 * (2water_r) * (2beam_r) * smallest(2beam_r, 2water_r)
	 * neglecting the curvature of the faces of the volume */
	if ( beam_r > water_r ) {
		width = 2.0 * water_r;
	} else {
		width = 2.0 * beam_r;
	}
	const double water_v = 2.0*beam_r * 2.0*water_r * width;

	/* Number of water molecules */
	const double n_water = water_v * WATER_DENSITY
	                                        * (AVOGADRO / WATER_MOLAR_MASS);

	/* s = sin(theta)/lambda = 1/2d = |q|/2 */
	modq = modulus(q.u, q.v, q.w);
	s = modq / 2.0;

	fH = get_sfac("H", s, en);
	fO = get_sfac("O", s, en);

	/* Four O-H cross terms */
	ifac = 4.0*fH*fO * sin(2.0*M_PI*modq*rOH)/(2.0*M_PI*modq*rOH);

	/* Three H-H cross terms */
	ifac += 3.0*fH*fH * sin(2.0*M_PI*modq*rHH)/(2.0*M_PI*modq*rHH);

	/* Three diagonal terms */
	ifac += 2.0*fH*fH + fO*fO;

	return cabs(ifac) * n_water;
}


struct rvec get_q(struct image *image, unsigned int xs, unsigned int ys,
                  unsigned int sampling, float *ttp, float k)
{
	struct rvec q;
	float twotheta, r, az;
	float rx = 0.0;
	float ry = 0.0;
	struct panel *p;

	const unsigned int x = xs / sampling;
	const unsigned int y = ys / sampling; /* Integer part only */

	p = find_panel(&image->det, x, y);

	rx = ((float)xs - (sampling*p->cx)) / (sampling * p->res);
	ry = ((float)ys - (sampling*p->cy)) / (sampling * p->res);

	/* Calculate q-vector for this sub-pixel */
	r = sqrt(pow(rx, 2.0) + pow(ry, 2.0));

	twotheta = atan2(r, p->clen);
	az = atan2(ry, rx);
	if ( ttp != NULL ) *ttp = twotheta;

	q.u = k * sin(twotheta)*cos(az);
	q.v = k * sin(twotheta)*sin(az);
	q.w = k - k * cos(twotheta);

	return quat_rot(q, image->orientation);
}


void get_diffraction(struct image *image, int na, int nb, int nc,
                     const double *intensities, const unsigned int *counts,
                     UnitCell *cell, int do_water, GradientMethod m)
{
	unsigned int xs, ys;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	float k, klow, bwstep;

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	/* Allocate (and zero) the "diffraction array" */
	image->data = calloc(image->width * image->height, sizeof(float));

	/* Needed later for Lorentz calculation */
	image->twotheta = malloc(image->width * image->height * sizeof(double));

	k = 1.0/image->lambda;  /* Centre value */
	klow = k - k*(BANDWIDTH/2.0);  /* Lower value */
	bwstep = k * BANDWIDTH / BWSAMPLING;

	for ( xs=0; xs<image->width*SAMPLING; xs++ ) {
	for ( ys=0; ys<image->height*SAMPLING; ys++ ) {

		struct rvec q;
		float twotheta;
		double sw = 1.0/(SAMPLING*SAMPLING);  /* Sample weight */

		const unsigned int x = xs / SAMPLING;
		const unsigned int y = ys / SAMPLING; /* Integer part only */

		int kstep;

		for ( kstep=0; kstep<BWSAMPLING; kstep++ ) {

			float k;
			double kw = 1.0/BWSAMPLING;
			double intensity;
			double f_lattice, I_lattice;
			double I_molecule;

			/* Calculate k this time round */
			k = klow + kstep * bwstep;

			q = get_q(image, xs, ys, SAMPLING, &twotheta, k);
			image->twotheta[x + image->width*y] = twotheta;

			f_lattice = lattice_factor(q, ax, ay, az,
			                              bx, by, bz,
			                              cx, cy, cz,
			                              na, nb, nc);

			if ( intensities == NULL ) {
				I_molecule = 1.0e10;
			} else {
				I_molecule = molecule_factor(intensities,
				                             counts, q,
			                                     ax,ay,az,
			                                     bx,by,bz,cx,cy,cz,
				                             m);
			}

			I_lattice = pow(f_lattice, 2.0);

			intensity = sw * kw * I_lattice * I_molecule;
			image->data[x + image->width*y] += intensity;

		}

		if ( do_water ) {

			/* Bandwidth not simulated for water */
			struct rvec q;

			q = get_q(image, x, y, 1, NULL, 1.0/image->lambda);

			/* Add intensity contribution from water */
			image->data[x + image->width*y] += water_diffraction(q,
			                        ph_lambda_to_en(image->lambda),
			                        BEAM_RADIUS, WATER_RADIUS) * sw;

		}


	}
	progress_bar(xs, SAMPLING*image->width-1, "Calculating diffraction");
	}
}
