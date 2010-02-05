/*
 * ewald.c
 *
 * Calculate q-vector arrays
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "image.h"
#include "utils.h"
#include "cell.h"
#include "ewald.h"
#include "detector.h"


#define SAMPLING (4)
#define BWSAMPLING (10)
#define BANDWIDTH (0.05)


static struct rvec quat_rot(struct rvec q, struct quaternion z)
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


static void add_sphere(struct image *image, double k, int soffs)
{
	int x, y;

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double rx = 0.0;
		double ry = 0.0;
		double r;
		double twothetax, twothetay, twotheta;
		double qx, qy, qz;
		struct rvec q1, q2, q3, q4;
		int p, sx, sy, i;

		/* Calculate q vectors for Ewald sphere */
		for ( p=0; p<image->det.n_panels; p++ ) {
			if ( (x >= image->det.panels[p].min_x)
			  && (x <= image->det.panels[p].max_x)
			  && (y >= image->det.panels[p].min_y)
			  && (y <= image->det.panels[p].max_y) ) {
				rx = ((double)x - image->det.panels[p].cx)
				                            / image->resolution;
				ry = ((double)y - image->det.panels[p].cy)
				                            / image->resolution;
				break;
			}
		}

		/* Bottom left corner */
		r = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
		twothetax = atan2(rx, image->camera_len);
		twothetay = atan2(ry, image->camera_len);
		twotheta = atan2(r, image->camera_len);
		qx = k * sin(twothetax);
		qy = k * sin(twothetay);
		qz = k - k * cos(twotheta);
		q1.u = qx;  q1.v = qy;  q1.w = qz;
		/* 2theta value is calculated at the bottom left only */
		image->twotheta[x + image->width*y] = twotheta;

		/* Bottom right corner (using the same panel configuration!) */
		rx = ((double)(x+1) - image->det.panels[p].cx)
		                            / image->resolution;
		ry = ((double)y - image->det.panels[p].cy)
		                            / image->resolution;
		twothetax = atan2(rx, image->camera_len);
		twothetay = atan2(ry, image->camera_len);
		twotheta = atan2(r, image->camera_len);
		qx = k * sin(twothetax);
		qy = k * sin(twothetay);
		qz = k - k * cos(twotheta);
		q2.u = qx;  q2.v = qy;  q2.w = qz;

		/* Top left corner (using the same panel configuration!) */
		rx = ((double)x - image->det.panels[p].cx)
		                            / image->resolution;
		ry = ((double)(y+1) - image->det.panels[p].cy)
		                            / image->resolution;
		twothetax = atan2(rx, image->camera_len);
		twothetay = atan2(ry, image->camera_len);
		twotheta = atan2(r, image->camera_len);
		qx = k * sin(twothetax);
		qy = k * sin(twothetay);
		qz = k - k * cos(twotheta);
		q3.u = qx;  q3.v = qy;  q3.w = qz;

		/* Top right corner (using the same panel configuration!) */
		rx = ((double)(x+1) - image->det.panels[p].cx)
		                            / image->resolution;
		ry = ((double)(y+1) - image->det.panels[p].cy)
		                            / image->resolution;
		twothetax = atan2(rx, image->camera_len);
		twothetay = atan2(ry, image->camera_len);
		twotheta = atan2(r, image->camera_len);
		qx = k * sin(twothetax);
		qy = k * sin(twothetay);
		qz = k - k * cos(twotheta);
		q4.u = qx;  q4.v = qy;  q4.w = qz;

		/* Now interpolate between the values to get
		 * the sampling points */
		i = soffs;
		for ( sx=0; sx<SAMPLING; sx++ ) {
		for ( sy=0; sy<SAMPLING; sy++ ) {

			struct rvec q;

			q.u = q1.u + ((q2.u - q1.u)/SAMPLING)*sx
			           + ((q3.u - q1.u)/SAMPLING)*sy;
			q.v = q1.v + ((q2.v - q1.v)/SAMPLING)*sx
			           + ((q3.v - q1.v)/SAMPLING)*sy;
			q.w = q1.w + ((q2.w - q1.w)/SAMPLING)*sx
			           + ((q3.w - q1.w)/SAMPLING)*sy;
			image->qvecs[i++][x + image->width*y] = quat_rot(q,
		                                            image->orientation);

		}
		}

		if ( (x==0) && (y==(int)image->y_centre) ) {
			double s;
			s = 1.0e-9*modulus(qx, qy, qz)/2.0;
			STATUS("At left edge: 2theta = %5.3f deg,"
			       " sin(theta)/lambda = %5.3f nm^-1,"
			       " d = %5.3f nm\n",
			       rad2deg(twotheta), s, 1.0/(2.0*s));
		}
		if ( (x==0) && (y==0) ) {
			double s;
			s = 1.0e-9*modulus(qx, qy, qz)/2.0;
			STATUS("   At corner: 2theta = %5.3f deg,"
			       " sin(theta)/lambda = %5.3f nm^-1,"
			       " d = %5.3f nm\n",
			       rad2deg(twotheta), s, 1.0/(2.0*s));
		}

	}
	}

}


void get_ewald(struct image *image)
{
	double kc;  /* Wavenumber */
	int i, kstep;
	int mtotal = 0;

	kc = 1/image->lambda;  /* Centre */

	image->twotheta = malloc(image->width * image->height
	                       * sizeof(double));

	/* Create the spheres */
	image->nspheres = SAMPLING*SAMPLING*BWSAMPLING;
	image->qvecs = malloc(image->nspheres * sizeof(struct rvec *));

	for ( i=0; i<image->nspheres; i++ ) {
		mtotal += image->width * image->height * sizeof(struct rvec);
		image->qvecs[i] = malloc(image->width * image->height
                                                 * sizeof(struct rvec));
	}
	STATUS("%i spheres, %i Mbytes\n", image->nspheres, mtotal/(1024*1024));

	for ( kstep=0; kstep<BWSAMPLING; kstep++ ) {

		double k;

		k = kc + (kstep-(BWSAMPLING/2))*kc*(BANDWIDTH/BWSAMPLING);
		add_sphere(image, k, kstep*SAMPLING);

	}
}
