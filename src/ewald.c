/*
 * ewald.c
 *
 * Calculate q-vector arrays
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "image.h"
#include "utils.h"
#include "cell.h"
#include "ewald.h"


static struct threevec quat_rot(struct threevec q, struct quaternion z)
{
	struct threevec res;
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


void get_ewald(struct image *image)
{
	int x, y;
	double k;  /* Wavenumber */

	k = 1/image->lambda;

	image->qvecs = malloc(image->width * image->height
	                       * sizeof(struct threevec));

	image->twotheta = malloc(image->width * image->height
	                       * sizeof(double));

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double rx, ry, r;
		double twothetax, twothetay, twotheta;
		double qx, qy, qz;
		struct threevec q;

		/* Calculate q vectors for Ewald sphere */
		rx = ((double)x - image->x_centre) / image->resolution;
		ry = ((double)y - image->y_centre) / image->resolution;
		r = sqrt(pow(rx, 2.0) + pow(ry, 2.0));

		twothetax = atan2(rx, image->camera_len);
		twothetay = atan2(ry, image->camera_len);
		twotheta = atan2(r, image->camera_len);

		qx = k * sin(twothetax);
		qy = k * sin(twothetay);
		qz = k - k * cos(twotheta);

		q.u = qx;  q.v = qy;  q.w = qz;
		image->qvecs[x + image->width*y] = quat_rot(q,
		                                            image->orientation);

		image->twotheta[x + image->width*y] = twotheta;

		if ( (x==0) && (y==(int)image->y_centre) ) {
			double s;
			s = 1.0e-9*modulus(qx, qy, qz)/2.0;
			printf("At left edge: 2theta = %5.3f deg,"
			       " sin(theta)/lambda = %5.3f nm^-1,"
			       " d = %f nm\n",
			       rad2deg(twotheta), s, 1.0/(2.0*s));
		}
		if ( (x==0) && (y==0) ) {
			double s;
			s = 1.0e-9*modulus(qx, qy, qz)/2.0;
			printf("At corner: 2theta = %5.3f deg,"
			       " sin(theta)/lambda = %5.3f nm^-1,"
			       " d = %f nm\n",
			       rad2deg(twotheta), s, 1.0/(2.0*s));
		}

	}
	}
}
