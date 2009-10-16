/*
 * relrod.c
 *
 * Calculate reflection positions via line-sphere intersection test
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * template_index - Indexing diffraction patterns by template matching
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "image.h"
#include "utils.h"
#include "cell.h"


static void mapping_rotate(double x, double y, double z,
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


void get_reflections(struct image *image, UnitCell *cell, double smax)
{
	ImageFeatureList *flist;
	double tilt, omega, wavenumber;
	double nx, ny, nz; /* "normal" vector */
	double kx, ky, kz; /* Electron wavevector ("normal" times 1/lambda) */
	double ux, uy, uz; /* "up" vector */
	double rx, ry, rz; /* "right" vector */
	double asx, asy, asz; /* "a*" lattice parameter */
	double bsx, bsy, bsz; /* "b*" lattice parameter */
	double csx, csy, csz; /* "c*" lattice parameter */
	signed int h, k, l;
	double res_max;

	/* Get the reciprocal unit cell */
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	/* Prepare list and some parameters */
	flist = image_feature_list_new();
	tilt = image->tilt;
	omega = image->omega;
	wavenumber = 1.0/image->lambda;

	/* Calculate (roughly) the maximum resolution */
	if ( image->fmode == FORMULATION_CLEN ) {
		double w2, h2;
		w2 = image->width/2;  h2 = image->height/2;
		w2 = pow(w2, 2.0);  h2 = pow(h2, 2.0);
		res_max = sqrt(w2 + h2) / image->resolution;
		res_max *= (wavenumber / image->camera_len);
	} else {
		fprintf(stderr,
			"Unrecognised formulation mode in get_reflections"
			" (resolution cutoff calculation)\n");
		return;
	}
	printf("Resolution cutoff is %5.2f nm^-1\n", res_max/1e9);
	res_max = pow(res_max, 2.0);

	/* Calculate the (normalised) incident electron wavevector */
	mapping_rotate(0.0, 0.0, 1.0, &nx, &ny, &nz, omega, tilt);
	kx = nx / image->lambda;
	ky = ny / image->lambda;
	kz = nz / image->lambda;	/* This is the centre of the Ewald sphere */

	/* Determine where "up" is */
	mapping_rotate(0.0, 1.0, 0.0, &ux, &uy, &uz, omega, tilt);

	/* Determine where "right" is */
	mapping_rotate(1.0, 0.0, 0.0, &rx, &ry, &rz, omega, tilt);

	for ( h=-50; h<50; h++ ) {
	for ( k=-50; k<50; k++ ) {
	for ( l=-50; l<50; l++ ) {

		double xl, yl, zl;
		double a, b, c;
		double s1, s2, s, t;
		double g_sq, gn;

		/* Get the coordinates of the reciprocal lattice point */
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;
		zl = h*asz + k*bsz + l*csz;
		g_sq = modulus_squared(xl, yl, zl);
		gn = xl*nx + yl*ny + zl*nz;

		/* Early bailout if resolution is clearly too high */
		if ( g_sq > res_max ) continue;

		/* Next, solve the relrod equation to calculate
		 * the excitation error */
		a = 1.0;
		b = 2.0*(wavenumber + gn);
		c = -2.0*gn*wavenumber + g_sq;
		t = -0.5*(b + sign(b)*sqrt(b*b - 4.0*a*c));
		s1 = t/a;
		s2 = c/t;
		if ( fabs(s1) < fabs(s2) ) s = s1; else s = s2;

		/* Skip this reflection if s is large */
		if ( fabs(s) <= smax ) {

			double xi, yi, zi;
			double gx, gy, gz;
			double theta;
			double x, y;
			double dx, dy, psi;

			/* Determine the intersection point */
			xi = xl + s*nx;  yi = yl + s*ny;  zi = zl + s*nz;

			/* Calculate Bragg angle */
			gx = xi - kx;
			gy = yi - ky;
			gz = zi - kz;	/* This is the vector from the centre of
					*  the sphere to the intersection */
			theta = angle_between(-kx, -ky, -kz, gx, gy, gz);

			/* Calculate azimuth of point in image
			 * (anticlockwise from +x) */
			dx = xi*rx + yi*ry + zi*rz;
			dy = xi*ux + yi*uy + zi*uz;
			psi = atan2(dy, dx);

			/* Get image coordinates from polar
			 * representation */
			if ( image->fmode == FORMULATION_CLEN ) {
				x = image->camera_len*tan(theta)*cos(psi);
				y = image->camera_len*tan(theta)*sin(psi);
				x *= image->resolution;
				y *= image->resolution;
			} else if ( image->fmode==FORMULATION_PIXELSIZE ) {
				x = tan(theta)*cos(psi) / image->lambda;
				y = tan(theta)*sin(psi) / image->lambda;
				x /= image->pixel_size;
				y /= image->pixel_size;
			} else {
				fprintf(stderr,
					"Unrecognised formulation mode "
					"in get_reflections\n");
				return;
			}

			x += image->x_centre;
			y += image->y_centre;

			/* Sanity check */
			if ( (x>=0) && (x<image->width)
			  && (y>=0) && (y<image->height) ) {

				/* Record the reflection.
				 * Intensity should be multiplied by relrod
				 * spike function, except reprojected
				 * reflections aren't used quantitatively for
				 * anything. */
				image_add_feature(flist, x, y, image, 1.0);

			} /* else it's outside the picture somewhere */

		} /* else reflection is not excited in this orientation */

	}
	}
	}

	image->rflist = flist;
}
