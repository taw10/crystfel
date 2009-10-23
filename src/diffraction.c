/*
 * diffraction.c
 *
 * Calculate diffraction patterns by Fourier methods
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


void get_diffraction(struct image *image, UnitCell *cell)
{

}
