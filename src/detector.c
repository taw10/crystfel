/*
 * detector.c
 *
 * Detector properties
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

#include "image.h"
#include "utils.h"
#include "diffraction.h"
#include "detector.h"
#include "parameters-lcls.tmp"


int atob(const char *a)
{
	if ( strcasecmp(a, "true") == 0 ) return 1;
	if ( strcasecmp(a, "false") == 0 ) return 1;
	return atoi(a);
}


/* x,y in pixels relative to image origin */
int map_position(struct image *image, double dx, double dy,
                 double *rx, double *ry, double *rz)
{
	double d;
	double twotheta, psi;
	const double k = 1.0 / image->lambda;
	struct panel *p;
	double x = 0.0;
	double y = 0.0;

	p = find_panel(image->det, dx, dy);
	if ( p == NULL ) return 1;
	if ( p->no_index ) return 1;

	x = ((double)dx - p->cx);
	y = ((double)dy - p->cy);

	/* Convert pixels to metres */
	x /= p->res;
	y /= p->res;	/* Convert pixels to metres */
	d = sqrt((x*x) + (y*y));
	twotheta = atan2(d, p->clen);

	psi = atan2(y, x);

	*rx = k*sin(twotheta)*cos(psi);
	*ry = k*sin(twotheta)*sin(psi);
	*rz = k - k*cos(twotheta);

	return 0;
}


void record_image(struct image *image, int do_poisson)
{
	int x, y;
	double total_energy, energy_density;
	double ph_per_e;
	double area;
	double max_tt = 0.0;

	/* How many photons are scattered per electron? */
	area = M_PI*pow(BEAM_RADIUS, 2.0);
	total_energy = FLUENCE * ph_lambda_to_en(image->lambda);
	energy_density = total_energy / area;
	ph_per_e = (FLUENCE/area) * pow(THOMSON_LENGTH, 2.0);
	STATUS("Fluence = %8.2e photons, "
	       "Energy density = %5.3f kJ/cm^2, "
	       "Total energy = %5.3f microJ\n",
	       FLUENCE, energy_density/1e7, total_energy*1e6);

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double counts;
		double cf;
		double intensity, sa;
		double pix_area, Lsq;
		double dsq, proj_area;
		struct panel *p;

		intensity = (double)image->data[x + image->width*y];
		if ( isinf(intensity) ) {
			ERROR("Infinity at %i,%i\n", x, y);
		}
		if ( intensity < 0.0 ) {
			ERROR("Negative at %i,%i\n", x, y);
		}
		if ( isnan(intensity) ) {
			ERROR("NaN at %i,%i\n", x, y);
		}

		p = find_panel(image->det, x, y);

		/* Area of one pixel */
		pix_area = pow(1.0/p->res, 2.0);
		Lsq = pow(p->clen, 2.0);

		/* Area of pixel as seen from crystal (approximate) */
		proj_area = pix_area * cos(image->twotheta[x + image->width*y]);

		/* Calculate distance from crystal to pixel */
		dsq = pow(((double)x - p->cx) / p->res, 2.0);
		dsq += pow(((double)y - p->cy) / p->res, 2.0);

		/* Projected area of pixel divided by distance squared */
		sa = proj_area / (dsq + Lsq);

		if ( do_poisson ) {
			counts = poisson_noise(intensity * ph_per_e * sa * DQE);
		} else {
			cf = intensity * ph_per_e * sa * DQE;
			counts = cf;
		}

		image->data[x + image->width*y] = counts * DETECTOR_GAIN;
		if ( isinf(image->data[x+image->width*y]) ) {
			ERROR("Processed infinity at %i,%i\n", x, y);
		}
		if ( isnan(image->data[x+image->width*y]) ) {
			ERROR("Processed NaN at %i,%i\n", x, y);
		}
		if ( image->data[x+image->width*y] < 0.0 ) {
			ERROR("Processed negative at %i,%i %f\n", x, y, counts);
		}

		if ( image->twotheta[x + image->width*y] > max_tt ) {
			max_tt = image->twotheta[x + image->width*y];
		}

	}
	progress_bar(x, image->width-1, "Post-processing");
	}

	STATUS("Max 2theta = %.2f deg, min d = %.2f nm\n",
	        rad2deg(max_tt), (image->lambda/(2.0*sin(max_tt/2.0)))/1e-9);

	double tt_side = image->twotheta[512+image->width*0];
	STATUS("At 512,0: %.2f deg, min d = %.2f nm\n",
	        rad2deg(tt_side), (image->lambda/(2.0*sin(tt_side/2.0)))/1e-9);

	tt_side = image->twotheta[0+image->width*512];
	STATUS("At 0,512: %.2f deg, min d = %.2f nm\n",
	        rad2deg(tt_side), (image->lambda/(2.0*sin(tt_side/2.0)))/1e-9);

	STATUS("Halve the d values to get the voxel size for a synthesis.\n");
}


struct panel *find_panel(struct detector *det, int x, int y)
{
	int p;

	for ( p=0; p<det->n_panels; p++ ) {
		if ( (x >= det->panels[p].min_x)
		  && (x <= det->panels[p].max_x)
		  && (y >= det->panels[p].min_y)
		  && (y <= det->panels[p].max_y) ) {
			return &det->panels[p];
		}
	}
	ERROR("No mapping found for %i,%i\n", x, y);

	return NULL;
}


struct detector *get_detector_geometry(const char *filename)
{
	FILE *fh;
	struct detector *det;
	char *rval;
	char **bits;
	int i;

	fh = fopen(filename, "r");
	if ( fh == NULL ) return NULL;

	det = malloc(sizeof(struct detector));
	if ( det == NULL ) {
		fclose(fh);
		return NULL;
	}
	det->n_panels = -1;

	do {

		int n1, n2;
		char **path;
		char line[1024];
		int np;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) break;
		chomp(line);

		n1 = assplode(line, " \t", &bits, ASSPLODE_NONE);
		if ( n1 < 3 ) {
			for ( i=0; i<n1; i++ ) free(bits[i]);
			free(bits);
			continue;
		}

		if ( bits[1][0] != '=' ) {
			for ( i=0; i<n1; i++ ) free(bits[i]);
			free(bits);
			continue;
		}

		if ( strcmp(bits[0], "n_panels") == 0 ) {
			if ( det->n_panels != -1 ) {
				ERROR("Duplicate n_panels statement.\n");
				fclose(fh);
				free(det);
				for ( i=0; i<n1; i++ ) free(bits[i]);
				free(bits);
				return NULL;
			}
			det->n_panels = atoi(bits[2]);
			det->panels = malloc(det->n_panels
			                      * sizeof(struct panel));
			for ( i=0; i<n1; i++ ) free(bits[i]);
			free(bits);
			continue;
		}

		n2 = assplode(bits[0], "/\\.", &path, ASSPLODE_NONE);
		if ( n2 < 2 ) {
			/* This was a top-level option, but not handled above. */
			for ( i=0; i<n1; i++ ) free(bits[i]);
			free(bits);
			for ( i=0; i<n2; i++ ) free(path[i]);
			free(path);
			continue;
		}

		np = atoi(path[0]);
		if ( det->n_panels == -1 ) {
			ERROR("n_panels statement must come first in "
			      "detector geometry file.\n");
			return NULL;
		}

		if ( np > det->n_panels ) {
			ERROR("The detector geometry file said there were %i "
			      "panels, but then tried to specify number %i\n",
			      det->n_panels, np);
			ERROR("Note: panel indices are counted from zero.\n");
			return NULL;
		}

		if ( strcmp(path[1], "min_x") == 0 ) {
			det->panels[np].min_x = atof(bits[2]);
		} else if ( strcmp(path[1], "max_x") == 0 ) {
			det->panels[np].max_x = atof(bits[2]);
		} else if ( strcmp(path[1], "min_y") == 0 ) {
			det->panels[np].min_y = atof(bits[2]);
		} else if ( strcmp(path[1], "max_y") == 0 ) {
			det->panels[np].max_y = atof(bits[2]);
		} else if ( strcmp(path[1], "cx") == 0 ) {
			det->panels[np].cx = atof(bits[2]);
		} else if ( strcmp(path[1], "cy") == 0 ) {
			det->panels[np].cy = atof(bits[2]);
		} else if ( strcmp(path[1], "clen") == 0 ) {
			det->panels[np].clen = atof(bits[2]);
		} else if ( strcmp(path[1], "res") == 0 ) {
			det->panels[np].res = atof(bits[2]);
		} else if ( strcmp(path[1], "badrow_direction") == 0 ) {
			det->panels[np].badrow = bits[2][0];
			if ( (det->panels[np].badrow != 'x')
			  && (det->panels[np].badrow != 'y') ) {
				ERROR("badrow_direction must be 'x' or 'y'\n");
				ERROR("Assuming 'x'\n.");
				det->panels[np].badrow = 'x';
			}
		} else if ( strcmp(path[1], "no_index") == 0 ) {
			det->panels[np].no_index = atob(bits[2]);
		} else {
			ERROR("Unrecognised field '%s'\n", path[1]);
		}

		for ( i=0; i<n1; i++ ) free(bits[i]);
		for ( i=0; i<n2; i++ ) free(path[i]);
		free(bits);
		free(path);

	} while ( rval != NULL );

	if ( det->n_panels == -1 ) {
		ERROR("No panel descriptions in geometry file.\n");
		fclose(fh);
		free(det->panels);
		free(det);
		return NULL;
	}

	for ( i=0; i<det->n_panels; i++ ) {
		STATUS("Panel %i, min_x = %i\n", i, det->panels[i].min_x);
		STATUS("Panel %i, max_x = %i\n", i, det->panels[i].max_x);
		STATUS("Panel %i, min_y = %i\n", i, det->panels[i].min_y);
		STATUS("Panel %i, max_y = %i\n", i, det->panels[i].max_y);
		STATUS("Panel %i, cx = %f\n", i, det->panels[i].cx);
		STATUS("Panel %i, cy = %f\n", i, det->panels[i].cy);
		STATUS("Panel %i, clen = %f\n", i, det->panels[i].clen);
		STATUS("Panel %i, res = %f\n", i, det->panels[i].res);
	}

	fclose(fh);

	return det;
}
