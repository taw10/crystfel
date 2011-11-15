/*
 * detector.c
 *
 * Detector properties
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 * (c) 2011 Rick Kirian <rkirian@asu.edu>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#define _ISOC99_SOURCE
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "image.h"
#include "utils.h"
#include "detector.h"
#include "beam-parameters.h"
#include "hdf5-file.h"


/**
 * SECTION:detector
 * @short_description: Detector geometry
 * @title: Detector
 * @section_id:
 * @see_also:
 * @include: "detector.h"
 * @Image:
 *
 * This structure represents the detector geometry
 */


static int atob(const char *a)
{
	if ( strcasecmp(a, "true") == 0 ) return 1;
	if ( strcasecmp(a, "false") == 0 ) return 0;
	return atoi(a);
}


static int assplode_algebraic(const char *a_orig, char ***pbits)
{
	int len, i;
	int nexp;
	char **bits;
	char *a;
	int idx, istr;

	len = strlen(a_orig);

	/* Add plus at start if no sign there already */
	if ( (a_orig[0] != '+') && (a_orig[0] != '-') ) {
		len += 1;
		a = malloc(len+1);
		snprintf(a, len+1, "+%s", a_orig);
		a[len] = '\0';

	} else {
		a = strdup(a_orig);
	}

	/* Count the expressions */
	nexp = 0;
	for ( i=0; i<len; i++ ) {
		if ( (a[i] == '+') || (a[i] == '-') ) nexp++;
	}

	bits = calloc(nexp, sizeof(char *));

	/* Break the string up */
	idx = -1;
	istr = 0;
	assert((a[0] == '+') || (a[0] == '-'));
	for ( i=0; i<len; i++ ) {

		char ch;

		ch = a[i];

		if ( (ch == '+') || (ch == '-') ) {
			if ( idx >= 0 ) bits[idx][istr] = '\0';
			idx++;
			bits[idx] = malloc(len+1);
			istr = 0;
		}

		if ( !isdigit(ch) && (ch != '.') && (ch != 'x') && (ch != 'y')
		  && (ch != '+') && (ch != '-') )
		{
			ERROR("Invalid character '%C' found.\n", ch);
			return 0;
		}

		assert(idx >= 0);
		bits[idx][istr++] = ch;

	}
	if ( idx >= 0 ) bits[idx][istr] = '\0';

	*pbits = bits;
	free(a);

	return nexp;
}


/* Parses the scan directions (accounting for possible rotation)
 * Assumes all white spaces have been already removed */
static int dir_conv(const char *a, double *sx, double *sy)
{
	int n;
	char **bits;
	int i;

	*sx = 0.0;  *sy = 0.0;

	n = assplode_algebraic(a, &bits);

	if ( n == 0 ) {
		ERROR("Invalid direction '%s'\n", a);
		return 1;
	}

	for ( i=0; i<n; i++ ) {

		int len;
		double val;
		char axis;
		int j;

		len = strlen(bits[i]);
		assert(len != 0);
		axis = bits[i][len-1];
		if ( (axis != 'x') && (axis != 'y') ) {
			ERROR("Invalid symbol '%C' - must be x or y.\n", axis);
			return 1;
		}

		/* Chop off the symbol now it's dealt with */
		bits[i][len-1] = '\0';

		/* Check for anything that isn't part of a number */
		for ( j=0; j<strlen(bits[i]); j++ ) {
			if ( isdigit(bits[i][j]) ) continue;
			if ( bits[i][j] == '+' ) continue;
			if ( bits[i][j] == '-' ) continue;
			if ( bits[i][j] == '.' ) continue;
			ERROR("Invalid coefficient '%s'\n", bits[i]);
		}

		if ( strlen(bits[i]) == 0 ) {
			val = 1.0;
		} else {
			val = atof(bits[i]);
		}
		if ( strlen(bits[i]) == 1 ) {
			if ( bits[i][0] == '+' ) val = 1.0;
			if ( bits[i][0] == '-' ) val = -1.0;
		}
		if ( axis == 'x' ) {
			*sx += val;
		} else if ( axis == 'y' ) {
			*sy += val;
		}

		free(bits[i]);

	}
	free(bits);

	//STATUS("'%s' -> %5.2fx + %5.2fy\n", a, *sx, *sy);

	return 0;
}


struct rvec get_q_for_panel(struct panel *p, double fs, double ss,
                            double *ttp, double k)
{
	struct rvec q;
	double twotheta, r, az;
	double rx, ry;
	double xs, ys;

	/* Convert xs and ys, which are in fast scan/slow scan coordinates,
	 * to x and y */
	xs = fs*p->fsx + ss*p->ssx;
	ys = fs*p->fsy + ss*p->ssy;

	rx = (xs + p->cnx) / p->res;
	ry = (ys + p->cny) / p->res;

	/* Calculate q-vector for this sub-pixel */
	r = sqrt(pow(rx, 2.0) + pow(ry, 2.0));

	twotheta = atan2(r, p->clen);
	az = atan2(ry, rx);
	if ( ttp != NULL ) *ttp = twotheta;

	q.u = k * sin(twotheta)*cos(az);
	q.v = k * sin(twotheta)*sin(az);
	q.w = k * (cos(twotheta) - 1.0);

	return q;
}


struct rvec get_q(struct image *image, double fs, double ss,
                  double *ttp, double k)
{
	struct panel *p;
	const unsigned int fsi = fs;
	const unsigned int ssi = ss;  /* Explicit rounding */

	/* Determine which panel to use */
	p = find_panel(image->det, fsi, ssi);
	assert(p != NULL);

	return get_q_for_panel(p, fs-(double)p->min_fs, ss-(double)p->min_ss,
	                       ttp, k);
}


int in_bad_region(struct detector *det, double fs, double ss)
{
	double rx, ry;
	struct panel *p;
	double xs, ys;
	int i;

	/* Determine which panel to use */
	p = find_panel(det, fs, ss);

	/* No panel found -> definitely bad! */
	if ( p == NULL ) return 1;

	/* Convert xs and ys, which are in fast scan/slow scan coordinates,
	 * to x and y */
	xs = (fs-(double)p->min_fs)*p->fsx + (ss-(double)p->min_ss)*p->ssx;
	ys = (fs-(double)p->min_fs)*p->fsy + (ss-(double)p->min_ss)*p->ssy;

	rx = xs + p->cnx;
	ry = ys + p->cny;

	for ( i=0; i<det->n_bad; i++ ) {
		struct badregion *b = &det->bad[i];
		if ( rx < b->min_x ) continue;
		if ( rx > b->max_x ) continue;
		if ( ry < b->min_y ) continue;
		if ( ry > b->max_y ) continue;
		return 1;
	}

	return 0;
}


double get_tt(struct image *image, double fs, double ss)
{
	double r, rx, ry;
	struct panel *p;
	double xs, ys;

	p = find_panel(image->det, fs, ss);

	/* Convert xs and ys, which are in fast scan/slow scan coordinates,
	 * to x and y */
	xs = (fs-p->min_fs)*p->fsx + (ss-p->min_ss)*p->ssx;
	ys = (fs-p->min_fs)*p->fsy + (ss-p->min_ss)*p->ssy;

	rx = (xs + p->cnx) / p->res;
	ry = (ys + p->cny) / p->res;

	r = sqrt(pow(rx, 2.0) + pow(ry, 2.0));

	return atan2(r, p->clen);
}


void record_image(struct image *image, int do_poisson)
{
	int x, y;
	double total_energy, energy_density;
	double ph_per_e;
	double area;
	double max_tt = 0.0;
	int n_inf1 = 0;
	int n_neg1 = 0;
	int n_nan1 = 0;
	int n_inf2 = 0;
	int n_neg2 = 0;
	int n_nan2 = 0;

	/* How many photons are scattered per electron? */
	area = M_PI*pow(image->beam->beam_radius, 2.0);
	total_energy = image->beam->fluence * ph_lambda_to_en(image->lambda);
	energy_density = total_energy / area;
	ph_per_e = (image->beam->fluence /area) * pow(THOMSON_LENGTH, 2.0);
	STATUS("Fluence = %8.2e photons, "
	       "Energy density = %5.3f kJ/cm^2, "
	       "Total energy = %5.3f microJ\n",
	       image->beam->fluence, energy_density/1e7, total_energy*1e6);

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double counts;
		double cf;
		double intensity, sa;
		double pix_area, Lsq;
		double xs, ys, rx, ry;
		double dsq, proj_area;
		struct panel *p;

		intensity = (double)image->data[x + image->width*y];
		if ( isinf(intensity) ) n_inf1++;
		if ( intensity < 0.0 ) n_neg1++;
		if ( isnan(intensity) ) n_nan1++;

		p = find_panel(image->det, x, y);

		/* Area of one pixel */
		pix_area = pow(1.0/p->res, 2.0);
		Lsq = pow(p->clen, 2.0);

		/* Area of pixel as seen from crystal (approximate) */
		proj_area = pix_area * cos(image->twotheta[x + image->width*y]);

		/* Calculate distance from crystal to pixel */
		xs = (x-p->min_fs)*p->fsx + (y-p->min_ss)*p->ssx;
		ys = (x-p->min_fs)*p->fsy + (y-p->min_ss)*p->ssy;
		rx = (xs + p->cnx) / p->res;
		ry = (ys + p->cny) / p->res;
		dsq = sqrt(pow(rx, 2.0) + pow(ry, 2.0));

		/* Projected area of pixel divided by distance squared */
		sa = proj_area / (dsq + Lsq);

		if ( do_poisson ) {
			counts = poisson_noise(intensity * ph_per_e
			                              * sa * image->beam->dqe );
		} else {
			cf = intensity * ph_per_e * sa * image->beam->dqe;
			counts = cf;
		}

		image->data[x + image->width*y] = counts
		                                  * image->beam->adu_per_photon;

		/* Sanity checks */
		if ( isinf(image->data[x+image->width*y]) ) n_inf2++;
		if ( isnan(image->data[x+image->width*y]) ) n_nan2++;
		if ( image->data[x+image->width*y] < 0.0 ) n_neg2++;

		if ( image->twotheta[x + image->width*y] > max_tt ) {
			max_tt = image->twotheta[x + image->width*y];
		}

	}
	progress_bar(x, image->width-1, "Post-processing");
	}

	STATUS("Max 2theta = %.2f deg, min d = %.2f nm\n",
	        rad2deg(max_tt), (image->lambda/(2.0*sin(max_tt/2.0)))/1e-9);

	double tt_side = image->twotheta[(image->width/2)+image->width*0];
	STATUS("At middle of bottom edge: %.2f deg, min d = %.2f nm\n",
	        rad2deg(tt_side), (image->lambda/(2.0*sin(tt_side/2.0)))/1e-9);

	tt_side = image->twotheta[0+image->width*(image->height/2)];
	STATUS("At middle of left edge: %.2f deg, min d = %.2f nm\n",
	        rad2deg(tt_side), (image->lambda/(2.0*sin(tt_side/2.0)))/1e-9);

	STATUS("Halve the d values to get the voxel size for a synthesis.\n");

	if ( n_neg1 + n_inf1 + n_nan1 + n_neg2 + n_inf2 + n_nan2 ) {
		ERROR("WARNING: The raw calculation produced %i negative"
		      " values, %i infinities and %i NaNs.\n",
		      n_neg1, n_inf1, n_nan1);
		ERROR("WARNING: After processing, there were %i negative"
		      " values, %i infinities and %i NaNs.\n",
		      n_neg2, n_inf2, n_nan2);
	}
}


struct panel *find_panel(struct detector *det, double fs, double ss)
{
	int p;

	for ( p=0; p<det->n_panels; p++ ) {
		if ( (fs >= det->panels[p].min_fs)
		  && (fs <= det->panels[p].max_fs)
		  && (ss >= det->panels[p].min_ss)
		  && (ss <= det->panels[p].max_ss) ) {
			return &det->panels[p];
		}
	}

	return NULL;
}


int find_panel_number(struct detector *det, int fs, int ss)
{
	int p;

	for ( p=0; p<det->n_panels; p++ ) {
		if ( (fs >= det->panels[p].min_fs)
		  && (fs <= det->panels[p].max_fs)
		  && (ss >= det->panels[p].min_ss)
		  && (ss <= det->panels[p].max_ss) ) return p;
	}

	return -1;
}


void fill_in_values(struct detector *det, struct hdfile *f)
{
	int i;

	for ( i=0; i<det->n_panels; i++ ) {

		struct panel *p = &det->panels[i];

		if ( p->clen_from != NULL ) {
			p->clen = get_value(f, p->clen_from) * 1.0e-3;
			free(p->clen_from);
			p->clen_from = NULL;
		}

		p->clen += p->coffset;

	}
}


static struct panel *new_panel(struct detector *det, const char *name)
{
	struct panel *new;

	det->n_panels++;
	det->panels = realloc(det->panels, det->n_panels*sizeof(struct panel));

	new = &det->panels[det->n_panels-1];
	memcpy(new, &det->defaults, sizeof(struct panel));

	/* Create a new copy of the camera length location if needed */
	if ( new->clen_from != NULL ) {
		new->clen_from = strdup(new->clen_from);
	}
	strcpy(new->name, name);

	return new;
}


static struct badregion *new_bad_region(struct detector *det, const char *name)
{
	struct badregion *new;

	det->n_bad++;
	det->bad = realloc(det->bad, det->n_bad*sizeof(struct badregion));

	new = &det->bad[det->n_bad-1];
	new->min_x = NAN;
	new->max_x = NAN;
	new->min_y = NAN;
	new->max_y = NAN;
	strcpy(new->name, name);

	return new;
}


struct panel *find_panel_by_name(struct detector *det, const char *name)
{
	int i;

	for ( i=0; i<det->n_panels; i++ ) {
		if ( strcmp(det->panels[i].name, name) == 0 ) {
			return &det->panels[i];
		}
	}

	return NULL;
}


static struct badregion *find_bad_region_by_name(struct detector *det,
                                                 const char *name)
{
	int i;

	for ( i=0; i<det->n_bad; i++ ) {
		if ( strcmp(det->bad[i].name, name) == 0 ) {
			return &det->bad[i];
		}
	}

	return NULL;
}


static char *find_or_add_rg(struct detector *det, const char *name)
{
	int i;
	char **new;
	char *tmp;

	for ( i=0; i<det->num_rigid_groups; i++ ) {

		if ( strcmp(det->rigid_groups[i], name) == 0 ) {
			return det->rigid_groups[i];
		}

	}

	new = realloc(det->rigid_groups,
	              (1+det->num_rigid_groups)*sizeof(char *));
	if ( new == NULL ) return NULL;

	det->rigid_groups = new;

	tmp = strdup(name);
	det->rigid_groups[det->num_rigid_groups] = tmp;

	det->num_rigid_groups++;

	return tmp;
}


static int parse_field_for_panel(struct panel *panel, const char *key,
                                 const char *val, struct detector *det)
{
	int reject = 0;

	if ( strcmp(key, "min_fs") == 0 ) {
		panel->min_fs = atof(val);
	} else if ( strcmp(key, "max_fs") == 0 ) {
		panel->max_fs = atof(val);
	} else if ( strcmp(key, "min_ss") == 0 ) {
		panel->min_ss = atof(val);
	} else if ( strcmp(key, "max_ss") == 0 ) {
		panel->max_ss = atof(val);
	} else if ( strcmp(key, "corner_x") == 0 ) {
		panel->cnx = atof(val);
	} else if ( strcmp(key, "corner_y") == 0 ) {
		panel->cny = atof(val);
	} else if ( strcmp(key, "rigid_group") == 0 ) {
		panel->rigid_group = find_or_add_rg(det, val);
	} else if ( strcmp(key, "clen") == 0 ) {

		char *end;
		double v = strtod(val, &end);
		if ( end == val ) {
			/* This means "fill in later" */
			panel->clen = -1.0;
			panel->clen_from = strdup(val);
		} else {
			panel->clen = v;
			panel->clen_from = NULL;
		}

	} else if ( strcmp(key, "coffset") == 0) {
		panel->coffset = atof(val);
	} else if ( strcmp(key, "res") == 0 ) {
		panel->res = atof(val);
	} else if ( strcmp(key, "peak_sep") == 0 ) {
		panel->peak_sep = atof(val);
	} else if ( strcmp(key, "integr_radius") == 0 ) {
		panel->integr_radius = atof(val);
	} else if ( strcmp(key, "badrow_direction") == 0 ) {
		panel->badrow = val[0]; /* First character only */
		if ( (panel->badrow != 'x') && (panel->badrow != 'y')
		  && (panel->badrow != 'f') && (panel->badrow != 's')
		  && (panel->badrow != '-') ) {
			ERROR("badrow_direction must be x, y, f, s or '-'\n");
			ERROR("Assuming '-'\n.");
			panel->badrow = '-';
		}
		if ( panel->badrow == 'x' ) panel->badrow = 'f';
		if ( panel->badrow == 'y' ) panel->badrow = 's';
	} else if ( strcmp(key, "no_index") == 0 ) {
		panel->no_index = atob(val);
	} else if ( strcmp(key, "fs") == 0 ) {
		if ( dir_conv(val, &panel->fsx, &panel->fsy) != 0 ) {
			ERROR("Invalid fast scan direction '%s'\n", val);
			reject = 1;
		}
	} else if ( strcmp(key, "ss") == 0 ) {
		if ( dir_conv(val, &panel->ssx, &panel->ssy) != 0 ) {
			ERROR("Invalid slow scan direction '%s'\n", val);
			reject = 1;
		}
	} else {
		ERROR("Unrecognised field '%s'\n", key);
	}

	return reject;
}


static int parse_field_bad(struct badregion *panel, const char *key,
                                  const char *val)
{
	int reject = 0;

	if ( strcmp(key, "min_x") == 0 ) {
		panel->min_x = atof(val);
	} else if ( strcmp(key, "max_x") == 0 ) {
		panel->max_x = atof(val);
	} else if ( strcmp(key, "min_y") == 0 ) {
		panel->min_y = atof(val);
	} else if ( strcmp(key, "max_y") == 0 ) {
		panel->max_y = atof(val);
	} else {
		ERROR("Unrecognised field '%s'\n", key);
	}

	return reject;
}


static void parse_toplevel(struct detector *det, const char *key,
                           const char *val)
{
	if ( strcmp(key, "mask") == 0 ) {

		det->mask = strdup(val);

	} else if ( strcmp(key, "mask_bad") == 0 ) {

		char *end;
		double v = strtod(val, &end);

		if ( end != val ) {
			det->mask_bad = v;
		}

	} else if ( strcmp(key, "mask_good") == 0 ) {

		char *end;
		double v = strtod(val, &end);

		if ( end != val ) {
			det->mask_good = v;
		}

	} else if ( strcmp(key, "peak_sep") == 0 ) {
		det->defaults.peak_sep = atof(val);
	} else if ( strcmp(key, "integr_radius") == 0 ) {
		det->defaults.integr_radius = atof(val);
	} else if ( parse_field_for_panel(&det->defaults, key, val, det) ) {
		ERROR("Unrecognised top level field '%s'\n", key);
	}
}


struct detector *get_detector_geometry(const char *filename)
{
	FILE *fh;
	struct detector *det;
	char *rval;
	char **bits;
	int i;
	int reject = 0;
	int x, y, max_fs, max_ss;

	fh = fopen(filename, "r");
	if ( fh == NULL ) return NULL;

	det = calloc(1, sizeof(struct detector));
	if ( det == NULL ) {
		fclose(fh);
		return NULL;
	}

	det->n_panels = 0;
	det->panels = NULL;
	det->n_bad = 0;
	det->bad = NULL;
	det->mask_good = 0;
	det->mask_bad = 0;
	det->mask = NULL;
	det->num_rigid_groups = 0;
	det->rigid_groups = NULL;

	/* The default defaults... */
	det->defaults.min_fs = -1;
	det->defaults.min_ss = -1;
	det->defaults.max_fs = -1;
	det->defaults.max_ss = -1;
	det->defaults.cnx = NAN;
	det->defaults.cny = NAN;
	det->defaults.clen = -1.0;
	det->defaults.coffset = 0.0;
	det->defaults.res = -1.0;
	det->defaults.badrow = '-';
	det->defaults.no_index = 0;
	det->defaults.peak_sep = 50.0;
	det->defaults.integr_radius = 3.0;
	det->defaults.fsx = 1.0;
	det->defaults.fsy = 0.0;
	det->defaults.ssx = 0.0;
	det->defaults.ssy = 1.0;
	det->defaults.rigid_group = NULL;
	strncpy(det->defaults.name, "", 1023);

	do {

		int n1, n2;
		char **path;
		char line[1024];
		struct badregion *badregion = NULL;
		struct panel *panel = NULL;
		char *key;
		char wholeval[1024];

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) break;
		chomp(line);

		if ( line[0] == ';' ) continue;

		n1 = assplode(line, " \t", &bits, ASSPLODE_NONE);
		if ( n1 < 3 ) {
			for ( i=0; i<n1; i++ ) free(bits[i]);
			free(bits);
			continue;
		}

		/* Stitch the pieces of the "value" back together */
		wholeval[0] = '\0';  /* Empty string */
		for ( i=2; i<n1; i++ ) {
			if ( bits[i][0] == ';' ) break;  /* Stop on comment */
			strncat(wholeval, bits[i], 1023);
		}

		if ( bits[1][0] != '=' ) {
			for ( i=0; i<n1; i++ ) free(bits[i]);
			free(bits);
			continue;
		}

		n2 = assplode(bits[0], "/\\.", &path, ASSPLODE_NONE);
		if ( n2 < 2 ) {
			/* This was a top-level option, not handled above. */
			parse_toplevel(det, bits[0], bits[2]);
			for ( i=0; i<n1; i++ ) free(bits[i]);
			free(bits);
			for ( i=0; i<n2; i++ ) free(path[i]);
			free(path);
			continue;
		}

		if ( strncmp(path[0], "bad", 3) == 0 ) {
			badregion = find_bad_region_by_name(det, path[0]);
			if ( badregion == NULL ) {
				badregion = new_bad_region(det, path[0]);
			}
		} else {
			panel = find_panel_by_name(det, path[0]);
			if ( panel == NULL ) {
				panel = new_panel(det, path[0]);
			}
		}

		key = path[1];

		if ( panel != NULL ) {
			if ( parse_field_for_panel(panel, path[1],
			                           wholeval, det) )
			{
				reject = 1;
			}
		} else {
			if ( parse_field_bad(badregion, path[1], wholeval) ) {
				reject = 1;
			}
		}

		for ( i=0; i<n1; i++ ) free(bits[i]);
		for ( i=0; i<n2; i++ ) free(path[i]);
		free(bits);
		free(path);

	} while ( rval != NULL );

	if ( det->n_panels == -1 ) {
		ERROR("No panel descriptions in geometry file.\n");
		fclose(fh);
		free(det);
		return NULL;
	}

	max_fs = 0;
	max_ss = 0;
	for ( i=0; i<det->n_panels; i++ ) {

		if ( det->panels[i].min_fs < 0 ) {
			ERROR("Please specify the minimum FS coordinate for"
			      " panel %s\n", det->panels[i].name);
			reject = 1;
		}
		if ( det->panels[i].max_fs < 0 ) {
			ERROR("Please specify the maximum FS coordinate for"
			      " panel %s\n", det->panels[i].name);
			reject = 1;
		}
		if ( det->panels[i].min_ss < 0 ) {
			ERROR("Please specify the minimum SS coordinate for"
			      " panel %s\n", det->panels[i].name);
			reject = 1;
		}
		if ( det->panels[i].max_ss < 0 ) {
			ERROR("Please specify the maximum SS coordinate for"
			      " panel %s\n", det->panels[i].name);
			reject = 1;
		}
		if ( isnan(det->panels[i].cnx)  ) {
			ERROR("Please specify the corner X coordinate for"
			      " panel %s\n", det->panels[i].name);
			reject = 1;
		}
		if ( isnan(det->panels[i].cny) ) {
			ERROR("Please specify the corner Y coordinate for"
			      " panel %s\n", det->panels[i].name);
			reject = 1;
		}
		if ( (det->panels[i].clen < 0.0)
		  && (det->panels[i].clen_from == NULL) ) {
			ERROR("Please specify the camera length for"
			      " panel %s\n", det->panels[i].name);
			reject = 1;
		}
		if ( det->panels[i].res < 0 ) {
			ERROR("Please specify the resolution for"
			      " panel %s\n", det->panels[i].name);
			reject = 1;
		}
		/* It's OK if the badrow direction is '0' */
		/* It's not a problem if "no_index" is still zero */
		/* The default peak_sep is OK (maybe) */
		/* The default transformation matrix is at least valid */

		if ( det->panels[i].max_fs > max_fs ) {
			max_fs = det->panels[i].max_fs;
		}
		if ( det->panels[i].max_ss > max_ss ) {
			max_ss = det->panels[i].max_ss;
		}

	}

	for ( i=0; i<det->n_bad; i++ ) {

		if ( isnan(det->bad[i].min_x) ) {
			ERROR("Please specify the minimum x coordinate for"
			      " bad region %s\n", det->bad[i].name);
			reject = 1;
		}
		if ( isnan(det->bad[i].min_y) ) {
			ERROR("Please specify the minimum y coordinate for"
			      " bad region %s\n", det->bad[i].name);
			reject = 1;
		}
		if ( isnan(det->bad[i].max_x) ) {
			ERROR("Please specify the maximum x coordinate for"
			      " bad region %s\n", det->bad[i].name);
			reject = 1;
		}
		if ( isnan(det->bad[i].max_y) ) {
			ERROR("Please specify the maximum y coordinate for"
			      " bad region %s\n", det->bad[i].name);
			reject = 1;
		}
	}

	for ( x=0; x<=max_fs; x++ ) {
	for ( y=0; y<=max_ss; y++ ) {
		if ( find_panel(det, x, y) == NULL ) {
			ERROR("Detector geometry invalid: contains gaps.\n");
			reject = 1;
			goto out;
		}
	}
	}

out:
	det->max_fs = max_fs;
	det->max_ss = max_ss;

	/* Calculate matrix inverse */
	for ( i=0; i<det->n_panels; i++ ) {

		struct panel *p;
		double d;

		p = &det->panels[i];

		if ( p->fsx*p->ssy == p->ssx*p->fsy ) {
			ERROR("Panel %i transformation singular.\n", i);
			reject = 1;
		}

		d = (double)p->fsx*p->ssy - p->ssx*p->fsy;
		p->xfs = p->ssy / d;
		p->yfs = -p->ssx / d;
		p->xss = -p->fsy / d;
		p->yss = p->fsx / d;

	}

	if ( reject ) return NULL;

	fclose(fh);

	return det;
}


void free_detector_geometry(struct detector *det)
{
	int i;

	for ( i=0; i<det->num_rigid_groups; i++ ) {
		free(det->rigid_groups[i]);
	}
	free(det->rigid_groups);

	free(det->panels);
	free(det->bad);
	free(det->mask);
	free(det);
}


struct detector *copy_geom(const struct detector *in)
{
	struct detector *out;
	int i;

	out = malloc(sizeof(struct detector));
	memcpy(out, in, sizeof(struct detector));

	if ( in->mask != NULL ) {
		out->mask = strdup(in->mask);
	} else {
		out->mask = NULL;  /* = in->mask */
	}

	out->panels = malloc(out->n_panels * sizeof(struct panel));
	memcpy(out->panels, in->panels, out->n_panels * sizeof(struct panel));

	out->bad = malloc(out->n_bad * sizeof(struct badregion));
	memcpy(out->bad, in->bad, out->n_bad * sizeof(struct badregion));

	if ( in->rigid_groups != NULL ) {

		out->rigid_groups = malloc(out->num_rigid_groups*sizeof(char *));
		memcpy(out->rigid_groups, in->rigid_groups,
		       out->num_rigid_groups*sizeof(char *));

		for ( i=0; i<in->num_rigid_groups; i++ ) {
			out->rigid_groups[i] = strdup(in->rigid_groups[i]);
		}

	}

	for ( i=0; i<out->n_panels; i++ ) {

		struct panel *p;

		p = &out->panels[i];

		if ( p->clen_from != NULL ) {
			/* Make a copy of the clen_from fields unique to this
			 * copy of the structure. */
			p->clen_from = strdup(p->clen_from);
		}

	}

	for ( i=0; i<in->num_rigid_groups; i++ ) {

		int j;
		char *rg = in->rigid_groups[i];
		char *rgn = out->rigid_groups[i];

		for ( j=0; j<in->n_panels; j++ ) {

			if ( in->panels[j].rigid_group == rg ) {
				out->panels[j].rigid_group = rgn;
			}

		}

	}

	return out;
}


struct detector *simple_geometry(const struct image *image)
{
	struct detector *geom;

	geom = calloc(1, sizeof(struct detector));

	geom->n_panels = 1;
	geom->panels = calloc(1, sizeof(struct panel));

	geom->panels[0].min_fs = 0;
	geom->panels[0].max_fs = image->width-1;
	geom->panels[0].min_ss = 0;
	geom->panels[0].max_ss = image->height-1;
	geom->panels[0].cnx = -image->width / 2.0;
	geom->panels[0].cny = -image->height / 2.0;
	geom->panels[0].rigid_group = NULL;

	geom->panels[0].fsx = 1;
	geom->panels[0].fsy = 0;
	geom->panels[0].ssx = 0;
	geom->panels[0].ssy = 1;

	geom->panels[0].xfs = 1;
	geom->panels[0].xss = 0;
	geom->panels[0].yfs = 0;
	geom->panels[0].yss = 1;

	return geom;
}


int reverse_2d_mapping(double x, double y, double *pfs, double *pss,
                       struct detector *det)
{
	int i;

	for ( i=0; i<det->n_panels; i++ ) {

		struct panel *p = &det->panels[i];
		double cx, cy, fs, ss;

		/* Get position relative to corner */
		cx = x - p->cnx;
		cy = y - p->cny;

		/* Reverse the transformation matrix */
		fs = cx*p->xfs + cy*p->yfs;
		ss = cx*p->xss + cy*p->yss;

		/* In range? */
		if ( fs < 0 ) continue;
		if ( ss < 0 ) continue;
		if ( fs > (p->max_fs-p->min_fs+1) ) continue;
		if ( ss > (p->max_ss-p->min_ss+1) ) continue;

		*pfs = fs + p->min_fs;
		*pss = ss + p->min_ss;
		return 0;

	}

	return 1;
}


static void check_extents(struct panel p, double *min_x, double *min_y,
                          double *max_x, double *max_y, double fs, double ss)
{
	double xs, ys, rx, ry;

	xs = fs*p.fsx + ss*p.ssx;
	ys = fs*p.fsy + ss*p.ssy;

	rx = xs + p.cnx;
	ry = ys + p.cny;

	if ( rx > *max_x ) *max_x = rx;
	if ( ry > *max_y ) *max_y = ry;
	if ( rx < *min_x ) *min_x = rx;
	if ( ry < *min_y ) *min_y = ry;
}


double largest_q(struct image *image)
{
	int fs, ss;
	double ttm = 0.0;
	double qmax = 0.0;

	for ( fs=0; fs<image->width; fs++ ) {
	for ( ss=0; ss<image->height; ss++ ) {

		struct rvec q;
		double tt;

		q = get_q(image, fs, ss, &tt, 1.0/image->lambda);

		if ( tt > ttm ) {
			qmax = modulus(q.u, q.v, q.w);
			ttm = tt;
		}

	}
	}

	return qmax;
}


double smallest_q(struct image *image)
{
	int fs, ss;
	double ttm  = +INFINITY;
	double qmin = +INFINITY;
	for ( fs=0; fs<image->width; fs++ ) {
	for ( ss=0; ss<image->height; ss++ ) {

		struct rvec q;
		double tt;

		q = get_q(image, fs, ss, &tt, 1.0/image->lambda);

		if ( tt < ttm ) {
			qmin = modulus(q.u, q.v, q.w);
			ttm = tt;
		}

	}
	}

	return qmin;
}


void get_pixel_extents(struct detector *det,
                       double *min_x, double *min_y,
                       double *max_x, double *max_y)
{
	int i;

	*min_x = 0.0;
	*max_x = 0.0;
	*min_y = 0.0;
	*max_y = 0.0;

	/* To determine the maximum extents of the detector, put all four
	 * corners of each panel through the transformations and watch for the
	 * biggest */

	for ( i=0; i<det->n_panels; i++ ) {

		check_extents(det->panels[i], min_x, min_y, max_x, max_y,
		              0.0,
		              0.0);

		check_extents(det->panels[i], min_x, min_y, max_x, max_y,
		              0.0,
		              det->panels[i].max_ss-det->panels[i].min_ss+1);

		check_extents(det->panels[i], min_x, min_y, max_x, max_y,
		              det->panels[i].max_fs-det->panels[i].min_fs+1,
		              0.0);

		check_extents(det->panels[i], min_x, min_y, max_x, max_y,
		              det->panels[i].max_fs-det->panels[i].min_fs+1,
		              det->panels[i].max_ss-det->panels[i].min_ss+1);


	}
}


int write_detector_geometry(const char *filename, struct detector *det)
{
	struct panel *p;
	int pi;
	FILE *fh;

	if ( filename == NULL ) return 2;
	if ( det->n_panels < 1 ) return 3;

	fh = fopen(filename, "w");
	if ( fh == NULL ) return 1;

	for ( pi=0; pi<det->n_panels; pi++) {

		p = &(det->panels[pi]);

		if ( p == NULL ) return 4;

		if ( pi > 0 ) fprintf(fh, "\n");

		fprintf(fh, "%s/min_fs = %d\n", p->name, p->min_fs);
		fprintf(fh, "%s/min_ss = %d\n", p->name, p->min_ss);
		fprintf(fh, "%s/max_fs = %d\n", p->name, p->max_fs);
		fprintf(fh, "%s/max_ss = %d\n", p->name, p->max_ss);
		fprintf(fh, "%s/badrow_direction = %C\n", p->name, p->badrow);
		fprintf(fh, "%s/res = %g\n", p->name, p->res);
		fprintf(fh, "%s/peak_sep = %g\n", p->name, p->peak_sep);
		fprintf(fh, "%s/clen = %s\n", p->name, p->clen_from);
		fprintf(fh, "%s/fs = %+fx %+fy\n", p->name, p->fsx, p->fsy);
		fprintf(fh, "%s/ss = %+fx %+fy\n", p->name, p->ssx, p->ssy);
		fprintf(fh, "%s/corner_x = %g\n", p->name, p->cnx);
		fprintf(fh, "%s/corner_y = %g\n", p->name, p->cny);

		if ( p->no_index ) {
			fprintf(fh, "%s/no_index = 1\n", p->name);
		} /* else don't clutter up the file */

		if ( p->rigid_group != NULL ) {
			fprintf(fh, "%s/rigid_group = %s\n",
			        p->name, p->rigid_group);
		}

	}
	fclose(fh);

	return 0;
}
