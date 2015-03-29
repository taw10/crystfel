/*
 * detector.c
 *
 * Detector properties
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2009-2015 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
 *   2014      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
 *   2011      Andrew Aquila
 *   2011      Richard Kirian <rkirian@asu.edu>
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

#define _ISOC99_SOURCE
#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "image.h"
#include "utils.h"
#include "detector.h"
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



struct rg_definition {
	char *name;
	char *pns;
};


struct rgc_definition {
	char *name;
	char *rgs;
};


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


static int count_trailing_spaces(const char *string) {

	int i;

	for ( i=0; i<strlen(string); i++ ) {
		if ( !isspace(string[i]) ) {
			return i;
		}
	}

	return -1;
}


static void build_output_line(const char *line, char *new_line,
                              const char *string_to_write)
{
	int nsp, i, w, ww;
	int trailsp;
	int n_bits;
	char **bits;

	n_bits = assplode(line, "=", &bits, ASSPLODE_NONE);
	trailsp = count_trailing_spaces(bits[1]);

	strcat(new_line, bits[0]);
	strcat(new_line, "=");

	for ( nsp=0; nsp < trailsp; nsp++ ) {
		strcat(new_line, " ");
	}
	strcat(new_line, string_to_write);

	for ( w=0; w<strlen(line); w++ ) {
		if ( strncmp(&line[w],";",1) == 0 ) {
			for ( ww=w-1; ww>=0; ww-- ) {
				if ( !isspace(line[ww])) {
					strcat(new_line, &line[ww]);
					break;
				}
			}
			break;
		}
	}
	if ( w == strlen(line) ) strcat(new_line, "\n");

	for ( i=0; i<n_bits; i++) free(bits[i]);
	free(bits);
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

		if ( (b->panel != NULL)
		  && (strcmp(b->panel, p->name) != 0) ) continue;

		if ( b->is_fsss ) {

			int nfs, nss;

			/* fs/ss bad regions are specified according to the
			 * original coordinates */
			nfs = (fs-p->min_fs) + p->orig_min_fs;
			nss = (ss-p->min_ss) + p->orig_min_ss;

			if ( nfs < b->min_fs ) continue;
			if ( nfs > b->max_fs ) continue;
			if ( nss < b->min_ss ) continue;
			if ( nss > b->max_ss ) continue;

		} else {

			if ( rx < b->min_x ) continue;
			if ( rx > b->max_x ) continue;
			if ( ry < b->min_y ) continue;
			if ( ry > b->max_y ) continue;

		}

		return 1;
	}

	return 0;
}


double get_tt(struct image *image, double fs, double ss, int *err)
{
	double r, rx, ry;
	struct panel *p;
	double xs, ys;

	*err = 0;

	p = find_panel(image->det, fs, ss);
	if ( p == NULL ) {
		*err = 1;
		return 0.0;
	}

	/* Convert xs and ys, which are in fast scan/slow scan coordinates,
	 * to x and y */
	xs = (fs-p->min_fs)*p->fsx + (ss-p->min_ss)*p->ssx;
	ys = (fs-p->min_fs)*p->fsy + (ss-p->min_ss)*p->ssy;

	rx = (xs + p->cnx) / p->res;
	ry = (ys + p->cny) / p->res;

	r = sqrt(pow(rx, 2.0) + pow(ry, 2.0));

	return atan2(r, p->clen);
}


void record_image(struct image *image, int do_poisson, double background,
                  gsl_rng *rng, double beam_radius, double nphotons)
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
	area = M_PI*pow(beam_radius, 2.0);
	total_energy = nphotons * ph_lambda_to_en(image->lambda);
	energy_density = total_energy / area;
	ph_per_e = (nphotons /area) * pow(THOMSON_LENGTH, 2.0);
	STATUS("Fluence = %8.2e photons, "
	       "Energy density = %5.3f kJ/cm^2, "
	       "Total energy = %5.3f microJ\n",
	       nphotons, energy_density/1e7, total_energy*1e6);

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double counts;
		double cf;
		double intensity, sa;
		double pix_area, Lsq;
		double xs, ys, rx, ry;
		double dsq, proj_area;
		float dval;
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
		dsq = pow(rx, 2.0) + pow(ry, 2.0);

		/* Projected area of pixel divided by distance squared */
		sa = proj_area / (dsq + Lsq);

		if ( do_poisson ) {
			counts = poisson_noise(rng, intensity * ph_per_e * sa);
		} else {
			cf = intensity * ph_per_e * sa;
			counts = cf;
		}

		/* Number of photons in pixel */
		dval = counts + poisson_noise(rng, background);

		/* Convert to ADU */
		dval *= p->adu_per_eV * ph_lambda_to_eV(image->lambda);

		/* Saturation */
		if ( dval > p->max_adu ) dval = p->max_adu;

		image->data[x + image->width*y] = dval;

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


signed int find_panel_number(struct detector *det, double fs, double ss)
{
	int p;

	/* Fractional pixel coordinates are allowed to be a little further along
	 * than "== max_{f,s}s" for an integer. */

	for ( p=0; p<det->n_panels; p++ ) {
		if ( (fs >= det->panels[p].min_fs)
		  && (fs < det->panels[p].max_fs+1)
		  && (ss >= det->panels[p].min_ss)
		  && (ss < det->panels[p].max_ss+1) ) return p;
	}

	return -1;
}


struct panel *find_panel(struct detector *det, double fs, double ss)
{
	signed int pn = find_panel_number(det, fs, ss);
	if ( pn == -1 ) return NULL;
	return &det->panels[pn];
}


/* Like find_panel(), but uses the original panel bounds, i.e. referring to
 * what's in the HDF5 file */
struct panel *find_orig_panel(struct detector *det, double fs, double ss)
{
	int p;

	for ( p=0; p<det->n_panels; p++ ) {
		if ( (fs >= det->panels[p].orig_min_fs)
		  && (fs < det->panels[p].orig_max_fs+1)
		  && (ss >= det->panels[p].orig_min_ss)
		  && (ss < det->panels[p].orig_max_ss+1) )
		{
			return &det->panels[p];
		}
	}

	return NULL;
}


void fill_in_values(struct detector *det, struct hdfile *f, struct event* ev)
{
	int i;

	for ( i=0; i<det->n_panels; i++ ) {

		struct panel *p = &det->panels[i];

		if ( p->clen_from != NULL ) {

			double val;
			int r;

			r = hdfile_get_value(f, p->clen_from, ev, &val,
			                     H5T_NATIVE_DOUBLE);
			if ( r ) {
				ERROR("Failed to read '%s'\n", p->clen_from);
			} else {
				p->clen = val * 1.0e-3;
			}

		}

		p->clen += p->coffset;

	}
}


int panel_is_in_rigid_group(const struct rigid_group *rg, struct panel *p)
{
	int i;

	for ( i=0; i<rg->n_panels; i++ ) {
		if ( rg->panels[i] == p ) {
			return 1;
		}
	}

	return 0;
}


int rigid_group_is_in_collection(struct rg_collection *c,
                                 struct rigid_group *rg)
{
	int i;

	for ( i=0; i<c->n_rigid_groups; i++ ) {
		if ( c->rigid_groups[i] == rg  ) {
			return 1;
		}
	}

	return 0;
}


struct rg_collection *find_rigid_group_collection_by_name(struct detector *det,
                                                          const char *name)
{
	int i;

	for ( i=0; i<det->n_rg_collections; i++ ) {
		if ( strcmp(det->rigid_group_collections[i]->name,
		            name) == 0 ) {
			return det->rigid_group_collections[i];
		}
	}

	return NULL;
}


static struct panel *new_panel(struct detector *det, const char *name)
{
	struct panel *new;

	det->n_panels++;
	det->panels = realloc(det->panels, det->n_panels*sizeof(struct panel));

	new = &det->panels[det->n_panels-1];
	memcpy(new, &det->defaults, sizeof(struct panel));

	strcpy(new->name, name);

	/* Create a new copy of the camera length location if needed */
	if ( new->clen_from != NULL ) {
		new->clen_from = strdup(new->clen_from);
	}

	/* Create a new copy of the data location if needed */
	if ( new->data != NULL ) {
		new->data = strdup(new->data);
	}

	/* Create a new copy of the dim_structure if needed */
	if ( new->dim_structure != NULL ) {

		struct dim_structure *dim_copy;
		int di;

		dim_copy = initialize_dim_structure();
		dim_copy->num_dims = new->dim_structure->num_dims;
		dim_copy->dims = malloc(dim_copy->num_dims*sizeof(int));
		for ( di=0; di<dim_copy->num_dims; di++ ) {
			dim_copy->dims[di] = new->dim_structure->dims[di];
		}

		new->dim_structure = dim_copy;
	}

	/* Create a new copy of the bad pixel mask location */
	if ( new->mask != NULL ) {
		new->mask = strdup(new->mask);
	}

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
	new->min_fs = 0;
	new->max_fs = 0;
	new->min_ss = 0;
	new->max_ss = 0;
	new->is_fsss = 99; /* Slightly nasty: means "unassigned" */
	new->panel = NULL;
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


static struct rigid_group *find_or_add_rg(struct detector *det,
                                          const char *name)
{
	int i;
	struct rigid_group **new;
	struct rigid_group *rg;

	for ( i=0; i<det->n_rigid_groups; i++ ) {

		if ( strcmp(det->rigid_groups[i]->name, name) == 0 ) {
			return det->rigid_groups[i];
		}

	}

	new = realloc(det->rigid_groups,
	              (1+det->n_rigid_groups)*sizeof(struct rigid_group *));
	if ( new == NULL ) return NULL;

	det->rigid_groups = new;

	rg = malloc(sizeof(struct rigid_group));
	if ( rg == NULL ) return NULL;

	det->rigid_groups[det->n_rigid_groups++] = rg;

	rg->name = strdup(name);
	rg->panels = NULL;
	rg->n_panels = 0;
	rg->have_deltas = 0;

	return rg;
}


static struct rg_collection *find_or_add_rg_coll(struct detector *det,
                                                 const char *name)
{
	int i;
	struct rg_collection **new;
	struct rg_collection *rgc;

	for ( i=0; i<det->n_rg_collections; i++ ) {
		if ( strcmp(det->rigid_group_collections[i]->name, name) == 0 )
		{
			return det->rigid_group_collections[i];
		}
	}

	new = realloc(det->rigid_group_collections,
	              (1+det->n_rg_collections)*sizeof(struct rg_collection *));
	if ( new == NULL ) return NULL;

	det->rigid_group_collections = new;

	rgc = malloc(sizeof(struct rg_collection));
	if ( rgc == NULL ) return NULL;

	det->rigid_group_collections[det->n_rg_collections++] = rgc;

	rgc->name = strdup(name);
	rgc->rigid_groups = NULL;
	rgc->n_rigid_groups = 0;

	return rgc;
}


static void add_to_rigid_group(struct rigid_group *rg, struct panel *p)
{
	struct panel **pn;

	pn = realloc(rg->panels, (1+rg->n_panels)*sizeof(struct panel *));
	if ( pn == NULL ) {
		ERROR("Couldn't add panel to rigid group.\n");
		return;
	}

	rg->panels = pn;
	rg->panels[rg->n_panels++] = p;
}


static void add_to_rigid_group_coll(struct rg_collection *rgc,
                                    struct rigid_group *rg)
{
	struct rigid_group **r;

	r = realloc(rgc->rigid_groups, (1+rgc->n_rigid_groups)*
	            sizeof(struct rigid_group *));
	if ( r == NULL ) {
		ERROR("Couldn't add rigid group to collection.\n");
		return;
	}

	rgc->rigid_groups = r;
	rgc->rigid_groups[rgc->n_rigid_groups++] = rg;
}


/* Free all rigid groups in detector */
static void free_all_rigid_groups(struct detector *det)
{
	int i;

	if ( det->rigid_groups == NULL ) return;
	for ( i=0; i<det->n_rigid_groups; i++ ) {
		free(det->rigid_groups[i]->name);
		free(det->rigid_groups[i]->panels);
		free(det->rigid_groups[i]);
	}
	free(det->rigid_groups);
}


/* Free all rigid groups in detector */
static void free_all_rigid_group_collections(struct detector *det)
{
	int i;

	if ( det->rigid_group_collections == NULL ) return;
	for ( i=0; i<det->n_rg_collections; i++ ) {
		free(det->rigid_group_collections[i]->name);
		free(det->rigid_group_collections[i]->rigid_groups);
		free(det->rigid_group_collections[i]);
	}
	free(det->rigid_group_collections);
}


static struct rigid_group *find_rigid_group_by_name(struct detector *det,
                                                    char *name)
{
	int i;

	for ( i=0; i<det->n_rigid_groups; i++ ) {
		if ( strcmp(det->rigid_groups[i]->name, name) == 0 ) {
			return det->rigid_groups[i];
		}
	}

	return NULL;
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
	} else if ( strcmp(key, "adu_per_eV") == 0 ) {
		panel->adu_per_eV = atof(val);
	} else if ( strcmp(key, "rigid_group") == 0 ) {
		add_to_rigid_group(find_or_add_rg(det, val), panel);
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

	} else if ( strcmp(key, "data") == 0 ) {
		if ( strncmp(val,"/",1) != 0 ) {
			ERROR("Invalid data location '%s'\n", val);
			reject = -1;
		}
		panel->data = strdup(val);

	} else if ( strcmp(key, "mask") == 0 ) {
		if ( strncmp(val,"/",1) != 0 ) {
			ERROR("Invalid mask location '%s'\n", val);
			reject = -1;
		}
		panel->mask = strdup(val);

	} else if ( strcmp(key, "coffset") == 0) {
		panel->coffset = atof(val);
	} else if ( strcmp(key, "res") == 0 ) {
		panel->res = atof(val);
	} else if ( strcmp(key, "max_adu") == 0 ) {
		panel->max_adu = atof(val);
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
	} else if ( strncmp(key, "dim", 3) == 0) {
		if  ( panel->dim_structure == NULL ) {
			panel->dim_structure = initialize_dim_structure();
		}
		set_dim_structure_entry(panel->dim_structure, key, val);
	} else {
		ERROR("Unrecognised field '%s'\n", key);
	}

	return reject;
}


static int check_badr_fsss(struct badregion *badr, int is_fsss)
{
	/* First assignment? */
	if ( badr->is_fsss == 99 ) {
		badr->is_fsss = is_fsss;
		return 0;
	}

	if ( is_fsss != badr->is_fsss ) {
		ERROR("You can't mix x/y and fs/ss in a bad region.\n");
		return 1;
	}

	return 0;
}


static int parse_field_bad(struct badregion *badr, const char *key,
                           const char *val)
{
	int reject = 0;

	if ( strcmp(key, "min_x") == 0 ) {
		badr->min_x = atof(val);
		reject = check_badr_fsss(badr, 0);
	} else if ( strcmp(key, "max_x") == 0 ) {
		badr->max_x = atof(val);
		reject = check_badr_fsss(badr, 0);
	} else if ( strcmp(key, "min_y") == 0 ) {
		badr->min_y = atof(val);
		reject = check_badr_fsss(badr, 0);
	} else if ( strcmp(key, "max_y") == 0 ) {
		badr->max_y = atof(val);
		reject = check_badr_fsss(badr, 0);
	} else if ( strcmp(key, "min_fs") == 0 ) {
		badr->min_fs = atof(val);
		reject = check_badr_fsss(badr, 1);
	} else if ( strcmp(key, "max_fs") == 0 ) {
		badr->max_fs = atof(val);
		reject = check_badr_fsss(badr, 1);
	} else if ( strcmp(key, "min_ss") == 0 ) {
		badr->min_ss = atof(val);
		reject = check_badr_fsss(badr, 1);
	} else if ( strcmp(key, "max_ss") == 0 ) {
		badr->max_ss = atof(val);
		reject = check_badr_fsss(badr, 1);
	} else if ( strcmp(key, "panel") == 0 ) {
		badr->panel = strdup(val);
	} else {
		ERROR("Unrecognised field '%s'\n", key);
	}

	return reject;
}


static void parse_toplevel(struct detector *det, struct beam_params *beam,
                           const char *key, const char *val,
                           struct rg_definition ***rg_defl,
                           struct rgc_definition ***rgc_defl, int *n_rg_defs,
                           int *n_rgc_defs)
{

	if ( strcmp(key, "mask_bad") == 0 ) {

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

	} else if ( strcmp(key, "coffset") == 0 ) {
		det->defaults.coffset = atof(val);

	} else if ( strcmp(key, "photon_energy") == 0 ) {
		if ( beam != NULL ) {
			if ( strncmp(val, "/", 1) == 0 ) {
				beam->photon_energy = 0.0;
				beam->photon_energy_from = strdup(val);
			} else {
				beam->photon_energy = atof(val);
				beam->photon_energy_from = NULL;
			}
		}

	} else if ( strcmp(key, "photon_energy_scale") == 0 ) {
		if ( beam != NULL ) {
			beam->photon_energy_scale = atof(val);
		}

	} else if (strncmp(key, "rigid_group", 11) == 0
	        && strncmp(key, "rigid_group_collection", 22) != 0 ) {

		struct rg_definition **new;

		new = realloc(*rg_defl,
		             ((*n_rg_defs)+1)*sizeof(struct rg_definition*));
		*rg_defl = new;

		(*rg_defl)[*n_rg_defs] = malloc(sizeof(struct rg_definition));
		(*rg_defl)[*n_rg_defs]->name = strdup(key+12);
		(*rg_defl)[*n_rg_defs]->pns = strdup(val);
		*n_rg_defs = *n_rg_defs+1;

	} else if ( strncmp(key, "rigid_group_collection", 22) == 0 ) {

		struct rgc_definition **new;

		new = realloc(*rgc_defl, ((*n_rgc_defs)+1)*
		      sizeof(struct rgc_definition*));
		*rgc_defl = new;

		(*rgc_defl)[*n_rgc_defs] =
		                   malloc(sizeof(struct rgc_definition));
		(*rgc_defl)[*n_rgc_defs]->name = strdup(key+23);
		(*rgc_defl)[*n_rgc_defs]->rgs = strdup(val);
		*n_rgc_defs = *n_rgc_defs+1;

	} else if ( parse_field_for_panel(&det->defaults, key, val, det) ) {
		ERROR("Unrecognised top level field '%s'\n", key);
	}
}


/* Test if fs,ss in panel "p" is further {out,in} than {*p_max_d,*p_min_d}, and
 * if so update det->furthest_{out,in}_{panel,fs,ss}. */
static void check_point(struct panel *p, double fs, double ss,
                        double *p_min_d, double *p_max_d, struct detector *det)
{
	double xs, ys, rx, ry, d;

	xs = fs*p->fsx + ss*p->ssx;
	ys = fs*p->fsy + ss*p->ssy;

	rx = (xs + p->cnx) / p->res;
	ry = (ys + p->cny) / p->res;

	d = sqrt(pow(rx, 2.0) + pow(ry, 2.0));

	if ( d > *p_max_d ) {

		det->furthest_out_panel = p;
		det->furthest_out_fs = fs;
		det->furthest_out_ss = ss;
		*p_max_d = d;

	} else if ( d < *p_min_d ) {

		det->furthest_in_panel = p;
		det->furthest_in_fs = fs;
		det->furthest_in_ss = ss;
		*p_min_d = d;

	}
}


static void find_min_max_d(struct detector *det)
{
	double max_d, min_d;
	int i;

	min_d = +INFINITY;
	max_d = 0.0;
	for ( i=0; i<det->n_panels; i++ ) {

		struct panel *p;
		double w, h;

		p = &det->panels[i];
		w = p->max_fs - p->min_fs + 1;
		h = p->max_ss - p->min_ss + 1;

		check_point(p, 0, 0, &min_d, &max_d, det);
		check_point(p, w, 0, &min_d, &max_d, det);
		check_point(p, 0, h, &min_d, &max_d, det);
		check_point(p, w, h, &min_d, &max_d, det);

	}
}


struct detector *get_detector_geometry(const char *filename,
                                       struct beam_params *beam)
{
	FILE *fh;
	struct detector *det;
	char *rval;
	char **bits;
	int i;
	int rgi, rgci;
	int reject = 0;
	int path_dim;
	int dim_dim;
	int curr_ss;
	int x, y, max_fs, max_ss;
	int dim_reject = 0;
	int dim_dim_reject = 0;
	struct rg_definition **rg_defl = NULL;
	struct rgc_definition **rgc_defl = NULL;
	int n_rg_definitions = 0;
	int n_rgc_definitions = 0;

	fh = fopen(filename, "r");
	if ( fh == NULL ) return NULL;

	det = calloc(1, sizeof(struct detector));
	if ( det == NULL ) {
		fclose(fh);
		return NULL;
	}

	if ( beam != NULL ) {
		beam->photon_energy = 0.0;
		beam->photon_energy_from = NULL;
		beam->photon_energy_scale = 1.0;
	}

	det->n_panels = 0;
	det->panels = NULL;
	det->n_bad = 0;
	det->bad = NULL;
	det->mask_good = 0;
	det->mask_bad = 0;
	det->n_rigid_groups = 0;
	det->rigid_groups = NULL;
	det->path_dim = 0;
	det->dim_dim = 0;
	det->n_rg_collections = 0;
	det->rigid_group_collections = NULL;

	/* The default defaults... */
	det->defaults.min_fs = -1;
	det->defaults.min_ss = -1;
	det->defaults.max_fs = -1;
	det->defaults.max_ss = -1;
	det->defaults.orig_min_fs = -1;
	det->defaults.orig_min_ss = -1;
	det->defaults.orig_max_fs = -1;
	det->defaults.orig_max_ss = -1;
	det->defaults.cnx = NAN;
	det->defaults.cny = NAN;
	det->defaults.clen = NAN;
	det->defaults.coffset = 0.0;
	det->defaults.res = -1.0;
	det->defaults.badrow = '-';
	det->defaults.no_index = 0;
	det->defaults.fsx = 1.0;
	det->defaults.fsy = 0.0;
	det->defaults.ssx = 0.0;
	det->defaults.ssy = 1.0;
	det->defaults.adu_per_eV = NAN;
	det->defaults.max_adu = +INFINITY;
	det->defaults.mask = NULL;
	det->defaults.data = NULL;
	det->defaults.dim_structure = NULL;
	strncpy(det->defaults.name, "", 1023);

	do {

		int n1, n2;
		char **path;
		char line[1024];
		struct badregion *badregion = NULL;
		struct panel *panel = NULL;
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
			parse_toplevel(det, beam, bits[0], bits[2], &rg_defl,
			               &rgc_defl, &n_rg_definitions,
				       &n_rgc_definitions);
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

	path_dim = -1;
	dim_reject = 0;

	for ( i=0; i<det->n_panels; i++ ) {

		int panel_dim = 0;
		char *next_instance;

		next_instance = det->panels[i].data;

		while ( next_instance ) {
			next_instance = strstr(next_instance, "%");
			if ( next_instance != NULL ) {
				next_instance += 1*sizeof(char);
				panel_dim += 1;
			}
		}

		if ( path_dim == -1 ) {
			path_dim = panel_dim;
		} else {
			if ( panel_dim != path_dim ) {
				dim_reject = 1;
			}
		}

	}

	for ( i=0; i<det->n_panels; i++ ) {

		int panel_mask_dim = 0;
		char *next_instance;

		if ( det->panels[i].mask != NULL ) {

			next_instance = det->panels[i].mask;

			while(next_instance)
			{
				next_instance = strstr(next_instance, "%");
				if ( next_instance != NULL ) {
					next_instance += 1*sizeof(char);
					panel_mask_dim += 1;
				}
			}

			if ( panel_mask_dim != path_dim ) {
				dim_reject = 1;
			}
		}
	}

	if ( dim_reject ==  1) {
		ERROR("All panels' data and mask entries must have the same "
		      "number of placeholders\n");
		reject = 1;
	}

	det->path_dim = path_dim;

	dim_dim_reject = 0;
	dim_dim = -1;

	for ( i=0; i<det->n_panels; i++ ) {

		int di;
		int found_ss = 0;
		int found_fs = 0;
		int panel_dim_dim = 0;

		if ( det->panels[i].dim_structure == NULL ) {
			det->panels[i].dim_structure = default_dim_structure();
		}

		for ( di=0; di<det->panels[i].dim_structure->num_dims; di++ ) {

			if ( det->panels[i].dim_structure->dims[di] ==
			                                   HYSL_UNDEFINED  ) {
				dim_dim_reject = 1;
			}
			if ( det->panels[i].dim_structure->dims[di] ==
			                                   HYSL_PLACEHOLDER  ) {
				panel_dim_dim += 1;
			}
			if ( det->panels[i].dim_structure->dims[di] ==
			                                   HYSL_SS  ) {
				found_ss += 1;
			}
			if ( det->panels[i].dim_structure->dims[di] ==
			                                   HYSL_FS  ) {
				found_fs += 1;
			}

		}

		if ( found_ss != 1 ) {
			ERROR("Only one slow scan dim coordinate is allowed\n");
			dim_dim_reject = 1;
		}

		if ( found_fs != 1 ) {
			ERROR("Only one fast scan dim coordinate is allowed\n");
			dim_dim_reject = 1;
		}

		if ( panel_dim_dim > 1 ) {
			ERROR("Maximum one placeholder dim coordinate is "
			      "allowed\n");
			dim_dim_reject = 1;
		}

		if ( dim_dim == -1 ) {
			dim_dim = panel_dim_dim;
		} else {
			if ( panel_dim_dim != dim_dim ) {
				dim_dim_reject = 1;
			}
		}

	}

	if ( dim_dim_reject ==  1) {
		reject = 1;
	}

	det->dim_dim = dim_dim;

	curr_ss = 0;

	for ( i=0; i<det->n_panels; i++ ) {

		if ( det->panels[i].max_fs-det->panels[i].min_fs+1 !=
			det->panels[0].max_fs-det->panels[0].min_fs+1 ) {
			ERROR("All panels should have the same fs extent\n");
			reject = 1;
		}

		if ( det->panels[i ].min_fs < 0 ) {
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
		if ( isnan(det->panels[i].clen)
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
		if ( isnan(det->panels[i].adu_per_eV) ) {
			ERROR("Please specify the number of ADU per eV for"
			      " panel %s\n", det->panels[i].name);
			reject = 1;
		}

		/* It's OK if the badrow direction is '0' */
		/* It's not a problem if "no_index" is still zero */
		/* The default transformation matrix is at least valid */

		det->panels[i].orig_max_fs = det->panels[i].max_fs;
		det->panels[i].orig_min_fs = det->panels[i].min_fs;
		det->panels[i].orig_max_ss = det->panels[i].max_ss;
		det->panels[i].orig_min_ss = det->panels[i].min_ss;

		det->panels[i].w =
		                 det->panels[i].max_fs-det->panels[i].min_fs+1;
		det->panels[i].h =
		                 det->panels[i].max_ss-det->panels[i].min_ss+1;

		det->panels[i].min_fs = 0;
		det->panels[i].max_fs = det->panels[i].w-1;
		det->panels[i].min_ss = curr_ss;
		det->panels[i].max_ss = curr_ss+det->panels[i].h-1;
		curr_ss += det->panels[i].h;

		if ( det->panels[i].max_fs > max_fs ) {
			max_fs = det->panels[i].max_fs;
		}
		if ( det->panels[i].max_ss > max_ss ) {
			max_ss = det->panels[i].max_ss;
		}

	}

	for ( i=0; i<det->n_bad; i++ ) {
		if ( det->bad[i].is_fsss == 99 ) {
			ERROR("Please specify the coordinate ranges for"
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

	free(det->defaults.clen_from);
	free(det->defaults.data);
	free(det->defaults.mask);

	for ( rgi=0; rgi<n_rg_definitions; rgi++) {

		int pi, n1;
		struct rigid_group *rigidgroup = NULL;

		rigidgroup = find_or_add_rg(det, rg_defl[rgi]->name);

		n1 = assplode(rg_defl[rgi]->pns, ",", &bits, ASSPLODE_NONE);

		for ( pi=0; pi<n1; pi++ ) {

			struct panel *p;

			p = find_panel_by_name(det, bits[pi]);
			if ( p == NULL ) {
				ERROR("Cannot add panel to rigid group\n");
				ERROR("Panel not found: %s\n", bits[pi]);
				return NULL;
			}
			add_to_rigid_group(rigidgroup, p);
			free(bits[pi]);
		}
		free(bits);
		free(rg_defl[rgi]->name);
		free(rg_defl[rgi]->pns);
		free(rg_defl[rgi]);
	}
	free(rg_defl);

	for ( rgci=0; rgci<n_rgc_definitions; rgci++ ) {

		int rgi, n2;
		struct rg_collection *rgcollection = NULL;

		rgcollection = find_or_add_rg_coll(det, rgc_defl[rgci]->name);

		n2 = assplode(rgc_defl[rgci]->rgs, ",", &bits, ASSPLODE_NONE);

		for ( rgi=0; rgi<n2; rgi++ ) {

			struct rigid_group *r;

			r = find_rigid_group_by_name(det, bits[rgi]);
			if ( r == NULL ) {
				ERROR("Cannot add rigid group to collection\n");
				ERROR("Rigid group not found: %s\n", bits[rgi]);
				return NULL;
			}
			add_to_rigid_group_coll(rgcollection, r);
			free(bits[rgi]);
		}
		free(rgc_defl[rgci]->name);
		free(rgc_defl[rgci]->rgs);
		free(rgc_defl[rgci]);
	}
	free(rgc_defl);

	if ( n_rg_definitions == 0 ) {

		int pi;

		for ( pi=0; pi<det->n_panels; pi++ ) {

			struct rigid_group *rigidgroup = NULL;

			rigidgroup = find_or_add_rg(det, det->panels[pi].name);
			add_to_rigid_group(rigidgroup, &det->panels[pi]);

		}
	}

	if ( n_rgc_definitions == 0 ) {

		int rgi;
		struct rg_collection *rgcollection = NULL;

		rgcollection = find_or_add_rg_coll(det, "default");

		for ( rgi=0; rgi<det->n_rigid_groups; rgi++ ) {

			add_to_rigid_group_coll(rgcollection,
			                        det->rigid_groups[rgi]);

		}
	}

	/* Calculate matrix inverses and other stuff */
	for ( i=0; i<det->n_panels; i++ ) {

		struct panel *p;
		double d;

		p = &det->panels[i];

		if ( p->fsx*p->ssy == p->ssx*p->fsy ) {
			ERROR("Panel %i transformation singular.\n", i);
		}

		d = (double)p->fsx*p->ssy - p->ssx*p->fsy;
		p->xfs = p->ssy / d;
		p->yfs = -p->ssx / d;
		p->xss = -p->fsy / d;
		p->yss = p->fsx / d;

		p->w = p->max_fs - p->min_fs + 1;
		p->h = p->max_ss - p->min_ss + 1;

	}

	find_min_max_d(det);

	if ( reject ) return NULL;

	fclose(fh);

	return det;
}


void free_detector_geometry(struct detector *det)
{
	int i;

	free_all_rigid_groups(det);
	free_all_rigid_group_collections(det);

	for ( i=0; i<det->n_panels; i++ ) {
		free(det->panels[i].clen_from);
		free_dim_structure(det->panels[i].dim_structure);
	}

	free(det->panels);
	free(det->bad);
	free(det);
}


struct detector *copy_geom(const struct detector *in)
{
	struct detector *out;
	int i;

	out = malloc(sizeof(struct detector));
	memcpy(out, in, sizeof(struct detector));

	out->panels = malloc(out->n_panels * sizeof(struct panel));
	memcpy(out->panels, in->panels, out->n_panels * sizeof(struct panel));

	out->bad = malloc(out->n_bad * sizeof(struct badregion));
	memcpy(out->bad, in->bad, out->n_bad * sizeof(struct badregion));

	out->n_rigid_groups = 0;
	out->rigid_groups = NULL;
	out->n_rg_collections = 0;
	out->rigid_group_collections = NULL;

	for ( i=0; i<out->n_panels; i++ ) {

		struct panel *p;

		p = &out->panels[i];

		if ( p->clen_from != NULL ) {
			/* Make a copy of the clen_from fields unique to this
			 * copy of the structure. */
			p->clen_from = strdup(p->clen_from);
		}

		if ( p->data != NULL ) {
			/* Make a copy of the data fields unique to this
			 * copy of the structure. */
			p->data = strdup(p->data);
		}

		if ( p->dim_structure != NULL ) {
			/* Make a copy of the dim_structure fields unique to this
			 * copy of the structure. */

			struct dim_structure *dim_new;
			int di;

			dim_new = initialize_dim_structure();
			dim_new->num_dims = p->dim_structure->num_dims;
			dim_new->dims = malloc(dim_new->num_dims*sizeof(int));
			for ( di=0; di<dim_new->num_dims; di++ ) {
				dim_new->dims[di] = p->dim_structure->dims[di];
			}

			p->dim_structure = dim_new;

		}

		if ( p->clen_from != NULL ) {
			/* Make a copy of the mask fields unique to this
			 * copy of the structure. */
			p->clen_from = strdup(p->clen_from);
		}

	}

	for ( i=0; i<in->n_panels; i++ ) {

		int rgi;

		for ( rgi=0; rgi<in->n_rigid_groups; rgi++ ) {

			if ( panel_is_in_rigid_group(in->rigid_groups[rgi],
			                             &in->panels[i]) )
			{
				struct rigid_group *g;
				g = find_or_add_rg(out,
				                   in->rigid_groups[rgi]->name);
				add_to_rigid_group(g, &out->panels[i]);
			}
		}

		if ( &in->panels[i] == in->furthest_out_panel ) {
			out->furthest_out_panel = &out->panels[i];
		}
		if ( &in->panels[i] == in->furthest_in_panel ) {
			out->furthest_in_panel = &out->panels[i];
		}


	}

	for ( i=0; i<in->n_rigid_groups; i++ ) {

		int rgci;

		for ( rgci=0; rgci<in->n_rg_collections; rgci++ ) {

			const char *n = in->rigid_group_collections[rgci]->name;

			if ( rigid_group_is_in_collection(
			                      in->rigid_group_collections[rgci],
			                      in->rigid_groups[i]) )
			{
				struct rg_collection *rgcoll;
				rgcoll = find_or_add_rg_coll(out, n);
				add_to_rigid_group_coll(rgcoll,
				                        out->rigid_groups[i]);
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
	geom->panels[0].orig_min_fs = 0;
	geom->panels[0].orig_max_fs = image->width-1;
	geom->panels[0].orig_min_ss = 0;
	geom->panels[0].orig_max_ss = image->height-1;
	geom->panels[0].cnx = -image->width / 2.0;
	geom->panels[0].cny = -image->height / 2.0;
	geom->panels[0].max_adu = INFINITY;
	geom->panels[0].orig_min_fs = -1;
	geom->panels[0].orig_max_fs = -1;
	geom->panels[0].orig_min_ss = -1;
	geom->panels[0].orig_max_ss = -1;

	geom->panels[0].fsx = 1;
	geom->panels[0].fsy = 0;
	geom->panels[0].ssx = 0;
	geom->panels[0].ssy = 1;

	geom->panels[0].xfs = 1;
	geom->panels[0].xss = 0;
	geom->panels[0].yfs = 0;
	geom->panels[0].yss = 1;

	geom->panels[0].w = image->width;
	geom->panels[0].h = image->height;

	geom->panels[0].mask = NULL;
	geom->panels[0].data = NULL;

	find_min_max_d(geom);

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
		if ( fs > p->max_fs-p->min_fs ) continue;
		if ( ss > p->max_ss-p->min_ss ) continue;

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


static void process_panel_fields(const struct panel *p, char *line,
                                     FILE *fh, char **bits,
				     int write_panel_coffset)
{
	char new_line[1024];
	char string_to_write[512];

	strcpy(new_line,"\0");
	strcpy(string_to_write,"\0");

	if(strstr(bits[1], "fs") != NULL &&
	   strstr(bits[1], "min_fs") == NULL &&
	   strstr(bits[1], "max_fs") == NULL &&
	   strstr(bits[1], "offset") == NULL   ) {

		sprintf(string_to_write, "%+fx %+fy",
			p->fsx, p->fsy);
		build_output_line(line, new_line,
		                  string_to_write);
		fputs(new_line, fh);
		return;

	} else if ( strstr(bits[1], "ss") != NULL &&
		    strstr(bits[1], "min_ss") == NULL &&
		    strstr(bits[1], "max_ss") == NULL) {

		sprintf(string_to_write, "%+fx %+fy",
			p->ssx, p->ssy);
		build_output_line(line, new_line,
		                  string_to_write);
		fputs(new_line, fh);
		return;

	} else if ( strstr(bits[1], "corner_x") != NULL) {

		sprintf(string_to_write, "%g",
			p->cnx);
		build_output_line(line, new_line,
		                  string_to_write);
		fputs(new_line, fh);
		return;

	} else if ( strstr(bits[1], "corner_y") != NULL) {

		sprintf(string_to_write, "%g",
			p->cny);
		build_output_line(line, new_line,
		                  string_to_write);
		fputs(new_line, fh);
		return;

	} else if ( strstr(bits[1], "coffset") != NULL) {

		if ( write_panel_coffset ) {
			return;
		} else {
			fputs(line, fh);
			return;
		}

	} else {
		fputs(line, fh);
	}
}


double largest_q(struct image *image)
{
	struct rvec q;
	double tt;

	q = get_q_for_panel(image->det->furthest_out_panel,
	                    image->det->furthest_out_fs,
	                    image->det->furthest_out_ss,
	                    &tt, 1.0/image->lambda);

	return modulus(q.u, q.v, q.w);
}


double smallest_q(struct image *image)
{
	struct rvec q;
	double tt;

	q = get_q_for_panel(image->det->furthest_in_panel,
	                    image->det->furthest_in_fs,
	                    image->det->furthest_in_ss,
	                    &tt, 1.0/image->lambda);

	return modulus(q.u, q.v, q.w);
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


int write_detector_geometry_2(const char *geometry_filename,
                              const char *output_filename, struct detector *det,
                              const char *additional_comment,
                              int write_panel_coffset)
{
	FILE *ifh;
	FILE *fh;

	if ( geometry_filename == NULL ) return 2;
	if ( output_filename == NULL ) return 2;
	if ( det->n_panels < 1 ) return 3;

	ifh = fopen(geometry_filename, "r");
	if ( ifh == NULL ) return 1;

	fh = fopen(output_filename, "w");
	if ( fh == NULL ) {
		return 1;
	}

	if ( additional_comment != NULL ) {

		char additional_comment_tmp[1024];

		strcpy(additional_comment_tmp, "; ");
		strncat(additional_comment_tmp, additional_comment, 1021);
		strcat(additional_comment_tmp, "\n");
		if ( strlen(additional_comment) > 1021 ) {
			strncpy(&additional_comment_tmp[1023], "\0", 1);
		}
		fputs(additional_comment_tmp, fh);
	}
	if  ( write_panel_coffset ) {

		char new_line[1024];
		strcpy(new_line,"; Optimized panel offsets can be found at the "
		                "end of the file\n");
		fputs(new_line, fh);

	}

	do {

		char *rval;
		char line[1024];

		int n_bits;
		char **bits;

		int i;
		struct panel *p;

		rval = fgets(line, 1023, ifh);
		if ( rval == NULL ) break;

		n_bits = assplode(line, "/=", &bits, ASSPLODE_NONE);

		if ( n_bits != 3 ) {
			if ( strstr(bits[0], "coffset" ) != NULL &&
			     write_panel_coffset ) {
				continue;
			} else {
				fputs(line, fh);
			}
		} else {

			p = find_panel_by_name(det, bits[0]);

			if ( p != NULL ) {
				process_panel_fields(p, line, fh, bits,
						     write_panel_coffset);

			} else {
				fputs(line, fh);
			}
		}

		for ( i=0; i<n_bits; i++ ) free(bits[i]);

	} while ( 1 );

	if ( write_panel_coffset ) {

		int pi;
		char new_line[1024];

		strcpy(new_line, "\n");
		fputs(new_line, fh);

		for ( pi=0; pi<det->n_panels; pi++ ) {
			sprintf(new_line, "%s/coffset = %f\n",
			        det->panels[pi].name, det->panels[pi].coffset);
			fputs(new_line, fh);
		}
	}


	fclose(ifh);
	fclose(fh);

	return 0;
}


int write_detector_geometry(const char *geometry_filename,
                            const char *output_filename, struct detector *det)
{
	return write_detector_geometry_2(geometry_filename, output_filename,
	                                 det, NULL, 0);
}


/**
 * mark_resolution_range_as_bad:
 * @image: An image structure
 * @min: Minimum value of 1/d to be marked as bad
 * @max: Maximum value of 1/d to be marked as bad
 *
 * Flags, in the bad pixel mask for @image, every pixel whose resolution is
 * between @min and @max.
 *
 */

void mark_resolution_range_as_bad(struct image *image,
                                  double min, double max)
{
	int i;

	for ( i=0; i<image->det->n_panels; i++ ) {

		int fs, ss;
		struct panel *p = &image->det->panels[i];

		for ( ss=0; ss<p->h; ss++ ) {
		for ( fs=0; fs<p->w; fs++ ) {
			struct rvec q;
			double r;
			q = get_q_for_panel(p, fs, ss, NULL, 1.0/image->lambda);
			r = modulus(q.u, q.v, q.w);
			if ( (r >= min) && (r <= max) ) {
				image->bad[i][fs+p->w*ss] = 1;
			}
		}
		}

	}
}


extern int single_panel_data_source(struct detector *det, const char *element)
{
	int pi;
	char *first_datafrom = NULL;
	char *curr_datafrom = NULL;

	if ( det->panels[0].data == NULL ) {
		if ( element != NULL ) {
			first_datafrom = strdup(element);
		} else {
			first_datafrom = strdup("/data/data");
		}
	} else {
		first_datafrom = strdup(det->panels[0].data);
	}

	for ( pi=1;pi<det->n_panels;pi++ ) {

		if ( det->panels[pi].data == NULL ) {
			if ( element != NULL ) {
				curr_datafrom = strdup(element);
			} else {
				curr_datafrom = strdup("/data/data");
			}
		} else {
			curr_datafrom = strdup(det->panels[pi].data);
		}

		if ( strcmp(curr_datafrom, first_datafrom) != 0 ) {
			return 0;
		}
	}

	free(first_datafrom);
	free(curr_datafrom);

	return 1;
}
