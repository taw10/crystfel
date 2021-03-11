/*
 * image.c
 *
 * Handle images and image features
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2014      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
 *   2011-2021 Thomas White <taw@physics.org>
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <fenv.h>

#include "image.h"
#include "utils.h"
#include "detgeom.h"
#include "image-hdf5.h"
#include "image-cbf.h"

#include "datatemplate.h"
#include "datatemplate_priv.h"

/** \file image.h */

struct _imagefeaturelist
{
	struct imagefeature *features;
	int                  max_features;
	int                  n_features;
};


void image_add_feature(ImageFeatureList *flist, double fs, double ss,
                       int pn,
                       struct image *parent, double intensity, const char *name)
{
	if ( flist->n_features == flist->max_features ) {
		struct imagefeature *nf;
		int nmf = flist->max_features + 128;
		nf = realloc(flist->features, nmf*sizeof(struct imagefeature));
		if ( nf == NULL ) return;
		flist->features = nf;
		flist->max_features = nmf;
	}

	flist->features[flist->n_features].fs = fs;
	flist->features[flist->n_features].ss = ss;
	flist->features[flist->n_features].pn = pn;
	flist->features[flist->n_features].intensity = intensity;
	flist->features[flist->n_features].name = name;

	flist->n_features++;
}


ImageFeatureList *image_feature_list_new()
{
	ImageFeatureList *flist;

	flist = malloc(sizeof(ImageFeatureList));

	flist->n_features = 0;
	flist->max_features = 0;
	flist->features = NULL;

	return flist;
}


static int comp(const void *a, const void *b)
{
	const struct imagefeature *ap = a;
	const struct imagefeature *bp = b;

	return ap->intensity < bp->intensity;
}


ImageFeatureList *image_feature_list_copy(const ImageFeatureList *flist)
{
	ImageFeatureList *n;
	int nf, i;

	if ( flist == NULL ) return NULL;

	n = image_feature_list_new();
	if ( n == NULL ) return NULL;

	n->features = malloc(flist->n_features*sizeof(struct imagefeature));
	if ( n->features == NULL ) {
		free(n);
		return NULL;
	}

	nf = 0;
	for ( i=0; i<flist->n_features; i++ ) {
		const struct imagefeature *f;
		f = image_get_feature_const(flist, i);
		if ( f == NULL ) continue;
		n->features[nf++] = flist->features[i];
	}
	n->n_features = nf;

	return n;
}


/**
 * Strongest first.
 */
ImageFeatureList *sort_peaks(ImageFeatureList *flist)
{
	ImageFeatureList *n = image_feature_list_copy(flist);
	qsort(n->features, image_feature_count(n),
	      sizeof(struct imagefeature), comp);
	return n;
}


void image_feature_list_free(ImageFeatureList *flist)
{
	if ( flist == NULL ) return;
	free(flist->features);
	free(flist);
}


struct imagefeature *image_feature_closest(ImageFeatureList *flist,
                                           double fs, double ss,
                                           int pn, double *d, int *idx)
{
	int i;
	double dmin = +HUGE_VAL;
	int closest = 0;

	for ( i=0; i<flist->n_features; i++ ) {

		double ds;

		if ( pn != flist->features[i].pn ) continue;

		ds = distance(flist->features[i].fs, flist->features[i].ss,
		              fs, ss);

		if ( ds < dmin ) {
			dmin = ds;
			closest = i;
		}

	}

	if ( dmin < +HUGE_VAL ) {
		*d = dmin;
		*idx = closest;
		return &flist->features[closest];
	}

	*d = +INFINITY;
	return NULL;
}


int image_feature_count(ImageFeatureList *flist)
{
	if ( flist == NULL ) return 0;
	return flist->n_features;
}


const struct imagefeature *image_get_feature_const(const ImageFeatureList *flist,
                                                   int idx)
{
	/* Sanity check */
	if ( flist == NULL ) return NULL;
	if ( idx >= flist->n_features ) return NULL;

	return &flist->features[idx];
}


struct imagefeature *image_get_feature(ImageFeatureList *flist, int idx)
{
	/* Sanity check */
	if ( flist == NULL ) return NULL;
	if ( idx >= flist->n_features ) return NULL;

	return &flist->features[idx];
}


void image_remove_feature(ImageFeatureList *flist, int idx)
{
	memmove(&flist->features[idx], &flist->features[idx+1],
	        (flist->n_features-idx-1)*sizeof(struct imagefeature));
	flist->n_features--;
}


void image_add_crystal(struct image *image, Crystal *cryst)
{
	Crystal **crs;
	int n;

	n = image->n_crystals;
	crs = realloc(image->crystals, (n+1)*sizeof(Crystal *));
	if ( crs == NULL ) {
		ERROR("Failed to allocate memory for crystals.\n");
		return;
	}

	crs[n] = cryst;
	image->crystals = crs;
	image->n_crystals = n+1;
}


int remove_flagged_crystals(struct image *image)
{
	int i;
	int n_bad = 0;

	for ( i=0; i<image->n_crystals; i++ ) {
		if ( crystal_get_user_flag(image->crystals[i]) ) {
			int j;
			Crystal *deleteme = image->crystals[i];
			cell_free(crystal_get_cell(deleteme));
			crystal_free(deleteme);
			for ( j=i; j<image->n_crystals-1; j++ ) {
				image->crystals[j] = image->crystals[j+1];
			}
			image->n_crystals--;
			n_bad++;
			i--;
		}
	}

	return n_bad;
}


/* Free all crystals, including their RefLists and UnitCells */
void free_all_crystals(struct image *image)
{
	int i;
	if ( image->crystals == NULL ) return;
	for ( i=0; i<image->n_crystals; i++ ) {
		Crystal *cr = image->crystals[i];
		reflist_free(crystal_get_reflections(cr));
		cell_free(crystal_get_cell(cr));
		crystal_free(image->crystals[i]);
	}
	free(image->crystals);
	image->n_crystals = 0;
}


static double get_value(struct image *image, const char *from,
                        int *is_literal_number)
{
	double val;
	char *rval;

	if ( from == NULL ) return NAN;

	val = strtod(from, &rval);
	if ( (*rval == '\0') && (rval != from) ) {
		if ( is_literal_number != NULL ) {
			*is_literal_number = 1;
		}
		return val;
	}

	if ( image == NULL ) {
		ERROR("Attempt to retrieve a header value without an image\n");
		return NAN;
	}

	if ( is_hdf5_file(image->filename) ) {
		return image_hdf5_get_value(from,
		                            image->filename,
		                            image->ev);

	} else if ( is_cbf_file(image->filename) ) {
		/* FIXME: From headers */
		return NAN;

	} else if ( is_cbfgz_file(image->filename) ) {
		/* FIXME: From headers */
		return NAN;

	} else {
		ERROR("Unrecognised file type: %s\n", image->filename);
		return NAN;
	}
}


static char *get_value_and_units(struct image *image, const char *from,
                                 double *pvalue,
                                 int *is_literal_number)
{
	char *sp;
	char *fromcpy;
	char *unitscpy;

	if ( from == NULL ) {
		*pvalue = NAN;
		return NULL;
	}

	fromcpy = strdup(from);
	if ( fromcpy == NULL ) {
		*pvalue = NAN;
		return NULL;
	}

	sp = strchr(fromcpy, ' ');
	if ( sp == NULL ) {
		unitscpy = NULL;
	} else {
		unitscpy = strdup(sp+1);
		sp[0] = '\0';
	}

	*pvalue = get_value(image, fromcpy, is_literal_number);
	free(fromcpy);

	return unitscpy;
}


/* default_scale is a value to be used if both of the following
 * conditions are met:
 *
 *  1. The value is a reference to image headers/metadata,
 *      rather than a literal number.
 *  2. No units are specified in the number.
 *
 * This is totally horrible.  Sorry.  Blame history.
 */
static double im_get_length(struct image *image, const char *from,
                            double default_scale)
{
	char *units;
	double value;
	double scale;
	int is_literal_number = 0;

	if ( from == NULL ) return NAN;

	units = get_value_and_units(image, from,
	                            &value, &is_literal_number);
	if ( units == NULL ) {
		if ( is_literal_number ) {
			scale = 1.0;
		} else {
			scale = default_scale;
		}
	} else {
		if ( strcmp(units, "mm") == 0 ) {
			scale = 1e-3;
		} else if ( strcmp(units, "m") == 0 ) {
			scale = 1.0;
		} else {
			ERROR("Invalid length unit '%s'\n", units);
			free(units);
			return NAN;
		}
	}

	free(units);
	return value * scale;
}


static double convert_to_m(double val, int units)
{
	switch ( units ) {

		case WAVELENGTH_M :
		return val;

		case WAVELENGTH_A :
		return val * 1e-10;

		case WAVELENGTH_PHOTON_EV :
		return ph_eV_to_lambda(val);

		case WAVELENGTH_PHOTON_KEV :
		return ph_eV_to_lambda(val*1e3);

		case WAVELENGTH_ELECTRON_V :
		return el_V_to_lambda(val);

		case WAVELENGTH_ELECTRON_KV :
		return el_V_to_lambda(val*1e3);

	}

	return NAN;
}


int create_detgeom(struct image *image, const DataTemplate *dtempl)
{
	struct detgeom *detgeom;
	int i;

	if ( dtempl == NULL ) {
		ERROR("NULL data template!\n");
		return 1;
	}

	detgeom = malloc(sizeof(struct detgeom));
	if ( detgeom == NULL ) return 1;

	detgeom->panels = malloc(dtempl->n_panels*sizeof(struct detgeom_panel));
	if ( detgeom->panels == NULL ) return 1;

	detgeom->n_panels = dtempl->n_panels;

	for ( i=0; i<dtempl->n_panels; i++ ) {

		double shift_x, shift_y;

		detgeom->panels[i].name = safe_strdup(dtempl->panels[i].name);

		detgeom->panels[i].pixel_pitch = dtempl->panels[i].pixel_pitch;

		/* NB cnx,cny are in pixels, cnz is in m */
		detgeom->panels[i].cnx = dtempl->panels[i].cnx;
		detgeom->panels[i].cny = dtempl->panels[i].cny;
		detgeom->panels[i].cnz = im_get_length(image,
		                                       dtempl->panels[i].cnz_from,
		                                       1e-3);

		/* Apply offset (in m) and then convert cnz from
		 * m to pixels */
		detgeom->panels[i].cnz += dtempl->panels[i].cnz_offset;
		detgeom->panels[i].cnz /= detgeom->panels[i].pixel_pitch;

		/* Apply overall shift (already in m) */
		shift_x = im_get_length(image, dtempl->shift_x_from, 1.0);
		shift_y = im_get_length(image, dtempl->shift_y_from, 1.0);

		if ( !isnan(shift_x) ) {
			detgeom->panels[i].cnx += shift_x;
		}
		if ( !isnan(shift_y) ) {
			detgeom->panels[i].cny += shift_y;
		}

		detgeom->panels[i].max_adu = dtempl->panels[i].max_adu;

		switch ( dtempl->panels[i].adu_scale_unit ) {

			case ADU_PER_PHOTON:
			detgeom->panels[i].adu_per_photon = dtempl->panels[i].adu_scale;
			break;

			case ADU_PER_EV:
			detgeom->panels[i].adu_per_photon = dtempl->panels[i].adu_scale
				* ph_lambda_to_eV(image->lambda);
			break;

			default:
			detgeom->panels[i].adu_per_photon = 1.0;
			ERROR("Invalid ADU/ph scale unit (%i)\n",
			      dtempl->panels[i].adu_scale_unit);
			break;

		}

		detgeom->panels[i].w = dtempl->panels[i].orig_max_fs
		                        - dtempl->panels[i].orig_min_fs + 1;
		detgeom->panels[i].h = dtempl->panels[i].orig_max_ss
		                        - dtempl->panels[i].orig_min_ss + 1;

		detgeom->panels[i].fsx = dtempl->panels[i].fsx;
		detgeom->panels[i].fsy = dtempl->panels[i].fsy;
		detgeom->panels[i].fsz = dtempl->panels[i].fsz;
		detgeom->panels[i].ssx = dtempl->panels[i].ssx;
		detgeom->panels[i].ssy = dtempl->panels[i].ssy;
		detgeom->panels[i].ssz = dtempl->panels[i].ssz;

	}

	image->detgeom = detgeom;

	return 0;
}


int image_set_zero_data(struct image *image,
                        const DataTemplate *dtempl)
{
	int pi;

	image->dp = malloc(dtempl->n_panels*sizeof(float *));
	if ( image->dp == NULL ) return 1;

	for ( pi=0; pi<dtempl->n_panels; pi++ ) {

		struct panel_template *p;
		int p_w, p_h;

		p = &dtempl->panels[pi];
		p_w = p->orig_max_fs - p->orig_min_fs + 1;
		p_h = p->orig_max_ss - p->orig_min_ss + 1;

		image->dp[pi] = calloc(p_w*p_h, sizeof(float));
		if ( image->dp[pi] == NULL ) return 1;

	}

	return 0;
}


int image_set_zero_mask(struct image *image,
                        const DataTemplate *dtempl)
{
	int pi;

	image->bad = malloc(dtempl->n_panels*sizeof(int *));
	image->sat = malloc(dtempl->n_panels*sizeof(float *));
	if ( (image->bad == NULL) || (image->sat == NULL) ) return 1;

	for ( pi=0; pi<dtempl->n_panels; pi++ ) {

		struct panel_template *p;
		int p_w, p_h;
		long int i;

		p = &dtempl->panels[pi];
		p_w = p->orig_max_fs - p->orig_min_fs + 1;
		p_h = p->orig_max_ss - p->orig_min_ss + 1;

		image->bad[pi] = calloc(p_w*p_h, sizeof(int));
		image->sat[pi] = calloc(p_w*p_h, sizeof(float));
		if ( image->bad[pi] == NULL ) return 1;
		if ( image->sat[pi] == NULL ) return 1;

		for ( i=0; i<p_w*p_h; i++ ) {
			image->sat[pi][i] = INFINITY;
		}
	}

	return 0;
}


static int file_exists(const char *filename)
{
	struct stat statbuf;
	int r;

	r = stat(filename, &statbuf);
	if ( r != 0 ) {
		return 0;
	}

	return 1;
}


static int image_read_image_data(struct image *image,
                                 const DataTemplate *dtempl,
                                 const char *filename,
                                 const char *event)
{
	if ( !file_exists(filename) ) {
		ERROR("File not found: %s\n", filename);
		return image_set_zero_data(image, dtempl);
	}

	if ( is_hdf5_file(filename) ) {
		return image_hdf5_read(image, dtempl, filename, event);

	} else if ( is_cbf_file(filename) ) {
		return image_cbf_read(image, dtempl, filename, event, 0);

	} else if ( is_cbfgz_file(filename) ) {
		return image_cbf_read(image, dtempl, filename, event, 1);

	}

	ERROR("Unrecognised file type: %s\n", filename);
	return 1;
}


static void set_image_parameters(struct image *image,
                                 const DataTemplate *dtempl)
{
	/* Wavelength might be needed to create detgeom (adu_per_eV) */
	image->lambda = convert_to_m(get_value(image,
	                                       dtempl->wavelength_from,
	                                       NULL),
	                             dtempl->wavelength_unit);

	image->bw = dtempl->bandwidth;

	/* FIXME: Possibly load spectrum from file */
	image->spectrum = spectrum_generate_gaussian(image->lambda,
	                                             image->bw);

}


static void mark_flagged_pixels_lessthan(float *dp, int *bad,
                                         long int n, float val)
{
	long int i;
	for ( i=0; i<n; i++ ) {
		if ( dp[i] < val ) bad[i] = 1;
	}
}


static void mark_flagged_pixels_morethan(float *dp, int *bad,
                                         long int n, float val)
{
	long int i;
	for ( i=0; i<n; i++ ) {
		if ( dp[i] > val ) bad[i] = 1;
	}
}


static void mark_flagged_pixels_equal(float *dp, int *bad,
                                      long int n, float val)
{
	long int i;
	fenv_t envp;

	fegetenv(&envp);
	fesetround(1);  /* Round to nearest (for flag_value) */

	for ( i=0; i<n; i++ ) {
		if ( rint(dp[i]) == val ) bad[i] = 1;
	}

	fesetenv(&envp);
}


static void mark_flagged_pixels_naninf(float *dp, int *bad,
                                       long int n)
{
	long int i;
	for ( i=0; i<n; i++ ) {
		float val = dp[i];
		if ( isnan(val) || isinf(val) ) bad[i] = 1;
	}
}


static void mark_flagged_pixels(struct panel_template *p,
                                float *dp, int *bad)
{
	int p_w, p_h;
	long int n;
	int i;

	p_w = p->orig_max_fs - p->orig_min_fs + 1;
	p_h = p->orig_max_ss - p->orig_min_ss + 1;
	n = p_w * p_h;

	mark_flagged_pixels_naninf(dp, bad, n);

	for ( i=0; i<MAX_FLAG_VALUES; i++ ) {

		float fv = p->flag_values[i];

		switch ( p->flag_types[i] ) {

			case FLAG_NOTHING:
			break;

			case FLAG_LESSTHAN:
			mark_flagged_pixels_lessthan(dp, bad, n, fv);
			break;

			case FLAG_MORETHAN:
			mark_flagged_pixels_morethan(dp, bad, n, fv);
			break;

			case FLAG_EQUAL:
			mark_flagged_pixels_equal(dp, bad, n, fv);
			break;

		}
	}
}


static int region_within_panel(struct dt_badregion *region,
                               struct detgeom_panel *panel)
{
	assert(region->is_fsss);

	if ( region->min_fs < 0 ) return 0;
	if ( region->min_ss < 0 ) return 0;
	if ( region->max_fs >= panel->w ) return 0;
	if ( region->max_ss >= panel->h ) return 0;
	return 1;
}


static void draw_bad_region_fsss(struct dt_badregion *region,
                                 int **bad,
                                 struct detgeom *detgeom)
{
	struct detgeom_panel *panel;
	int fs, ss;

	panel = &detgeom->panels[region->panel_number];

	if ( !region_within_panel(region, panel) ) {
		ERROR("Bad pixel region %s is (partially) outside panel - ignoring it\n",
		      region->name);
		return;
	}

	for ( ss=region->min_ss; ss<=region->max_ss; ss++ ) {
		for ( fs=region->min_fs; fs<=region->max_fs; fs++ ) {
			bad[region->panel_number][fs+ss*panel->w] = 1;
		}
	}
}


static void draw_bad_region_xy(struct dt_badregion *region,
                               int **bad,
                               struct detgeom *detgeom)
{
	int i;

	for ( i=0; i<detgeom->n_panels; i++ ) {

		int fs, ss;

		struct detgeom_panel *p = &detgeom->panels[i];
		for ( ss=0; ss<p->h; ss++ ) {
			for ( fs=0; fs<p->w; fs++ ) {

				double x, y;

				x = fs*p->fsx + ss*p->ssx + p->cnx;
				y = fs*p->fsy + ss*p->ssy + p->cny;

				if ( (x > region->min_x )
				     && (x < region->max_x)
				     && (y > region->min_y)
				     && (y < region->max_y) )
				{
					bad[i][fs+ss*p->w] = 1;
				}

			}
		}
	}
}


static void mark_bad_regions(struct image *image,
                             const DataTemplate *dtempl)
{
	int i;

	for ( i=0; i<dtempl->n_bad; i++ ) {
		if ( dtempl->bad[i].is_fsss ) {
			draw_bad_region_fsss(&dtempl->bad[i],
			                     image->bad,
			                     image->detgeom);
		} else {
			draw_bad_region_xy(&dtempl->bad[i],
			                   image->bad,
			                   image->detgeom);
		}
	}
}


static int load_mask(struct panel_template *p,
                     const char *mask_fn,
                     const char *ev,
                     int *bad,
                     const char *mask_location,
                     unsigned int mask_good,
                     unsigned int mask_bad)
{
	if ( is_hdf5_file(mask_fn) ) {
		image_hdf5_read_mask(p, mask_fn, ev, bad, mask_location,
		                     mask_good, mask_bad);

	} else if ( is_cbf_file(mask_fn) ) {
		image_cbf_read_mask(p, mask_fn, ev, 0, bad,
		                    mask_good, mask_bad);

	} else if ( is_cbfgz_file(mask_fn) ) {
		image_cbf_read_mask(p, mask_fn, ev, 1, bad,
		                    mask_good, mask_bad);

	} else {
		ERROR("Unrecognised mask file type (%s)\n", mask_fn);
		return 1;
	}

	return 0;
}


static int create_badmap(struct image *image,
                         const DataTemplate *dtempl,
                         int no_mask_data)
{
	int i;

	image->bad = malloc(dtempl->n_panels * sizeof(int *));
	if ( image->bad == NULL ) {
		ERROR("Failed to allocate bad pixel mask\n");
		return 1;
	}

	for ( i=0; i<dtempl->n_panels; i++ ) {

		int p_w, p_h;
		struct panel_template *p = &dtempl->panels[i];

		p_w = p->orig_max_fs - p->orig_min_fs + 1;
		p_h = p->orig_max_ss - p->orig_min_ss + 1;

		image->bad[i] = calloc(p_w*p_h, sizeof(int));
		if ( image->bad[i] == NULL ) {
			ERROR("Failed to allocate bad pixel mask\n");
			return 1;
		}

		/* Panel marked as bad? */
		if ( p->bad ) {
			/* NB this sets every element to 0x1111,
			 * but that's OK - value is still 'true'. */
			memset(image->bad[i], 1, p_w*p_h);
		}

		/* Add bad regions (skip if panel is bad anyway) */
		if ( !p->bad ) {
			mark_flagged_pixels(p, image->dp[i],
			                    image->bad[i]);
		}

		/* Load masks (skip if panel is bad anyway) */
		if ( (!no_mask_data) && (!p->bad) ) {

			int j;

			for ( j=0; j<MAX_MASKS; j++ ) {

				const char *mask_fn;

				if ( p->masks[j].data_location == NULL ) {
					continue;
				}

				if ( p->masks[j].filename == NULL ) {
					mask_fn = image->filename;
				} else {
					mask_fn = p->masks[j].filename;
				}

				load_mask(p, mask_fn, image->ev, image->bad[i],
				          p->masks[j].data_location,
				          p->masks[j].good_bits,
				          p->masks[j].bad_bits);

			}
		}
	}

	mark_bad_regions(image, dtempl);

	return 0;
}


static int create_satmap(struct image *image,
                         const DataTemplate *dtempl)
{
	/* FIXME: Implementation */
	return 0;
}


/**
 * Create an image structure, suitable for simulation.
 *
 * WARNING: This is probably not the routine you are looking for!
 *   If you use this routine anywhere other than a simulation program, then
 *   you are abusing the API and can expect breakage.  In particular, your
 *   program will only work when the experiment is completely described by
 *   the DataTemplate, with no values whatsoever coming from image headers.
 *
 * \returns the new image structure.
 *
 */
struct image *image_create_for_simulation(const DataTemplate *dtempl)
{
	struct image *image;

	if ( dtempl == NULL ) {
		ERROR("NULL data template!\n");
		return NULL;
	}

	image = image_new();
	if ( image == NULL ) {
		ERROR("Couldn't allocate image structure.\n");
		return NULL;
	}

	if ( image_set_zero_data(image, dtempl) ) {
		image_free(image);
		return NULL;
	}

	set_image_parameters(image, dtempl);

	if ( create_detgeom(image, dtempl) ) {
		image_free(image);
		return NULL;
	}

	if ( create_badmap(image, dtempl, 1) ) {
		image_free(image);
		return NULL;
	}

	if ( create_satmap(image, dtempl) ) {
		image_free(image);
		return NULL;
	}

	return image;
}


struct image *image_read(const DataTemplate *dtempl,
                         const char *filename,
                         const char *event,
                         int no_image_data,
                         int no_mask_data)
{
	struct image *image;
	int r;

	if ( dtempl == NULL ) {
		ERROR("NULL data template!\n");
		return NULL;
	}

	image = image_new();
	if ( image == NULL ) {
		ERROR("Couldn't allocate image structure.\n");
		return NULL;
	}

	image->filename = strdup(filename);
	if ( event != NULL ) {
		image->ev = strdup(event);
	} else {
		image->ev = strdup("//");  /* Null event */
	}

	/* Load the image data */
	if ( !no_image_data ) {
		r = image_read_image_data(image, dtempl,
		                          filename, event);
	} else {
		r = image_set_zero_data(image, dtempl);
	}
	if ( r ) {
		image_free(image);
		return NULL;
	}

	set_image_parameters(image, dtempl);

	if ( create_detgeom(image, dtempl) ) {
		image_free(image);
		return NULL;
	}

	if ( create_badmap(image, dtempl, no_mask_data) ) {
		image_free(image);
		return NULL;
	}

	if ( create_satmap(image, dtempl) ) {
		image_free(image);
		return NULL;
	}

	return image;
}


void image_free(struct image *image)
{
	int i, np;

	if ( image == NULL ) return;
	image_feature_list_free(image->features);
	free_all_crystals(image);
	spectrum_free(image->spectrum);
	free(image->filename);
	free(image->ev);

	if ( image->detgeom != NULL ) {
		np = image->detgeom->n_panels;
		detgeom_free(image->detgeom);
	} else {
		np = 0;
	}

	for ( i=0; i<np; i++ ) {
		if ( image->dp != NULL ) free(image->dp[i]);
		if ( image->sat != NULL ) free(image->sat[i]);
		if ( image->bad != NULL ) free(image->bad[i]);
	}

	free(image->dp);
	free(image->sat);
	free(image->bad);

	free(image);
}


struct image *image_new()
{
       struct image *image;

       image = malloc(sizeof(struct image));
       if ( image == NULL ) return NULL;

       image->dp = NULL;
       image->bad = NULL;
       image->sat = NULL;
       image->hit = 0;
       image->crystals = NULL;
       image->n_crystals = 0;
       image->indexed_by = INDEXING_NONE;
       image->detgeom = NULL;
       image->filename = NULL;
       image->ev = NULL;
       image->copied_headers = NULL;
       image->id = 0;
       image->serial = 0;
       image->spectrum = NULL;
       image->lambda = -1.0;
       image->div = 0.0;
       image->bw = -1.0;
       image->peak_resolution = -1.0;
       image->features = NULL;

       return image;
}


ImageFeatureList *image_read_peaks(const DataTemplate *dtempl,
                                   const char *filename,
                                   const char *event,
                                   int half_pixel_shift)
{
	if ( is_hdf5_file(filename) ) {

		enum peak_layout layout;

		if ( dtempl->peak_list_type == PEAK_LIST_AUTO ) {

			const char *ext;
			ext = filename_extension(filename, NULL);

			if ( strcmp(ext, ".cxi") == 0 ) {
				layout = PEAK_LIST_CXI;
			} else if ( strcmp(ext, ".h5") == 0 ) {
				layout = PEAK_LIST_LIST3;
			} else {
				ERROR("Couldn't determine peak list layout.\n");
				ERROR("Specify peak_layout = cxi or list3n in geometry file.\n");
				return NULL;
			}

		} else {
			layout = dtempl->peak_list_type;
		}

		switch ( layout ) {

			case PEAK_LIST_CXI :
			return image_hdf5_read_peaks_cxi(dtempl,
			                                 filename,
			                                 event,
			                                 half_pixel_shift);

			case PEAK_LIST_LIST3 :
			return image_hdf5_read_peaks_hdf5(dtempl,
			                                  filename,
			                                  event,
			                                  half_pixel_shift);

			default :
			ERROR("Invalid peak list type %i\n", layout);
			return NULL;

		}

	} else  {
		ERROR("Peak lists can only be read from HDF5 files\n");
		return NULL;
	}
}


char **image_expand_frames(const DataTemplate *dtempl,
                           const char *filename, int *n_frames)
{
	if ( !file_exists(filename) ) {
		ERROR("File not found: %s\n", filename);
		return NULL;
	}

	if ( is_hdf5_file(filename) ) {
		return image_hdf5_expand_frames(dtempl, filename,
		                                n_frames);
	} else {
		char **list;
		list = malloc(sizeof(char *));
		if ( list == NULL ) return NULL;
		list[0] = strdup("//");
		if ( list[0] == NULL ) return NULL;
		*n_frames = 1;
		return list;
	}
}


void mark_resolution_range_as_bad(struct image *image,
                                  double min, double max)
{
	int i;

	if ( isinf(min) && isinf(max) ) return;  /* nothing to do */

	for ( i=0; i<image->detgeom->n_panels; i++ ) {

		int fs, ss;
		struct detgeom_panel *p = &image->detgeom->panels[i];

		for ( ss=0; ss<p->h; ss++ ) {
			for ( fs=0; fs<p->w; fs++ ) {
				double q[3];
				double r;
				detgeom_transform_coords(p, fs, ss,
				                         image->lambda,
				                         0.0, 0.0, q);
				r = modulus(q[0], q[1], q[2]);
				if ( (r >= min) && (r <= max) ) {
					image->bad[i][fs+p->w*ss] = 1;
				}
			}
		}

	}
}


int image_write(const struct image *image,
                const DataTemplate *dtempl,
                const char *filename)
{
	if ( is_hdf5_file(filename) ) {
		return image_hdf5_write(image, dtempl, filename);
	}

	ERROR("Can only write to HDF5 files.\n");
	return 1;
}
