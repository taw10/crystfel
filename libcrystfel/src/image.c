/*
 * image.c
 *
 * Handle images and image features
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2014      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
 *   2011-2017 Thomas White <taw@physics.org>
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

#include <config.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

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


static double get_value(struct image *image, const char *from)
{
	double val;
	char *rval;

	if ( from == NULL ) return NAN;

	val = strtod(from, &rval);
	if ( (*rval == '\0') && (rval != from) ) return val;

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
                                 double *pvalue)
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

	*pvalue = get_value(image, fromcpy);
	free(fromcpy);

	return unitscpy;
}


static double get_length(struct image *image, const char *from)
{
	char *units;
	double value;
	double scale;

	units = get_value_and_units(image, from, &value);
	if ( units == NULL ) {
		scale = 1.0e-3;
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


void create_detgeom(struct image *image, const DataTemplate *dtempl)
{
	struct detgeom *detgeom;
	int i;

	if ( dtempl == NULL ) {
		ERROR("NULL data template!\n");
		return;
	}

	detgeom = malloc(sizeof(struct detgeom));
	if ( detgeom == NULL ) return;

	detgeom->panels = malloc(dtempl->n_panels*sizeof(struct detgeom_panel));
	if ( detgeom->panels == NULL ) return;

	detgeom->n_panels = dtempl->n_panels;

	for ( i=0; i<dtempl->n_panels; i++ ) {

		detgeom->panels[i].name = safe_strdup(dtempl->panels[i].name);

		detgeom->panels[i].pixel_pitch = dtempl->panels[i].pixel_pitch;

		/* NB cnx,cny are in pixels, cnz is in m */
		detgeom->panels[i].cnx = dtempl->panels[i].cnx;
		detgeom->panels[i].cny = dtempl->panels[i].cny;
		detgeom->panels[i].cnz = get_length(image, dtempl->panels[i].cnz_from);

		/* Apply offset (in m) and then convert cnz from
		 * m to pixels */
		detgeom->panels[i].cnz += dtempl->panels[i].cnz_offset;
		detgeom->panels[i].cnz /= detgeom->panels[i].pixel_pitch;

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

	/* FIXME: spectrum */
}


struct image *image_read(DataTemplate *dtempl, const char *filename,
                         const char *event)
{
	struct image *image;
	int i;

	if ( dtempl == NULL ) {
		ERROR("NULL data template!\n");
		return NULL;
	}

	if ( is_hdf5_file(filename) ) {
		image = image_hdf5_read(dtempl, filename, event);

	} else if ( is_cbf_file(filename) ) {
		image = image_cbf_read(dtempl, filename, event, 0);

	} else if ( is_cbfgz_file(filename) ) {
		image = image_cbf_read(dtempl, filename, event, 1);

	} else {
		ERROR("Unrecognised file type: %s\n", filename);
		return NULL;
	}

	if ( image == NULL ) return NULL;

	/* Wavelength might be needed to create detgeom (adu_per_eV) */
	image->lambda = convert_to_m(get_value(image,
	                                       dtempl->wavelength_from),
	                             dtempl->wavelength_unit);

	create_detgeom(image, dtempl);

	image->bad = malloc(dtempl->n_panels * sizeof(int *));
	if ( image->bad == NULL ) {
		ERROR("Failed to allocate bad pixel mask\n");
		return NULL;
	}

	for ( i=0; i<dtempl->n_panels; i++ ) {

		const char *mask_fn;
		int p_w, p_h;
		struct panel_template *p = &dtempl->panels[i];

		p_w = p->orig_max_fs - p->orig_min_fs + 1;
		p_h = p->orig_max_ss - p->orig_min_ss + 1;

		image->bad[i] = calloc(p_w*p_h, sizeof(int));
		if ( image->bad[i] == NULL ) {
			ERROR("Failed to allocate bad pixel mask\n");
			return NULL;
		}

		/* Panel marked as bad? */
		if ( p->bad ) {
			/* NB this sets every element to 0x1111,
			 * but that's OK - value is still 'true'. */
			memset(image->bad[i], 1, p_w*p_h);
		}

		/* Add bad regions (skip if panel is bad anyway) */
		if ( !p->bad ) {
			int fs, ss;
			for ( fs=0; fs<p_w; fs++ ) {
			for ( ss=0; ss<p_h; ss++ ) {
				if ( data_template_in_bad_region(dtempl, i, fs, ss)
				     || isnan(image->dp[i][fs+ss*p_w])
				     || isinf(image->dp[i][fs+ss*p_w]) )
				{
					image->bad[i][fs+ss*p_w] = 1;
				}
			}
			}
		}

		/* Load mask (skip if panel is bad anyway) */
		if ( (!p->bad) && (p->mask != NULL) ) {
			if ( p->mask_file == NULL ) {
				mask_fn = filename;
			} else {
				mask_fn = p->mask_file;
			}
			if ( is_hdf5_file(mask_fn) ) {
				image_hdf5_read_mask(p, mask_fn, event,
				                     image->bad[i],
				                     dtempl->mask_good,
				                     dtempl->mask_bad);

			} else if ( is_cbf_file(filename) ) {
				image_cbf_read_mask(p, mask_fn, event,
				                    0, image->bad[i],
				                    dtempl->mask_good,
				                    dtempl->mask_bad);

			} else if ( is_cbfgz_file(filename) ) {
				image_cbf_read_mask(p, mask_fn, event,
				                    1, image->bad[i],
				                    dtempl->mask_good,
				                    dtempl->mask_bad);

			} else {
				ERROR("Unrecognised mask file type"
				      " (%s)\n", filename);
				return NULL;
			}
		}
	}

	/* FIXME: Load saturation map */

	return image;
}


void image_free(struct image *image)
{
	int i, np;

	if ( image == NULL ) return;
	image_feature_list_free(image->features);
	free_all_crystals(image);
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
       image->avg_clen = -1.0;
       image->id = 0;
       image->serial = 0;
       image->spectrum = NULL;
       image->lambda = -1.0;
       image->div = -1.0;
       image->bw = -1.0;
       image->peak_resolution = -1.0;
       image->features = NULL;

       return image;
}


int create_blank_arrays(struct image *image)
{
	int pn;
	int num_panels = image->detgeom->n_panels;

	image->dp = malloc(num_panels*sizeof(float *));
	image->bad = malloc(num_panels*sizeof(int *));
	image->sat = malloc(num_panels*sizeof(float *));

	if ( (image->dp == NULL) || (image->bad == NULL)
	  || (image->sat == NULL) ) return 1;

	for ( pn=0; pn<num_panels; pn++ ) {

		long int i;
		struct detgeom_panel *p = &image->detgeom->panels[pn];

		image->dp[pn] = malloc(p->w*p->h*sizeof(float));
		image->bad[pn] = malloc(p->w*p->h*sizeof(int));
		image->sat[pn] = malloc(p->w*p->h*sizeof(float));

		if ( (image->dp[pn] == NULL)
		  || (image->bad[pn] == NULL)
		  || (image->sat[pn] == NULL) )
		{
			return 1;
		}

		for ( i=0; i<p->w*p->h; i++ ) {
			image->dp[pn][i] = 0.0;
			image->bad[pn][i] = 0;
			image->sat[pn][i] = INFINITY;
		}

	}

	return 0;
}


ImageFeatureList *image_read_peaks(const DataTemplate *dtempl,
                                   const char *filename,
                                   const char *event,
                                   int half_pixel_shift)
{
	if ( is_hdf5_file(filename) ) {

		const char *ext;
		ext = filename_extension(filename, NULL);
		if ( strcmp(ext, ".cxi") == 0 ) {
			return image_hdf5_read_peaks_cxi(dtempl,
			                                 filename,
			                                 event,
			                                 half_pixel_shift);

		} else {
			return image_hdf5_read_peaks_hdf5(dtempl,
			                                  filename,
			                                  event,
			                                  half_pixel_shift);
		}

	} else  {
		ERROR("Peak lists can only be read from HDF5 files\n");
		return NULL;
	}
}


char **image_expand_frames(const DataTemplate *dtempl,
                           const char *filename, int *n_frames)
{
	if ( is_hdf5_file(filename) ) {
		return image_hdf5_expand_frames(dtempl, filename,
		                                n_frames);
	} else {
		ERROR("Can only expand HDF5 files\n");
		return NULL;
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
				                         q);
				r = modulus(q[0], q[1], q[2]);
				if ( (r >= min) && (r <= max) ) {
					image->bad[i][fs+p->w*ss] = 1;
				}
			}
		}

	}
}
