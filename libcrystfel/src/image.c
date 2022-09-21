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

#include <libcrystfel-config.h>

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
#include "image-msgpack.h"
#include "image-seedee.h"
#include "profile.h"

#include "datatemplate.h"
#include "datatemplate_priv.h"

/** \file image.h */

#ifndef HAVE_HDF5
int is_hdf5_file(const char *filename)
{
	return 0;
}
#endif

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


static struct header_cache_entry *find_cache_entry(struct image *image,
                                                   const char *name)
{
	int i;
	for ( i=0; i<image->n_cached_headers; i++ ) {
		if ( strcmp(name, image->header_cache[i]->header_name) == 0 ) {
			return image->header_cache[i];
		}
	}
	return NULL;
}


void image_cache_header_int(struct image *image,
                            const char *header_name,
                            int header_val)
{
	if ( image->n_cached_headers >= HEADER_CACHE_SIZE ) {
		ERROR("Too many headers to copy.\n");
	} else {

		struct header_cache_entry *ce;
		ce = malloc(sizeof(struct header_cache_entry));

		if ( ce != NULL ) {
			ce->header_name = strdup(header_name);
			ce->val_int = header_val;
			ce->type = HEADER_INT;
			image->header_cache[image->n_cached_headers++] = ce;
		} else {
			ERROR("Failed to add header cache entry.\n");
		}
	}
}


void image_cache_header_float(struct image *image,
                              const char *header_name,
                              double header_val)
{
	if ( image->n_cached_headers >= HEADER_CACHE_SIZE ) {
		ERROR("Too many headers to copy.\n");
	} else {

		struct header_cache_entry *ce;
		ce = malloc(sizeof(struct header_cache_entry));

		if ( ce != NULL ) {
			ce->header_name = strdup(header_name);
			ce->val_float = header_val;
			ce->type = HEADER_FLOAT;
			image->header_cache[image->n_cached_headers++] = ce;
		} else {
			ERROR("Failed to add header cache entry.\n");
		}
	}
}


void image_cache_header_str(struct image *image,
                            const char *header_name,
                            const char *header_val)
{
	if ( strchr(header_val, '\n') != NULL ) {
		ERROR("Header '%s' contains newline (not allowed!)\n");
		return;
	}

	if ( image->n_cached_headers >= HEADER_CACHE_SIZE ) {
		ERROR("Too many headers to copy.\n");
	} else {

		struct header_cache_entry *ce;
		ce = malloc(sizeof(struct header_cache_entry));

		if ( ce != NULL ) {
			ce->header_name = strdup(header_name);
			ce->val_str = strdup(header_val);
			ce->type = HEADER_STR;
			image->header_cache[image->n_cached_headers++] = ce;
		} else {
			ERROR("Failed to add header cache entry.\n");
		}
	}
}


static int read_header_to_cache(struct image *image, const char *from)
{
	switch ( image->data_source_type ) {

		case DATA_SOURCE_TYPE_NONE:
		ERROR("No data source for %s %s - not loading header\n",
		      image->filename, image->ev);
		return 1;

		case DATA_SOURCE_TYPE_HDF5:
		#ifdef HAVE_HDF5
		return image_hdf5_read_header_to_cache(image, from);
		#else
		return 1;
		#endif

		case DATA_SOURCE_TYPE_CBF:
		case DATA_SOURCE_TYPE_CBFGZ:
		return image_cbf_read_header_to_cache(image, from);

		case DATA_SOURCE_TYPE_MSGPACK:
		return image_msgpack_read_header_to_cache(image, from);

		default:
		ERROR("Unrecognised file type %i (read_header_to_cache)\n",
		      image->data_source_type);
		return 1;

	}
}


static struct header_cache_entry *cached_header(struct image *image, const char *from)
{
	struct header_cache_entry *ce;

	if ( image == NULL ) {
		ERROR("Attempt to retrieve a header value without an image\n");
		return NULL;
	}

	ce = find_cache_entry(image, from);
	if ( ce != NULL ) return ce;

	/* Try to get the value from the file */
	if ( read_header_to_cache(image, from) == 0 ) {
		return find_cache_entry(image, from);
	} else {
		ERROR("Couldn't find header value '%s'\n", from);
		return NULL;
	}
}


int image_read_header_float(struct image *image, const char *from, double *val)
{
	struct header_cache_entry *ce;

	ce = cached_header(image, from);
	if ( ce == NULL ) return 1;

	switch ( ce->type ) {

		case HEADER_FLOAT:
		*val = ce->val_float;
		return 0;

		case HEADER_INT:
		*val = ce->val_int;
		return 0;

		case HEADER_STR:
		if ( convert_float(ce->val_str, val) == 0 ) {
			return 0;
		} else {
			ERROR("Value '%s' (%s) can't be converted to float\n",
			      ce->val_str, from);
			return 1;
		}

		default:
		ERROR("Unrecognised header cache type %i\n", ce->type);
		return 1;
	}
}


static DataSourceType file_type(const char *filename)
{
	if ( !file_exists(filename) ) {
		ERROR("File not found: %s (file_type)\n", filename);
		return DATA_SOURCE_TYPE_NONE;
	}

	if ( is_hdf5_file(filename) ) {
		return DATA_SOURCE_TYPE_HDF5;

	} else if ( is_cbf_file(filename) ) {
		return DATA_SOURCE_TYPE_CBF;

	} else if ( is_cbfgz_file(filename) ) {
		return DATA_SOURCE_TYPE_CBFGZ;

	} else {
		ERROR("Unrecognised file type: %s (file_type)\n", filename);
		return DATA_SOURCE_TYPE_UNKNOWN;
	}
}


int image_set_zero_data(struct image *image,
                        const DataTemplate *dtempl)
{
	int pi;

	for ( pi=0; pi<dtempl->n_panels; pi++ ) {
		struct panel_template *p;
		long int i;
		p = &dtempl->panels[pi];
		for ( i=0; i<PANEL_WIDTH(p)*PANEL_HEIGHT(p); i++ ) {
			image->dp[pi][i] = 0.0;
		}
	}

	return 0;
}


int image_create_dp_bad_sat(struct image *image,
                            const DataTemplate *dtempl)
{
	int i;

	image->dp = malloc(dtempl->n_panels*sizeof(float *));
	if ( image->dp == NULL ) {
		ERROR("Failed to allocate data array.\n");
		return 1;
	}

	image->bad = malloc(dtempl->n_panels*sizeof(int *));
	if ( image->bad == NULL ) {
		ERROR("Failed to allocate bad pixel mask\n");
		free(image->dp);
		return 1;
	}

	image->sat = malloc(dtempl->n_panels*sizeof(float *));
	if ( image->sat == NULL ) {
		ERROR("Failed to allocate saturation map\n");
		free(image->dp);
		free(image->bad);
		return 1;
	}

	/* Set all pointers to NULL for easier clean-up */
	for ( i=0; i<dtempl->n_panels; i++ ) {
		image->dp[i] = NULL;
		image->bad[i] = NULL;
		image->sat[i] = NULL;
	}

	for ( i=0; i<dtempl->n_panels; i++ ) {

		int j;
		size_t nel = PANEL_WIDTH(&dtempl->panels[i]) * PANEL_HEIGHT(&dtempl->panels[i]);

		image->dp[i] = malloc(nel*sizeof(float));
		image->bad[i] = calloc(nel, sizeof(int));
		image->sat[i] = malloc(nel*sizeof(float));

		if ( (image->dp[i] == NULL)
		  || (image->bad[i] == NULL)
		  || (image->sat[i] == NULL) ) {
			ERROR("Failed to allocate panel data arrays\n");
			for ( i=0; i<dtempl->n_panels; i++ ) {
				free(image->dp[i]);
				free(image->bad[i]);
				free(image->sat[i]);
			}
			free(image->dp);
			free(image->bad);
			free(image->sat);
			return 1;
		}

		for ( j=0; j<nel; j++ ) {
			image->sat[i][j] = INFINITY;
		}

	}

	return 0;
}


static int image_read_image_data(struct image *image,
                                 const DataTemplate *dtempl)
{
	if ( (image->data_block == NULL)
	  && (!file_exists(image->filename)) )
	{
		ERROR("File not found: %s (read data)\n", image->filename);
		return image_set_zero_data(image, dtempl);
	}

	switch ( image->data_source_type ) {

		case DATA_SOURCE_TYPE_NONE:
		STATUS("No data source for %s %s - not loading.\n",
		       image->filename, image->ev);
		return 1;

		case DATA_SOURCE_TYPE_HDF5:
		#ifdef HAVE_HDF5
		return image_hdf5_read(image, dtempl);
		#else
		return 1;
		#endif

		case DATA_SOURCE_TYPE_CBF:
		return image_cbf_read(image, dtempl, 0);

		case DATA_SOURCE_TYPE_CBFGZ:
		return image_cbf_read(image, dtempl, 1);

		case DATA_SOURCE_TYPE_MSGPACK:
		return image_msgpack_read(image, dtempl, image->data_block,
		                          image->data_block_size);

		case DATA_SOURCE_TYPE_SEEDEE:
		return image_seedee_read(image, dtempl, image->data_block,
		                         image->data_block_size,
		                         image->meta_data);

		default:
		ERROR("Unrecognised file type %i (image_read_image_data)\n",
		      image->data_source_type);
		return 1;
	}
}


static int set_image_parameters(struct image *image,
                                const DataTemplate *dtempl)
{
	double wl_val;
	if ( convert_float(dtempl->wavelength_from, &wl_val) ) {

		/* Not a literal value - try header */
		if ( image_read_header_float(image,
		                             dtempl->wavelength_from,
		                             &wl_val) )
		{
			ERROR("Failed to read header value for wavelength (%s)\n",
			      dtempl->wavelength_from);
			return 1;
		}
	}

	image->lambda = convert_to_m(wl_val, dtempl->wavelength_unit);

	image->bw = dtempl->bandwidth;

	image->spectrum = spectrum_generate_gaussian(image->lambda,
	                                             image->bw);

	return 0;
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


static void mark_flagged_pixels(struct panel_template *p,
                                float *dp, int *bad)
{
	int p_w, p_h;
	long int n;
	int i;

	p_w = p->orig_max_fs - p->orig_min_fs + 1;
	p_h = p->orig_max_ss - p->orig_min_ss + 1;
	n = p_w * p_h;

	profile_start("flag-values");
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
	profile_end("flag-values");
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
		#ifdef HAVE_HDF5
		return image_hdf5_read_mask(p, mask_fn, ev, bad, mask_location,
		                            mask_good, mask_bad);
		#endif

	} else if ( is_cbf_file(mask_fn) ) {
		return image_cbf_read_mask(p, mask_fn, ev, 0, bad,
		                           mask_good, mask_bad);

	} else if ( is_cbfgz_file(mask_fn) ) {
		return image_cbf_read_mask(p, mask_fn, ev, 1, bad,
		                           mask_good, mask_bad);

	} else {
		ERROR("Unrecognised mask file type (%s)\n", mask_fn);
		return 1;
	}

	return 0;
}


static void mask_panel_edges(int *bad, int p_w, int p_h, int edgew)
{
	int i;

	/* Silly values should not cause memory errors */
	if ( edgew > p_w ) edgew = p_w/2 + 1;
	if ( edgew > p_h ) edgew = p_h/2 + 1;
	if ( edgew < 0 ) return;

	for ( i=0; i<edgew; i++ ) {
		int fs, ss;
		for ( fs=i; fs<p_w-i; fs++ ) {
			bad[fs+p_w*i] = 1;
			bad[fs+p_w*(p_h-i-1)] = 1;
		}
		for ( ss=i; ss<p_h-i; ss++ ) {
			bad[i+p_w*ss] = 1;
			bad[(p_w-i-1)+p_w*ss] = 1;
		}
	}
}


static int create_badmap(struct image *image,
                         const DataTemplate *dtempl,
                         int no_mask_data)
{
	int i;

	/* The bad pixel map array is already created (see image_create_dp_bad_sat),
	 * and a preliminary mask (with NaN/inf pixels marked) has already been
	 * created when the image data was loaded. */

	for ( i=0; i<dtempl->n_panels; i++ ) {

		int p_w, p_h;
		struct panel_template *p = &dtempl->panels[i];

		p_w = p->orig_max_fs - p->orig_min_fs + 1;
		p_h = p->orig_max_ss - p->orig_min_ss + 1;

		/* Panel marked as bad? */
		if ( p->bad ) {
			profile_start("whole-panel");
			/* NB this sets every element to 0x1111,
			 * but that's OK - value is still 'true'. */
			memset(image->bad[i], 1, p_w*p_h);
			profile_end("whole-panel");
		}

		/* Add bad regions (skip if panel is bad anyway) */
		if ( !p->bad ) {
			profile_start("flagged-pixels");
			mark_flagged_pixels(p, image->dp[i],
			                    image->bad[i]);
			profile_end("flagged-pixels");
		}

		/* Mask panel edges (skip if panel is bad anyway) */
		if ( (p->mask_edge_pixels > 0) && !p->bad ) {
			profile_start("panel-edges");
			mask_panel_edges(image->bad[i], p_w, p_h,
			                 p->mask_edge_pixels);
			profile_end("panel-edges");
		}

		/* Load masks (skip if panel is bad anyway) */
		if ( (!no_mask_data) && (!p->bad) ) {

			int j;
			profile_start("load-masks");

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

				if ( load_mask(p, mask_fn, image->ev, image->bad[i],
				               p->masks[j].data_location,
				               p->masks[j].good_bits,
				               p->masks[j].bad_bits) )
				{
					ERROR("Failed to load mask for %s\n",
					      p->name);
					profile_end("load-masks");
					return 1;
				}

			}
			profile_end("load-masks");
		}
	}

	profile_start("mark-regions");
	mark_bad_regions(image, dtempl);
	profile_end("mark-regions");

	return 0;
}


static int create_satmap(struct image *image,
                         const DataTemplate *dtempl)
{
	int i;
	int any;

	/* The panels will be treated separately, but we'll only bother at all
	 * if at least one of them has a saturation map. */
	any = 0;
	for ( i=0; i<dtempl->n_panels; i++ ) {
		if ( dtempl->panels[i].satmap != NULL ) {
			any = 1;
			break;
		}
	}

	if ( !any ) return 0;

	image->sat = malloc(dtempl->n_panels * sizeof(float *));
	if ( image->sat == NULL ) {
		ERROR("Failed to allocate saturation map\n");
		return 1;
	}

	for ( i=0; i<dtempl->n_panels; i++ ) {

		struct panel_template *p = &dtempl->panels[i];

		if ( p->satmap == NULL ) {

			/* At least one other panel has a saturation map,
			 * but it isn't this one.  Therefore make a fake
			 * saturation map */

			int p_w, p_h;

			p_w = p->orig_max_fs - p->orig_min_fs + 1;
			p_h = p->orig_max_ss - p->orig_min_ss + 1;

			image->sat[i] = malloc(p_w*p_h*sizeof(float));

			if ( image->sat[i] != NULL ) {
				long int j;
				for ( j=0; j<p_w*p_h; j++ ) {
					image->sat[i][j] = INFINITY;
				}
			}

		} else {

			const char *map_fn;

			if ( p->satmap_file == NULL ) {
				map_fn = image->filename;
			} else {
				map_fn = p->satmap_file;
			}

			if ( is_hdf5_file(map_fn) ) {
				#ifdef HAVE_HDF5
				image_hdf5_read_satmap(p, map_fn, image->ev, p->satmap);
				#endif

			} else {
				ERROR("Saturation map must be in HDF5 format\n");
				return 1;
			}
		}

		if ( image->sat[i] == NULL ) {
			ERROR("Failed to allocate saturation map (panel %s)\n",
			      p->name);
			return 1;
		}

	}

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

	if ( image_create_dp_bad_sat(image, dtempl) ) {
		image_free(image);
		return NULL;
	}

	if ( image_set_zero_data(image, dtempl) ) {
		image_free(image);
		return NULL;
	}

	if ( set_image_parameters(image, dtempl) ) {
		image_free(image);
		return NULL;
	}

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


static int do_image_read(struct image *image, const DataTemplate *dtempl,
                         int no_image_data, int no_mask_data)
{
	int i;
	int r;

	if ( image_create_dp_bad_sat(image, dtempl) ) return 1;

	/* Load the image data */
	if ( !no_image_data ) {
		int r;
		profile_start("load-image-data");
		r = image_read_image_data(image, dtempl);
		profile_end("load-image-data");
		if ( r ) return r;
	} else {
		int r;
		profile_start("set-zero-image-data");
		r = image_set_zero_data(image, dtempl);
		profile_end("set-zero-image-data");
		if ( r ) return 1;
	}

	profile_start("set-image-parameters");
	r = set_image_parameters(image, dtempl);
	profile_end("set-image-parameters");
	if ( r ) {
		ERROR("Failed to read image parameters\n");
		return 1;
	}

	profile_start("create-detgeom");
	r = create_detgeom(image, dtempl);
	profile_end("create-detgeom");
	if ( r ) {
		ERROR("Failed to read geometry information\n");
		return 1;
	}

	profile_start("create-badmap");
	r = create_badmap(image, dtempl, no_mask_data);
	profile_end("create-badmap");
	if ( r ) return 1;

	profile_start("create-satmap");
	r = create_satmap(image, dtempl);
	profile_end("create-satmap");
	if ( r ) return 1;

	profile_start("read-headers-to-cache");
	for ( i=0; i<dtempl->n_headers_to_copy; i++ ) {
		read_header_to_cache(image, dtempl->headers_to_copy[i]);
	}
	profile_end("read-headers-to-cache");

	return 0;
}


struct image *image_read(const DataTemplate *dtempl,
                         const char *filename,
                         const char *event,
                         int no_image_data,
                         int no_mask_data)
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

	image->filename = strdup(filename);
	if ( event != NULL ) {
		image->ev = strdup(event);
	} else {
		image->ev = strdup("//");  /* Null event */
	}
	image->data_block = NULL;
	image->data_block_size = 0;

	image->data_source_type = file_type(image->filename);

	if ( do_image_read(image, dtempl, no_image_data, no_mask_data) ) {
		image_free(image);
		return NULL;
	}

	return image;
}


struct image *image_read_data_block(const DataTemplate *dtempl,
                                    void *data_block,
                                    size_t data_block_size,
                                    char *meta_data,
                                    DataSourceType type,
                                    int serial,
                                    int no_image_data,
                                    int no_mask_data)
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

	image->filename = NULL;
	image->ev = NULL;
	image->data_block = data_block;
	image->data_block_size = data_block_size;
	image->meta_data = meta_data;

	image->data_source_type = type;

	if ( do_image_read(image, dtempl, no_image_data, no_mask_data) ) {
		image_free(image);
		ERROR("Failed to load image\n");
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
	free(image->data_block);
	free(image->meta_data);

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

	for ( i=0; i<image->n_cached_headers; i++ ) {
		free(image->header_cache[i]);
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
	image->data_block = NULL;
	image->data_block_size = 0;
	image->meta_data = NULL;
	image->data_source_type = DATA_SOURCE_TYPE_UNKNOWN;

	image->n_cached_headers = 0;
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

		#ifdef HAVE_HDF5
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
		#else
		ERROR("Can't read peak list - compiled without HDF5\n");
		return NULL;
		#endif

	} else  {
		ERROR("Peak lists can only be read from HDF5 files\n");
		return NULL;
	}
}


char **image_expand_frames(const DataTemplate *dtempl,
                           const char *filename, int *n_frames)
{
	if ( is_hdf5_file(filename) ) {
		#ifdef HAVE_HDF5
		return image_hdf5_expand_frames(dtempl, filename,
		                                n_frames);
		#else
		ERROR("Can't expand frames - compiled without HDF5\n");
		return NULL;
		#endif

	} else {
		char **list;
		list = malloc(sizeof(char *));
		if ( list == NULL ) return NULL;
		list[0] = strdup("//");
		if ( list[0] == NULL ) {
			free(list);
			return NULL;
		}
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
		#ifdef HAVE_HDF5
		return image_hdf5_write(image, dtempl, filename);
		#else
		ERROR("Can't write image - compiled without HDF5\n");
		return 1;
		#endif
	}

	ERROR("Can only write to HDF5 files.\n");
	return 1;
}
