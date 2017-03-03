/*
 * index.c
 *
 * Perform indexing (somehow)
 *
 * Copyright © 2012-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2017 Thomas White <taw@physics.org>
 *   2010-2011 Richard Kirian <rkirian@asu.edu>
 *   2012      Lorenzo Galli
 *   2013      Cornelius Gati <cornelius.gati@cfel.de>
 *   2015      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
 *   2014      Takanori Nakane <nakane.t@gmail.com>
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
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <fenv.h>

#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "dirax.h"
#include "asdf.h"
#include "mosflm.h"
#include "xds.h"
#include "detector.h"
#include "index.h"
#include "geometry.h"
#include "cell-utils.h"
#include "felix.h"
#include "predict-refine.h"


struct _indexingprivate
{
	UnitCell *target_cell;
	float tolerance[4];

	int n_methods;
	IndexingMethod *methods;
	void **engine_private;
};


static int debug_index(struct image *image)
{
	Crystal *cr = crystal_new();
	UnitCell *cell = cell_new();
	cell_set_reciprocal(cell, +0.0000e9, +0.0000e9, +0.0000e9,
	                          +0.0000e9, +0.0000e9, +0.0000e9,
	                          +0.0000e9, +0.0000e9, +0.0000e9);
	crystal_set_cell(cr, cell);
	image_add_crystal(image, cr);
	return 1;
}


static void *prepare_method(IndexingMethod *m, UnitCell *cell,
                            struct detector *det, float *ltl,
                            const char *options)
{
	char *str;
	IndexingMethod in = *m;
	void *priv = NULL;

	switch ( *m & INDEXING_METHOD_MASK ) {

		case INDEXING_DIRAX :
		priv = dirax_prepare(m, cell, det, ltl);
		break;

		case INDEXING_ASDF :
		priv = asdf_prepare(m, cell, det, ltl);
		break;

		case INDEXING_MOSFLM :
		priv = mosflm_prepare(m, cell, det, ltl);
		break;

		case INDEXING_XDS :
		priv = xds_prepare(m, cell, det, ltl);
		break;

		case INDEXING_DEBUG :
		priv = (IndexingPrivate *)strdup("Hello!");
		break;

		case INDEXING_FELIX :
		priv = felix_prepare(m, cell, det, ltl, options);
		break;

		default :
		ERROR("Don't know how to prepare indexing method %i\n", *m);
		break;

	}

	str = indexer_str(*m);

	if ( priv == NULL ) {
		ERROR("Failed to prepare indexing method %s\n", str);
		free(str);
		return NULL;
	}

	STATUS("Prepared indexing method %s\n", str);
	free(str);

	if ( in != *m ) {
		ERROR("Note: flags were altered to take into account "
		      "the information provided and/or the limitations "
		      "of the indexing method.\nPlease check the "
		      "methods listed above carefully.\n");
	}

	return priv;
}


IndexingPrivate *setup_indexing(const char *method_list, UnitCell *cell,
                                struct detector *det, float *ltl,
                                int no_refine, const char *options)
{
	int i, n;
	char **method_strings;
	IndexingPrivate *ipriv;

	/* Parse indexing methods */
	n = assplode(method_list, ",", &method_strings, ASSPLODE_NONE);

	IndexingMethod *methods = malloc(n * sizeof(IndexingMethod));
	if ( methods == NULL ) {
		ERROR("Failed to allocate indexing method list\n");
		return NULL;
	}

	for ( i=0; i<n; i++ ) {
		methods[i] = get_indm_from_string(method_strings[i]);
	}

	ipriv = malloc(sizeof(struct _indexingprivate));
	if ( ipriv == NULL ) {
		ERROR("Failed to allocate indexing data\n");
		return NULL;
	}

	ipriv->engine_private = malloc((n+1) * sizeof(void *));

	for ( i=0; i<n; i++ ) {

		int j;

		ipriv->engine_private[i] = prepare_method(&methods[i], cell,
		                                          det, ltl, options);
		for ( j=0; j<i; j++ ) {
			if ( methods[i] == methods[j] ) {
				ERROR("Duplicate indexing method.\n");
				ERROR("Have you specified some flags which "
				      "aren't accepted by one of your "
				      "chosen indexing methods?\n");
				return NULL;
			}
		}

	}

	ipriv->methods = methods;
	ipriv->n_methods = n;

	if ( cell != NULL ) {
		ipriv->target_cell = cell_new_from_cell(cell);
	} else {
		ipriv->target_cell = NULL;
	}
	for ( i=0; i<4; i++ ) ipriv->tolerance[i] = ltl[i];

	return ipriv;
}


void cleanup_indexing(IndexingPrivate *ipriv)
{
	int n;

	if ( ipriv == NULL ) return;  /* Nothing to do */

	for ( n=0; n<ipriv->n_methods; n++ ) {

		switch ( ipriv->methods[n] & INDEXING_METHOD_MASK ) {

			case INDEXING_NONE :
			break;

			case INDEXING_DIRAX :
			dirax_cleanup(ipriv->engine_private[n]);
			break;

			case INDEXING_ASDF :
			asdf_cleanup(ipriv->engine_private[n]);
			break;

			case INDEXING_MOSFLM :
			mosflm_cleanup(ipriv->engine_private[n]);
			break;

			case INDEXING_XDS :
			xds_cleanup(ipriv->engine_private[n]);
			break;

			case INDEXING_FELIX :
			felix_cleanup(ipriv->engine_private[n]);
			break;

			case INDEXING_DEBUG :
			free(ipriv->engine_private[n]);
			break;

			default :
			ERROR("Don't know how to clean up indexing method %i\n",
			      ipriv->methods[n]);
			break;

		}

		n++;

	}

	free(ipriv->methods);
	free(ipriv->engine_private);
	cell_free(ipriv->target_cell);
	free(ipriv);
}


void map_all_peaks(struct image *image)
{
	int i, n;

	n = image_feature_count(image->features);

	/* Map positions to 3D */
	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		struct rvec r;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		r = get_q_for_panel(f->p, f->fs, f->ss,
		                    NULL, 1.0/image->lambda);
		f->rx = r.u;  f->ry = r.v;  f->rz = r.w;

	}
}


/* Return non-zero for "success" */
static int try_indexer(struct image *image, IndexingMethod indm,
                       IndexingPrivate *ipriv, void *mpriv)
{
	int i, r;
	int n_bad = 0;

	switch ( indm & INDEXING_METHOD_MASK ) {

		case INDEXING_NONE :
		return 0;

		case INDEXING_DIRAX :
		r = run_dirax(image, mpriv);
		break;

		case INDEXING_ASDF :
		r = run_asdf(image, mpriv);
		break;

		case INDEXING_MOSFLM :
		r = run_mosflm(image, mpriv);
		break;

		case INDEXING_XDS :
		r = run_xds(image, mpriv);
		break;

		case INDEXING_DEBUG :
		r = debug_index(image);
		break;

		case INDEXING_FELIX :
		r = felix_index(image, mpriv);
		break;

		default :
		ERROR("Unrecognised indexing method: %i\n", indm);
		return 0;

	}

	for ( i=0; i<r; i++ ) {

		Crystal *cr = image->crystals[image->n_crystals-i-1];
		crystal_set_image(cr, image);
		crystal_set_user_flag(cr, 0);
		crystal_set_profile_radius(cr, 0.02e9);
		crystal_set_mosaicity(cr, 0.0);

		/* Prediction refinement if requested */
		if ( (indm & INDEXING_REFINE)
		   && (refine_prediction(image, cr) != 0) )
		{
			crystal_set_user_flag(cr, 1);
			n_bad++;
		}

	}

	remove_flagged_crystals(image);

	if ( n_bad == r ) return 0;

	return r;
}


static int delete_weakest_peaks(ImageFeatureList *peaks)
{
	int i;
	int np, ndel, n;

	np = image_feature_count(peaks);
	n = 0;
	for ( i=0; i<np; i++ ) {
		if ( image_get_feature(peaks, i) != NULL ) n++;
	}

	if ( n < 6 ) return 1;

	/* Delete 10% weakest peaks */
	ndel = n/10;
	if ( ndel == 0 ) return 1;

	for ( i=n-1; i>n-ndel-1; i-- ) {
		image_remove_feature(peaks, i);
	}

	return 0;
}


static int delete_explained_peaks(struct image *image, Crystal *cr)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	const double min_dist = 0.25;
	int i, nspots = 0, nindexed = 0;

	/* Round towards nearest */
	fesetround(1);

	/* Cell basis vectors for this image */
	cell_get_cartesian(crystal_get_cell(cr), &ax, &ay, &az,
	                   &bx, &by, &bz, &cx, &cy, &cz);

	/* Loop over peaks, checking proximity to nearest reflection */
	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		struct rvec q;
		double h, k, l, hd, kd, ld;
		double dsq;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;
		nspots++;

		/* Reciprocal space position of found peak */
		q = get_q_for_panel(f->p, f->fs, f->ss,
		                    NULL, 1.0/image->lambda);

		/* Decimal and fractional Miller indices of nearest
		 * reciprocal lattice point */
		hd = q.u * ax + q.v * ay + q.w * az;
		kd = q.u * bx + q.v * by + q.w * bz;
		ld = q.u * cx + q.v * cy + q.w * cz;
		h = lrint(hd);
		k = lrint(kd);
		l = lrint(ld);

		/* Check distance */
		dsq = pow(h-hd, 2.0) + pow(k-kd, 2.0) + pow(l-ld, 2.0);

		if ( sqrt(dsq) < min_dist ) {
			image_remove_feature(image->features, i);
			nindexed++;
		}
	}

	/* Return TRUE if not enough peaks to continue */
	return (nspots - nindexed) < 5;

}


/* indm = the current indexing method
 * r = the result from try_indexer() on this method just now
 * image = the image structure
 *
 * Returns false for "try again", true for "no, stop now"
 */
static int finished_retry(IndexingMethod indm, int r, struct image *image)
{
	if ( r == 0 ) {

		/* Indexing failed on the previous attempt.  Maybe try again
		 * after poking the peak list a bit */

		if ( indm & INDEXING_RETRY ) {
			/* Retry with fewer peaks */
			return delete_weakest_peaks(image->features);
		} else {
			/* Indexing failed, opted not to try again */
			return 1;
		}

	} else {

		/* Indexing succeeded on previous attempt.  Maybe try again
		 * after deleting the explained peaks */

		if ( indm & INDEXING_MULTI ) {
			/* Remove "used" spots and try for
			 * another lattice */
			Crystal *cr;
			cr = image->crystals[image->n_crystals-1];
			STATUS("WARNING: Multi-lattice indexing does not work"
			       " well in this version.\n");
			return delete_explained_peaks(image, cr);
		} else {
			return 1;
		}

	}
}

void index_pattern(struct image *image, IndexingPrivate *ipriv)
{
	index_pattern_2(image, ipriv, NULL);
}


void index_pattern_2(struct image *image, IndexingPrivate *ipriv, int *ping)
{
	int n = 0;
	ImageFeatureList *orig;

	if ( ipriv == NULL ) return;

	map_all_peaks(image);
	image->crystals = NULL;
	image->n_crystals = 0;

	orig = image->features;

	for ( n=0; n<ipriv->n_methods; n++ ) {

		int done = 0;
		int r;
		int ntry = 0;
		int success = 0;

		image->features = sort_peaks(orig);

		do {

			r = try_indexer(image, ipriv->methods[n],
			                ipriv, ipriv->engine_private[n]);
			success += r;
			ntry++;
			done = finished_retry(ipriv->methods[n], r, image);
			if ( ntry > 5 ) done = 1;
			if ( ping != NULL ) (*ping)++;

		} while ( !done );

		image_feature_list_free(image->features);

		/* Stop now if the pattern is indexed (don't try again for more
		 * crystals with a different indexing method) */
		if ( success ) break;

	}

	if ( n < ipriv->n_methods ) {
		image->indexed_by = ipriv->methods[n];
	} else {
		image->indexed_by = INDEXING_NONE;
	}

	image->features = orig;
}


/* Set the indexer flags for "raw mode" ("--cell-reduction=none") */
static IndexingMethod set_raw(IndexingMethod a)
{
	/* Disable all unit cell checks */
	a &= ~(INDEXING_CHECK_CELL_COMBINATIONS | INDEXING_CHECK_CELL_AXES);
	return a;
}


/* Set the indexer flags for "bad mode" ("--insane) */
static IndexingMethod set_bad(IndexingMethod a)
{
	/* Disable the peak check */
	return a & ~INDEXING_CHECK_PEAKS;
}


/* Set the indexer flags for "axes mode" ("--cell-reduction=compare") */
static IndexingMethod set_axes(IndexingMethod a)
{
	return (a & ~INDEXING_CHECK_CELL_COMBINATIONS)
	          | INDEXING_CHECK_CELL_AXES;
}


/* Set the indexer flags for "combination mode" ("--cell-reduction=reduce") */
static IndexingMethod set_comb(IndexingMethod a)
{
	return (a & ~INDEXING_CHECK_CELL_AXES)
	          | INDEXING_CHECK_CELL_COMBINATIONS;
}


/* Set the indexer flags for "use no lattice type information" */
static IndexingMethod set_nolattice(IndexingMethod a)
{
	return a & ~INDEXING_USE_LATTICE_TYPE;
}


/* Set the indexer flags for "use lattice type information" */
static IndexingMethod set_lattice(IndexingMethod a)
{
	return a | INDEXING_USE_LATTICE_TYPE;
}


/* Set the indexer flags for "use no unit cell parameters" */
static IndexingMethod set_nocellparams(IndexingMethod a)
{
	return a & ~INDEXING_USE_CELL_PARAMETERS;
}


/* Set the indexer flags for "use unit cell parameters" */
static IndexingMethod set_cellparams(IndexingMethod a)
{
	return a | INDEXING_USE_CELL_PARAMETERS;
}


char *indexer_str(IndexingMethod indm)
{
	char *str;

	str = malloc(256);
	if ( str == NULL ) {
		ERROR("Failed to allocate string.\n");
		return NULL;
	}
	str[0] = '\0';

	switch ( indm & INDEXING_METHOD_MASK ) {

		case INDEXING_NONE :
		strcpy(str, "none");
		return str;

		case INDEXING_DIRAX :
		strcpy(str, "dirax");
		break;

		case INDEXING_ASDF :
		strcpy(str, "asdf");
		break;

		case INDEXING_MOSFLM :
		strcpy(str, "mosflm");
		break;

		case INDEXING_FELIX :
		strcpy(str, "felix");
		break;

		case INDEXING_XDS :
		strcpy(str, "xds");
		break;

		case INDEXING_SIMULATION :
		strcpy(str, "simulation");
		break;

		case INDEXING_DEBUG :
		strcpy(str, "debug");
		break;

		default :
		ERROR("No test description for indexing method %i\n",
		      indm & INDEXING_METHOD_MASK);
		strcpy(str, "(unknown)");
		break;

	}

	if ( (indm & INDEXING_METHOD_MASK) == INDEXING_SIMULATION ) return str;

	if ( indm & INDEXING_CHECK_CELL_COMBINATIONS ) {
		strcat(str, "-comb");
	} else if ( indm & INDEXING_CHECK_CELL_AXES ) {
		strcat(str, "-axes");
	} else {
		strcat(str, "-raw");
	}

	if ( !(indm & INDEXING_CHECK_PEAKS) ) {
		strcat(str, "-bad");
	}

	if ( indm & INDEXING_USE_LATTICE_TYPE ) {
		strcat(str, "-latt");
	} else {
		strcat(str, "-nolatt");
	}

	if ( indm & INDEXING_USE_CELL_PARAMETERS ) {
		strcat(str, "-cell");
	} else {
		strcat(str, "-nocell");
	}

	if ( indm & INDEXING_RETRY ) {
		strcat(str, "-retry");
	} else {
		strcat(str, "-noretry");
	}

	if ( indm & INDEXING_MULTI ) {
		strcat(str, "-multi");
	} else {
		strcat(str, "-nomulti");
	}

	if ( indm & INDEXING_REFINE ) {
		strcat(str, "-refine");
	} else {
		strcat(str, "-norefine");
	}

	return str;
}


static IndexingMethod warn_method(const char *str)
{
	ERROR("Indexing method must contain exactly one engine name: '%s'\n",
	      str);
	return INDEXING_ERROR;
}


IndexingMethod get_indm_from_string(const char *str)
{
	int n, i;
	char **bits;
	IndexingMethod method = INDEXING_NONE;
	int have_method = 0;

	n = assplode(str, "-", &bits, ASSPLODE_NONE);

	for ( i=0; i<n; i++ ) {

		if ( strcmp(bits[i], "dirax") == 0) {
			if ( have_method ) return warn_method(str);
			method = INDEXING_DEFAULTS_DIRAX;
			have_method = 1;

		} else if ( strcmp(bits[i], "asdf") == 0) {
			if ( have_method ) return warn_method(str);
			method = INDEXING_DEFAULTS_ASDF;
			have_method = 1;

		} else if ( strcmp(bits[i], "mosflm") == 0) {
			if ( have_method ) return warn_method(str);
			method = INDEXING_DEFAULTS_MOSFLM;
			have_method = 1;

		} else if ( strcmp(bits[i], "xds") == 0) {
			if ( have_method ) return warn_method(str);
			method = INDEXING_DEFAULTS_XDS;
			have_method = 1;

		} else if ( strcmp(bits[i], "felix") == 0) {
			if ( have_method ) return warn_method(str);
			method = INDEXING_DEFAULTS_FELIX;
			have_method = 1;

		} else if ( strcmp(bits[i], "none") == 0) {
			if ( have_method ) return warn_method(str);
			method = INDEXING_NONE;
			have_method = 1;

		} else if ( strcmp(bits[i], "simulation") == 0) {
			if ( have_method ) return warn_method(str);
			method = INDEXING_SIMULATION;
			return method;

		} else if ( strcmp(bits[i], "debug") == 0) {
			if ( have_method ) return warn_method(str);
			method = INDEXING_DEBUG;
			return method;

		} else if ( strcmp(bits[i], "raw") == 0) {
			method = set_raw(method);

		} else if ( strcmp(bits[i], "bad") == 0) {
			method = set_bad(method);

		} else if ( strcmp(bits[i], "comb") == 0) {
			method = set_comb(method);  /* Default */

		} else if ( strcmp(bits[i], "axes") == 0) {
			method = set_axes(method);

		} else if ( strcmp(bits[i], "latt") == 0) {
			method = set_lattice(method);

		} else if ( strcmp(bits[i], "nolatt") == 0) {
			method = set_nolattice(method);

		} else if ( strcmp(bits[i], "cell") == 0) {
			method = set_cellparams(method);

		} else if ( strcmp(bits[i], "nocell") == 0) {
			method = set_nocellparams(method);

		} else if ( strcmp(bits[i], "retry") == 0) {
			method |= INDEXING_RETRY;

		} else if ( strcmp(bits[i], "noretry") == 0) {
			method &= ~INDEXING_RETRY;

		} else if ( strcmp(bits[i], "multi") == 0) {
			method |= INDEXING_MULTI;

		} else if ( strcmp(bits[i], "nomulti") == 0) {
			method &= ~INDEXING_MULTI;

		} else if ( strcmp(bits[i], "refine") == 0) {
			method |= INDEXING_REFINE;

		} else if ( strcmp(bits[i], "norefine") == 0) {
			method &= ~INDEXING_REFINE;

		} else {
			ERROR("Bad list of indexing methods: '%s'\n", str);
			return INDEXING_ERROR;
		}

		free(bits[i]);

	}
	free(bits);

	if ( !have_method ) return warn_method(str);

	return method;
}
