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
#include "taketwo.h"


struct _indexingprivate
{
	IndexingFlags flags;
	UnitCell *target_cell;
	float tolerance[4];

	struct taketwo_options *ttopts;

	int n_methods;
	IndexingMethod *methods;
	void **engine_private;
};


static const char *onoff(int a)
{
	if ( a ) return "on";
	return "off";
}


static void show_indexing_flags(IndexingFlags flags)
{
	char check[64];

	assert( !((flags & INDEXING_CHECK_CELL_COMBINATIONS)
	       && (flags & INDEXING_CHECK_CELL_AXES)) );
	STATUS("Indexing parameters:\n");
	strcpy(check, onoff(flags & (INDEXING_CHECK_CELL_COMBINATIONS | INDEXING_CHECK_CELL_AXES)));
	if ( flags & INDEXING_CHECK_CELL_AXES ) {
		strcat(check, " (axis permutations only)");
	}
	if ( flags & INDEXING_CHECK_CELL_COMBINATIONS ) {
		strcat(check, " (axis combinations)");
	}
	STATUS("                  Check unit cell parameters: %s\n", check);
	STATUS("                        Check peak alignment: %s\n",
	       onoff(flags & INDEXING_CHECK_PEAKS));
	STATUS("                   Refine indexing solutions: %s\n",
	       onoff(flags & INDEXING_REFINE));
	STATUS(" Multi-lattice indexing (\"delete and retry\"): %s\n",
	       onoff(flags & INDEXING_MULTI));
	STATUS("                              Retry indexing: %s\n",
	       onoff(flags & INDEXING_RETRY));
}


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


static char *base_indexer_str(IndexingMethod indm)
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

		case INDEXING_TAKETWO :
		strcpy(str, "taketwo");
		break;

		case INDEXING_SIMULATION :
		strcpy(str, "simulation");
		break;

		case INDEXING_DEBUG :
		strcpy(str, "debug");
		break;

		default :
		strcpy(str, "(unknown)");
		break;

	}

	return str;
}


static char *friendly_indexer_name(IndexingMethod m)
{
	char *base = base_indexer_str(m & INDEXING_METHOD_MASK);
	if ( (m & INDEXING_USE_CELL_PARAMETERS)
	  && (m & INDEXING_USE_LATTICE_TYPE) ) {
		strcat(base, " using cell parameters and Bravais lattice type "
		             "as prior information");
	} else if ( m & INDEXING_USE_CELL_PARAMETERS ) {
		strcat(base, " using cell parameters as prior information");
	} else if ( m & INDEXING_USE_LATTICE_TYPE ) {
		strcat(base, " using Bravais lattice type as prior information");
	} else {
		strcat(base, " - no prior information");
	}
	return base;
}


static void *prepare_method(IndexingMethod *m, UnitCell *cell,
                            const char *options)
{
	char *str;
	IndexingMethod in = *m;
	void *priv = NULL;

	switch ( *m & INDEXING_METHOD_MASK ) {

		case INDEXING_NONE :
		priv = "none";
		break;

		case INDEXING_DIRAX :
		priv = dirax_prepare(m, cell);
		break;

		case INDEXING_ASDF :
		priv = asdf_prepare(m, cell);
		break;

		case INDEXING_MOSFLM :
		priv = mosflm_prepare(m, cell);
		break;

		case INDEXING_XDS :
		priv = xds_prepare(m, cell);
		break;

		case INDEXING_DEBUG :
		priv = (IndexingPrivate *)strdup("Hello!");
		break;

		case INDEXING_FELIX :
		priv = felix_prepare(m, cell, options);
		break;

		case INDEXING_TAKETWO :
		priv = taketwo_prepare(m, cell);
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
                                IndexingFlags flags, const char *options,
                                struct taketwo_options *ttopts)
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
		int err = 0;
		methods[i] = get_indm_from_string_2(method_strings[i], &err);
		if ( err ) {
			ERROR("----- Notice -----\n");
			ERROR("The way indexing options are given has changed in this CrystFEL version.\n");
			ERROR("The indexing method should contain only the method itself and ");
			ERROR("prior information modifiers ('cell' or 'latt')\n");
			ERROR("To disable prediction refinement ('norefine'), use --no-refine.\n");
			ERROR("To check cell axes only ('axes'), use --no-cell-combinations.\n");
			ERROR("To disable all unit cell checks ('raw'), use --no-check-cell.\n");
			ERROR("To disable indexing retry ('noretry'), use --no-retry.\n");
			ERROR("Multi-lattice indexing ('multi') is now the default: "
			      "use --no-multi to disable it.\n");
			ERROR("------------------\n");
			free(methods);
			return NULL;
		}

	}

	/* No cell parameters -> no cell checking, no prior cell */
	if ( !cell_has_parameters(cell) ) {

		int warn = 0;

		if ( (flags & INDEXING_CHECK_CELL_COMBINATIONS)
		  || (flags & INDEXING_CHECK_CELL_COMBINATIONS) )
		{
			ERROR("WARNING: Forcing --no-check-cell because "
			      "reference unit cell parameters were not "
			      "given.\n");
			flags &= ~INDEXING_CHECK_CELL_COMBINATIONS;
			flags &= ~INDEXING_CHECK_CELL_AXES;
		}

		for ( i=0; i<n; i++ ) {
			if ( methods[i] & INDEXING_USE_CELL_PARAMETERS ) {
				methods[i] &= ~INDEXING_USE_CELL_PARAMETERS;
				warn = 1;
			}
		}
		if ( warn ) {
			ERROR("WARNING: Forcing all indexing methods to use "
			      "\"-nocell\", because reference cell parameters "
			      "were not given.\n");
		}
	}

	/* No cell at all -> no prior lattice type */
	if ( cell == NULL ) {

		int warn = 0;

		for ( i=0; i<n; i++ ) {
			if ( methods[i] & INDEXING_USE_LATTICE_TYPE ) {
				methods[i] &= ~INDEXING_USE_LATTICE_TYPE;
				warn = 1;
			}
		}
		if ( warn ) {
			ERROR("WARNING: Forcing all indexing methods to use "
			      "\"-nolatt\", because reference Bravais lattice "
			      "type was not given.\n");
		}

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
		                                          options);

		if ( ipriv->engine_private[i] == NULL ) return NULL;

		for ( j=0; j<i; j++ ) {
			if ( methods[i] == methods[j] ) {
				ERROR("Duplicate indexing method.\n");
				return NULL;
			}
		}

	}

	ipriv->methods = methods;
	ipriv->n_methods = n;
	ipriv->flags = flags;

	if ( cell != NULL ) {
		ipriv->target_cell = cell_new_from_cell(cell);
	} else {
		ipriv->target_cell = NULL;
	}
	for ( i=0; i<4; i++ ) ipriv->tolerance[i] = ltl[i];

	ipriv->ttopts = ttopts;

	STATUS("List of indexing methods:\n");
	for ( i=0; i<n; i++ ) {
		char *str = indexer_str(methods[i]);
		char *tmp = friendly_indexer_name(methods[i]);
		STATUS("  %2i: %-25s (%s)\n", i, str, tmp);
		free(str);
		free(tmp);
	}
	show_indexing_flags(flags);

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

			case INDEXING_TAKETWO :
			taketwo_cleanup(ipriv->engine_private[n]);
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


static int check_cell(IndexingFlags flags, Crystal *cr, UnitCell *target,
                      float *tolerance)
{
	if ( (flags & INDEXING_CHECK_CELL_COMBINATIONS)
	  || (flags & INDEXING_CHECK_CELL_AXES) )
	{
		UnitCell *out;
		int reduce;

		if ( flags & INDEXING_CHECK_CELL_COMBINATIONS )
		{
			reduce = 1;
		} else {
			reduce = 0;
		}

		out = match_cell(crystal_get_cell(cr),
		                 target, 0, tolerance, reduce);

		if ( out == NULL ) {
			return 1;
		}

		cell_free(crystal_get_cell(cr));
		crystal_set_cell(cr, out);
	}
	return 0;
}


/* Return non-zero for "success" */
static int try_indexer(struct image *image, IndexingMethod indm,
                       IndexingPrivate *ipriv, void *mpriv)
{
	int i, r;
	int n_bad = 0;
	int n_before = image->n_crystals;

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

		case INDEXING_TAKETWO :
		r = taketwo_index(image, ipriv->ttopts, mpriv);
		break;

		default :
		ERROR("Unrecognised indexing method: %i\n", indm);
		return 0;

	}

	/* Stop a really difficult to debug situation in its tracks */
	assert(image->n_crystals - n_before == r);

	/* For all the crystals found this time ... */
	for ( i=0; i<r; i++ ) {

		int j;
		int this_crystal = image->n_crystals - i - 1;

		/* ... starting at the end of the (complete) list ... */
		Crystal *cr = image->crystals[this_crystal];

		crystal_set_image(cr, image);
		crystal_set_profile_radius(cr, 0.02e9);
		crystal_set_mosaicity(cr, 0.0);

		assert( !((ipriv->flags & INDEXING_CHECK_CELL_COMBINATIONS)
		      && (ipriv->flags & INDEXING_CHECK_CELL_AXES)) );

		/* Pre-refinement unit cell check if requested */
		if ( check_cell(ipriv->flags, cr, ipriv->target_cell,
		                ipriv->tolerance) )
		{
			crystal_set_user_flag(cr, 1);
			continue;
		}

		/* Prediction refinement if requested */
		if ( ipriv->flags & INDEXING_REFINE )
		{
			if ( refine_prediction(image, cr) ) {
				crystal_set_user_flag(cr, 1);
				continue;
			}
		}

		/* After refinement unit cell check if requested */
		if ( check_cell(ipriv->flags, cr, ipriv->target_cell,
		                ipriv->tolerance) )
		{
			crystal_set_user_flag(cr, 1);
			continue;
		}

		/* Peak alignment check if requested */
		if ( ipriv->flags & INDEXING_CHECK_PEAKS )
		{
			if ( !peak_sanity_check(image, &cr, 1) ) {
				crystal_set_user_flag(cr, 1);
				continue;
			}
		}

		/* Don't do similarity check if this crystal is bad */
		if ( crystal_get_user_flag(cr) ) continue;

		/* Check if cell is too similar to existing ones */
		for ( j=0; j<this_crystal; j++ ) {

			Crystal *that_cr = image->crystals[j];

			/* Don't do similarity check against bad crystals */
			if ( crystal_get_user_flag(that_cr) ) continue;

			if ( compare_cells(crystal_get_cell(cr),
			                   crystal_get_cell(that_cr),
			                   0.1, deg2rad(5.0), NULL) )
			{
				crystal_set_user_flag(cr, 1);
			}
		}

	}

	n_bad = remove_flagged_crystals(image);
	assert(r >= n_bad);

	return r - n_bad;
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
static int finished_retry(IndexingMethod indm, IndexingFlags flags,
                          int r, struct image *image)
{
	if ( r == 0 ) {

		/* Indexing failed on the previous attempt.  Maybe try again
		 * after poking the peak list a bit */

		if ( flags & INDEXING_RETRY ) {
			/* Retry with fewer peaks */
			return delete_weakest_peaks(image->features);
		} else {
			/* Indexing failed, opted not to try again */
			return 1;
		}

	} else {

		/* Indexing succeeded on previous attempt.  Maybe try again
		 * after deleting the explained peaks */

		if ( flags & INDEXING_MULTI ) {
			/* Remove "used" spots and try for another lattice */
			Crystal *cr;
			cr = image->crystals[image->n_crystals-1];
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
			done = finished_retry(ipriv->methods[n], ipriv->flags,
			                      r, image);
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
	char *str = base_indexer_str(indm);

	if ( (indm & INDEXING_METHOD_MASK) == INDEXING_SIMULATION ) return str;
	if ( (indm & INDEXING_METHOD_MASK) == INDEXING_NONE ) return str;

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

	return str;
}


static IndexingMethod warn_method(const char *str)
{
	ERROR("Indexing method must contain exactly one engine name: '%s'\n",
	      str);
	return INDEXING_ERROR;
}


IndexingMethod get_indm_from_string_2(const char *str, int *err)
{
	int n, i;
	char **bits;
	IndexingMethod method = INDEXING_NONE;
	int have_method = 0;

	if ( err != NULL ) *err = 0;

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

		} else if ( strcmp(bits[i], "taketwo") == 0) {
			if ( have_method ) return warn_method(str);
			method = INDEXING_DEFAULTS_TAKETWO;
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

		} else if ( strcmp(bits[i], "latt") == 0) {
			method = set_lattice(method);

		} else if ( strcmp(bits[i], "nolatt") == 0) {
			method = set_nolattice(method);

		} else if ( strcmp(bits[i], "cell") == 0) {
			method = set_cellparams(method);

		} else if ( strcmp(bits[i], "nocell") == 0) {
			method = set_nocellparams(method);

		/* Deprecated options */
		} else if ( strcmp(bits[i], "retry") == 0) {
			if ( err != NULL ) *err = 1;
		} else if ( strcmp(bits[i], "noretry") == 0) {
			if ( err != NULL ) *err = 1;
		} else if ( strcmp(bits[i], "multi") == 0) {
			if ( err != NULL ) *err = 1;
		} else if ( strcmp(bits[i], "nomulti") == 0) {
			if ( err != NULL ) *err = 1;
		} else if ( strcmp(bits[i], "refine") == 0) {
			if ( err != NULL ) *err = 1;
		} else if ( strcmp(bits[i], "norefine") == 0) {
			if ( err != NULL ) *err = 1;
		} else if ( strcmp(bits[i], "raw") == 0) {
			if ( err != NULL ) *err = 1;
		} else if ( strcmp(bits[i], "bad") == 0) {
			if ( err != NULL ) *err = 1;
		} else if ( strcmp(bits[i], "comb") == 0) {
			if ( err != NULL ) *err = 1;
		} else if ( strcmp(bits[i], "axes") == 0) {
			if ( err != NULL ) *err = 1;

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


IndexingMethod get_indm_from_string(const char *str)
{
	return get_indm_from_string_2(str, NULL);
}


static void do_probe(const char *(*func)(UnitCell *cell),
                     UnitCell *cell, char *methods)
{
	const char *probe;
	probe = func(cell);
	if ( probe != NULL ) {
		if ( methods[0] != '\0' ) {
			strcat(methods, ",");
		}
		strcat(methods, probe);
	}
}


char *detect_indexing_methods(UnitCell *cell)
{
	char *methods;

	methods = malloc(1024);
	if ( methods == NULL ) return NULL;
	methods[0] = '\0';

	do_probe(mosflm_probe, cell, methods);
	do_probe(dirax_probe, cell, methods);
	do_probe(asdf_probe, cell, methods);
	do_probe(xds_probe, cell, methods);
	do_probe(taketwo_probe, cell, methods);
	/* Don't automatically use Felix (yet) */
	//do_probe(felix_probe, cell, methods);

	if ( strlen(methods) == 0 ) {
		free(methods);
		return NULL;
	}

	return methods;
}
