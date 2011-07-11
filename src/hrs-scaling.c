/*
 * hrs-scaling.c
 *
 * Intensity scaling using generalised HRS target function
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include "image.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"
#include "utils.h"
#include "reflist.h"


/* Maximum number of iterations of scaling per macrocycle. */
#define MAX_CYCLES (50)

/* ESD of restraint driving scale factors to unity */
#define SCALING_RESTRAINT (1.0)


static double iterate_scale(struct image *images, int n, RefList *reference)
{
	double max_shift = 0.0;
	int frame;

	assert(reference != NULL);

	for ( frame=0; frame<n; frame++ ) {

		struct image *image = &images[frame];
		Reflection *refl;
		RefListIterator *iter;
		double num = 0.0;
		double den = 0.0;
		double corr;

		for ( refl = first_refl(image->reflections, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			signed int h, k, l;
			double Ih, Ihl, esd;
			Reflection *r;

			if ( !get_scalable(refl) ) continue;

			/* Look up by asymmetric indices */
			get_indices(refl, &h, &k, &l);
			r = find_refl(reference, h, k, l);
			if ( r == NULL ) {
				ERROR("%3i %3i %3i isn't in the "
				      "reference list, so why is it "
				      "marked as scalable?\n", h, k, l);
				Ih = 0.0;
			} else {
				if ( get_redundancy(r) < 2 ) continue;
				Ih = get_intensity(r);
			}

			Ihl = get_intensity(refl) / get_partiality(refl);
			esd = get_esd_intensity(refl) / get_partiality(refl);

			num += Ih * (Ihl/image->osf) / pow(esd/image->osf, 2.0);
			den += pow(Ih, 2.0)/pow(esd/image->osf, 2.0);

		}

		//num += image->osf / pow(SCALING_RESTRAINT, 2.0);
		//den += pow(image->osf, 2.0)/pow(SCALING_RESTRAINT, 2.0);

		corr = num / den;
		if ( !isnan(corr) && !isinf(corr) ) {
			image->osf *= corr;
		}
		if ( fabs(corr-1.0) > max_shift ) max_shift = fabs(corr-1.0);

	}

	return max_shift;
}


static RefList *lsq_intensities(struct image *images, int n)
{
	RefList *full;
	int frame;
	Reflection *refl;
	RefListIterator *iter;

	full = reflist_new();

	for ( frame=0; frame<n; frame++ ) {

		double G;

		if ( images[frame].pr_dud ) continue;
		G = images[frame].osf;

		for ( refl = first_refl(images[frame].reflections, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			Reflection *f;
			signed int h, k, l;
			double num, den;
			int red;
			double Ihl, esd;

			if ( !get_scalable(refl) ) continue;

			get_indices(refl, &h, &k, &l);
			f = find_refl(full, h, k, l);
			if ( f == NULL ) {
				f = add_refl(full, h, k, l);
				num = 0.0;
				den = 0.0;
				red = 0;
			} else {
				num = get_temp1(f);
				den = get_temp2(f);
				red = get_redundancy(f);
			}

			Ihl = get_intensity(refl) / get_partiality(refl);
			esd = get_esd_intensity(refl) / get_partiality(refl);

			num += (Ihl/G) / pow(esd/G, 2.0);
			den += 1.0 / pow(esd/G, 2.0);
			red++;

			set_temp1(f, num);
			set_temp2(f, den);
			set_redundancy(f, red);
		}

	}

	for ( refl = first_refl(full, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double Ih;

		Ih = get_temp1(refl) / get_temp2(refl);
		set_int(refl, Ih);

	}

	return full;
}


static UNUSED void show_scale_factors(struct image *images, int n)
{
	int i;
	for ( i=0; i<n; i++ ) {
		STATUS("Image %4i: scale factor %5.2f\n", i, images[i].osf);
	}
}


static UNUSED double total_dev(struct image *image, const RefList *full)
{
	double dev = 0.0;

	/* For each reflection */
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(image->reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double G, p;
		signed int h, k, l;
		Reflection *full_version;
		double I_full, I_partial;

		if ( !get_scalable(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		assert((h!=0) || (k!=0) || (l!=0));

		full_version = find_refl(full, h, k, l);
		if ( full_version == NULL ) continue;
		/* Some reflections may have recently become scalable, but
		 * scale_intensities() might not yet have been called, so the
		 * full version may not have been calculated yet. */

		G = image->osf;
		p = get_partiality(refl);
		I_partial = get_intensity(refl);
		I_full = get_intensity(full_version);
		//STATUS("%3i %3i %3i  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f\n",
		//       h, k, l, G, p, I_partial, I_full,
		//       I_partial - p*G*I_full);

		dev += pow(I_partial - p*G*I_full, 2.0);

	}

	return dev;
}


static UNUSED void plot_graph(struct image *image, RefList *reference)
{
	double sc;

	for ( sc=0.0; sc<3.0; sc+=0.1 ) {

		image->osf = sc;
		STATUS("%5.2f: %e\n", sc, total_dev(image, reference));

	}
}


/* Scale the stack of images */
RefList *scale_intensities(struct image *images, int n, RefList *gref)
{
	int i;
	double max_corr;
	RefList *full = NULL;

	for ( i=0; i<n; i++ ) {
		images[i].osf = 1.0;
	}

	/* No reference -> create an initial list to refine against */
	if ( gref == NULL ) {
		full = lsq_intensities(images, n);
	}

	/* Iterate */
	i = 0;
	do {

		RefList *reference;

		/* Refine against reference or current "full" estimates */
		if ( gref != NULL ) {
			reference = gref;
		} else {
			reference = full;
		}

		max_corr = iterate_scale(images, n, reference);
		STATUS("Scaling iteration %2i: max correction = %5.2f\n",
		       i+1, max_corr);

		/* No reference -> generate list for next iteration */
		if ( gref == NULL ) {
			reflist_free(full);
			full = lsq_intensities(images, n);
		}

		//show_scale_factors(images, n);

		i++;

	} while ( (max_corr > 0.01) && (i < MAX_CYCLES) );

	if ( gref != NULL ) {
		full = lsq_intensities(images, n);
	} /* else we already did it */

	return full;
}
