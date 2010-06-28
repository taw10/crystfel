/*
 * likelihood.c
 *
 * Likelihood maximisation
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "statistics.h"
#include "utils.h"

void detwin_intensities(const double *model, double *new_pattern,
                        const unsigned int *model_counts,
                        ReflItemList *items)
{
	/* Placeholder... */
}

void scale_intensities(const double *model, double *new_pattern,
                       const unsigned int *model_counts,
                       ReflItemList *items, double f0, int f0_valid)
{
	double s;
	unsigned int i;
	unsigned int *new_counts;

	new_counts = items_to_counts(items);

	s = stat_scale_intensity(model, model_counts, new_pattern, new_counts);
	if ( f0_valid ) printf("%f %f\n", s, f0);

	/* NaN -> abort */
	if ( isnan(s) ) return;

	/* Multiply the new pattern up by "s" */
	for ( i=0; i<LIST_SIZE; i++ ) {
		new_counts[i] *= s;
	}

	free(new_counts);
}
