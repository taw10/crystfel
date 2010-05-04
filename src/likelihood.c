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
                        unsigned int *new_counts)
{
	/* Placeholder... */
}

void scale_intensities(const double *model, double *new_pattern,
                       const unsigned int *model_counts,
                       unsigned int *new_counts)
{
	double s;
	unsigned int i;

	s = stat_scale_intensity(model, model_counts, new_pattern, new_counts);
	printf("%f\n", s);

	/* NaN -> abort */
	if ( isnan(s) ) return;

	/* Multiply the new pattern up by "s" */
	for ( i=0; i<LIST_SIZE; i++ ) {
		new_counts[i] *= s;
	}
}
