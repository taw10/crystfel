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


void scale_intensities(const double *model, ReflItemList *model_items,
                       double *new_pattern, ReflItemList *new_items,
                       double f0, int f0_valid)
{
	double s;
	unsigned int i;
	ReflItemList *items;

	items = intersection_items(model_items, new_items);
	s = stat_scale_intensity(model, new_pattern, items);
	delete_items(items);
	if ( f0_valid ) printf("%f %f\n", s, f0);

	/* NaN -> abort */
	if ( isnan(s) ) return;

	/* Multiply the new pattern up by "s" */
	for ( i=0; i<LIST_SIZE; i++ ) {
		new_pattern[i] *= s;
	}
}
