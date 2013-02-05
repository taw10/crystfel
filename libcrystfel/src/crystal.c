/*
 * crystal.c
 *
 * A class representing a single crystal
 *
 * Copyright © 2013 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2013 Thomas White <taw@physics.org>
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

#include "crystal.h"
#include "utils.h"


/**
 * SECTION:crystal
 * @short_description: Crystal
 * @title: Crystal
 * @section_id:
 * @see_also:
 * @include: "crystal.h"
 * @Image:
 *
 * This structure represents a single crystal.
 */


struct _crystal
{
	/* Information about the crystal */
	UnitCell                *cell;
	double                  m;     /* Mosaicity in radians */
	double                  osf;
	double                  profile_radius;
	int                     pr_dud;
	double                  resolution_limit;

	/* Integrated (or about-to-be-integrated) reflections */
	RefList                 *reflections;
	long long int           n_saturated;  /* Number of overloads */
};


/************************** Constructor / Destructor **************************/


/**
 * crystal_new:
 *
 * Create a new %Crystal.
 *
 * Returns: the new unit cell, or NULL on failure.
 *
 */
Crystal *crystal_new()
{
	Crystal *cryst;

	cryst = malloc(sizeof(Crystal));
	if ( cryst == NULL ) return NULL;

	cryst->cell = NULL;
	cryst->reflections = NULL;
	cryst->resolution_limit = 0.0;
	cryst->n_saturated = 0;

	return cryst;
}


/**
 * crystal_free:
 * @cryst: A %Crystal to free.
 *
 * Frees a %Crystal, and all internal resources concerning that crystal.
 *
 */
void crystal_free(Crystal *cryst)
{
	if ( cryst == NULL ) return;
	if ( cryst->cell != NULL ) cell_free(cryst->cell);
	if ( cryst->reflections != NULL ) reflist_free(cryst->reflections);
	free(cryst);
}


/********************************** Getters ***********************************/


UnitCell *crystal_get_cell(Crystal *cryst)
{
	return cryst->cell;
}


double crystal_get_profile_radius(Crystal *cryst)
{
	return cryst->profile_radius;
}


RefList *crystal_get_reflections(Crystal *cryst)
{
	return cryst->reflections;
}


double crystal_get_resolution_limit(Crystal *cryst)
{
	return cryst->resolution_limit;
}


long long int crystal_get_num_saturated_reflections(Crystal *cryst)
{
	return cryst->n_saturated;
}


/********************************** Setters ***********************************/


void crystal_set_cell(Crystal *cryst, UnitCell *cell)
{
	if ( cryst->cell != NULL ) cell_free(cryst->cell);
	cryst->cell = cell_new_from_cell(cell);
}


void crystal_set_profile_radius(Crystal *cryst, double r)
{
	cryst->profile_radius = r;
}


void crystal_set_reflections(Crystal *cryst, RefList *reflist)
{
	if ( cryst->reflections != NULL ) reflist_free(reflist);
	cryst->reflections = reflist;
}


void crystal_set_resolution_limit(Crystal *cryst, double res)
{
	cryst->resolution_limit = res;
}


void crystal_set_num_saturated_reflections(Crystal *cryst, long long int n)
{
	cryst->n_saturated = n;
}
