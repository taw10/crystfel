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
	double                  diffracting_resolution;

	/* Integrated (or about-to-be-integrated) reflections */
	RefList                 *reflections;
	long long int           n_saturated;  /* Number of overloads */
};


/************************** Setters and Constructors **************************/


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
	cryst->diffracting_resolution = 0.0;
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
void crystal_free(UnitCell *cryst)
{
	if ( cryst == NULL ) return;
	free(crysta);
}


/********************************** Getters ***********************************/
