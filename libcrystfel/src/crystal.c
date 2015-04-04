/*
 * crystal.c
 *
 * A class representing a single crystal
 *
 * Copyright © 2013-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2013-2015 Thomas White <taw@physics.org>
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
	/* The image containing the crystal */
	struct image            *image;

	/* Information about the crystal */
	UnitCell                *cell;
	double                  m;     /* Mosaicity in radians */
	double                  osf;
	double                  Bfac;
	double                  profile_radius;
	int                     pr_dud;
	double                  resolution_limit;

	/* Integrated (or about-to-be-integrated) reflections */
	RefList                 *reflections;
	long long int           n_saturated;  /* Number of overloads */
	long long int           n_implausible;  /* Number of implausibly
	                                         * negative reflectionss */

	/* User flag, e.g. for "this is a bad crystal". */
	int                     user_flag;

	/* Text notes, which go in the stream */
	char                    *notes;
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
	cryst->n_implausible = 0;
	cryst->notes = NULL;
	cryst->user_flag = 0;

	return cryst;
}


/**
 * crystal_copy:
 * @cryst: A %Crystal to copy.
 *
 * Creates a new %Crystal which is a copy of @cryst.  The copy is a "shallow
 * copy", which means that copies are NOT made of the data structures which
 * @cryst contains references to, for example its %RefList.
 *
 * Returns: a (shallow) copy of @cryst, or NULL on failure.
 *
 */
Crystal *crystal_copy(Crystal *cryst)
{
	Crystal *c;

	c = crystal_new();
	if ( c == NULL ) return NULL;

	memcpy(c, cryst, sizeof(Crystal));

	return c;
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


long long int crystal_get_num_implausible_reflections(Crystal *cryst)
{
	return cryst->n_implausible;
}


struct image *crystal_get_image(Crystal *cryst)
{
	return cryst->image;
}


double crystal_get_osf(Crystal *cryst)
{
	return cryst->osf;
}


double crystal_get_Bfac(Crystal *cryst)
{
	return cryst->Bfac;
}


int crystal_get_user_flag(Crystal *cryst)
{
	return cryst->user_flag;
}


double crystal_get_mosaicity(Crystal *cryst)
{
	return cryst->m;
}


const char *crystal_get_notes(Crystal *cryst)
{
	return cryst->notes;
}


/********************************** Setters ***********************************/


void crystal_set_cell(Crystal *cryst, UnitCell *cell)
{
	cryst->cell = cell;
}


void crystal_set_profile_radius(Crystal *cryst, double r)
{
	cryst->profile_radius = r;
}


void crystal_set_reflections(Crystal *cryst, RefList *reflist)
{
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


void crystal_set_num_implausible_reflections(Crystal *cryst, long long int n)
{
	cryst->n_implausible = n;
}


void crystal_set_image(Crystal *cryst, struct image *image)
{
	cryst->image = image;
}


void crystal_set_osf(Crystal *cryst, double osf)
{
	cryst->osf = osf;
}


void crystal_set_Bfac(Crystal *cryst, double Bfac)
{
	cryst->Bfac = Bfac;
}


void crystal_set_user_flag(Crystal *cryst, int user_flag)
{
	cryst->user_flag = user_flag;
}


void crystal_set_mosaicity(Crystal *cryst, double m)
{
	cryst->m = m;
}


void crystal_set_notes(Crystal *cryst, const char *notes)
{
	free(cryst->notes);  /* free(NULL) is OK */
	cryst->notes = strdup(notes);
}


void crystal_add_notes(Crystal *cryst, const char *notes_add)
{
	size_t len;
	char *nnotes;

	if ( cryst->notes == NULL ) {
		crystal_set_notes(cryst, notes_add);
		return;
	}

	len = strlen(notes_add) + strlen(cryst->notes) + 2;
	nnotes = malloc(len);
	if ( nnotes == NULL ) {
		ERROR("Failed to add notes to crystal.\n");
		return;
	}

	strcpy(nnotes, cryst->notes);
	strcat(nnotes, "\n");
	strcat(nnotes, notes_add);
	free(cryst->notes);
	cryst->notes = nnotes;
}
