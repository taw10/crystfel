/*
 * reflist-utils.c
 *
 * Utilities to complement the core reflist.c
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2011-2014 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
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

#define _ISOC99_SOURCE
#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <assert.h>


#include "reflist.h"
#include "cell.h"
#include "cell-utils.h"
#include "utils.h"
#include "reflist-utils.h"
#include "symmetry.h"


/**
 * SECTION:reflist-utils
 * @short_description: Reflection list utilities
 * @title: RefList utilities
 * @section_id:
 * @see_also:
 * @include: "reflist-utils.h"
 * @Image:
 *
 * There are some utility functions associated with the core %RefList.
 **/


int check_list_symmetry(RefList *list, const SymOpList *sym)
{
	Reflection *refl;
	RefListIterator *iter;
	SymOpMask *mask;

	mask = new_symopmask(sym);
	if ( mask == NULL ) {
		ERROR("Couldn't create mask for list symmetry check.\n");
		return 1;
	}

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		int j;
		int found = 0;
		signed int h, k, l;
		int n;

		get_indices(refl, &h, &k, &l);

		special_position(sym, mask, h, k, l);
		n = num_equivs(sym, mask);

		for ( j=0; j<n; j++ ) {

			signed int he, ke, le;
			Reflection *f;

			get_equiv(sym, mask, j, h, k, l, &he, &ke, &le);

			f = find_refl(list, he, ke, le);
			if ( f != NULL ) found++;

		}

		assert(found != 0);  /* That'd just be silly */
		if ( found > 1 ) {

			STATUS("Found %i %i %i: %i times:\n", h, k, l, found);

			for ( j=0; j<n; j++ ) {

				signed int he, ke, le;
				Reflection *f;

				get_equiv(sym, mask, j, h, k, l, &he, &ke, &le);

				f = find_refl(list, he, ke, le);
				if ( f != NULL ) {
					STATUS("%3i %3i %3i\n", he, ke, le);
				}

			}
			free_symopmask(mask);

			return 1;  /* Symmetry is wrong! */
		}

	}

	free_symopmask(mask);

	return 0;
}


int find_equiv_in_list(RefList *list, signed int h, signed int k,
                       signed int l, const SymOpList *sym, signed int *hu,
                       signed int *ku, signed int *lu)
{
	int i;
	int found = 0;

	for ( i=0; i<num_equivs(sym, NULL); i++ ) {

		signed int he, ke, le;
		Reflection *f;
		get_equiv(sym, NULL, i, h, k, l, &he, &ke, &le);
		f = find_refl(list, he, ke, le);

		/* There must only be one equivalent.  If there are more, it
		 * indicates that the user lied about the input symmetry.
		 * This situation should have been checked for earlier by
		 * calling check_symmetry() with 'items' and 'mero'. */

		if ( (f != NULL) && !found ) {
			*hu = he;  *ku = ke;  *lu = le;
			return 1;
		}

	}

	return 0;
}


/**
 * write_reflections_to_file:
 * @fh: File handle to write to
 * @list: The reflection list to write
 *
 * This function writes the contents of @list to @fh,
 *
 * Reflections which have a redundancy of zero will not be written.
 *
 * The resulting list can be read back with read_reflections_from_file().
 **/
void write_reflections_to_file(FILE *fh, RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;

	fprintf(fh, "   h    k    l          I    phase   sigma(I)   nmeas\n");

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{

		signed int h, k, l;
		double intensity, esd_i, ph;
		int red;
		char phs[16];
		int have_phase;

		get_indices(refl, &h, &k, &l);
		intensity = get_intensity(refl);
		esd_i = get_esd_intensity(refl);
		red = get_redundancy(refl);
		ph = get_phase(refl, &have_phase);

		/* Reflections with redundancy = 0 are not written */
		if ( red == 0 ) continue;

		if ( have_phase ) {
			snprintf(phs, 16, "%8.2f", rad2deg(ph));
		} else {
			strncpy(phs, "       -", 15);
		}

		fprintf(fh,
		       "%4i %4i %4i %10.2f %s %10.2f %7i\n",
		       h, k, l, intensity, phs, esd_i, red);

	}
}


/**
 * write_reflist_2:
 * @filename: Filename
 * @list: The reflection list to write
 * @sym: A %SymOpList describing the symmetry of the list
 *
 * This function writes the contents of @list to @file,
 *
 * Reflections which have a redundancy of zero will not be written.
 *
 * The resulting list can be read back with read_reflections_from_file() or
 * read_reflections().
 *
 * Returns: zero on success, non-zero on failure.
 **/
int write_reflist_2(const char *filename, RefList *list, SymOpList *sym)
{
	FILE *fh;
	const char *ssym;

	if ( filename == NULL ) {
		fh = stdout;
	} else {
		fh = fopen(filename, "w");
	}

	if ( fh == NULL ) {
		ERROR("Couldn't open output file '%s'.\n", filename);
		return 1;
	}

	fprintf(fh, "CrystFEL reflection list version 2.0\n");

	if ( sym == NULL ) {
		ssym = "unknown";
	} else {
		ssym = symmetry_name(sym);
	}
	fprintf(fh, "Symmetry: %s\n", ssym);

	write_reflections_to_file(fh, list);
	fprintf(fh, REFLECTION_END_MARKER"\n");

	fclose(fh);

	return 0;
}


/**
 * write_reflist:
 * @filename: Filename
 * @list: The reflection list to write
 *
 * This function writes the contents of @list to @file,
 *
 * Reflections which have a redundancy of zero will not be written.
 *
 * The resulting list can be read back with read_reflections_from_file() or
 * read_reflections().
 *
 * This is a convenience function which simply opens @filename and then calls
 * write_reflections_to_file.
 *
 * Deprecated: use write_reflist_2() instead.
 *
 * Returns: zero on success, non-zero on failure.
 **/
int write_reflist(const char *filename, RefList *list)
{
	return write_reflist_2(filename, list, NULL);
}


#define HEADER_1_0 "  h   k   l          I    phase   sigma(I)  counts  " \
	                  "fs/px  ss/px"

#define HEADER_2_0 "CrystFEL reflection list version 2.0"


/**
 * read_reflections_from_file:
 * @fh: File handle to read from
 *
 * This function reads a reflection list from @fh.
 *
 * Returns: a %RefList read from the file, or NULL on error
 **/
RefList *read_reflections_from_file(FILE *fh)
{
	char *rval = NULL;
	RefList *out;
	int major_version;  /* Minor version as well, but not used yet */
	char line[1024];

	rval = fgets(line, 1023, fh);
	if ( rval == NULL ) return NULL;
	chomp(line);
	if ( strcmp(line, HEADER_1_0) == 0 ) {
		major_version = 1;
	} else if ( strcmp(line, HEADER_2_0) == 0 ) {
		major_version = 2;
	} else {
		fprintf(stderr, "Unrecognised header '%s'\n", line);
		return NULL;
	}

	if ( major_version >= 2 ) {

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) return NULL;
		chomp(line);
		if ( strncmp(line, "Symmetry: ", 10) != 0 ) return NULL;

		/* FIXME: Do something with the symmetry */

		/* Read (and ignore) the header */
		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) return NULL;

	}

	out = reflist_new();

	do {

		Reflection *refl;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		if ( strcmp(line, REFLECTION_END_MARKER) == 0 ) return out;

		if ( major_version >= 2 ) {

			double ph;
			char *v;
			signed int h, k, l;
			float intensity, sigma;
			char phs[1024];
			int cts;
			int r;

			r = sscanf(line, "%i %i %i %f %s %f %i",
				   &h, &k, &l, &intensity, phs, &sigma, &cts);

			if ( r != 7 ) {
				reflist_free(out);
				printf("Bad line '%s'\n", line);
				return NULL;
			}

			refl = add_refl(out, h, k, l);
			set_intensity(refl, intensity);
			set_esd_intensity(refl, sigma);
			set_redundancy(refl, cts);

			ph = strtod(phs, &v);
			if ( v != phs ) set_phase(refl, deg2rad(ph));

		} else {

			/* Deprecated reflection format */

			double ph;
			char *v;
			signed int h, k, l;
			float intensity, sigma, fs, ss;
			char phs[1024];
			int cts;
			int r;

			r = sscanf(line, "%i %i %i %f %s %f %i %f %f",
				   &h, &k, &l, &intensity, phs, &sigma,
				   &cts, &fs, &ss);

			if ( r != 9 ) {
				reflist_free(out);
				return NULL;
			}

			refl = add_refl(out, h, k, l);
			set_intensity(refl, intensity);
			set_detector_pos(refl, fs, ss);

			set_esd_intensity(refl, sigma);
			set_redundancy(refl, cts);

			ph = strtod(phs, &v);
			if ( v != phs ) set_phase(refl, deg2rad(ph));

		}

	} while ( rval != NULL );

	/* Got read error of some kind before finding PEAK_LIST_END_MARKER */
	return NULL;
}


RefList *read_reflections(const char *filename)
{
	FILE *fh;
	RefList *out;

	if ( filename == NULL ) {
		fh = stdout;
	} else {
		fh = fopen(filename, "r");
	}

	if ( fh == NULL ) {
		ERROR("Couldn't open input file '%s'.\n", filename);
		return NULL;
	}

	out = read_reflections_from_file(fh);

	fclose(fh);

	return out;
}


/**
 * asymmetric_indices:
 * @in: A %RefList
 * @sym: A %SymOpList
 *
 * This function creates a newly allocated copy of @in, but indexed using the
 * asymmetric indices according to @sym instead of the original indices.  The
 * original indices are stored and can be retrieved using
 * get_symmetric_indices() if required.
 *
 * Returns: the new %RefList, or NULL on failure.
 **/
RefList *asymmetric_indices(RefList *in, const SymOpList *sym)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *new;

	new = reflist_new();
	if ( new == NULL ) return NULL;

	for ( refl = first_refl(in, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		signed int ha, ka, la;
		Reflection *cr;

		get_indices(refl, &h, &k, &l);

		get_asymm(sym, h, k, l, &ha, &ka, &la);

		cr = add_refl(new, ha, ka, la);
		assert(cr != NULL);

		copy_data(cr, refl);
		set_symmetric_indices(cr, h, k, l);

	}

	return new;
}


/**
 * resolution_limits:
 * @list: A %RefList
 * @cell: A %UnitCell
 * @rmin: Place to store the minimum 1/d value
 * @rmax: Place to store the maximum 1/d value
 *
 * This function calculates the minimum and maximum values of 1/d, where
 * 2dsin(theta) = wavelength.  The answers are in m^-1.
 **/
void resolution_limits(RefList *list, UnitCell *cell,
                       double *rmin, double *rmax)
{
	Reflection *refl;
	RefListIterator *iter;

	*rmin = INFINITY;
	*rmax = 0.0;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double r;
		signed int h, k, l;

		get_indices(refl, &h, &k, &l);
		r = 2.0 * resolution(cell, h, k, l);

		if ( r > *rmax ) *rmax = r;
		if ( r < *rmin ) *rmin = r;
	}
}


/**
 * max_intensity:
 * @list: A %RefList
 *
 * Returns: The maximum intensity in @list.
 **/
double max_intensity(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;
	double max;

	max = -INFINITY;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double val = get_intensity(refl);
		if ( val > max ) max = val;
	}

	return max;
}


/**
 * res_cutoff:
 * @list: A %RefList
 * @cell: A %UnitCell with which to calculate 1/d values for @list
 * @min: Minimum acceptable value of 1/d
 * @max: Maximum acceptable value of 1/d
 *
 * Applies a resolution cutoff to @list, returning the new version and freeing
 * the old version.
 *
 * Returns: A new %RefList with resolution cutoff applied
 **/
RefList *res_cutoff(RefList *list, UnitCell *cell, double min, double max)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *new;

	new = reflist_new();

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double one_over_d;
		signed int h, k, l;
		Reflection *n;

		get_indices(refl, &h, &k, &l);

		one_over_d = 2.0 * resolution(cell, h, k, l);
		if ( one_over_d < min ) continue;
		if ( one_over_d > max ) continue;

		n = add_refl(new, h, k, l);
		copy_data(n, refl);
	}

	reflist_free(list);
	return new;
}


/**
 * copy_reflist:
 * @list: A %RefList
 *
 * Returns: A copy of %RefList.
 **/
RefList *copy_reflist(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *new;

	new = reflist_new();

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		Reflection *n;

		get_indices(refl, &h, &k, &l);

		n = add_refl(new, h, k, l);
		copy_data(n, refl);
	}

	return new;
}
