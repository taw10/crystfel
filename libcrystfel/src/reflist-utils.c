/*
 * reflist-utils.c
 *
 * Utilities to complement the core reflist.c
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2011-2020 Thomas White <taw@physics.org>
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "reflist.h"
#include "cell.h"
#include "cell-utils.h"
#include "utils.h"
#include "reflist-utils.h"
#include "symmetry.h"
#include "libcrystfel-version.h"


/** \file reflist-utils.h
 */


/**
 * Checks that the symmetry of \p list is indeed \p sym.
 *
 * \param list A list of reflections
 * \param sym Symmetry of the reflection list
 *
 * \returns 0 if the symmetry is correct, otherwise 1
 */
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


/*
 * Write the actual reflections to the file, no headers etc.
 * Reflections which have a redundancy of zero will not be written.
 * The resulting list can be read back with read_reflections_from_file().
 **/
static void write_reflections_to_file(FILE *fh, RefList *list)
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
 * This function writes the contents of list to file,
 *
 * Reflections which have a redundancy of zero will not be written.
 *
 * The resulting list can be read back with read_reflections_from_file() or
 * read_reflections().
 *
 * \param filename Filename
 * \param list The reflection list to write
 * \param sym A %SymOpList describing the symmetry of the list
 *
 * \returns zero on success, non-zero on failure.
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

	if ( reflist_get_notes(list) != NULL ) {
		fprintf(fh, "%s\n", reflist_get_notes(list));
	}

	fclose(fh);

	return 0;
}


/**
 * This function writes the contents of list to file,
 *
 * Reflections which have a redundancy of zero will not be written.
 *
 * The resulting list can be read back with read_reflections_from_file() or
 * read_reflections().
 *
 * This is a convenience function which simply opens filename and then calls
 * write_reflections_to_file.
 *
 * \deprecated Use write_reflist_2() instead.
 *
 * \param filename Filename
 * \param list The reflection list to write
 *
 * \returns Zero on success, non-zero on failure.
 **/
int write_reflist(const char *filename, RefList *list)
{
	return write_reflist_2(filename, list, NULL);
}


#define HEADER_1_0 "  h   k   l          I    phase   sigma(I)  counts  " \
	                  "fs/px  ss/px"

#define HEADER_2_0 "CrystFEL reflection list version 2.0"


/* fh: File handle to read from
 * sym: Location at which to store pointer to symmetry, or NULL if you don't
 *  need it
 *
 * Returns: a %RefList read from the file, or NULL on error
 */
static RefList *read_reflections_from_file(FILE *fh, char **sym)
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

		if ( sym != NULL ) {
			*sym = strdup(line+10);
		}

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

		if ( strcmp(line, REFLECTION_END_MARKER) == 0 ) break;

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

	if ( strcmp(line, REFLECTION_END_MARKER) != 0 ) {
		/* Got read error of some kind before finding
		 * PEAK_LIST_END_MARKER */
		return NULL;
	}

	/* We are now in the notes region */
	do {
		char line[1024];
		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);
		reflist_add_notes(out, line);
	} while ( rval != NULL );

	return out;
}


/**
 * read_reflections_2:
 * \param filename: Filename to read from
 * \param sym: Pointer to a "char *" at which to store the symmetry
 *
 * This function reads a reflection list from a file, including the
 * symmetry from the header (e.g. "Symmetry: 4/mmm").
 *
 * Returns: A %RefList read from the file, or NULL on error
 */
RefList *read_reflections_2(const char *filename, char **sym)
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

	out = read_reflections_from_file(fh, sym);

	fclose(fh);

	return out;
}


/**
 * read_reflections:
 * \param filename: Filename to read from
 *
 * This function reads a reflection list from a file.
 *
 * Returns: A %RefList read from the file, or NULL on error
 */
RefList *read_reflections(const char *filename)
{
	return read_reflections_2(filename, NULL);
}


/**
 * asymmetric_indices:
 * \param in: A %RefList
 * \param sym: A %SymOpList
 *
 * This function creates a newly allocated copy of in, but indexed using the
 * asymmetric indices according to sym instead of the original indices.  The
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
 * \param list: A %RefList
 * \param cell: A %UnitCell
 * \param rmin: Place to store the minimum 1/d value
 * \param rmax: Place to store the maximum 1/d value
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
 * \param list: A %RefList
 *
 * Returns: The maximum intensity in \p list.
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
 * \param list: A %RefList
 * \param cell: A %UnitCell with which to calculate 1/d values for \p list
 * \param min: Minimum acceptable value of 1/d
 * \param max: Maximum acceptable value of 1/d
 *
 * Applies a resolution cutoff to \p list, returning the new version and freeing
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
 * \param list: A %RefList
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


/**
 * free_contribs:
 * \param list: A %RefList
 *
 * Goes through \p list and frees all the reflection contribution structures.
 **/
void free_contribs(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		struct reflection_contributions *c;
		c = get_contributions(refl);
		free(c->contribs);
		free(c->contrib_crystals);
		free(c);
	}
}


static char *full_command_line(int argc, char *argv[])
{
	int i;
	size_t len = 1;
	char *cl;

	if ( argc == 0 ) return strdup("");
	for ( i=0; i<argc; i++ ) {
		len += strlen(argv[i]) + 1;
	}

	cl = malloc(len);
	if ( cl == NULL ) return strdup("");

	cl[0] = '\0';
	for ( i=0; i<argc; i++ ) {
		if ( i > 0 ) strcat(cl, " ");
		strcat(cl, argv[i]);
	}

	return cl;
}


void reflist_add_command_and_version(RefList *list, int argc, char *argv[])
{
	char *tmp;
	char vers[128];

	snprintf(vers, 128, "Generated by CrystFEL %s",
	         libcrystfel_version_string());
	reflist_add_notes(list, vers);

	tmp = full_command_line(argc, argv);
	reflist_add_notes(list, tmp);
	free(tmp);
}
