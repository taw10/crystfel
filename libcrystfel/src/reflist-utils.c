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

#include <libcrystfel-config.h>

#include <math.h>
#include <stdio.h>
#include <assert.h>

#ifdef HAVE_LIBCCP4
#include <ccp4/cmtzlib.h>
#include <ccp4/csymlib.h>
#include <ccp4/ccp4_parser.h>
#endif

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


struct point_group_conversion
{
	char centering;
	const char *crystfel;
	int friedel;

	int xds_spgnum;

	const char *ccp4;
};


/* Table for converting CrystFEL's knowledge about centering, point group
 * and unique axis into something that can be recognised by external programs.
 * If xds_spgnum=0, ccp4=NULL, or something is missing form this table entirely,
 * it means that there is no way to represent the situation to that program
 * without re-indexing the dataset. */
struct point_group_conversion pg_conversions[] = {

	/* Triclinic  */
	{'P', "1",       0,      1,     "P 1"},
	{'P', "-1",      1,      1,     "P 1"},

	/* Monoclinic */
	{'P', "2_uaa",   0,      0,     "P211"},
	{'P', "m_uaa",   0,      0,     "Pm11"},
	{'P', "2/m_uaa", 1,      0,     "P211"},
	{'P', "2_uab",   0,      3,     "P121"},
	{'P', "m_uab",   0,      0,     "P1m1"},
	{'P', "2/m_uab", 1,      3,     "P121"},
	{'P', "2_uac",   0,      0,     "P112"},
	{'P', "m_uac",   0,      0,     "P11m"},
	{'P', "2/m_uac", 1,      0,     "P112"},
	{'P', "2",       0,      0,     "P121"}, /* unique axis c */
	{'P', "m",       0,      0,     "P11m"}, /* unique axis c */
	{'P', "2/m",     1,      0,     "P121"}, /* unique axis c */

	{'A', "2_uab",   0,      0,     "A121"},
	{'A', "m_uab",   0,      0,     "A1m1"},
	{'A', "2/m_uab", 1,      0,     "A121"},
	{'A', "2_uac",   0,      0,     "A112"},
	{'A', "m_uac",   0,      0,     "A11m"},
	{'A', "2/m_uac", 1,      0,     "A112"},
	{'A', "2",       0,      0,     "A121"}, /* unique axis c */
	{'A', "m",       0,      0,     "A11m"}, /* unique axis c */
	{'A', "2/m",     1,      0,     "A121"}, /* unique axis c */

	{'B', "2_uaa",   0,      0,     "B211"},
	{'B', "m_uaa",   0,      0,     "Bm11"},
	{'B', "2/m_uaa", 1,      0,     "B211"},
	{'B', "2_uac",   0,      0,     "B112"},
	{'B', "m_uac",   0,      0,     "B11m"},
	{'B', "2/m_uac", 1,      0,     "B112"},
	{'B', "2",       0,      0,     "B112"}, /* unique axis c */
	{'B', "m",       0,      0,     "B11m"}, /* unique axis c */
	{'B', "2/m",     1,      0,     "B112"}, /* unique axis c */

	{'C', "2_uaa",   0,      0,     "C211"},
	{'C', "m_uaa",   0,      0,     "Cm11"},
	{'C', "2/m_uaa", 1,      0,     "C211"},
	{'C', "2_uab",   0,      5,     "C121"},
	{'C', "m_uab",   0,      0,     "C1m1"},
	{'C', "2/m_uab", 1,      5,     "C121"},

	{'I', "2_uaa",   0,      0,     "I211"},
	{'I', "m_uaa",   0,      0,     "Im11"},
	{'I', "2/m_uaa", 1,      0,     "I211"},
	{'I', "2_uab",   0,      0,     "I121"},
	{'I', "m_uab",   0,      0,     "I1m1"},
	{'I', "2/m_uab", 1,      0,     "I121"},
	{'I', "2_uac",   0,      0,     "I112"},
	{'I', "m_uac",   0,      0,     "I11m"},
	{'I', "2/m_uac", 1,      0,     "I112"},
	{'I', "2",       0,      0,     "I121"}, /* unique axis c */
	{'I', "m",       0,      0,     "I11m"}, /* unique axis c */
	{'I', "2/m",     1,      0,     "I121"}, /* unique axis c */

	/* Orthorhombic */
	{'P', "222",       0,     16,     "P222"},
	{'P', "mmm",       1,     16,     "P222"},
	{'P', "mm2",       0,     25,     "Pmm2"},
	{'A', "222",       0,      0,     "A222"},
	{'A', "mmm",       1,      0,     "A222"},
	{'A', "mm2",       0,     38,     "Amm2"},
	{'B', "222",       0,      0,     "B222"},
	{'B', "mmm",       1,      0,     "B222"},
	{'B', "mm2",       0,      0,     "Bmm2"},
	{'C', "222",       0,     21,     "C222"},
	{'C', "mmm",       1,     21,     "C222"},
	{'C', "mm2",       0,     35,     "Cmm2"},
	{'F', "222",       0,     22,     "F222"},
	{'F', "mmm",       1,     22,     "F222"},
	{'F', "mm2",       0,     42,     "Fmm2"},
	{'I', "222",       0,     23,     "I222"},
	{'I', "mmm",       1,     23,     "I222"},
	{'I', "mm2",       0,     45,     "Imm2"},

	/* Tetragonal */
	{'P', "4",         0,     75,     "P4"},    /* unique axis c */
	{'P', "4/m",       1,     75,     "P4"},    /* unique axis c */
	{'P', "422",       0,     89,     "P422"},  /* unique axis c */
	{'P', "4/mmm",     1,     89,     "P422"},  /* unique axis c */
	{'P', "4mm",       0,     99,     "P4mm"},  /* unique axis c */
	{'P', "-4",        0,     81,     "P-4"},   /* unique axis c */
	{'P', "-42m",      0,    111,     "P-42m"}, /* unique axis c */
	{'P', "-4m2",      0,    115,     "P-4m2"}, /* unique axis c */
	{'P', "4_uac",     0,     75,     "P4"},
	{'P', "4/m_uac",   1,     75,     "P4"},
	{'P', "422_uac",   0,     89,     "P422"},
	{'P', "4/mmm_uac", 1,     89,     "P422"},
	{'P', "4mm_uac",   0,     99,     "P4mm"},
	{'P', "-4_uac",    0,     81,     "P-4"},
	{'P', "-42m_uac",  0,    111,     "P-42m"},
	{'P', "-4m2_uac",  0,    115,     "P-4m2"},
	{'I', "4",         0,     79,     "I4"},    /* unique axis c */
	{'I', "4/m",       1,     79,     "I4"},    /* unique axis c */
	{'I', "422",       0,     97,     "I422"},  /* unique axis c */
	{'I', "4/mmm",     1,     97,     "I422"},  /* unique axis c */
	{'I', "4mm",       0,    107,     "I4mm"},  /* unique axis c */
	{'I', "-4",        0,     82,     "I-4"},   /* unique axis c */
	{'I', "-42m",      0,    121,     "I-42m"}, /* unique axis c */
	{'I', "-4m2",      0,    119,     "I-4m2"}, /* unique axis c */
	{'I', "4_uac",     0,     79,     "I4"},
	{'I', "4/m_uac",   1,     79,     "I4"},
	{'I', "422_uac",   0,     97,     "I422"},
	{'I', "4/mmm_uac", 1,     97,     "I422"},
	{'I', "4mm_uac",   0,    107,     "I4mm"},
	{'I', "-4_uac",    0,     82,     "I-4"},
	{'I', "-42m_uac",  0,    121,     "I-42m"},
	{'I', "-4m2_uac",  0,    119,     "I-4m2"},

	/* Trigonal (rhombohedral) */
	{'R', "3_R",       0,      0,     "R3:R"},
	{'R', "-3_R",      1,      0,     "R3:R"},
	{'R', "32_R",      0,      0,     "R32:R"},
	{'R', "-3m_R",     1,      0,     "R32:R"},
	{'R', "3m_R",      0,      0,     "R3m:R"},

	/* Trigonal (rhombohedral on hexagonal axes) */
	{'H', "3_H",       0,    146,     "R3:H"},
	{'H', "-3_H",      1,    146,     "R3:H"},
	{'H', "32_H",      0,    155,     "R3:H"},
	{'H', "-3m_H",     1,    155,     "R3:H"},
	{'H', "3m_H",      0,      0,     "R3m:H"},

	/* Trigonal (hexagonal) */
	{'P', "3_H",       0,    143,     "P3"},
	{'P', "-3_H",      1,    143,     "P3"},
	{'P', "312_H",     0,    149,     "P312"},
	{'P', "-31m_H",    1,    149,     "P312"},
	{'P', "321_H",     0,    150,     "P321"},
	{'P', "-3m1_H",    1,    150,     "P321"},
	{'P', "3m1_H",     0,    156,     "P3m1"},
	{'P', "31m_H",     0,    157,     "P31m"},

	/* Hexagonal */
	{'P', "6",         0,    168,     "P6"},
	{'P', "6/m",       1,    168,     "P6"},
	{'P', "622",       0,    177,     "P622"},
	{'P', "6/mmm",     1,    177,     "P622"},
	{'P', "6mm",       0,    177,     "P6mm"},
	{'P', "-6m2",      0,    187,     "P-6m2"},
	{'P', "-62m",      0,    189,     "P-62m"},

	/* Cubic */
	{'P', "23",        0,    195,     "P23"},
	{'P', "m-3",       1,    195,     "P23"},
	{'P', "432",       0,    207,     "P432"},
	{'P', "m-3m",      1,    207,     "P432"},
	{'P', "-43m",      0,    215,     "P -4 3 m"},
	{'I', "23",        0,    197,     "I23"},
	{'I', "m-3",       1,    197,     "I23"},
	{'I', "432",       0,    211,     "I432"},
	{'I', "m-3m",      1,    211,     "I432"},
	{'I', "-43m",      0,    217,     "I -4 3 m"},
	{'F', "23",        0,    196,     "F23"},
	{'F', "m-3",       1,    196,     "F23"},
	{'F', "432",       0,    209,     "F432"},
	{'F', "m-3m",      1,    209,     "F432"},
	{'F', "-43m",      0,    216,     "F -4 3 m"},

	{'*', NULL,  0, 0, NULL}
};


static int space_group_for_xds(const char *sym_str, char cen)
{
	int i = 0;
	do {
		if ( (pg_conversions[i].centering == cen)
		  && (strcmp(sym_str, pg_conversions[i].crystfel) == 0) )
		{
			return pg_conversions[i].xds_spgnum;
		}
		i++;
	} while (pg_conversions[i].centering != '*');

	ERROR("Couldn't derive XDS representation of symmetry.\n");
	return 0;
}


#ifdef HAVE_LIBCCP4
static const char *space_group_for_mtz(const char *sym_str, char cen)
{
	int i = 0;
	do {
		if ( (pg_conversions[i].centering == cen)
		  && (strcmp(sym_str, pg_conversions[i].crystfel) == 0) )
		{
			return pg_conversions[i].ccp4;
		}
		i++;
	} while (pg_conversions[i].centering != '*');

	ERROR("Couldn't derive CCP4 representation of symmetry.\n");
	return NULL;
}
#endif


int write_to_xds(RefList *reflist,
                 SymOpList *sym,
                 UnitCell *cell,
                 double min_res,
                 double max_res,
                 const char *filename)
{
	FILE *fh;
	RefListIterator *iter;
	Reflection *refl;
	double a, b, c, al, be,ga;
	int spg;

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);

	spg = space_group_for_xds(symmetry_name(sym),
	                          cell_get_centering(cell));
	if ( spg == 0 ) return 1;

	fh = fopen(filename, "w");
	if ( fh == NULL ) return 1;

	fprintf(fh, "!FORMAT=XDS_ASCII MERGE=TRUE FRIEDEL'S_LAW=%s\n",
	        is_centrosymmetric(sym) ? "TRUE" : "FALSE");
	fprintf(fh, "!SPACE_GROUP_NUMBER=%i\n", spg);
	fprintf(fh, "!UNIT_CELL_CONSTANT= %.2f %.2f %.2f %.2f %.2f %.2f\n",
	        a*1e10, b*1e10, c*1e10, rad2deg(al), rad2deg(be), rad2deg(ga));
	fprintf(fh, "!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=5\n");
	fprintf(fh, "!ITEM_H=1\n");
	fprintf(fh, "!ITEM_K=2\n");
	fprintf(fh, "!ITEM_L=3\n");
	fprintf(fh, "!ITEM_IOBS=4\n");
	fprintf(fh, "!ITEM_SIGMA(IOBS)=5\n");
	fprintf(fh, "!END_OF_HEADER\n");

	for ( refl = first_refl(reflist, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double one_over_d;

		get_indices(refl, &h, &k, &l);

		one_over_d = 2.0*resolution(cell, h, k, l);
		if ( (one_over_d > min_res) && (one_over_d < max_res) ) {

			fprintf(fh, "%6i %6i %6i %9.2f %9.2f\n",
			        h, k, l,
			        get_intensity(refl),
			        get_esd_intensity(refl));

		}
	}

	fprintf(fh, "!END_OF_DATA\n");

	fclose(fh);
	return 0;
}


#ifdef HAVE_LIBCCP4
static CCP4SPG *add_mtz_symmetry_header(MTZ *mtz, const char *spg_name)
{
	CCP4SPG *spg;
	float rsymx[192][4][4];
	char ltypex[2];
	int i;

	spg = ccp4spg_load_by_spgname(spg_name);
	if ( spg == NULL ) {
		ERROR("Couldn't look up CCP4 space group '%s'\n", spg_name);
		return NULL;
	}

	for ( i=0; i<spg->nsymop; i++ ) {
		rotandtrn_to_mat4(rsymx[i], spg->symop[i]);
	}
	ltypex[0] = spg->symbol_old[0];
	ltypex[1] = '\0';

	ccp4_lwsymm(mtz, spg->nsymop, spg->nsymop_prim,
	            rsymx, ltypex, spg->spg_ccp4_num, spg->symbol_old,
	            spg->point_group);

	return spg;
}
#endif


int write_to_mtz(RefList *reflist,
                 SymOpList *sym,
                 UnitCell *cell,
                 double min_res,
                 double max_res,
                 const char *filename,
                 const char *dataset_name)
{
#ifdef HAVE_LIBCCP4
	MTZ *mtz;
	MTZXTAL *cr;
	MTZSET *ds;
	MTZCOL *columns[7];
	double a, b, c, al, be, ga;
	int r;
	char tmp[128];
	float cellp[6];
	int refl_i;
	Reflection *refl;
	RefListIterator *iter;
	CCP4SPG *spg;
	const char *spg_name;

	spg_name = space_group_for_mtz(symmetry_name(sym),
	                               cell_get_centering(cell));
	if ( spg_name == NULL ) {
		reflist_free(reflist);
		return 1;
	}

	mtz = MtzMalloc(0, 0);

	snprintf(tmp, 128, "Data exported via CrystFEL version %s",
	         libcrystfel_version_string());
	ccp4_lwtitl(mtz, tmp, 0);

	mtz->refs_in_memory = 0;
	mtz->fileout = MtzOpenForWrite(filename);

	spg = add_mtz_symmetry_header(mtz, spg_name);
	if ( spg == NULL ) {
		return 1;
	}

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
	cellp[0] = a*1e10;
	cellp[1] = b*1e10;
	cellp[2] = c*1e10;
	cellp[3] = rad2deg(al);
	cellp[4] = rad2deg(be);
	cellp[5] = rad2deg(ga);

	/* FIXME: Proposed labelling:
	 *  title = as above
	 *  project = basename of folder containing crystfel.project
	 *  crystal = name of indexing results run
	 *  dataset = name of merge results run */
	cr = MtzAddXtal(mtz, "Crystal_name", "Project_name", cellp);
	ds = MtzAddDataset(mtz, cr, dataset_name, 0.0);
	columns[0] = MtzAddColumn(mtz, ds, "H", "H");
	columns[1] = MtzAddColumn(mtz, ds, "K", "H");
	columns[2] = MtzAddColumn(mtz, ds, "L", "H");
	columns[3] = MtzAddColumn(mtz, ds, "I+", "J");
	columns[4] = MtzAddColumn(mtz, ds, "SIGI+", "Q");
	columns[5] = MtzAddColumn(mtz, ds, "I-", "J");
	columns[6] = MtzAddColumn(mtz, ds, "SIGI-", "Q");

	refl_i = 1;
	for ( refl = first_refl(reflist, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double one_over_d;
		int isym;

		get_indices(refl, &h, &k, &l);

		one_over_d = 2.0*resolution(cell, h, k, l);
		if ( (one_over_d > min_res) && (one_over_d < max_res) ) {

			float refldata[7];
			signed int nh, nk, nl;
			signed int fh, fk, fl;
			Reflection *friedel;
			Reflection *refl_plus;
			Reflection *refl_minus;

			/* Look for Friedel partner */
			if ( find_equiv_in_list(reflist, -h, -k, -l,
			                        sym, &fh, &fk, &fl) )
			{
				friedel = find_refl(reflist, fh, fk, fl);
			} else {
				friedel = NULL;
			}

			/* Move to CCP4's idea of the ASU */
			isym = ccp4spg_put_in_asu(spg, h, k, l, &nh, &nk, &nl);

			/* Ok, do we have an I+ or an I- ? */
			if ( is_odd(isym) ) {
				/* I+ */
				refl_plus = refl;
				refl_minus = friedel;
			} else {
				/* I- */
				refl_minus = refl;
				refl_plus = friedel;
			}

			/* If we are looking at an I-, only write it out now
			 * if the corresponding I+ if not in 'reflist'.
			 * If I+ is present, then this I- will get written when
			 * the Friedel pair is processed. */
			if ( !is_odd(isym) && (refl_plus != NULL) ) continue;

			refldata[0] = nh;
			refldata[1] = nk;
			refldata[2] = nl;
			if ( refl_plus != NULL ) {
				refldata[3] = get_intensity(refl_plus);
				refldata[4] = get_esd_intensity(refl_plus);
			} else {
				refldata[3] = NAN;
				refldata[4] = NAN;
			}
			if ( refl_minus != NULL ) {
				refldata[5] = get_intensity(refl_minus);
				refldata[6] = get_esd_intensity(refl_minus);
			} else {
				refldata[5] = NAN;
				refldata[6] = NAN;
			}

			ccp4_lwrefl(mtz, refldata, columns, 7, refl_i++);

		}
	}

	r = MtzPut(mtz, " ");
	ccp4spg_free(&spg);
	MtzFree(mtz);
	return 1-r; /* Yes, really.  MtzPut return values are backwards */
#else
	return 1;
#endif
}


int libcrystfel_can_write_mtz()
{
#ifdef HAVE_LIBCCP4
	return 1;
#else
	return 0;
#endif
}
