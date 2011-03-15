/*
 * reflist-utils.c
 *
 * Utilities to complement the core reflist.c
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <stdio.h>


#include "reflist.h"
#include "cell.h"


void write_reflections_to_file(FILE *fh, RefList *list, UnitCell *cell)
{
	Reflection *refl;
	RefListIterator *iter;

	fprintf(fh, "  h   k   l          I    phase   sigma(I) "
		     " 1/d(nm^-1)  counts  fs/px  ss/px\n");

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		double intensity, esd_i, s;
		double fs, ss;

		get_indices(refl, &h, &k, &l);
		get_detector_pos(refl, &fs, &ss);
		intensity = get_intensity(refl);
		esd_i = 0.0; /* FIXME! */
		s = 0.0; /* FIXME! */

		/* h, k, l, I, sigma(I), s */
		fprintf(fh,
		       "%3i %3i %3i %10.2f %s %10.2f  %10.2f %7i %6.1f %6.1f\n",
		       h, k, l, intensity, "       -", esd_i, s/1.0e9, 1,
		       fs, ss);

	}
}


int write_reflist(const char *filename, RefList *list, UnitCell *cell)
{
	FILE *fh;

	if ( filename == NULL ) {
		fh = stdout;
	} else {
		fh = fopen(filename, "w");
	}

	if ( fh == NULL ) {
		ERROR("Couldn't open output file '%s'.\n", filename);
		return 1;
	}

	write_reflections_to_file(fh, list, cell);

	fclose(fh);

	return 0;
}
