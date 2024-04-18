/*
 * fromfile.c
 *
 * Perform indexing from solution file
 *
 * Copyright © 2020-2021 Max-Planck-Gesellschaft
 *                       zur Förderung der Wissenschaften e.V.
 * Copyright © 2021 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Pascal Hogan-Lamarre <pascal.hogan.lamarre@mail.utoronto.ca>
 *   2021 Thomas White <thomas.white@desy.de>
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <fenv.h>
#include <unistd.h>
#include <argp.h>

#include "image.h"
#include "uthash.h"
#include "index.h"

/** \file fromfile.h */

#define MAX_KEY_LEN (256)
#define MAX_CRYSTALS (16)

struct fromfile_key
{
	char filename[MAX_KEY_LEN];
	char event[MAX_KEY_LEN];
};


struct fromfile_entry
{
	struct fromfile_key key_field;
	Crystal *crystals[MAX_CRYSTALS];
	int n_crystals;
	UT_hash_handle hh;
};


struct fromfile_private
{
	struct fromfile_entry *sol_hash;
};


static int make_key(struct fromfile_key *key,
                    const char *filename, const char *ev)
{
	if ( (strlen(filename) > MAX_KEY_LEN) || (strlen(ev) > MAX_KEY_LEN) ) {
		ERROR("Filename/event too long: %s %s\n", filename, ev);
		return 1;
	}

	/* The entire structure is used as a key, not just the pre-terminator
	 * parts of the strings.  Therefore it must be initialised to zero */
	memset(key, 0, sizeof(struct fromfile_key));

	strcpy(key->filename, filename);
	strcpy(key->event, ev);

	return 0;
}


struct fromfile_entry *add_unique(struct fromfile_entry **phead,
                                  struct fromfile_key key)
{
	struct fromfile_entry *p;
	struct fromfile_entry *head = *phead;

	HASH_FIND(hh, head, &key, sizeof(struct fromfile_key), p);
	if ( p == NULL ) {

		struct fromfile_entry *item;

		item = cfmalloc(sizeof(struct fromfile_entry));
		if ( item == NULL ) return NULL;

		item->n_crystals = 0;
		item->key_field = key;

		HASH_ADD(hh, head, key_field, sizeof(struct fromfile_key), item);
		*phead = head;
		return item;

	} else {
		return p;
	}
}


static int set_ua(UnitCell *cell, const char *ltsym)
{
	if ( strlen(ltsym) != 3 ) return 1;
	cell_set_unique_axis(cell, ltsym[2]);
	return 0;
}


static int set_lattice(UnitCell *cell, const char *ltsym)
{
	if ( (strlen(ltsym) != 2) && (strlen(ltsym) != 3) ) return 1;

	switch ( ltsym[1] ) {
		case 'P':
		case 'A':
		case 'B':
		case 'C':
		case 'I':
		case 'F':
		case 'R':
		case 'H':
		break;

		default:
		return 1;
	}
	cell_set_centering(cell, ltsym[1]);

	switch ( ltsym[0] ) {

		case 'a' :
		cell_set_lattice_type(cell, L_TRICLINIC);
		break;

		case 'm' :
		cell_set_lattice_type(cell, L_MONOCLINIC);
		return set_ua(cell, ltsym);

		case 'o' :
		cell_set_lattice_type(cell, L_ORTHORHOMBIC);
		break;

		case 't' :
		cell_set_lattice_type(cell, L_TETRAGONAL);
		return set_ua(cell, ltsym);

		case 'c' :
		cell_set_lattice_type(cell, L_CUBIC);
		break;

		case 'r' :
		cell_set_lattice_type(cell, L_RHOMBOHEDRAL);
		break;

		case 'h' :
		cell_set_lattice_type(cell, L_HEXAGONAL);
		return set_ua(cell, ltsym);

		default :
		return 1;
	}

	return 0;
}


void *fromfile_prepare(IndexingMethod *indm, struct fromfile_options *opts)
{
	FILE *fh;
	struct fromfile_private *dp;

	if ( opts->filename == NULL ) {
		ERROR("Please try again with --fromfile-input-file\n");
		return NULL;
	}

	/* If filename is not absolute, jump out of working directory */
	if ( opts->filename[0] == '/' ) {
		fh = fopen(opts->filename, "r");
	} else {
		char *prefixed_fn = cfmalloc(4+strlen(opts->filename));
		if ( prefixed_fn == NULL ) return NULL;
		strcpy(prefixed_fn, "../");
		strcat(prefixed_fn, opts->filename);
		fh = fopen(prefixed_fn, "r");
		cffree(prefixed_fn);
	}

	if ( fh == NULL ) {
		ERROR("Couldn't find solution file '%s'\n", opts->filename);
		return NULL;
	}

	dp = cfmalloc(sizeof(struct fromfile_private));
	if ( dp == NULL ) {
		fclose(fh);
		return NULL;
	}

	dp->sol_hash = NULL;

	/* Read indexing solutions */
	do {

		char *rval;
		char line[1024];
		int i, n;
		char **bits;
		float vals[11];
		struct fromfile_key key;
		Crystal *cr;
		size_t len;
		int n_sp;
		struct fromfile_entry *item = NULL;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) break;

		chomp(line);
		notrail(line);

		len = strlen(line);
		n_sp = 0;
		for ( i=len-1; i>0; i-- ) {
			if ( line[i] == ' ' ) {
				n_sp++;
				if ( n_sp == 13 ) {
					line[i] = '\0';
					break;
				}
			}
		}

		n = assplode(line+i+1, " \t,", &bits, ASSPLODE_NONE);
		if ( n < 13 ) {
			ERROR("Badly formatted line '%s'\n", line);
			return NULL;
		}

		/* filename, event, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz,
		 * det_shift_x, det_shift_y, latticetype+centering */
		for ( i=1; i<12; i++ ) {
			if (sscanf(bits[i], "%f", &vals[i-1]) != 1)
			{
				ERROR("Invalid value for number %i\n", i);
				return NULL;
			}
		}

		if ( make_key(&key, line, bits[0]) ) {
			ERROR("Failed to make key for %s %s\n",
			      line, bits[0]);
			continue;
		}

		item = add_unique(&dp->sol_hash, key);
		if ( item == NULL ) {
			ERROR("Failed to add/find entry for %s %s\n",
			      line, bits[0]);
			continue;
		}

		if ( item->n_crystals == MAX_CRYSTALS ) {

			ERROR("Too many crystals for %s %s\n", line, bits[0]);

		} else {

			UnitCell *cell;

			cr = crystal_new();

			/* mm -> m */
			crystal_set_det_shift(cr, vals[9]*1e-3, vals[10]*1e-3);

			cell = cell_new();
			cell_set_reciprocal(cell, vals[0]*1e9, vals[1]*1e9, vals[2]*1e9,
			                          vals[3]*1e9, vals[4]*1e9, vals[5]*1e9,
			                          vals[6]*1e9, vals[7]*1e9, vals[8]*1e9);
			if ( set_lattice(cell, bits[12]) ) {
				ERROR("Invalid lattice type '%s'\n", bits[12]);
			} else {
				crystal_set_cell(cr, cell);
				item->crystals[item->n_crystals++] = cr;
			}

		}

		for ( i=0; i<n; i++ ) cffree(bits[i]);
		cffree(bits);

	} while ( 1 );

	fclose(fh);

	STATUS("Read %i crystals from %s\n",
	       HASH_CNT(hh, dp->sol_hash), opts->filename);

	return dp;
}


int fromfile_index(struct image *image, void *mpriv)
{
	struct fromfile_entry *p;
	struct fromfile_private *dp = mpriv;
	struct fromfile_key key;
	int i;

	make_key(&key, image->filename, image->ev);

	HASH_FIND(hh, dp->sol_hash, &key, sizeof(struct fromfile_key), p);
	if ( p == NULL ) {
		STATUS("WARNING: No solution for %s %s\n",
		       image->filename, image->ev);
		return 0;
	}

	for ( i=0; i<p->n_crystals; i++ ) {
		Crystal *cr;
		cr = crystal_copy(p->crystals[i]);
		image_add_crystal(image, cr);
	}

	return p->n_crystals;
}


void fromfile_cleanup(void *mpriv)
{
	struct fromfile_private *dp = mpriv;
	struct fromfile_entry *item, *tmp;

	HASH_ITER(hh, dp->sol_hash, item, tmp) {
		int i;
		HASH_DEL(dp->sol_hash, item);
		for ( i=0; i<item->n_crystals; i++ ) {
			Crystal *cr = item->crystals[i];
			cell_free(crystal_get_cell(cr));
			crystal_free(cr);
		}
	}

	cffree(dp);
}


static void fromfile_show_help()
{
	printf("Parameters for 'fromfile' indexing:\n"
"     --fromfile-input-file\n"
"                           Filename of indexing solution file\n"
);
}


int fromfile_default_options(struct fromfile_options **opts_ptr)
{
	struct fromfile_options *opts;
	opts = cfmalloc(sizeof(struct fromfile_options));
	if ( opts == NULL ) return ENOMEM;
	opts->filename = NULL;
	*opts_ptr = opts;
	return 0;
}


static error_t fromfile_parse_arg(int key, char *arg,
                                 struct argp_state *state)
{
	struct fromfile_options **opts_ptr = state->input;
	int r;

	switch ( key ) {

		case ARGP_KEY_INIT :
		r = fromfile_default_options(opts_ptr);
		if ( r ) return r;
		break;

		case 1 :
		fromfile_show_help();
		return EINVAL;

		case 2 :
		(*opts_ptr)->filename = cfstrdup(arg);
		break;

		default :
		return ARGP_ERR_UNKNOWN;

	}

	return 0;
}


static struct argp_option fromfile_options[] = {

	{"help-fromfile", 1, NULL, OPTION_NO_USAGE,
	 "Show options for 'from file' indexing", 99},

	{"fromfile-input-file", 2, "filename", OPTION_HIDDEN, NULL},
	{0}
};


struct argp fromfile_argp = { fromfile_options, fromfile_parse_arg,
                              NULL, NULL, NULL, NULL, NULL };
