/*
 * fromfile.c
 *
 * Perform indexing from solution file
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

/** \file fromfile.h */

/* There are 9 vector components,
 * 2 detector shifts, 1 profile radius,
 *  1 resolution limit */
#define NPARAMS_PER_LINE 11
/* The keys read from file
 * are the filename, event */
#define NKEYS_PER_LINE 2


struct fromfile_options
{
	char *filename;
};


struct fromfile_keys
{
	char filename[100];
	char event[100];
	int crystal_number;
};


struct fromfile_entries
{
	struct fromfile_keys key;
	float solution[NPARAMS_PER_LINE];
	UT_hash_handle hh;
};


struct fromfile_private
{
	UnitCell *cellTemplate;
	struct fromfile_entries *sol_hash;
};


void print_struct(struct fromfile_entries *sol_hash)
{
	struct fromfile_entries *s;

	s = calloc(1, sizeof(struct fromfile_entries));
	if ( s == NULL ) return;

	for( s=sol_hash; s != NULL; s=s->hh.next ) {
		printf("File %s, event %s, and crystal_number %d \n",
		       s->key.filename, s->key.event, s->key.crystal_number);
	}
}


void full_print_struct(struct fromfile_entries *sol_hash)
{
	struct fromfile_entries *s;
	s = calloc(1, sizeof(struct fromfile_entries));
	if ( s == NULL ) return;

	for( s=sol_hash; s != NULL; s=s->hh.next ) {
		printf("File %s, event %s, and crystal_number %d \n",
		       s->key.filename, s->key.event,
		       s->key.crystal_number);

		printf("Solution parameters:\n");
		for( int i = 0; i < NPARAMS_PER_LINE; i++ ){
			printf("%e", s->solution[i]);
		}
		printf("\n");
	}
}


int ncrystals_in_sol(char *path)
{
	FILE *fh;
	int count = 0;  /* Line counter (result) */
	char c;         /* To store a character read from file */

	fh = fopen(path, "r");

	if ( fh == NULL ) {
		ERROR("%s not found by ncrystals_in_sol\n", path);
		return 0;
	}

	for ( c = getc(fh); c != EOF; c = getc(fh) ) {
		if ( c == '\n' ) {
			count = count + 1;
		}
	}

	/* For the last line, which has no \n at the end*/
	count = count + 1;

	fclose(fh);

	return count;
}


char *read_unknown_string(FILE *fh)
{
	/* Source: "https://stackoverflow.com/questions/16870485/
	 * how-can-i-read-an-input-string-of-unknown-length" */

	char *str = NULL;
	int ch;
	size_t len = 0;
	size_t size = 1;

	str = realloc(NULL, sizeof(char)*size); //size is start size
	if ( !str ) {
		ERROR("Can't reallocate string size");
	}

	while( ( ch = fgetc(fh) ) != ' ' && ch != EOF ){
		if (ch != '\n'){
			str[len++]=ch;
		}
		if(len==size){
			size+=64;
			str = realloc(str, sizeof(char)*(size));
			if ( !str ) {
				ERROR("Can't reallocate string size");
			}
		}
	}

	return realloc(str, sizeof(char)*len);
}


void *fromfile_prepare(IndexingMethod *indm, struct fromfile_options *opts)
{
	FILE *fh;
	int nlines;
	int nparams_in_solution;
	int nentries;
	char *filename;
	char *event;
	int crystal_number;
	int current_line;
	int position_in_current_line;
	struct fromfile_entries *sol_hash = NULL;
	struct fromfile_entries *item = NULL;
	float params[NPARAMS_PER_LINE];

	fh = fopen(opts->filename, "r");
	if ( fh == NULL ) {
		ERROR("%s not found by fromfile_prepare\n", opts->filename);
		return NULL;
	} else {
		STATUS("Found solution file %s\n", opts->filename);
	}

	nlines = ncrystals_in_sol(opts->filename);
	/* Total crystal parameters in solution file */
	nparams_in_solution = nlines*NPARAMS_PER_LINE;
	/* Total entries in solution file */
	nentries = nlines*(NPARAMS_PER_LINE+NKEYS_PER_LINE);

	STATUS("Parsing solution file containing %d lines...\n", nlines);

	/* Reads indexing solutions */
	int j = 0; /* follows the solution parameter [0,NPARAMS_PER_LINE] */
	for(int i = 0; i < nentries; i++) {

		crystal_number = 0;

		current_line = i/(NPARAMS_PER_LINE+NKEYS_PER_LINE);

		position_in_current_line = (i)%(NPARAMS_PER_LINE+NKEYS_PER_LINE);

		if ( position_in_current_line == 0 ) {

			filename = read_unknown_string(fh);

			if ( !filename ){
				if ( current_line == nlines-1 ) break;
				printf("Failed to read a filename\n");
				return 0;
			}

		}

		if ( position_in_current_line == 1 ) {
			event = read_unknown_string(fh);
			if ( !event ){
				printf("Failed to read a event\n");
				return 0;
			}

		}

		if ( position_in_current_line > 1 ) {
			if ( fscanf(fh, "%e", &params[j]) != 1 ) {
				printf("Failed to read a parameter\n");
				return 0;
			}
			j+=1;
		}

		if ( j == (NPARAMS_PER_LINE) ) {

			/* Prepare to add to the hash table */
			item = (struct fromfile_entries *)malloc(sizeof *item);
			memset(item, 0, sizeof *item);
			strcpy(item->key.filename, filename);
			strcpy(item->key.event, event);
			item->key.crystal_number = crystal_number;
			for ( int k = 0; k < NPARAMS_PER_LINE; k++){
				item->solution[k] = params[k];
			}

			/* Verify the uniqueness of the key */
			struct fromfile_entries *uniqueness_test;
			HASH_FIND(hh, sol_hash, &item->key,
			          sizeof(struct fromfile_keys), uniqueness_test);

			if ( uniqueness_test == NULL ) {
				HASH_ADD(hh, sol_hash, key,
				         sizeof(struct fromfile_keys), item);
			} else {
				/* Look for the next available set of keys */
				do {
					uniqueness_test = NULL;
					crystal_number += 1;
					item->key.crystal_number = crystal_number;
					HASH_FIND(hh, sol_hash, &item->key,
					          sizeof(struct fromfile_keys),
					          uniqueness_test);
				} while ( uniqueness_test != NULL );

				HASH_ADD(hh, sol_hash, key,
				         sizeof(struct fromfile_keys), item);
			}

			j=0;

		}
	}

	fclose(fh);

	STATUS("Solution parsing done. Have %d parameters and %d total entries.\n",
	       nparams_in_solution, nentries);

	struct fromfile_private *dp;
	dp = (struct fromfile_private *) malloc( sizeof(struct fromfile_private));

	if ( dp == NULL ){
		return NULL;
	}

	dp->cellTemplate = cell;
	dp->sol_hash = sol_hash;

	STATUS("Solution lookup table initialized!\n");

	return (void *)dp;
}


int fromfile_index(struct image *image, void *mpriv, int crystal_number)
{
	Crystal *cr;
	UnitCell *cell;
	float asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;
	float xshift, yshift;
	struct fromfile_entries *item, *p, *pprime;
	int ncryst = 0;
	float *sol;
	struct fromfile_private *dp = mpriv;

	/* Look up the hash table */
	item = calloc(1, sizeof(struct fromfile_entries));
	strcpy(item->key.filename, image->filename);
	strcpy(item->key.event, image->ev);
	item->key.crystal_number = crystal_number;

	/* key already in the hash? */
	HASH_FIND(hh, dp->sol_hash, &item->key, sizeof(struct fromfile_keys), p);
	if ( p == NULL ) return 0;

	sol = &(p->solution)[0];

	asx = sol[0] * 1e9;
	asy = sol[1] * 1e9;
	asz = sol[2] * 1e9;
	bsx = sol[3] * 1e9;
	bsy = sol[4] * 1e9;
	bsz = sol[5] * 1e9;
	csx = sol[6] * 1e9;
	csy = sol[7] * 1e9;
	csz = sol[8] * 1e9;
	xshift = sol[9] * 1e-3;
	yshift = sol[10] * 1e-3;

	cell = cell_new();
	cell_set_reciprocal(cell, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);
	cell_set_lattice_type(cell, cell_get_lattice_type(dp->cellTemplate));
	cell_set_centering(cell, cell_get_centering(dp->cellTemplate));
	cell_set_unique_axis(cell, cell_get_unique_axis(dp->cellTemplate));

	cr = crystal_new();
	ncryst += 1;
	crystal_set_cell(cr, cell);
	crystal_set_det_shift(cr, xshift , yshift);
	image_add_crystal(image, cr);

	/* Look for additional crystals */
	item->key.crystal_number = crystal_number+1;
	HASH_FIND(hh,  dp->sol_hash, &item->key,
	          sizeof(struct fromfile_keys), pprime);

	/* If a similar tag exist,
	 * recursive call increasing the crystal_number by 1 */
	if ( pprime != NULL ) {
		ncryst += fromfile_index(image, mpriv, crystal_number+1);
	}

	return ncryst;
}


void fromfile_cleanup(void *mpriv)
{
	struct fromfile_private *dp = mpriv;

	/* FIXME: Implementation */

	free(dp);
}


static void fromfile_show_help()
{
	printf("Parameters for 'fromfile' indexing:\n"
"     --fromfile-input-file\n"
"                           Filename of indexing solution file\n"
);
}


int fromfile_default_options(FromFileOptions **opts_ptr)
{
	FromFileOptions *opts;
	opts = malloc(sizeof(struct fromfile_options));
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
		(*opts_ptr)->filename = strdup(arg);
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
