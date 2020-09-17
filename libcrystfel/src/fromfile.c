/*
 * fromfile.c
 *
 * Perform indexing from solution file
 *
 *
 * Authors:
 *   2020 Pascal Hogan-Lamarre <pascal.hogan.lamarre@mail.utoronto.ca>
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

#include "image.h"
#include "detector.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <fenv.h>
#include <unistd.h>

#include "uthash.h"

/** \file fromfile.h */

/* There are 9 vector components, 
 * 2 detector shifts, 1 profile radius,
 *  1 resolution limit */
#define NPARAMS_PER_LINE 13  
/* The keys are the filename, 
 * event, and crystal number */
#define NKEYS_PER_LINE 3

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
    s = (struct fromfile_entries *)malloc(sizeof *s);
    memset(s, 0, sizeof *s);

    for( s=sol_hash; s != NULL; s=(struct fromfile_entries*)(s->hh.next) ) {
        printf("File %s, event %s, and crystal_number %d \n",
		       s->key.filename, s->key.event, 
			   s->key.crystal_number);
    }
}

void full_print_struct(struct fromfile_entries *sol_hash)
{
    struct fromfile_entries *s;
    s = (struct fromfile_entries *)malloc(sizeof *s);
    memset(s, 0, sizeof *s);

    for( s=sol_hash; s != NULL; s=(struct fromfile_entries*)(s->hh.next) ) {
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
		ERROR("%s not found by ncrystals_in_sol\n",path);
		return 0;
	}

	for ( c = getc(fh); c != EOF; c = getc(fh) ){
        if ( c == '\n' ){
            count = count + 1; 
		}
	}
	
	/* For the last line, which has no \n at the end*/
	count = count + 1;
  
    fclose(fh); 

	return count;

}

void *fromfile_prepare(char *solution_filename, UnitCell *cell)
{	
	FILE *fh;
	int nlines;
	int nparams_in_solution;
	int nentries;
	char filename[100];   					
	char event[100];                 
	int crystal_number;
	int current_line;
	int position_in_current_line;
	struct fromfile_entries *sol_hash = NULL;
	struct fromfile_entries *item = NULL;
	float params[NPARAMS_PER_LINE];
	char path_to_sol[50],extension[10];
	char cwd[PATH_MAX];

	/* Assembling solution file name from input file name*/
	strcpy(path_to_sol, "../");
	strcpy(extension, ".sol");
	strcat(path_to_sol, strtok(solution_filename, "."));
	strcat(path_to_sol, extension);

	if (getcwd(cwd, sizeof(cwd)) != NULL) {
		ERROR("Cannot identify current directory\n");
	}

	fh = fopen(path_to_sol, "r");

	if ( fh == NULL ) {
		ERROR("%s not found by fromfile_prepare in %s\n", path_to_sol, cwd);
		return 0;
	}
	else {
		STATUS("Found solution file %s at %s\n", path_to_sol, cwd);
	}

	nlines = ncrystals_in_sol(path_to_sol);
	/* Total crystal parameters in solution file */
	nparams_in_solution = nlines*NPARAMS_PER_LINE; 
	/* Total entries in solution file */
	nentries = nlines*(NPARAMS_PER_LINE+NKEYS_PER_LINE); 
  
	STATUS("Parsing solution file containing %d lines...\n", nlines);

	/* Reads indexing solutions */
	int j = 0; /* follows the solution parameter [0,NPARAMS_PER_LINE] */
	for(int i = 0; i < nentries; i++)
	{	

		current_line = i/(NPARAMS_PER_LINE+NKEYS_PER_LINE);
		
		position_in_current_line = (i)%(NPARAMS_PER_LINE+NKEYS_PER_LINE);

		if ( position_in_current_line == 0 ){
			if ( fscanf(fh, "%s", filename) != 1 ) {
				if ( current_line == (nlines-1) ){
					break;
				}
				else{
				printf("Failed to read a filename\n");
				return 0;
				}
			}
		}

		if ( position_in_current_line == 1 ){

			if ( fscanf(fh, "%s", event) != 1 ) {
				printf("Failed to read an event\n");
				return 0;
			}
		}
			
		if ( position_in_current_line == 2 ){
			if ( fscanf(fh, "%d", &crystal_number) != 1 ) {
				printf("Failed to read a crystal number\n");
				return 0;
			}			
		}

		if ( position_in_current_line > 2 ){
			if ( fscanf(fh, "%e", &params[j]) != 1 ) {
				printf("Failed to read a parameter\n");
				return 0;
			}
			j+=1;			
		}

		if ( j == (NPARAMS_PER_LINE) ){

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
    		if (uniqueness_test==NULL) {
				HASH_ADD(hh, sol_hash, key, sizeof(struct fromfile_keys), item);
			}
			else{
				printf("Keys must be unique! Verify the combinations");
				return 0;
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

static void update_detector(struct detector *det, double xoffs, double yoffs)
{
	int i;

	for ( i = 0; i < det->n_panels; i++ ) {
		struct panel *p = &det->panels[i];
		p->cnx += xoffs * p->res;
		p->cny += yoffs * p->res;
	}
}

int fromfile_index(struct image *image, void *mpriv, int crystal_number)
{
	Crystal *cr;
	UnitCell *cell;
	float asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;
	float xshift, yshift, profile_radius, resolution_limit;
	struct fromfile_entries *item, *p, *pprime;
	float *sol;
	
	struct fromfile_private *dp = (struct fromfile_private *)mpriv;

	/* Look up the hash table */
	item = (struct fromfile_entries *)malloc(sizeof *item);
	memset(item, 0, sizeof *item);
	strcpy(item->key.filename, image->filename);
	strcpy(item->key.event, get_event_string(image->event));
	item->key.crystal_number = crystal_number;

	/* key already in the hash? */
    HASH_FIND(hh,  dp->sol_hash, &item->key, sizeof(struct fromfile_keys), p);
    if ( p == NULL ) {
		return 0;
	}

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
	profile_radius = sol[11] * 1e9;
	resolution_limit = sol[12] * 1e9;

	cell = cell_new();
	cell_set_reciprocal(cell, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);
    cell_set_lattice_type(cell, cell_get_lattice_type(dp->cellTemplate));
	cell_set_centering(cell, cell_get_centering(dp->cellTemplate));
	cell_set_unique_axis(cell, cell_get_unique_axis(dp->cellTemplate));
	
	cr = crystal_new();
	crystal_set_cell(cr, cell);
	crystal_set_profile_radius(cr, profile_radius);
	crystal_set_resolution_limit(cr, resolution_limit);
	crystal_set_det_shift(cr, xshift , yshift);
	update_detector(image->det, xshift , yshift);
	image_add_crystal(image, cr);

	/*Look for additional crystals*/
	item->key.crystal_number = crystal_number+1;
	HASH_FIND(hh,  dp->sol_hash, &item->key, 
	          sizeof(struct fromfile_keys), pprime);

    if ( pprime == NULL ) {
		/* If no more crystal, done */
		return 1;
	}
	else{
		/* If more crystals, recursive call for next crystal in line */
		fromfile_index(image, mpriv, crystal_number+1);
	}
	
	dp->sol_hash = NULL; /* Clean up local copy */
	
	return 1;
	
}