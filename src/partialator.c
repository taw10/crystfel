/*
 * partialator.c
 *
 * Scaling and post refinement for coherent nanocrystallography
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2021 Thomas White <taw@physics.org>
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

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <pthread.h>
#include <gsl/gsl_errno.h>
#include <sys/stat.h>

#include <image.h>
#include <utils.h>
#include <symmetry.h>
#include <stream.h>
#include <geometry.h>
#include <peaks.h>
#include <thread-pool.h>
#include <reflist.h>
#include <reflist-utils.h>
#include <cell.h>
#include <cell-utils.h>

#include "scaling.h"
#include "post-refinement.h"
#include "merge.h"
#include "rejection.h"
#include "version.h"
#include "json-utils.h"


struct csplit_hash_entry
{
	int    n_events;
	char **events;
	int   *datasets;
};

#define CSPLIT_HASH_MAX (65521)

struct custom_split
{
	int    n_events_total;
	int    n_datasets;
	char **dataset_names;
	struct csplit_hash_entry hashes[CSPLIT_HASH_MAX];
};


static int csplit_hash(const char *id)
{
	int i;
	size_t len = strlen(id);
	int h = 0;

	for ( i=0; i<len; i++ ) {
		h = (h*31 + id[i]) % CSPLIT_HASH_MAX;
	}
	assert(h < CSPLIT_HASH_MAX);

	return h;
}


static void add_to_hash_entry(struct csplit_hash_entry *he, const char *id,
                              int dsn)
{
	he->events = realloc(he->events, (1+he->n_events)*sizeof(char *));
	he->datasets = realloc(he->datasets, (1+he->n_events)*sizeof(int));
	if ( (he->events == NULL) || (he->datasets == NULL) ) {
		ERROR("Failed to grow csplit hash entry.\n");
		abort();
	}

	he->events[he->n_events] = strdup(id);
	he->datasets[he->n_events] = dsn;
	he->n_events++;
}


static signed int find_dsn_for_id(struct custom_split *csplit, const char *id)
{
	int hash = csplit_hash(id);
	int i;
	struct csplit_hash_entry *he = &csplit->hashes[hash];

	for ( i=0; i<he->n_events; i++ ) {
		if ( strcmp(he->events[i], id) == 0 ) {
			return he->datasets[i];
		}
	}

	return -1;
}


/* Find dataset number */
static int find_dsn(struct custom_split *csplit, const char *ds)
{
	int i;

	for ( i=0; i<csplit->n_datasets; i++ ) {
		if ( strcmp(csplit->dataset_names[i], ds) == 0 ) {
			return i;
		}
	}

	csplit->dataset_names = realloc(csplit->dataset_names,
	                                (1+csplit->n_datasets)*sizeof(char *));
	if ( csplit->dataset_names == NULL ) {
		ERROR("Failed to grow list of dataset names\n");
		abort();
	}

	csplit->n_datasets++;
	csplit->dataset_names[csplit->n_datasets-1] = strdup(ds);
	return csplit->n_datasets-1;
}


/* Add arbitrary ID 'id' to dataset table with name 'ds' */
static void add_to_csplit(struct custom_split *csplit, const char *id,
                          const char *ds)
{
	int dsn;
	int hash;
	struct csplit_hash_entry *he;

	dsn = find_dsn(csplit, ds);

	hash = csplit_hash(id);
	he = &csplit->hashes[hash];
	add_to_hash_entry(he, id, dsn);
	csplit->n_events_total++;
}


/* Write two-way split results (i.e. for CC1/2 etc) for this list of crystals */
static void write_split(Crystal **crystals, int n_crystals, const char *outfile,
                        int nthreads, PartialityModel pmodel,
                        int min_measurements, SymOpList *sym, double push_res)
{
	char tmp[1024];
	RefList *split;
	Crystal **crystals1;
	Crystal **crystals2;
	int n_crystals1 = 0;
	int n_crystals2 = 0;
	int i;

	if ( n_crystals == 0 ) {
		ERROR("No crystals for split!\n");
		return;
	}

	crystals1 = malloc(n_crystals * sizeof(Crystal *));
	if ( crystals1 == NULL ) return;

	crystals2 = malloc(n_crystals * sizeof(Crystal *));
	if ( crystals2 == NULL ) return;


	for ( i=0; i<n_crystals; i++ ) {
		if ( i % 2 ) {
			crystals1[n_crystals1] = crystals[i];
			n_crystals1++;
		} else {
			crystals2[n_crystals2] = crystals[i];
			n_crystals2++;
		}
	}
	snprintf(tmp, 1024, "%s1", outfile);
	split = merge_intensities(crystals1, n_crystals1, nthreads,
		                  min_measurements, push_res, 1, 0);

	if ( split == NULL ) {
		ERROR("Not enough crystals for two way split!\n");
		free(crystals1);
		free(crystals2);
		return;
	}

	STATUS("Writing two-way split to %s ", tmp);
	write_reflist_2(tmp, split, sym);
	free_contribs(split);
	reflist_free(split);
	snprintf(tmp, 1024, "%s2", outfile);
	split = merge_intensities(crystals2, n_crystals2, nthreads,
		                  min_measurements, push_res, 1, 0);
	STATUS("and %s\n", tmp);
	write_reflist_2(tmp, split, sym);
	free_contribs(split);
	reflist_free(split);

	free(crystals1);
	free(crystals2);
}


static char *insert_into_filename(const char *fn, const char *add)
{
	int i;
	char *out;

	out = malloc(strlen(fn) + strlen(add) + 2);
	if ( out == NULL ) return NULL;

	for ( i=strlen(fn); i>0; i-- ) {
		if ( fn[i] == '.' ) {
			strncpy(out, fn, i);
			out[i] = '\0';
			strcat(out, "-");
			strcat(out, add);
			strcat(out, fn+i);
			return out;
		}
	}

	/* Fallback if fn does not contain a dot */
	strcpy(out, fn);
	strcat(out, "-");
	strcat(out, add);
	return out;
}


/* Write custom split results (including a two-way split) */
static void write_custom_split(struct custom_split *csplit, int dsn,
                               Crystal **crystals, int n_crystals,
                               PartialityModel pmodel, int min_measurements,
                               double push_res, SymOpList *sym, int nthreads,
                               const char *outfile)
{
	char *tmp;
	RefList *split;
	Crystal *crystalsn[n_crystals];
	int n_crystalsn = 0;
	int i;

	for ( i=0; i<n_crystals; i++ ) {

		const char *fn;
		char *evs;
		char *id;
		int dsn_crystal;

		fn = crystal_get_image(crystals[i])->filename;
		evs = crystal_get_image(crystals[i])->ev;

		id = malloc(strlen(evs)+strlen(fn)+2);
		if ( id == NULL ) {
			ERROR("Failed to allocate ID\n");
			return;
		}
		strcpy(id, fn);
		strcat(id, " ");
		strcat(id, evs);
		dsn_crystal = find_dsn_for_id(csplit, id);
		free(id);
		if ( dsn == dsn_crystal ) {
			crystalsn[n_crystalsn] = crystals[i];
			n_crystalsn++;
		}

	}

	tmp = insert_into_filename(outfile, csplit->dataset_names[dsn]);

	if ( n_crystalsn == 0 ) {
		ERROR("Not writing dataset '%s' because it contains no "
		      "crystals\n", csplit->dataset_names[dsn]);
		return;
	}

	STATUS("Writing dataset '%s' to %s (%i crystals)\n",
	       csplit->dataset_names[dsn], tmp, n_crystalsn);
	split = merge_intensities(crystalsn, n_crystalsn, nthreads,
		                  min_measurements, push_res, 1, 0);
	write_reflist_2(tmp, split, sym);
	free_contribs(split);
	reflist_free(split);

	write_split(crystalsn, n_crystalsn, tmp, nthreads, pmodel,
	            min_measurements, sym, push_res);
	free(tmp);
}


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Scaling and post refinement for coherent nanocrystallography.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"      --version              Print CrystFEL version number and exit.\n"
"\n"
"  -i, --input=<filename>     Specify the name of the input 'stream'.\n"
"  -o, --output=<filename>    Output filename.  Default: partialator.hkl.\n"
"      --output-every-cycle   Write .hkl* and .params files in every cycle.\n"
"  -y, --symmetry=<sym>       Merge according to symmetry <sym>.\n"
"      --start-after=<n>      Skip <n> crystals at the start of the stream.\n"
"      --stop-after=<n>       Stop after merging <n> crystals.\n"
"  -n, --iterations=<n>       Run <n> cycles of scaling and post-refinement.\n"
"      --no-scale             Disable scale factor (G, B) refinement.\n"
"      --no-Bscale            Disable B factor scaling.\n"
"      --no-pr                Disable orientation/physics refinement.\n"
"      --no-deltacchalf       Disable rejection based on deltaCChalf.\n"
"  -m, --model=<model>        Specify partiality model.\n"
"      --min-measurements=<n> Minimum number of measurements to require.\n"
"      --no-polarisation      Disable polarisation correction.\n"
"      --polarisation=<p>     Specify type of polarisation correction.\n"
"      --max-adu=<n>          Saturation value of detector.\n"
"      --min-res=<n>          Merge only crystals which diffract above <n> A.\n"
"      --push-res=<n>         Merge higher than apparent resolution cutoff.\n"
"  -j <n>                     Run <n> analyses in parallel.\n"
"      --no-free              Disable cross-validation (testing only).\n"
"      --custom-split         List of files for custom dataset splitting.\n"
"      --max-rel-B            Maximum allowable relative |B| factor.\n"
"      --no-logs              Do not write extensive log files.\n"
"      --log-folder=<fn>      Location for log folder.\n"
"  -w <pg>                    Apparent point group for resolving ambiguities.\n"
"      --operator=<op>        Indexing ambiguity operator for resolving.\n"
"      --force-bandwidth=<n>  Set all bandwidths to <n> (fraction).\n"
"      --force-radius=<n>     Set all profile radii to <n> nm^-1.\n"
"      --force-lambda=<n>     Set all wavelengths to <n> A.\n");
}


static signed int find_first_crystal(Crystal **crystals, int n_crystals,
                                     struct custom_split *csplit, int dsn)
{
	int i;

	for ( i=0; i<n_crystals; i++ ) {
		const char *fn;
		char *evs;
		char *id;
		int dsn_crystal;

		fn = crystal_get_image(crystals[i])->filename;
		evs = crystal_get_image(crystals[i])->ev;

		id = malloc(strlen(evs)+strlen(fn)+2);
		if ( id == NULL ) {
			ERROR("Failed to allocate ID\n");
			return -1;
		}
		strcpy(id, fn);
		strcat(id, " ");
		strcat(id, evs);
		dsn_crystal = find_dsn_for_id(csplit, id);
		free(id);

		if ( dsn == dsn_crystal ) return i;
	}
	return -1;
}


static void check_csplit(Crystal **crystals, int n_crystals,
                         struct custom_split *csplit)
{
	int i;
	int n_nosplit = 0;
	int n_split = 0;
	int n_cry = 0;
	int n_nocry = 0;

	STATUS("Checking your custom split datasets...\n");

	for ( i=0; i<n_crystals; i++ ) {

		const char *fn;
		char *evs;
		char *id;
		int dsn_crystal;

		fn = crystal_get_image(crystals[i])->filename;
		evs = crystal_get_image(crystals[i])->ev;

		id = malloc(strlen(evs)+strlen(fn)+2);
		if ( id == NULL ) {
			ERROR("Failed to allocate ID\n");
			return;
		}
		strcpy(id, fn);
		strcat(id, " ");
		strcat(id, evs);
		dsn_crystal = find_dsn_for_id(csplit, id);
		free(id);
		if ( dsn_crystal == -1 ) {
			n_nosplit++;
		} else {
			n_split++;
		}

	}

	for ( i=0; i<csplit->n_datasets; i++ ) {

		/* Try to find a crystal with dsn = i */
		if ( find_first_crystal(crystals, n_crystals, csplit, i) != -1 )
		{
			n_cry++;
		} else {
			n_nocry++;
			STATUS("Dataset %s has no crystals.\n",
			       csplit->dataset_names[i]);
		}
	}

	STATUS("Please check that these numbers match your expectations:\n");
	STATUS("    Number of crystals assigned to a dataset: %i\n", n_split);
	STATUS("Number of crystals with no dataset asssigned: %i\n", n_nosplit);
	STATUS("Number of datasets with at least one crystal: %i\n", n_cry);
	STATUS("         Number of datasets with no crystals: %i\n", n_nocry);
}


static int looks_like_event(const char *str)
{
	if ( strstr(str, "//") == NULL ) {
		return 0;
	} else {
		return 1;
	}
}


static struct custom_split *load_custom_split(const char *filename)
{
	struct custom_split *csplit;
	FILE *fh;
	int i;
	int lno = 0;

	csplit = malloc(sizeof(struct custom_split));
	if ( csplit == NULL ) return NULL;
	csplit->n_datasets = 0;
	csplit->n_events_total = 0;
	csplit->dataset_names = NULL;
	for ( i=0; i<CSPLIT_HASH_MAX; i++ ) {
		csplit->hashes[i].n_events = 0;
		csplit->hashes[i].events = NULL;
		csplit->hashes[i].datasets = NULL;
	}

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		ERROR("Failed to open '%s'\n", filename);
		free(csplit);
		return NULL;
	}

	do {

		char *rval;
		char line[1024];
		char *fn;
		char *evs;
		char *ds;
		char *id;
		size_t n, ev_start, ds_start;

		lno++;
		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) break;

		chomp(line);
		notrail(line);

		/* Look for start of dataset */
		n = strlen(line);
		while ( line[n] != ' ' && n > 0 ) n--;
		if ( n == 0 ) {
			ERROR("Custom split file line %i has too few (only 1) "
			      "fields.\n", lno);
			free(csplit);
			return NULL;
		}
		ds_start = n+1;
		ds = strdup(&line[ds_start]);

		n--;
		while ( line[n] != ' ' && n > 0 ) n--;
		if ( n == 0 ) {
			ev_start = 0;
		} else {
			ev_start = n+1;
		}

		evs = strndup(&line[ev_start], ds_start-ev_start-1);
		if ( !looks_like_event(evs) || (ev_start == 0) ) {
			/* It doesn't look like an event ID - assume it's part
			 * of the filename (which contains spaces) */
			ev_start = 0;
		}

		if ( ev_start > 0 ) {
			evs = strndup(&line[ev_start], ds_start-ev_start-1);
			fn = strndup(line, ev_start-1);
		} else {
			evs = strdup("//");
			fn = strndup(line, ds_start-1);
		}

		id = malloc(strlen(fn) + strlen(evs) + 2);
		strcpy(id, fn);
		strcat(id, " ");
		strcat(id, evs);
		add_to_csplit(csplit, id, ds);
		free(id);
		free(fn);
		free(evs);
		free(ds);

	} while ( 1 );

	fclose(fh);

	int max = 0;
	for ( i=0; i<CSPLIT_HASH_MAX; i++ ) {
		if ( csplit->hashes[i].n_events > max ) {
			max = csplit->hashes[i].n_events;
		}
	}

	STATUS("Hash table load factor = %.2f (max %i)\n",
	       (double)csplit->n_events_total / CSPLIT_HASH_MAX, max);

	return csplit;
}


static void display_progress(int n_images, int n_crystals)
{
	if ( !isatty(STDERR_FILENO) ) return;
	if ( tcgetpgrp(STDERR_FILENO) != getpgrp() ) return;

	pthread_mutex_lock(&stderr_lock);
	fprintf(stderr, "\r%i images loaded, %i crystals.",
	        n_images, n_crystals);
	pthread_mutex_unlock(&stderr_lock);

	fflush(stdout);
}


static RefList *apply_max_adu(RefList *list, double max_adu)
{
	RefList *nlist;
	Reflection *refl;
	RefListIterator *iter;

	nlist = reflist_new();
	if ( nlist == NULL ) return NULL;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		if ( get_peak(refl) < max_adu ) {
			signed int h, k, l;
			get_indices(refl, &h, &k, &l);
			Reflection *nrefl = add_refl(nlist, h, k, l);
			if ( nrefl == NULL ) {
				ERROR("Failed to add reflection\n");
				return NULL;
			}
			copy_data(nrefl, refl);
		}
	}
	return nlist;
}


static void skip_to_end(FILE *fh)
{
	int c;
	do {
		c = fgetc(fh);
	} while ( (c != '\n') && (c != EOF) );
}


static int set_initial_params(Crystal *cr, FILE *fh, double force_bandwidth,
                              double force_radius, double force_lambda)
{
	struct image *image = crystal_get_image(cr);

	if ( fh != NULL ) {

		int err;
		int n;
		float osf, B;

		err = fscanf(fh, "%i %f %f", &n, &osf, &B);
		if ( err != 3 ) {
			ERROR("Failed to read parameters.\n");
			return 1;
		}

		crystal_set_osf(cr, osf);
		crystal_set_Bfac(cr, B*1e-20);

		skip_to_end(fh);

	} else {

		crystal_set_osf(cr, 1.0);
		crystal_set_Bfac(cr, 0.0);

	}

	if ( force_bandwidth > 0.0 ) {
		image->bw = force_bandwidth;
	}
	if ( force_radius > 0.0 ) {
		crystal_set_profile_radius(cr, force_radius);
	}
	if ( force_lambda > 0.0 ) {
		image->lambda = force_lambda;
	}
	spectrum_free(image->spectrum);
	image->spectrum = spectrum_generate_gaussian(image->lambda, image->bw);

	return 0;
}


/* Flag a random 5% of reflections */
static void select_free_reflections(RefList *list, gsl_rng *rng)
{
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		set_flag(refl, random_flat(rng, 1.0) > 0.95);
	}
}


static void write_to_pgraph(FILE *fh, RefList *list, RefList *full, Crystal *cr,
                            int fr, signed int inum)
{
	Reflection *refl;
	RefListIterator *iter;
	double G = crystal_get_osf(cr);
	double B = crystal_get_Bfac(cr);
	UnitCell *cell = crystal_get_cell(cr);
	char ins[16];

	if ( inum >= 0 ) {
		snprintf(ins, 12, "%i", inum);
	} else {
		ins[0] = 'F';
		ins[1] = '\0';
	}


	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double pobs, pcalc;
		double res, Ipart;
		Reflection *match;

		if ( !get_flag(refl) ) continue;  /* Not free-flagged */

		/* Strong reflections only */
		if ( get_intensity(refl) < 3.0*get_esd_intensity(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		res = resolution(cell, h, k, l);

		match = find_refl(full, h, k, l);
		if ( match == NULL ) continue;

		/* Don't calculate pobs if reference reflection is weak */
		if ( fabs(get_intensity(match)) / get_esd_intensity(match) < 3.0 ) continue;

		/* Calculated partiality */
		pcalc = get_partiality(refl);

		/* Observed partiality */
		Ipart = correct_reflection_nopart(get_intensity(refl), refl, G, B, res);
		pobs = Ipart / get_intensity(match);

		fprintf(fh, "%5i %4i %4i %4i %e %e %8.3f %8.3f %s\n",
		        fr, h, k, l, 2*res, Ipart, pcalc, pobs, ins);

	}
}


static void write_pgraph(RefList *full, Crystal **crystals, int n_crystals,
                         signed int iter, const char *suff,
                         const char *log_folder)
{
	FILE *fh;
	char tmp[256];
	int i;

	snprintf(tmp, 256, "%s/pgraph%s.dat", log_folder, suff);

	if ( iter == 0 ) {
		fh = fopen(tmp, "w");
	} else {
		fh = fopen(tmp, "a");
	}

	if ( fh == NULL ) {
		ERROR("Failed to open '%s'\n", tmp);
		return;
	}

	if ( iter == 0 ) {
		fprintf(fh, "  Crystal    h    k    l  1/d(m)   Ipart    pcalc   pobs   iteration\n");
	}

	for ( i=0; i<n_crystals; i++ ) {
		if ( crystal_get_user_flag(crystals[i]) != 0 ) continue;
		write_to_pgraph(fh, crystal_get_reflections(crystals[i]), full,
		                crystals[i], i, iter);
	}

	fclose(fh);
}


static void all_residuals(Crystal **crystals, int n_crystals, RefList *full,
                          int no_free,
                          double *presidual, double *pfree_residual,
                          double *plog_residual, double *pfree_log_residual,
                          int *pn_used)
{
	int i;
	int n_used = 0;
	int n_nan_linear = 0;
	int n_nan_linear_free = 0;
	int n_nan_log = 0;
	int n_nan_log_free = 0;
	int n_non_linear = 0;
	int n_non_linear_free = 0;
	int n_non_log = 0;
	int n_non_log_free = 0;

	*presidual = 0.0;
	*pfree_residual = 0.0;
	*plog_residual = 0.0;
	*pfree_log_residual = 0.0;

	for ( i=0; i<n_crystals; i++ ) {

		double r, free_r, log_r, free_log_r;
		int n;

		if ( crystal_get_user_flag(crystals[i]) ) continue;

		/* Scaling should have been done right before calling this */
		r = residual(crystals[i], full, 0, &n, NULL);
		if ( n == 0 ) {
			n_non_linear++;
		} else if ( isnan(r) ) {
			n_nan_linear++;
		}
		free_r = residual(crystals[i], full, 1, &n, NULL);
		if ( n == 0 ) {
			n_non_linear_free++;
		} else if ( isnan(free_r) ) {
			n_nan_linear_free++;
		}
		log_r = log_residual(crystals[i], full, 0, &n, NULL);
		if ( n == 0 ) {
			n_non_log++;
		} else if ( isnan(log_r) ) {
			n_nan_log++;
		}
		free_log_r = log_residual(crystals[i], full, 1, &n, NULL);
		if ( n == 0 ) {
			n_non_log_free++;
		} else if ( isnan(free_log_r) ) {
			n_nan_log_free++;
		}

		if ( isnan(r) || isnan(log_r) ) continue;

		if ( !no_free && (isnan(free_r) || isnan(free_log_r)) ) continue;

		*presidual += r;
		*pfree_residual += free_r;
		*plog_residual += log_r;
		*pfree_log_residual += free_log_r;

		n_used++;
	}

	if ( n_non_linear ) {
		ERROR("WARNING: %i crystals had no reflections in linear "
		      "residual calculation\n", n_non_linear);
	}
	if ( n_non_linear_free ) {
		ERROR("WARNING: %i crystals had no reflections in linear free "
		      "residual calculation\n", n_non_linear_free);
	}
	if ( n_non_log ) {
		ERROR("WARNING: %i crystals had no reflections in log "
		      "residual calculation\n", n_non_log);
	}
	if ( n_non_log_free ) {
		ERROR("WARNING: %i crystals had no reflections in log free "
		      "residual calculation\n", n_non_log_free);
	}

	if ( n_nan_linear ) {
		ERROR("WARNING: %i crystals had NaN linear residuals\n",
		      n_nan_linear);
	}
	if ( n_nan_linear_free ) {
		ERROR("WARNING: %i crystals had NaN linear free residuals\n",
		      n_nan_linear_free);
	}
	if ( n_nan_log ) {
		ERROR("WARNING: %i crystals had NaN log residuals\n",
		      n_nan_log);
	}
	if ( n_nan_log_free ) {
		ERROR("WARNING: %i crystals had NaN log free residuals\n",
		      n_nan_log_free);
	}

	*pn_used = n_used;
}


static void show_all_residuals(Crystal **crystals, int n_crystals,
                               RefList *full, int no_free)
{
	double dev, free_dev, log_dev, free_log_dev;
	int n;

	all_residuals(crystals, n_crystals, full, no_free,
	              &dev, &free_dev, &log_dev, &free_log_dev, &n);
	STATUS("Residuals:"
	       "  linear          linear free     log             log free\n");
	STATUS("         ");
	STATUS("%15e %15e %15e %15e    (%i crystals)\n",
	       dev, free_dev, log_dev, free_log_dev, n);
}


struct log_qargs
{
	int iter;
	int next;
	Crystal **crystals;
	int n_crystals;
	RefList *full;
	int scaleflags;
	PartialityModel pmodel;
	int n_done;
	const char *log_folder;
};


struct log_args
{
	Crystal *cr;
	RefList *full;
	int scaleflags;
	PartialityModel pmodel;
	int iter;
	int cnum;
	const char *log_folder;
};


static void *get_log_task(void *vp)
{
	struct log_qargs *qargs = vp;
	struct log_args *task;

	if ( qargs->next >= qargs->n_crystals ) return NULL;

	task = malloc(sizeof(struct log_args));
	if ( task == NULL ) return NULL;

	task->cr = qargs->crystals[qargs->next];
	task->full = qargs->full;
	task->iter = qargs->iter;
	task->cnum = qargs->next;
	task->scaleflags = qargs->scaleflags;
	task->pmodel = qargs->pmodel;
	task->log_folder = qargs->log_folder;

	qargs->next += 20;
	return task;
}


static void write_logs(void *vp, int cookie)
{
	struct log_args *args = vp;
	write_specgraph(args->cr, args->full, args->iter, args->cnum,
	                args->log_folder);
	write_gridscan(args->cr, args->full, args->iter, args->cnum,
	               args->scaleflags, args->pmodel, args->log_folder);
	write_test_logs(args->cr, args->full, args->iter, args->cnum,
	                args->log_folder);
}


static void done_log(void *vqargs, void *vp)
{
	struct log_args *task = vp;
	struct log_qargs *qargs = vqargs;
	qargs->n_done++;
	progress_bar(qargs->n_done, qargs->n_crystals/20, "Writing logs/grid scans");
	free(task);
}


static void write_logs_parallel(Crystal **crystals, int n_crystals,
                                RefList *full, int iter, int n_threads,
                                int scaleflags, PartialityModel pmodel,
                                const char *log_folder)
{
	struct log_qargs qargs;

	qargs.iter = iter;
	qargs.next = 0;
	qargs.full = full;
	qargs.crystals = crystals;
	qargs.n_done = 0;
	qargs.n_crystals = n_crystals;
	qargs.scaleflags = scaleflags;
	qargs.pmodel = pmodel;
	qargs.log_folder = log_folder;

	run_threads(n_threads, write_logs, get_log_task, done_log, &qargs,
	            n_crystals/20, 0, 0, 0);
}


struct stream_list
{
	int n;
	int max_n;
	const char **filenames;
	Stream **streams;
};


static int add_stream(const char *filename, struct stream_list *list)
{
	if ( list->n == list->max_n ) {
		const char **new_filenames = realloc(list->filenames,
		                                     (list->n+16)*sizeof(const char *));
		Stream **new_streams = realloc(list->streams,
		                                     (list->n+16)*sizeof(Stream *));
		if ( (new_filenames == NULL) || (new_streams == NULL) ) return 1;
		list->max_n += 16;
		list->filenames = new_filenames;
		list->streams = new_streams;
	}

	list->filenames[list->n] = filename;
	list->streams[list->n] = NULL;
	list->n++;
	return 0;
}


static void write_polarisation(FILE *fh, const char *name,
                               struct polarisation p)
{
	fprintf(fh, "    \"%s\": {\n", name);
	fprintf(fh, "      \"angle_from_horizontal_deg\": %f,\n", rad2deg(p.angle));
	fprintf(fh, "      \"fraction\": %f\n", p.fraction);
	fprintf(fh, "    },\n");
}

static void write_harvest_file(const char *filename,
                               const char *model, const char *symmetry,
                               int scale, int bscale, int postref,
                               int niter, int deltacchalf,
                               int min_measurements, double max_adu,
                               double min_res, double push_res,
                               struct polarisation p)
{
	FILE *fh;

	fh = fopen(filename, "w");
	if ( fh == NULL ) {
		ERROR("Unable to write parameter harvesting file.\n");
		return;
	}

	fprintf(fh, "{\n");
	fprintf(fh, "  \"merging\": {\n");
	write_str(fh, 1, "partiality_model", model);
	write_str(fh, 1, "symmetry", symmetry);
	write_bool(fh, 1, "scale", scale);
	write_bool(fh, 1, "Bscale", bscale);
	write_bool(fh, 1, "post_refine", postref);
	write_int(fh, 1, "num_iterations", niter);
	write_polarisation(fh, "polarisation_correction", p);
	write_bool(fh, 1, "deltaCChalf", deltacchalf);
	write_int(fh, 1, "min_measurements_per_unique_reflection", min_measurements);
	write_float(fh, 1, "max_adu", max_adu);
	write_float(fh, 1, "min_resolution_invm", min_res);
	write_float(fh, 0, "push_res_invm", push_res);
	fprintf(fh, "  }\n");
	fprintf(fh, "}\n");

	fclose(fh);
}


int main(int argc, char *argv[])
{
	int c;
	struct stream_list stream_list = {.n = 0,
	                                  .max_n = 0,
	                                  .filenames = NULL,
	                                  .streams = NULL};
	char *outfile = NULL;
	char *sym_str = NULL;
	SymOpList *sym;
	SymOpList *amb;
	SymOpList *w_sym;
	int nthreads = 1;
	int istream, icmd, icryst, itn;
	int n_iter = 10;
	RefList *full;
	int n_images = 0;
	int n_crystals = 0;
	int n_crystals_seen = 0;
	char cmdline[1024];
	int no_scale = 0;
	int no_Bscale = 0;
	int no_pr = 0;
	Crystal **crystals;
	char *pmodel_str = NULL;
	PartialityModel pmodel = PMODEL_XSPHERE;
	int min_measurements = 2;
	char *rval;
	struct polarisation polarisation = {.fraction = 1.0,
	                                    .angle = 0.0,
	                                    .disable = 0};
	int start_after = 0;
	int stop_after = 0;
	double max_adu = +INFINITY;
	char *sparams_fn = NULL;
	FILE *sparams_fh;
	double push_res = +INFINITY;
	gsl_rng *rng;
	int no_free = 0;
	int output_everycycle = 0;
	char *csplit_fn = NULL;
	struct custom_split *csplit = NULL;
	double max_B = 1e-18;
	char *rfile = NULL;
	RefList *reference = NULL;
	int no_logs = 0;
	char *w_sym_str = NULL;
	char *operator = NULL;
	double force_bandwidth = -1.0;
	double force_radius = -1.0;
	double force_lambda = -1.0;
	char *audit_info;
	int scaleflags = 0;
	double min_res = 0.0;
	int do_write_logs = 0;
	int no_deltacchalf = 0;
	char *harvest_file = NULL;
	char *log_folder = "pr-logs";

	/* Long options */
	const struct option longopts[] = {

		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,               'v'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"start-after",        1, NULL,               's'},
		{"stop-after",         1, NULL,               'f'},
		{"symmetry",           1, NULL,               'y'},
		{"iterations",         1, NULL,               'n'},
		{"model",              1, NULL,               'm'},

		{"min-measurements",   1, NULL,                2},
		{"max-adu",            1, NULL,                3},
		{"start-params",       1, NULL,                4},
		{"push-res",           1, NULL,                5},
		{"res-push",           1, NULL,                5}, /* compat */
		{"custom-split",       1, NULL,                6},
		{"max-rel-B",          1, NULL,                7},
		{"max-rel-b",          1, NULL,                7}, /* compat */
		{"reference",          1, NULL,                8}, /* ssshhh! */
		{"operator",           1, NULL,                9},
		{"force-bandwidth",    1, NULL,               10},
		{"force-radius",       1, NULL,               11},
		{"min-res",            1, NULL,               12},
		{"force-lambda",       1, NULL,               13},
		{"polarisation",       1, NULL,               14},
		{"polarization",       1, NULL,               14}, /* compat */
		{"no-polarisation",    0, NULL,               15},
		{"no-polarization",    0, NULL,               15}, /* compat */
		{"harvest-file",       1, NULL,               16},
		{"log-folder",         1, NULL,               17},

		{"no-scale",           0, &no_scale,           1},
		{"no-Bscale",          0, &no_Bscale,          1},
		{"no-pr",              0, &no_pr,              1},
		{"no-free",            0, &no_free,            1},
		{"output-every-cycle", 0, &output_everycycle,  1},
		{"no-logs",            0, &no_logs,            1},
		{"no-deltacchalf",     0, &no_deltacchalf,     1},

		{0, 0, NULL, 0}
	};

	cmdline[0] = '\0';
	for ( icmd=1; icmd<argc; icmd++ ) {
		strncat(cmdline, argv[icmd], 1023-strlen(cmdline));
		strncat(cmdline, " ", 1023-strlen(cmdline));
	}

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:o:g:b:y:n:j:m:w:",
	                        longopts, NULL)) != -1)
	{

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'v' :
			printf("CrystFEL: %s\n",
			       crystfel_version_string());
			printf("%s\n",
			       crystfel_licence_string());
			return 0;

			case 'i' :
			add_stream(optarg, &stream_list);
			break;

			case 'j' :
			nthreads = atoi(optarg);
			break;

			case 's' :
			errno = 0;
			start_after = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --start-after (%s)\n",
				      optarg);
				return 1;
			}
			break;

			case 'f' :
			errno = 0;
			stop_after = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --stop-after (%s)\n",
				      optarg);
				return 1;
			}
			break;

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 'o' :
			outfile = strdup(optarg);
			break;

			case 'n' :
			n_iter = atoi(optarg);
			break;

			case 'm' :
			pmodel_str = strdup(optarg);
			break;

			case 'w' :
			w_sym_str = strdup(optarg);
			break;

			case 2 :
			errno = 0;
			min_measurements = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-measurements.\n");
				return 1;
			}
			break;

			case 3 :
			errno = 0;
			max_adu = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --max-adu.\n");
				return 1;
			}
			break;

			case 4 :
			sparams_fn = strdup(optarg);
			break;

			case 5 :
			errno = 0;
			push_res = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --push-res.\n");
				return 1;
			}
			push_res = push_res*1e9;
			break;

			case 6 :
			csplit_fn = strdup(optarg);
			break;

			case 7 :
			errno = 0;
			max_B = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --max-rel-B.\n");
				return 1;
			}
			max_B = max_B * 1e-20;
			break;

			case 8 :
			rfile = strdup(optarg);
			break;

			case 9 :
			operator = strdup(optarg);
			break;

			case 10 :
			errno = 0;
			force_bandwidth = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --force-bandwidth.\n");
				return 1;
			}
			break;

			case 11 :
			errno = 0;
			force_radius = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --force-radius.\n");
				return 1;
			}
			force_radius *= 1e9;
			break;

			case 12  :
			errno = 0;
			min_res = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-res.\n");
				return 1;
			}
			min_res = 1e10/min_res;
			break;

			case 13 :
			errno = 0;
			force_lambda = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --force-lambda.\n");
				return 1;
			}
			force_lambda *= 1e-10;
			break;

			case 14 :
			polarisation = parse_polarisation(optarg);
			break;

			case 15 :
			polarisation = parse_polarisation("none");
			break;

			case 16 :
			harvest_file = strdup(optarg);
			break;

			case 17 :
			log_folder = strdup(optarg);
			break;

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	while ( optind < argc ) {
		add_stream(argv[optind++], &stream_list);
	}

	if ( nthreads < 1 ) {
		ERROR("Invalid number of threads.\n");
		return 1;
	}

	if ( stream_list.n == 0 ) {
		ERROR("Please give at least one input filename\n");
		return 1;
	}

	if ( outfile == NULL ) {
		outfile = strdup("partialator.hkl");
	}

	if ( sym_str == NULL ) sym_str = strdup("1");
	pointgroup_warning(sym_str);
	sym = get_pointgroup(sym_str);
	free(sym_str);
	if ( sym == NULL ) return 1;

	if ( (w_sym_str != NULL) && (operator != NULL) ) {
		ERROR("Specify the apparent symmetry (-w) or the operator, "
		      "not both.\n");
		return 1;
	}

	if ( w_sym_str == NULL ) {
		w_sym = NULL;
		amb = NULL;
	} else {
		pointgroup_warning(w_sym_str);
		w_sym = get_pointgroup(w_sym_str);
		free(w_sym_str);
		if ( w_sym == NULL ) return 1;
		amb = get_ambiguities(w_sym, sym);
		if ( amb == NULL ) {
			ERROR("Couldn't find ambiguity operator.\n");
			ERROR("Check that your values for -y and -w are "
			      "correct.\n");
			return 1;
		}

	}

	if ( operator ) {
		amb = parse_symmetry_operations(operator);
		if ( amb == NULL ) return 1;
		set_symmetry_name(amb, "Ambiguity");
	}

	if ( amb != NULL ) {
		STATUS("Indexing ambiguity resolution enabled.  "
		       "The ambiguity operation(s) are:\n");
		describe_symmetry(amb);
		/* In contrast to ambigator, partialator can deal with multiple
		 * ambiguities at once */
	}

	if ( pmodel_str != NULL ) {
		if ( strcmp(pmodel_str, "unity") == 0 ) {
			pmodel = PMODEL_UNITY;
		} else if ( strcmp(pmodel_str, "xsphere") == 0 ) {
			pmodel = PMODEL_XSPHERE;
		} else if ( strcmp(pmodel_str, "offset") == 0 ) {
			pmodel = PMODEL_OFFSET;
		} else if ( strcmp(pmodel_str, "random") == 0 ) {
			pmodel = PMODEL_RANDOM;
		} else if ( strcmp(pmodel_str, "ggpm") == 0 ) {
			pmodel = PMODEL_GGPM;
		} else {
			ERROR("Unknown partiality model '%s'.\n", pmodel_str);
			return 1;
		}
	}

	if ( (pmodel == PMODEL_UNITY) && !no_pr ) {
		no_pr = 1;
		STATUS("Setting --no-pr because we are not modelling "
		       "partialities (--model=unity).\n");
	}

	if ( no_Bscale ) {
		scaleflags |= SCALE_NO_B;
	}

	/* Decide whether or not to create stuff in the pr-logs folder */
	if ( !(no_logs || (no_pr && pmodel == PMODEL_UNITY)) ) {
		do_write_logs = 1;
	} else {
		do_write_logs = 0;
	}

	if ( do_write_logs ) {
		int r = mkdir(log_folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if ( r ) {
			if ( errno == EEXIST ) {
				ERROR("WARNING: Log folder (%s) already exists. "
				      "Beware of mixing old and new log files!\n",
				      log_folder);
			} else {
				ERROR("Failed to create log folder (%s).\n",
				      log_folder);
				return 1;
			}
		}
	} else {
		struct stat s;
		if ( stat(log_folder, &s) != -1 ) {
			ERROR("WARNING: Log folder (%s) exists, but I will not "
			      "write anything in it with these settings.\n",
			      log_folder);
		}
	}

	/* Read the custom split list (if applicable) */
	if ( csplit_fn != NULL ) {
		csplit = load_custom_split(csplit_fn);
		if ( csplit == NULL ) {
			ERROR("Failed to load custom split list.\n");
			return 1;
		}
		free(csplit_fn);
	}

	if ( rfile != NULL ) {
		RefList *rread;
		rread = read_reflections(rfile);
		if ( rread == NULL ) {
			ERROR("Failed to read reference reflections\n");
			return 1;
		}
		reference = asymmetric_indices(rread, sym);
		reflist_free(rread);
		ERROR("WARNING: Using an external reference.\n");
		ERROR("WARNING: If you publish a structure based on the result,"
		      " expect to have to retract your paper!\n");
	}

	if ( harvest_file != NULL ) {
		write_harvest_file(harvest_file,
		                   pmodel_str,
		                   symmetry_name(sym),
		                   1-no_scale, 1-no_Bscale, 1-no_pr,
		                   n_iter, 1-no_deltacchalf, min_measurements,
		                   max_adu, min_res, push_res, polarisation);
	}

	free(pmodel_str);

	gsl_set_error_handler_off();
	rng = gsl_rng_alloc(gsl_rng_mt19937);

	/* Fill in what we know about the images so far */
	n_images = 0;
	n_crystals = 0;
	n_crystals_seen = 0;
	crystals = NULL;
	if ( sparams_fn != NULL ) {
		char line[1024];
		sparams_fh = fopen(sparams_fn, "r");
		if ( sparams_fh == NULL ) {
			ERROR("Failed to open '%s'\n", sparams_fn);
			return 1;
		}
		if ( fgets(line, 1024, sparams_fh) == NULL ) {
			ERROR("Failed to read header from %s\n", sparams_fn);
			return 1;
		}
		STATUS("Reading initial scaling factors (G,B) from '%s'\n",
		       sparams_fn);
		free(sparams_fn);
	} else {
		sparams_fh = NULL;
	}

	audit_info = NULL;
	for ( istream=0; istream<stream_list.n; istream++ ) {

		Stream *st = stream_open_for_read(stream_list.filenames[istream]);
		if ( st == NULL ) {
			ERROR("Couldn't open %s\n", stream_list.filenames[istream]);
			return 1;
		}

		if ( audit_info == NULL ) {
			audit_info = stream_audit_info(st);
		}

		do {

			struct image *image;
			RefList *as;
			int i;

			image = stream_read_chunk(st, STREAM_REFLECTIONS);
			if ( image == NULL ) break;

			if ( isnan(image->div) || isnan(image->bw) ) {
				ERROR("Chunk doesn't contain beam parameters.\n");
				return 1;
			}

			for ( i=0; i<image->n_crystals; i++ ) {

				Crystal *cr;
				Crystal **crystals_new;
				RefList *cr_refl;
				RefList *cr_refl_raw;
				struct image *image_for_crystal;
				double lowest_r;

				n_crystals_seen++;
				if ( n_crystals_seen <= start_after ) continue;

				if ( crystal_get_resolution_limit(image->crystals[i]) < min_res ) continue;

				lowest_r = lowest_reflection(crystal_get_cell(image->crystals[i]));
				if ( crystal_get_profile_radius(image->crystals[i]) > 0.5*lowest_r ) {
					ERROR("Rejecting %s %s crystal %i because "
					      "profile radius is obviously too big (%e %e).\n",
					      image->filename, image->ev, i,
					      crystal_get_profile_radius(image->crystals[i]),
					      lowest_r);
					continue;
				}

				crystals_new = realloc(crystals,
				                       (n_crystals+1)*sizeof(Crystal *));
				if ( crystals_new == NULL ) {
					ERROR("Failed to allocate memory for crystal "
					      "list.\n");
					return 1;
				}
				crystals = crystals_new;
				crystals[n_crystals] = crystal_copy_deep(image->crystals[i]);
				cr = crystals[n_crystals];

				/* Create a completely new, separate image
				 * structure for this crystal. */
				image_for_crystal = image_new();
				if ( image_for_crystal == NULL ) {
					ERROR("Failed to allocate memory for image.\n");
					return 1;
				}

				crystal_set_image(cr, image_for_crystal);
				*image_for_crystal = *image;
				image_for_crystal->n_crystals = 1;
				image_for_crystal->crystals = malloc(sizeof(Crystal *));
				image_for_crystal->crystals[0] = cr;
				image_for_crystal->filename = strdup(image->filename);
				image_for_crystal->ev = safe_strdup(image->ev);
				image_for_crystal->detgeom = NULL;
				image_for_crystal->features = NULL;
				image_for_crystal->spectrum = NULL;
				image_for_crystal->n_cached_headers = 0;
				image_for_crystal->dp = NULL;
				image_for_crystal->bad = NULL;
				image_for_crystal->sat = NULL;

				/* This is the raw list of reflections */
				cr_refl_raw = crystal_get_reflections(cr);

				cr_refl = apply_max_adu(cr_refl_raw, max_adu);
				reflist_free(cr_refl_raw);

				if ( !no_free ) select_free_reflections(cr_refl, rng);

				as = asymmetric_indices(cr_refl, sym);
				crystal_set_reflections(cr, as);
				crystal_set_user_flag(cr, PRFLAG_OK);
				reflist_free(cr_refl);

				if ( set_initial_params(cr, sparams_fh, force_bandwidth,
				                        force_radius, force_lambda) )
				{
					ERROR("Failed to set initial parameters\n");
					return 1;
				}

				n_crystals++;

				if ( n_crystals == stop_after ) break;

			}

			image_free(image);

			n_images++;

			if ( n_images % 100 == 0 ) {
				display_progress(n_images, n_crystals);
			}

			if ( (stop_after>0) && (n_crystals == stop_after) ) break;

		} while ( 1 );

		stream_close(st);

	}

	free(stream_list.filenames);
	free(stream_list.streams);

	display_progress(n_images, n_crystals);
	fprintf(stderr, "\n");
	if ( sparams_fh != NULL ) fclose(sparams_fh);

	STATUS("Initial partiality calculation...\n");
	for ( icryst=0; icryst<n_crystals; icryst++ ) {

		Crystal *cr = crystals[icryst];
		update_predictions(cr);

		/* Polarisation correction requires kpred values */
		polarisation_correction(crystal_get_reflections(cr),
		                        crystal_get_cell(cr), polarisation);

		calculate_partialities(cr, pmodel);
	}

	if (csplit != NULL) check_csplit(crystals, n_crystals, csplit);

	/* Make a first pass at cutting out crap */
	//STATUS("Early rejection...\n");
	//early_rejection(crystals, n_crystals);

	/* Create reference data set if we don't already have one */
	if ( reference == NULL ) {
		if ( !no_scale ) {
			STATUS("Initial scaling...\n");
			scale_all(crystals, n_crystals, nthreads, scaleflags);
		}
		full = merge_intensities(crystals, n_crystals, nthreads,
		                         min_measurements, push_res, 1, 0);
	} else {
		full = reference;
	}

	/* Check rejection and write figures of merit */
	check_rejection(crystals, n_crystals, full, max_B, no_deltacchalf,
	                nthreads);
	show_all_residuals(crystals, n_crystals, full, no_free);

	if ( do_write_logs ) {
		write_pgraph(full, crystals, n_crystals, 0, "", log_folder);
		write_logs_parallel(crystals, n_crystals, full, 0, nthreads,
		                    scaleflags, pmodel, log_folder);
	}

	/* Iterate */
	for ( itn=0; itn<n_iter; itn++ ) {

		STATUS("Scaling and refinement cycle %i of %i\n", itn+1, n_iter);

		if ( !no_pr ) {
			refine_all(crystals, n_crystals, full, nthreads, pmodel,
			           itn+1, no_logs, sym, amb, scaleflags,
			           log_folder);
		}

		/* Create new reference if needed */
		if ( reference == NULL ) {
			free_contribs(full);
			reflist_free(full);
			if ( !no_scale ) {
				scale_all(crystals, n_crystals, nthreads,
				          scaleflags);
			}
			full = merge_intensities(crystals, n_crystals, nthreads,
			                         min_measurements,
			                         push_res, 1, 0);
		} /* else full still equals reference */

		check_rejection(crystals, n_crystals, full, max_B,
		                no_deltacchalf, nthreads);
		show_all_residuals(crystals, n_crystals, full, no_free);

		if ( do_write_logs ) {
			write_pgraph(full, crystals, n_crystals, itn+1, "",
			             log_folder);
		}

		if ( output_everycycle ) {

			char tmp[1024];
			snprintf(tmp, 1024, "iter%.2d_%s", itn+1, outfile);

			/* Output results */
			STATUS("Writing overall results to %s\n", tmp);
			write_reflist_2(tmp, full, sym);

			/* Output split results */
			write_split(crystals, n_crystals, tmp, nthreads, pmodel,
			            min_measurements, sym, push_res);

			/* Output custom split results */
			if ( csplit != NULL ) {
				int j;
				for ( j=0; j<csplit->n_datasets; j++ ) {
					write_custom_split(csplit, j, crystals,
					                   n_crystals, pmodel,
					                   min_measurements,
							   push_res, sym,
					                   nthreads, tmp);
				}
			}

		}
	}

	/* Final merge */
	STATUS("Final merge...\n");
	if ( reference == NULL ) {
		free_contribs(full);
		reflist_free(full);
		if ( !no_scale ) {
			scale_all(crystals, n_crystals, nthreads, scaleflags);
		}
		full = merge_intensities(crystals, n_crystals, nthreads,
		                         min_measurements,
		                         push_res, 1, 0);
	} else {
		full = merge_intensities(crystals, n_crystals, nthreads,
		                         min_measurements, push_res, 1, 0);
	}

	/* Write final figures of merit (no rejection any more) */
	show_all_residuals(crystals, n_crystals, full, no_free);
	if ( do_write_logs ) {
		write_pgraph(full, crystals, n_crystals, -1, "", log_folder);
		write_logs_parallel(crystals, n_crystals, full, -1, nthreads,
		                    scaleflags, pmodel, log_folder);
	}

	/* Output results */
	STATUS("Writing overall results to %s\n", outfile);
	reflist_add_command_and_version(full, argc, argv);
	if ( audit_info != NULL ) {
		reflist_add_notes(full, "Audit information from stream:");
		reflist_add_notes(full, audit_info);
		free(audit_info);
	}
	write_reflist_2(outfile, full, sym);

	/* Output split results */
	write_split(crystals, n_crystals, outfile, nthreads, pmodel,
	            min_measurements, sym, push_res);

	/* Output custom split results */
	if ( csplit != NULL ) {
		int i;
		for ( i=0; i<csplit->n_datasets; i++ ) {
			write_custom_split(csplit, i, crystals, n_crystals,
			                   pmodel, min_measurements, push_res,
			                   sym, nthreads, outfile);
		}
	}

	/* Clean up */
	gsl_rng_free(rng);
	for ( icryst=0; icryst<n_crystals; icryst++ ) {
		struct image *image = crystal_get_image(crystals[icryst]);
		image_free(image);
	}
	free_contribs(full);
	reflist_free(full);
	free_symoplist(sym);
	free(outfile);
	free(crystals);

	return 0;
}
