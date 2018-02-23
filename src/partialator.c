/*
 * partialator.c
 *
 * Scaling and post refinement for coherent nanocrystallography
 *
 * Copyright © 2012-2018 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2018 Thomas White <taw@physics.org>
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

#include "version.h"
#include "scaling.h"
#include "post-refinement.h"
#include "merge.h"
#include "rejection.h"


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
	Crystal *crystals1[n_crystals];
	Crystal *crystals2[n_crystals];
	int n_crystals1 = 0;
	int n_crystals2 = 0;
	int i;

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
		                  pmodel, min_measurements, push_res, 1);

	if ( split == NULL ) {
		ERROR("Not enough crystals for two way split!\n");
		return;
	}

	STATUS("Writing two-way split to %s ", tmp);
	write_reflist_2(tmp, split, sym);
	reflist_free(split);
	snprintf(tmp, 1024, "%s2", outfile);
	split = merge_intensities(crystals2, n_crystals2, nthreads,
		                  pmodel, min_measurements, push_res, 1);
	STATUS("and %s\n", tmp);
	write_reflist_2(tmp, split, sym);
	reflist_free(split);
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
		struct event *ev;
		char *evs;
		char *id;
		int dsn_crystal;

		fn = crystal_get_image(crystals[i])->filename;
		ev = crystal_get_image(crystals[i])->event;
		evs = get_event_string(ev);

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
		                  pmodel, min_measurements, push_res, 1);
	write_reflist_2(tmp, split, sym);
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
"      --no-pr                Disable orientation/physics refinement.\n"
"  -m, --model=<model>        Specify partiality model.\n"
"      --min-measurements=<n> Minimum number of measurements to require.\n"
"      --no-polarisation      Disable polarisation correction.\n"
"      --max-adu=<n>          Saturation value of detector.\n"
"      --push-res=<n>         Merge higher than apparent resolution cutoff.\n"
"  -j <n>                     Run <n> analyses in parallel.\n"
"      --no-free              Disable cross-validation (testing only).\n"
"      --custom-split         List of files for custom dataset splitting.\n"
"      --max-rel-B            Maximum allowable relative |B| factor.\n"
"      --no-logs              Do not write extensive log files.\n"
"  -w <pg>                    Apparent point group for resolving ambiguities.\n"
"      --operator=<op>        Indexing ambiguity operator for resolving.\n");
}


static signed int find_first_crystal(Crystal **crystals, int n_crystals,
                                     struct custom_split *csplit, int dsn)
{
	int i;

	for ( i=0; i<n_crystals; i++ ) {
		const char *fn;
		struct event *ev;
		char *evs;
		char *id;
		int dsn_crystal;

		fn = crystal_get_image(crystals[i])->filename;
		ev = crystal_get_image(crystals[i])->event;
		evs = get_event_string(ev);

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
		struct event *ev;
		char *evs;
		char *id;
		int dsn_crystal;

		fn = crystal_get_image(crystals[i])->filename;
		ev = crystal_get_image(crystals[i])->event;
		evs = get_event_string(ev);

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


static struct custom_split *load_custom_split(const char *filename)
{
	struct custom_split *csplit;
	FILE *fh;
	int i;

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
		int n;
		char **bits;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) break;

		chomp(line);
		notrail(line);
		n = assplode(line, " \t,", &bits, ASSPLODE_NONE);
		if ( n < 2 ) {
			ERROR("Badly formatted line '%s'\n", line);
			return NULL;
		}

		if ( n == 3 ) {
			/* Filename, event, dataset */
			fn = bits[0];
			evs = bits[1];
			ds = bits[2];
		} else {
			fn = bits[0];
			evs = get_event_string(NULL);
			ds = bits[1];
		}
		free(bits);

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
			copy_data(nrefl, refl);
		}
	}
	reflist_free(list);
	return nlist;
}


static void skip_to_end(FILE *fh)
{
	int c;
	do {
		c = fgetc(fh);
	} while ( (c != '\n') && (c != EOF) );
}


static int set_initial_params(Crystal *cr, FILE *fh)
{
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
	char ins[5];

	if ( inum >= 0 ) {
		snprintf(ins, 4, "%i", inum);
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
		double res, corr, Ipart;
		Reflection *match;

		if ( !get_flag(refl) ) continue;  /* Not free-flagged */

		/* Strong reflections only */
		if ( get_intensity(refl) < 3.0*get_esd_intensity(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		res = resolution(cell, h, k, l);
		if ( 2.0*res > crystal_get_resolution_limit(cr) ) continue;

		match = find_refl(full, h, k, l);
		if ( match == NULL ) continue;

		/* Calculated partiality */
		pcalc = get_partiality(refl);

		/* Observed partiality */
		corr = G * exp(B*res*res) * get_lorentz(refl);
		Ipart = get_intensity(refl) * corr;
		pobs = Ipart / get_intensity(match);

		fprintf(fh, "%5i %4i %4i %4i %e %e %8.3f %8.3f %s\n",
		        fr, h, k, l, 2*res, Ipart, pcalc, pobs, ins);

	}
}


static void write_pgraph(RefList *full, Crystal **crystals, int n_crystals,
                         signed int iter, const char *suff)
{
	FILE *fh;
	char tmp[256];
	int i;

	snprintf(tmp, 256, "pgraph%s.dat", suff);

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
                          double *presidual, double *pfree_residual,
                          double *plog_residual, double *pfree_log_residual)
{
	int i;

	*presidual = 0.0;
	*pfree_residual = 0.0;
	*plog_residual = 0.0;
	*pfree_log_residual = 0.0;

	for ( i=0; i<n_crystals; i++ ) {

		double r, free_r, log_r, free_log_r;

		if ( crystal_get_user_flag(crystals[i]) ) continue;

		r = residual(crystals[i], full, 0, NULL, NULL, 1);
		free_r = residual(crystals[i], full, 1, NULL, NULL, 1);
		log_r = log_residual(crystals[i], full, 0, NULL, NULL);
		free_log_r = log_residual(crystals[i], full, 1, NULL, NULL);

		if ( isnan(r) || isnan(free_r)
		  || isnan(log_r) || isnan(free_log_r) ) continue;

		*presidual += r;
		*pfree_residual += free_r;
		*plog_residual += log_r;
		*pfree_log_residual += free_log_r;
	}
}


static void show_all_residuals(Crystal **crystals, int n_crystals,
                               RefList *full)
{
	double dev, free_dev, log_dev, free_log_dev;

	all_residuals(crystals, n_crystals, full,
	              &dev, &free_dev, &log_dev, &free_log_dev);
	STATUS("Residuals:"
	       "  linear          linear free     log             log free\n");
	STATUS("         ");
	STATUS("%15e %15e %15e %15e\n", dev, free_dev, log_dev, free_log_dev);
}


struct log_qargs
{
	int iter;
	int next;
	Crystal **crystals;
	int n_crystals;
	RefList *full;
};


struct log_args
{
	Crystal *cr;
	RefList *full;
	int iter;
	int cnum;
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

	qargs->next += 20;
	return task;
}


static void write_logs(void *vp, int cookie)
{
	struct log_args *args = vp;
	write_specgraph(args->cr, args->full, args->iter, args->cnum);
	write_gridscan(args->cr, args->full, args->iter, args->cnum);
}


static void done_log(void *qargs, void *vp)
{
	struct log_args *task = vp;
	free(task);
}


static void write_logs_parallel(Crystal **crystals, int n_crystals,
                                RefList *full, int iter, int n_threads)
{
	struct log_qargs qargs;


	qargs.iter = iter;
	qargs.next = 0;
	qargs.full = full;
	qargs.crystals = crystals;
	qargs.n_crystals = n_crystals;

	STATUS("Writing logs...\n");
	run_threads(n_threads, write_logs, get_log_task, done_log, &qargs,
	            n_crystals/20, 0, 0, 0);
	STATUS("Done.\n");
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	char *outfile = NULL;
	char *sym_str = NULL;
	SymOpList *sym;
	SymOpList *amb;
	SymOpList *w_sym;
	int nthreads = 1;
	int i;
	int n_iter = 10;
	RefList *full;
	int n_images = 0;
	int n_crystals = 0;
	int n_crystals_seen = 0;
	char cmdline[1024];
	int no_scale = 0;
	int no_pr = 0;
	Stream *st;
	Crystal **crystals;
	char *pmodel_str = NULL;
	PartialityModel pmodel = PMODEL_XSPHERE;
	int min_measurements = 2;
	char *rval;
	int polarisation = 1;
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

		{"no-scale",           0, &no_scale,           1},
		{"no-pr",              0, &no_pr,              1},
		{"no-polarisation",    0, &polarisation,       0},
		{"no-polarization",    0, &polarisation,       0}, /* compat */
		{"polarisation",       0, &polarisation,       1},
		{"polarization",       0, &polarisation,       1}, /* compat */
		{"no-free",            0, &no_free,            1},
		{"output-every-cycle", 0, &output_everycycle,  1},
		{"no-logs",            0, &no_logs,            1},

		{0, 0, NULL, 0}
	};

	cmdline[0] = '\0';
	for ( i=1; i<argc; i++ ) {
		strncat(cmdline, argv[i], 1023-strlen(cmdline));
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
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			case 'i' :
			infile = strdup(optarg);
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

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( nthreads < 1 ) {
		ERROR("Invalid number of threads.\n");
		return 1;
	}

	if ( infile == NULL ) {
		ERROR("Please give the input filename (with -i)\n");
		return 1;
	}
	st = open_stream_for_read(infile);
	if ( st == NULL ) {
		ERROR("Failed to open input stream '%s'\n", infile);
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
		} else if ( strcmp(pmodel_str, "random") == 0 ) {
			pmodel = PMODEL_RANDOM;
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

	if ( !no_logs ) {
		int r = mkdir("pr-logs", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if ( r ) {
			if ( errno == EEXIST ) {
				ERROR("A folder called 'pr-logs' exists in the "
				      "working directory.\n");
				ERROR("Please delete or move it first.\n");
			} else {
				ERROR("Failed to create pr-logs folder.\n");
			}
			return 1;
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

	do {

		RefList *as;
		int i;
		struct image cur;

		cur.div = NAN;
		cur.bw = NAN;
		cur.det = NULL;
		if ( read_chunk_2(st, &cur, STREAM_READ_REFLECTIONS
		                            | STREAM_READ_UNITCELL) != 0 ) {
			break;
		}

		if ( isnan(cur.div) || isnan(cur.bw) ) {
			ERROR("Chunk doesn't contain beam parameters.\n");
			return 1;
		}

		for ( i=0; i<cur.n_crystals; i++ ) {

			Crystal *cr;
			Crystal **crystals_new;
			RefList *cr_refl;
			struct image *image;

			n_crystals_seen++;
			if ( n_crystals_seen <= start_after ) continue;

			crystals_new = realloc(crystals,
			                      (n_crystals+1)*sizeof(Crystal *));
			if ( crystals_new == NULL ) {
				ERROR("Failed to allocate memory for crystal "
				      "list.\n");
				return 1;
			}
			crystals = crystals_new;
			crystals[n_crystals] = cur.crystals[i];
			cr = crystals[n_crystals];

			image = malloc(sizeof(struct image));
			if ( image == NULL ) {
				ERROR("Failed to allocatea memory for image.\n");
				return 1;
			}

			crystal_set_image(cr, image);
			*image = cur;
			image->n_crystals = 1;
			image->crystals = &crystals[n_crystals];

			/* This is the raw list of reflections */
			cr_refl = crystal_get_reflections(cr);

			cr_refl = apply_max_adu(cr_refl, max_adu);

			if ( polarisation ) {
				polarisation_correction(cr_refl,
						        crystal_get_cell(cr),
						        image);
			}

			if ( !no_free ) select_free_reflections(cr_refl, rng);

			as = asymmetric_indices(cr_refl, sym);
			crystal_set_reflections(cr, as);
			crystal_set_user_flag(cr, PRFLAG_OK);
			reflist_free(cr_refl);

			if ( set_initial_params(cr, sparams_fh) ) {
				ERROR("Failed to set initial parameters\n");
				return 1;
			}

			n_crystals++;

			if ( n_crystals == stop_after ) break;

		}

		n_images++;

		if ( n_images % 100 == 0 ) {
			display_progress(n_images, n_crystals);
		}

		if ( (stop_after>0) && (n_crystals == stop_after) ) break;

	} while ( 1 );
	display_progress(n_images, n_crystals);
	fprintf(stderr, "\n");
	if ( sparams_fh != NULL ) fclose(sparams_fh);

	close_stream(st);

	STATUS("Initial partiality calculation...\n");
	for ( i=0; i<n_crystals; i++ ) {
		Crystal *cr = crystals[i];
		update_predictions(cr);
		calculate_partialities(cr, pmodel);
	}

	if (csplit != NULL) check_csplit(crystals, n_crystals, csplit);

	/* Make a first pass at cutting out crap */
	STATUS("Early rejection...\n");
	//early_rejection(crystals, n_crystals);

	/* Create reference data set if we don't already have one */
	if ( reference == NULL ) {
		full = merge_intensities(crystals, n_crystals, nthreads, pmodel,
		                         min_measurements, push_res, 1);
	} else {
		full = reference;
	}

	/* Scale everything to the reference */
	STATUS("Initial scaling...\n");
	if ( !no_scale ) {
		scale_all_to_reference(crystals, n_crystals, full, nthreads);
	}

	/* Check rejection and write figures of merit */
	check_rejection(crystals, n_crystals, full, max_B);
	show_all_residuals(crystals, n_crystals, full);
	if ( !no_pr && !no_logs ) {
		write_pgraph(full, crystals, n_crystals, 0, "");
		write_logs_parallel(crystals, n_crystals, full, 0, nthreads);
	}

	/* Iterate */
	for ( i=0; i<n_iter; i++ ) {

		STATUS("Scaling and refinement cycle %i of %i\n", i+1, n_iter);

		if ( !no_pr ) {
			refine_all(crystals, n_crystals, full, nthreads, pmodel,
			           0, i+1, no_logs, sym, amb);
		} else if ( !no_scale ) {
			scale_all_to_reference(crystals, n_crystals, full, nthreads);
		}

		/* Create new reference if needed */
		if ( reference == NULL ) {
			reflist_free(full);
			full = merge_intensities(crystals, n_crystals, nthreads,
			                         pmodel, min_measurements,
			                         push_res, 1);
		} /* else full still equals reference */

		check_rejection(crystals, n_crystals, full, max_B);
		show_all_residuals(crystals, n_crystals, full);

		if ( !no_logs ) {
			write_pgraph(full, crystals, n_crystals, i+1, "");
		}

		if ( output_everycycle ) {

			char tmp[1024];
			snprintf(tmp, 1024, "iter%.2d_%s", i+1, outfile);

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
		reflist_free(full);
		full = merge_intensities(crystals, n_crystals, nthreads,
		                         pmodel, min_measurements,
		                         push_res, 1);
	} else {
		scale_all_to_reference(crystals, n_crystals, reference, nthreads);
		full = merge_intensities(crystals, n_crystals, nthreads,
		                         pmodel, min_measurements, push_res, 1);
	}

	/* Write final figures of merit (no rejection any more) */
	show_all_residuals(crystals, n_crystals, full);
	if ( !no_pr && !no_logs ) {
		write_pgraph(full, crystals, n_crystals, -1, "");
		write_logs_parallel(crystals, n_crystals, full, -1, nthreads);
	}

	/* Output results */
	STATUS("Writing overall results to %s\n", outfile);
	reflist_add_command_and_version(full, argc, argv);
	write_reflist_2(outfile, full, sym);

	/* Output split results */
	write_split(crystals, n_crystals, outfile, nthreads, pmodel,
	            min_measurements, sym, push_res);

	/* Output custom split results */
	if ( csplit != NULL ) {
		for ( i=0; i<csplit->n_datasets; i++ ) {
			write_custom_split(csplit, i, crystals, n_crystals,
			                   pmodel, min_measurements, push_res,
			                   sym, nthreads, outfile);
		}
	}

	/* Clean up */
	gsl_rng_free(rng);
	for ( i=0; i<n_crystals; i++ ) {
		reflist_free(crystal_get_reflections(crystals[i]));
		crystal_free(crystals[i]);
	}
	reflist_free(full);
	free_symoplist(sym);
	free(outfile);
	free(crystals);
	free(infile);

	return 0;
}
