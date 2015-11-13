/*
 * geoptimiser.c
 *
 * Refine detector geometry
 *
 * Copyright © 2014-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2014-2015 Oleksandr Yefanov
 *   2014-2015 Valerio Mariani
 *   2014-2015 Thomas White <taw@physics.org>
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <float.h>

#ifdef HAVE_CAIRO
#ifdef HAVE_GTK
#define HAVE_SAVE_TO_PNG 1
#include <cairo.h>
#include <gdk/gdk.h>
#endif /* HAVE_GTK */
#endif /* HAVE_CAIRO */

#include <detector.h>
#include <stream.h>
#include <version.h>
#include <crystal.h>
#include <image.h>
#include <utils.h>
#include <render.h>

#include "hdfsee-render.h"

struct imagefeature;

static void show_help(const char *s)
{
	printf("Syntax: %s -i input.stream -g input.geom -o refined.geom "
	       "-c connected_rgcollection -q quadrant_rgcollection [options]\n",
	       s);
	printf(
"Refines detector geometry.\n"
"\n"
"  -h, --help                                   Display this help message.\n"
"\n"
"      --version                                Print CrystFEL version number and\n"
"                                                exit.\n"
"  -i, --input=<filename>                       Specify stream file to be used for \n"
"                                                geometry optimization.\n"
"  -g. --geometry=<file>                        Get detector geometry from file.\n"
"  -o, --output=<filename>                      Output stream.\n"
"  -q, --quadrants=<rg_coll>                    Rigid group collection for quadrants.\n"
"  -c, --connected=<rg_coll>                    Rigid group collection for connected\n"
"                                                ASICs.\n"
"      --no-error-maps                          Do not generate error map PNGs.\n"
"  -x, --min-num-peaks-per-pixel=<num>          Minimum number of peaks per pixel.\n"
"	                                         Default: 3. \n"
"      --min-num-peaks-per-panel=<num>          DEPRECATED. This option has been\n"
"                                                renamed to  --min-num-pixels-per-conn-group.\n"
"  -p, --min-num-pixels-per-conn-group=<num>    Minimum number of useful pixels pern"
"                                                connected group.\n"
"                         f                        Default: 100.\n"
"  -l, --most-freq-clen                         Use only the most frequent camera\n"
"                                                length.\n"
"  -s, --individual-dist-offset                 Use a distance offset for each panel.\n"
"                                                Default: whole-detector offset.\n"
"      --no-stretch                             Do not optimize distance offset.\n"
"                                                Default: distance offset is optimized\n"
"  -m  --max-peak-dist=<num>                    Maximum distance between predicted and\n"
"                                                detected peaks (in pixels)\n"
"                                                 Default: 4.0 pixels.\n"
);
}

struct geoptimiser_params
{
	char *infile;
	char *outfile;
	char *geometry_filename;
	int min_num_peaks_per_pix;
	int min_num_pix_per_conn_group;
	int only_best_distance;
	int nostretch;
	int individual_coffset;
	int error_maps;
	int enforce_cspad_layout;
	int no_cspad;
	double max_peak_dist;
	const char *command_line;
};


struct connected_data
{
	double sh_x;
	double sh_y;
	double cang;
	double cstr;
	int num_quad;
	int num_peaks_per_pixel;
	unsigned int n_peaks_in_conn;
	char *name;
};


struct pattern {
	ImageFeatureList *im_list;
	RefList *ref_list;
	double clen;
	UnitCell **unit_cells;
	int n_unit_cells;
	double lambda;
	char *filename;
};


struct pattern_list {
	struct pattern **patterns;
	int n_patterns;
};


struct single_pixel_displ
{
	double dfs;
	double dss;
	struct single_pixel_displ* ne;
};


struct pixel_displ_list
{
	struct single_pixel_displ *pix_displ_list;
	struct single_pixel_displ **curr_pix_displ;
	int *num_pix_displ;
};


struct connected_stretch_and_angles
{
	double *stretch_coeff;
	long *num_angles;
	long num_coeff;
};

struct cell_params
{
	double a;
	double b;
	double c;
	double alpha;
	double beta;
	double gamma;
};


struct avg_displacements
{
	double *displ_x;
	double *displ_y;
	double *displ_abs;
};


struct avg_rot_and_stretch
{
	double *aver_ang;
	double *aver_str;
	int *aver_num_ang;
};


struct avg_shift
{
	double *aver_x;
	double *aver_y;
	int *aver_num_sh;
};

struct pixel_maps
{
	double *pix_to_x;
	double *pix_to_y;
};


struct enhanced_det
{
	struct detector *det;
	int width;
	int height;
	int num_pix;
};


static void compute_x_y(struct detector *det, double fs, double ss,
                        double * x, double *y)
{
	struct panel *p;
	double xs, ys;
	double dfs, dss;

	p = find_panel(det, fs, ss);

	dss = ss-p->min_ss;
	dfs = fs-p->min_fs;

	xs = dfs*p->fsx + dss*p->ssx;
	ys = dfs*p->fsy + dss*p->ssy;

	*x = xs + p->cnx;
	*y = ys + p->cny;
}


static Reflection *find_closest_reflection(RefList *rlist,
                                         double fx, double fy,
                                         struct detector *det,
                                         double *d)
{

	double dmin = HUGE_VAL;
	Reflection *closest = NULL;
	Reflection *refl;
	RefListIterator *iter;


	for ( refl = first_refl(rlist, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {
		double ds;
		double rfs, rss;
		double rx, ry;

		get_detector_pos(refl, &rfs, &rss);

		compute_x_y(det, rfs, rss, &rx, &ry);

		ds = distance(rx, ry, fx, fy);

		if ( ds < dmin ) {
			dmin = ds;
			closest = refl;
		}

	}

	if ( dmin < HUGE_VAL ) {
		*d = dmin;
		return closest;
	}

	*d = +INFINITY;
	return NULL;
}


static double compute_average_clen (struct detector *det, char **clen_from,
                                    double *offset)
{

	int np, num_pan;
	double sum_clen;

	sum_clen = 0;
	num_pan = 0;

	for ( np=0; np<det->n_panels; np++ ) {

		struct panel p = det->panels[np];

		if ( p.clen_from != NULL ) {
			*clen_from = strdup(p.clen_from);
			*offset = p.coffset;
			return -1;
		} else {
			sum_clen += p.clen+p.coffset;
			num_pan += 1;
		}

	}

	return sum_clen/num_pan;

}


static double extract_f_from_stuff(const char *field_name,
                                   struct stuff_from_stream* stuff)
{
	int i;

	char field_name_plus_equal[256];
	snprintf(field_name_plus_equal, 256, "hdf5%s = ", field_name);

	for ( i=0; i<stuff->n_fields; i++ ) {

		if ( strncmp(stuff->fields[i], field_name_plus_equal,
		     strlen(field_name_plus_equal)) == 0 ) {
			return atof(stuff->fields[i]+
			       strlen(field_name_plus_equal));
		}
	}

	ERROR("Failed to recover camera length from stream file\n");
	return -1;
}


static struct pattern_list *read_patterns_from_steam_file(const char *infile,
                                                          struct detector *det)
{
	Stream *st;

	struct pattern_list *pattern_list;

	int max_patterns, n_chunks;

	n_chunks = 0;
	max_patterns = 0;

	pattern_list = malloc(sizeof(struct pattern_list));
	if ( pattern_list == NULL ) {
		ERROR("Failed to allocate memory for loaded patterns.\n");
		return NULL;
	}
	pattern_list->n_patterns =0;

	pattern_list->patterns = malloc(1024*sizeof(struct pattern*));
	if ( pattern_list->patterns == NULL ) {
		ERROR("Failed to allocate memory for loaded patterns.\n");
		free(pattern_list);
		return NULL;
	}
	pattern_list->n_patterns = 0;
	max_patterns = 1024;

	st = open_stream_for_read(infile);
	if ( st == NULL ) {
		ERROR("Failed to open input stream '%s'\n", infile);
		free(pattern_list->patterns);
		free(pattern_list);
		return NULL;
	}

	do {

		struct image cur;
		int i;

		cur.det = det;
		cur.stuff_from_stream = NULL;

		if ( read_chunk_2(st, &cur, STREAM_READ_REFLECTIONS
		     | STREAM_READ_PEAKS | STREAM_READ_UNITCELL) != 0 ) {
			break;
		}

		n_chunks +=1;

		if ( cur.n_crystals !=0 ) {

			struct pattern *patn;
			double avg_clen = 0.0;
			double offset = 0.0;
			char *clen_from;

			if ( pattern_list->n_patterns == max_patterns ) {

				struct pattern **patterns_new;

				patterns_new = realloc(pattern_list->patterns,
				               (max_patterns+1024)*
					       sizeof(struct pattern *));
				if ( patterns_new == NULL ) {
					ERROR("Failed to allocate "
					      "memory for loaded patterns.\n");
					free(pattern_list->patterns);
					free(pattern_list);
					return NULL;
				}

				max_patterns += 1024;
				pattern_list->patterns = patterns_new;
			}

			patn = malloc(sizeof(struct pattern));
			if ( patn == NULL ) {
				ERROR("Failed to allocate memory for loaded "
				      "patterns.\n");
				free(pattern_list->patterns);
				free(pattern_list);
				return NULL;
			}
			patn->filename = cur.filename;
			patn->unit_cells = NULL;
			patn->n_unit_cells = 0;
			patn->im_list = cur.features;
			patn->ref_list = reflist_new();

			clen_from = NULL;
			avg_clen = compute_average_clen(det, &clen_from, &offset);
			if ( avg_clen == -1 ) {
				avg_clen = extract_f_from_stuff(clen_from,
				           cur.stuff_from_stream)*1e-3;
				avg_clen += offset;
			}

			patn->clen = avg_clen;
			free(clen_from);

			patn->lambda = cur.lambda;

			for ( i=0; i<cur.n_crystals; i++ ) {

				RefList *crystal_reflist;
				Reflection *refl;
				RefListIterator *iter;
				UnitCell *cell;
				UnitCell **new_unit_cells;

				cell = crystal_get_cell(cur.crystals[i]);

				new_unit_cells = realloc(patn->unit_cells,
				                 (patn->n_unit_cells+1)*
						 sizeof(UnitCell *));
				if ( new_unit_cells == NULL ) {
					ERROR("Failed to allocate memory for "
					      "loaded patterns.\n");
					free(pattern_list->patterns);
					free(pattern_list);
					free(patn);
					return NULL;
				}

				new_unit_cells[patn->n_unit_cells] = cell;
				patn->n_unit_cells++;
				patn->unit_cells = new_unit_cells;

				crystal_reflist = crystal_get_reflections(
				                  cur.crystals[i]);

				for ( refl = first_refl(crystal_reflist, &iter);
				      refl != NULL;
				      refl = next_refl(refl, iter) )
				{
					Reflection *n;
					int h, k, l;

					get_indices(refl, &h, &k, &l);
					n = add_refl(patn->ref_list, h, k, l);
					copy_data(n, refl);
				}

			}

			pattern_list->patterns[pattern_list->n_patterns] = patn;
			pattern_list->n_patterns++;

			if ( pattern_list->n_patterns%1000 == 0 ) {
				STATUS("Loaded %i indexed patterns from %i "
				       "total patterns.\n",
				       pattern_list->n_patterns, ++n_chunks);
			}

		}

	} while ( 1 );

	close_stream(st);

	STATUS("Found %d indexed patterns in file %s (from a total of %d).\n",
	       pattern_list->n_patterns, infile, n_chunks );

	return pattern_list;
}


static struct rvec get_q_from_xyz(double rx, double ry, double dist, double l)
{

	struct rvec q;
	double r = sqrt(rx*rx + ry*ry);
	double twotheta = atan2(r, dist);
	double az = atan2(ry, rx);

	q.u = 1.0/(l*1e9) * sin(twotheta)*cos(az);
	q.v = 1.0/(l*1e9) * sin(twotheta)*sin(az);
	q.w = 1.0/(l*1e9) * (cos(twotheta) - 1.0);

	return q;
}


static struct cell_params *compute_avg_cell_parameters(struct pattern_list
						      *pattern_list)
{
	int numavc;
	int j, i;
	double minc[6];
	double maxc[6];
	double avg_cpar[6];

	for (j=0; j<6; j++) {
		minc[j] = 1e10;
		maxc[j] = -1e10;
	}
	numavc = 0;
	for (i=0; i<pattern_list->n_patterns; i++) {

		struct pattern *ptn;
		double cpar[6];
		int j, cri;

		ptn = pattern_list->patterns[i];

		for ( cri=0; cri<ptn->n_unit_cells; cri++ ) {

			cell_get_parameters(ptn->unit_cells[cri],
			                    &cpar[0],  // a
			                    &cpar[1],  // b
			                    &cpar[2],  // c
			                    &cpar[3],  // alpha
			                    &cpar[4],  // beta
			                    &cpar[5]); // gamma

			cpar[0] *= 1e9;
			cpar[1] *= 1e9;
			cpar[2] *= 1e9;
			cpar[3] = rad2deg(cpar[3]);
			cpar[4] = rad2deg(cpar[4]);
			cpar[5] = rad2deg(cpar[5]);

			for ( j=0; j<6; j++ ) {
				avg_cpar[j] += cpar[j];
				if (cpar[j]<minc[j]) minc[j] = cpar[j];
				if (cpar[j]>maxc[j]) maxc[j] = cpar[j];
			}
			numavc++;

		}

	}

	if ( numavc>0 ) {
		for ( j=0; j<6; j++ ) avg_cpar[j] /= numavc;
	}

	struct cell_params *cparams = malloc(sizeof(struct cell_params));

	cparams->a = avg_cpar[0];
	cparams->b = avg_cpar[1];
	cparams->c = avg_cpar[2];
	cparams->alpha = avg_cpar[3];
	cparams->beta = avg_cpar[4];
	cparams->gamma = avg_cpar[5];

	STATUS("Average cell coordinates:\n");
	STATUS("Average a, b, c (in nm): %6.3f, %6.3f, %6.3f\n",
	       cparams->a, cparams->b, cparams->c);
	STATUS("Minimum -Maximum a, b, c:\n"
	       "\t%6.3f - %6.3f,\n"
	       "\t%6.3f - %6.3f,\n"
	       "\t%6.3f - %6.3f\n",
	       minc[0], maxc[0], minc[1], maxc[1], minc[2], maxc[2]);
	STATUS("Average alpha,beta,gamma in degrees: %6.3f, %6.3f, %6.3f\n",
	       cparams->alpha, cparams->beta, cparams->gamma);
	STATUS("Minimum - Maximum alpha,beta,gamma in degrees:\n"
	       "\t%5.2f - %5.2f,\n"
	       "\t%5.2f - %5.2f,\n"
	       "\t%5.2f - %5.2f\n",
	       minc[3], maxc[3], minc[4], maxc[4], minc[5], maxc[5]);

	return cparams;
}


static double pick_clen_to_use(struct geoptimiser_params *gparams,
                                struct pattern_list *pattern_list,
                                double avg_res, struct cell_params *cparams)
{
	int cp, i, u;
	int num_clens;
	int max_clens;
	int best_clen;
	int *clens_population;
	double *clens;
	double *lambdas;
	double min_braggp_dist;
	double clen_to_use;
	struct rvec cqu;

	max_clens = 1024;

	clens = calloc(max_clens,sizeof(double));
	if ( clens == NULL ) {
		ERROR("Failed to allocate memory for clen calculation.\n");
		return -1.0;
	}
	clens_population = calloc(max_clens,sizeof(int));
	if ( clens_population == NULL ) {
		ERROR("Failed to allocate memory for clen calculation.\n");
		free(clens);
		return -1.0;
	}
	lambdas = calloc(max_clens,sizeof(double));
	if ( lambdas == NULL ) {
		ERROR("Failed to allocate memory for clen calculation.\n");
		free(clens);
		free(clens_population);
		return -1.0;
	}

	num_clens = 0;

	for ( cp=0; cp<pattern_list->n_patterns; cp++ )	{

		int i;
		int found = 0;

		for ( i=0; i<num_clens; i++ ) {
			if ( fabs(pattern_list->patterns[cp]->clen-clens[i])
			     <0.0001 ) {
				clens_population[i]++;
				lambdas[i] += pattern_list->patterns[cp]->lambda;
				found = 1;
				break;
			}
		}

		if ( found == 1) continue;

		if ( num_clens == max_clens ) {

			int *clens_population_new;
			double *clens_new;
			double *lambdas_new;

			clens_population_new = realloc(clens_population,
			                       (max_clens+1024)*sizeof(int));
			clens_new = realloc(clens_population,
			            (max_clens+1024)*sizeof(double));
			lambdas_new = realloc(clens_population,
				     (max_clens+1024)*sizeof(double));

			if ( clens_new == NULL ||  clens_population_new == NULL
			     || lambdas_new == NULL) {
				ERROR("Failed to allocate memory for "
				      "camera length list\n");
				free(clens);
				free(clens_population);
				free(lambdas);
				return -1.0;
			}

			max_clens += 1024;
			clens_population_new = clens_population;
			clens = clens_new;
			lambdas = lambdas_new;
		}

		clens[num_clens] = pattern_list->patterns[cp]->clen;
		clens_population[num_clens] = 1;
		lambdas[num_clens] = pattern_list->patterns[cp]->lambda;
		num_clens++;

	}


	for ( u=0; u<num_clens; u++ ) {
		lambdas[u] /= clens_population[u];
	}

	if ( num_clens == 1 ) {
		STATUS("All patterns have the same camera length: %f m.\n",
		       clens[0]);
	} else {
		STATUS("%i different camera lengths were found for the input "
		       "patterns:\n", num_clens);
	}

	best_clen = 0;
	clen_to_use = clens[0];
	for ( i=0; i<num_clens; i++) {
		if ( clens_population[i] >0) {
			cqu = get_q_from_xyz(1.0/avg_res, 0, clens[i], lambdas[i]);

			min_braggp_dist = fmin(fmin((1.0/cqu.u)/cparams->a,
					       (1.0/cqu.u)/cparams->b),
					       (1.0/cqu.u)/cparams->c);
			STATUS("Camera length %0.4f m was found %i times.\n"
			       "Minimum inter-bragg peak distance (based on "
			       "average cell parameters): %0.1f pixels.\n",
			       clens[i], clens_population[i],
			       min_braggp_dist);
			if ( min_braggp_dist<1.2*gparams->max_peak_dist ) {
				STATUS("WARNING: The distance between Bragg "
				       "peaks is too small: "
				       "%0.1f < 1.2*%0.1f pixels.\n",
			               min_braggp_dist,
				       gparams->max_peak_dist);
			}
			if ( clens_population[i] > clens_population[best_clen] ) {
				best_clen = i;
				clen_to_use = clens[best_clen];
			}
		}
	}

	if ( gparams->only_best_distance ) {
		STATUS("Only %i patterns with  CLEN=%0.4f m will be used.\n",
		       clens_population[best_clen], clen_to_use);
	}

	free(clens);
	free(lambdas);
	free(clens_population);

	return clen_to_use;
}


static double comp_median(double *arr, long n)
{

	long low, high, median, middle, ll, hh;
	double A;

	if (n<1) return 0.0;

	low = 0;
	high = n-1 ;
	median = (low + high) / 2;
	while (1) {
		if (high <= low) return arr[median] ;

		if (high == low + 1) {
			if (arr[low] > arr[high]) {
				A = arr[low];
				arr[low] = arr[high];
				arr[high] = A;
			}
			return arr[median] ;
		}

		// Find median of low, middle and high items; swap into position
		// low
		middle = (low + high) / 2;
		if ( arr[middle]>arr[high] ) {
			A = arr[middle];
			arr[middle] = arr[high];
			arr[high] = A;
		}
		if ( arr[low]>arr[high] ) {
			A = arr[low];
			arr[low] = arr[high];
			arr[high] = A;
		}
		if ( arr[middle]>arr[low] ) {
			A = arr[middle];
			arr[middle] = arr[low];
			arr[low] = A;
		}

		// Swap low item (now in position middle) into position
		// (low+1)
		A = arr[middle];
		arr[middle] = arr[low+1];
		arr[low+1] = A;

		// Nibble from each end towards middle, swapping items when
		// stuck
		ll = low + 1;
		hh = high;
		while (1) {
			do ll++; while (arr[low] > arr[ll]);
			do hh--; while (arr[hh]  > arr[low]);

			if (hh < ll) break;

			A = arr[ll];
			arr[ll] = arr[hh];
			arr[hh] = A;
		}

		A = arr[low];
		arr[low] = arr[hh];
		arr[hh] = A;

		/* Re-set active partition */
		if ( hh<=median ) low = ll;
		if ( hh>=median ) high = hh-1;
	}

	return 0.0;
}


static int find_quad_for_connected(struct rigid_group *rg,
                                   struct rg_collection *quadrants)
{
	struct panel *p;
	int qi;

	// The quadrant for a group of connected panels is the quadrant to which
	// the first panel in the connected set belong
	p = rg->panels[0];

	for ( qi=0; qi<quadrants->n_rigid_groups; qi++ ) {
		if ( panel_is_in_rigid_group(quadrants->rigid_groups[qi], p) ) {
			return qi;
		}
	}

	// Hopefully never reached
	return -1;
}


static int fill_avg_pixel_displ(struct pixel_displ_list *pix_displ_list,
				int pix_index,
				struct avg_displacements *avg_displ)
{
	double *list_fs_displ;
	double *list_ss_displ;
	int count = 0;
	int ei;

	list_fs_displ = calloc(pix_displ_list->num_pix_displ[pix_index],
			sizeof(double));
	if ( list_fs_displ == NULL ) {
		ERROR("Failed to allocate memory for pixel statistics.\n");
		return 1;
	}
	list_ss_displ = calloc(pix_displ_list->num_pix_displ[pix_index],
			sizeof(double));
	if ( list_ss_displ == NULL ) {
		ERROR("Failed to allocate memory for pixel statistics.\n");
		free(list_fs_displ);
		return 1;
	}

	pix_displ_list->curr_pix_displ[pix_index] =
	                &pix_displ_list->pix_displ_list[pix_index];

	for ( ei=0; ei<pix_displ_list->num_pix_displ[pix_index]; ei++ ) {

		struct single_pixel_displ * pixel_displ_entry;
		pixel_displ_entry = pix_displ_list->curr_pix_displ[pix_index];

		if ( pixel_displ_entry->dfs == -10000.0 ) break;
		list_fs_displ[count] = pixel_displ_entry->dfs;
		list_ss_displ[count] = pixel_displ_entry->dss;
		count++;
		if ( pixel_displ_entry->ne == NULL ) {
			break;
		} else {
			pix_displ_list->curr_pix_displ[pix_index] =
				pix_displ_list->curr_pix_displ[pix_index]->ne;
		}
	}

	if ( count < 1 ) {
		free(list_fs_displ);
		free(list_ss_displ);
		return 0;
	}

	avg_displ->displ_x[pix_index] = comp_median(list_fs_displ, count);
	avg_displ->displ_y[pix_index] = comp_median(list_ss_displ, count);
	avg_displ->displ_abs[pix_index] = modulus2d(avg_displ->displ_x[pix_index],
						    avg_displ->displ_y[pix_index]);
	free(list_fs_displ);
	free(list_ss_displ);

	return 0;
}


static int allocate_next_element(struct single_pixel_displ** curr_pix_displ,
				 int pix_index)
{
	curr_pix_displ[pix_index]->ne = malloc(sizeof(struct single_pixel_displ));
	if ( curr_pix_displ[pix_index]->ne == NULL ) {
		ERROR("Failed to allocate memory for pixel statistics.\n");
		return 1;
	}

	curr_pix_displ[pix_index] = curr_pix_displ[pix_index]->ne;

	return 0;
}


static int add_distance_to_list(struct enhanced_det *edet,
				struct imagefeature *imfe,
				struct pixel_displ_list *pix_displ_list,
				Reflection *refl, double fx, double fy)
{
	int pix_index;
	double rfs, rss;
	double crx, cry;

	pix_index = ((int)rint(imfe->fs) + edet->width*(int)rint(imfe->ss));

	if ( pix_displ_list->num_pix_displ[pix_index]>0 ) {

		int ret;

		ret = allocate_next_element(pix_displ_list->curr_pix_displ,
		                            pix_index);

		if ( ret != 0) return ret;

	}

	get_detector_pos(refl, &rfs, &rss);
	compute_x_y(edet->det, rfs, rss, &crx, &cry);
	pix_displ_list->curr_pix_displ[pix_index]->dfs = (fx-crx);
	pix_displ_list->curr_pix_displ[pix_index]->dss = (fy-cry);
	pix_displ_list->curr_pix_displ[pix_index]->ne = NULL;
	pix_displ_list->num_pix_displ[pix_index]++;

	return 0;
}


static int count_pixels_with_min_peaks( struct panel *p, int min_num_peaks,
                                        struct pixel_displ_list *pix_displ_list,
                                        struct enhanced_det *edet)
{
	int pix_index;
	int pixel_count = 0;
	int ifs, iss;

	for ( ifs=p->min_fs; ifs<p->max_fs+1; ifs++ ) {
		for ( iss=p->min_ss; iss<p->max_ss+1; iss++ ) {

			pix_index = ifs+edet->width*iss;

			if ( pix_displ_list->num_pix_displ[pix_index] >=
			     min_num_peaks ) {
				pixel_count += 1;
			}

		}
	}
	return pixel_count;
}


static void adjust_min_peaks_per_conn(struct enhanced_det *edet,
                                      struct rg_collection *connected,
                                      struct geoptimiser_params *gparams,
                                      struct connected_data *conn_data,
                                      struct pixel_displ_list *pix_displ_list)
{

	int min_num_peaks, di, ip;
	struct panel *p;


	STATUS("Adjusting the minimum number of measurements per pixel in "
	       "order to have enough measurements for each connected group.\n");
	for ( di=0; di<connected->n_rigid_groups; di++ ) {

		for ( min_num_peaks=gparams->min_num_peaks_per_pix;
		      min_num_peaks>0; min_num_peaks-- ) {

			int di_count = 0;

			for (ip=0; ip<connected->rigid_groups[di]->n_panels;
			     ip++) {

				int pix_count;

				p = connected->rigid_groups[di]->panels[ip];


				pix_count = count_pixels_with_min_peaks(p,
				                                min_num_peaks,
				                                pix_displ_list,
				                                edet);

				di_count += pix_count;
			}
			conn_data[di].n_peaks_in_conn = di_count;
			if ( di_count >=
			     gparams->min_num_pix_per_conn_group ) {
				conn_data[di].num_peaks_per_pixel =
				              min_num_peaks;
				break;
			}
		}

		STATUS("Minimum number of measurement "
		       "per pixel for connected group "
		       "%s has been set to %i\n",
		       conn_data[di].name,
		       conn_data[di].num_peaks_per_pixel);
	}
}


static int compute_pixel_displacements(struct pattern_list *pattern_list,
                                       struct enhanced_det *edet,
                                       struct rg_collection *connected,
                                       struct geoptimiser_params *gparams,
                                       double clen_to_use,
                                       struct connected_data *conn_data,
                                       struct avg_displacements *avg_displ,
                                       struct pixel_displ_list *pix_displ_list)
{
	int cp, ich;


	STATUS("Computing pixel displacements.\n");
	for ( cp=0; cp<pattern_list->n_patterns; cp++ ) {

		ImageFeatureList *flist = pattern_list->patterns[cp]->im_list;

		if ( gparams->only_best_distance ) {
			if ( fabs(pattern_list->patterns[cp]->clen-clen_to_use) >
			     0.0001 ) {
				continue;
			}
		}

		for ( ich=0;
		      ich<image_feature_count(pattern_list->patterns[cp]->im_list);
		      ich++ ) {

			double min_dist;
			double fx, fy;
			Reflection *refl;
			int ret;

			RefList *rlist = pattern_list->patterns[cp]->ref_list;

			struct imagefeature *imfe = image_get_feature(flist, ich);
			compute_x_y(edet->det, imfe->fs, imfe->ss, &fx, &fy);

			refl = find_closest_reflection(rlist, fx, fy, edet->det,
			                               &min_dist);

			if ( refl == NULL ) continue;

			if ( min_dist < gparams->max_peak_dist ) {


				ret = add_distance_to_list(edet, imfe,
				                           pix_displ_list,
				                           refl, fx, fy);

				if ( ret != 0 ) return ret;

			}
		}
	}

	return 0;
}


static int compute_avg_pix_displ(struct pixel_displ_list *pix_displ_list,
                                 int num_peaks_per_pixel,int pix_index,
                                 struct avg_displacements *avg_displ)
{
	int ret;

	if ( pix_displ_list->num_pix_displ[pix_index] >=
	     num_peaks_per_pixel ) {

		ret = fill_avg_pixel_displ(pix_displ_list,
		                           pix_index,
		                           avg_displ);
		if ( ret !=0 ) return ret;
	} else {
		avg_displ->displ_x[pix_index] = -10000.0;
		avg_displ->displ_y[pix_index] = -10000.0;
		avg_displ->displ_abs[pix_index] = -10000.0;
	}

	return 0;

}


static int compute_avg_displacements(struct enhanced_det *edet,
                                     struct rg_collection *connected,
                                     struct pixel_displ_list *pix_displ_list,
                                     struct connected_data *conn_data,
                                     struct avg_displacements *avg_displ)
{
	int di, ip, ifs, iss;
	int pix_index, ret;
	struct panel *p;

	for ( di=0; di<connected->n_rigid_groups; di++ ) {

		for (ip=0; ip<connected->rigid_groups[di]->n_panels;
		     ip++) {

			p = connected->rigid_groups[di]->panels[ip];

			for ( ifs=p->min_fs; ifs<p->max_fs+1; ifs++ ) {
				for ( iss=p->min_ss; iss<p->max_ss+1; iss++ ) {

					pix_index = ifs+edet->width*iss;

					ret = compute_avg_pix_displ(pix_displ_list,
					      conn_data[di].num_peaks_per_pixel,
					      pix_index, avg_displ);

					if ( ret != 0 ) return ret;

				}
			}
		}

	}
	return 0;
}


static double compute_error(struct rg_collection *connected,
                            struct enhanced_det* edet,
                            struct connected_data *conn_data,
                            int *num_pix_displ,
                            double *displ_abs)
{
	double total_error = 0;
	int num_total_error = 0;
	int di, ip;

	for ( di=0;di<connected->n_rigid_groups;di++ ) {

		struct panel *p;
		double connected_error = 0;
		int num_connected_error = 0;
		int ifs, iss;

		for (ip=0; ip<connected->rigid_groups[di]->n_panels; ip++) {

			p = connected->rigid_groups[di]->panels[ip];

			for (ifs=p->min_fs; ifs<p->max_fs+1; ifs++) {
				for (iss=p->min_ss; iss<p->max_ss+1; iss++) {

					int pix_index;
					pix_index = ifs+edet->width*iss;

					if ( num_pix_displ[pix_index]>=
					     conn_data[di].num_peaks_per_pixel ) {

						double cer;

						cer = displ_abs[pix_index]*
						      displ_abs[pix_index];
						connected_error += cer;
						num_connected_error++;
						total_error += cer;
						num_total_error++;
					}
				}
			}
		}

		if ( num_connected_error > 0 ) {

			connected_error /= (double)num_connected_error;
			connected_error = sqrt(connected_error);

			STATUS("Error for connected group %s: %d pixels with "
			       "more than %d peaks: <delta^2> = %0.4f pixels.\n",
			       conn_data[di].name, num_connected_error,
			       conn_data[di].num_peaks_per_pixel,
			       connected_error);
		}
	}

	if ( num_total_error>0 ) {
		total_error /= (double)num_total_error;
		total_error = sqrt(total_error);
	} else {
		total_error = -1;
	}

	return total_error;
}


static struct pixel_maps *initialize_pixel_maps(struct enhanced_det *edet)
{

	int pi;
	struct pixel_maps *pixel_maps;

	pixel_maps = malloc(sizeof(struct pixel_maps));
	pixel_maps->pix_to_x = malloc(edet->num_pix*sizeof(double));
	if ( pixel_maps->pix_to_x == NULL ) {
		free(pixel_maps);
		return NULL;
	}

	pixel_maps->pix_to_y = malloc(edet->num_pix*sizeof(double));
	if ( pixel_maps->pix_to_x == NULL ) {
		ERROR("Failed to allocate memory for pixel maps.\n");
		free(pixel_maps->pix_to_x);
		free(pixel_maps);
		return NULL;
	}

	for ( pi=0; pi<edet->det->n_panels; pi++ ) {

		struct panel *p;
		int iss, ifs;

		p = &edet->det->panels[pi];

		for (iss=p->min_ss; iss < p->max_ss+1; iss++) {
			for (ifs=p->min_fs; ifs < p->max_fs+1; ifs++) {

				double xs, ys;
				int pix_index;

				pix_index = iss*edet->width+ifs;

				compute_x_y(edet->det, ifs, iss, &xs, &ys);
				pixel_maps->pix_to_x[pix_index] = xs;
				pixel_maps->pix_to_y[pix_index] = ys;
			}
		}
	}

	return pixel_maps;
}


void free_pixel_maps(struct pixel_maps* pixel_maps)
{
	free(pixel_maps->pix_to_x);
	free(pixel_maps->pix_to_y);
	free(pixel_maps);
}



struct avg_rot_and_stretch* initialize_avg_rot_stretch(
				 int num_rigid_groups)
{
	int i;

	struct avg_rot_and_stretch *avg_rot_str;

	avg_rot_str = malloc(sizeof(struct avg_rot_and_stretch));
	if ( avg_rot_str == NULL ) {
		ERROR("Failed to allocate memory to correct empty panels.|n");
		return NULL;
	}

	avg_rot_str->aver_ang = malloc(num_rigid_groups*sizeof(double));
	if ( avg_rot_str->aver_ang == NULL ) {
		ERROR("Failed to allocate memory to correct empty panels.|n");
		free(avg_rot_str);
		return NULL;
	}

	avg_rot_str->aver_str = malloc(num_rigid_groups*sizeof(double));
	if ( avg_rot_str->aver_str == NULL ) {
		ERROR("Failed to allocate memory to correct empty panels.|n");
		free(avg_rot_str->aver_ang);
		free(avg_rot_str);
		return NULL;
	}

	avg_rot_str->aver_num_ang = malloc(num_rigid_groups*sizeof(int));
	if ( avg_rot_str->aver_num_ang == NULL ) {
		ERROR("Failed to allocate memory to correct empty panels.|n");
		free(avg_rot_str->aver_ang);
		free(avg_rot_str->aver_str);
		free(avg_rot_str);
		return NULL;
	}

	for (i=0; i<num_rigid_groups; i++) {
		avg_rot_str->aver_ang[i] = 0.0;
		avg_rot_str->aver_str[i] = 0.0;
		avg_rot_str->aver_num_ang[i] = 0;
	}

	return avg_rot_str;
}


void free_avg_angle_and_stretches(struct avg_rot_and_stretch* avg_rot_str)
{
	free(avg_rot_str->aver_ang);
	free(avg_rot_str->aver_str);
	free(avg_rot_str->aver_num_ang);
}


static int compute_rot_stretch_for_empty_panels(
                                            struct rg_collection *quadrants,
                                            struct rg_collection *connected,
                                            struct geoptimiser_params *gparams,
                                            struct connected_data *conn_data)
{
	struct avg_rot_and_stretch *avg_rot_str;
	int di,i;

	STATUS("Computing rotation and elongation corrections for groups "
	       "without the required number of measurements.\n");

	avg_rot_str = initialize_avg_rot_stretch(quadrants->n_rigid_groups);
	if ( avg_rot_str == NULL ) {
		return 1;
	}

	for (di=0; di<connected->n_rigid_groups; di++) {
		if ( conn_data[di].n_peaks_in_conn >=
		     gparams->min_num_pix_per_conn_group ) {
			avg_rot_str->aver_ang[conn_data[di].num_quad] +=
			                      conn_data[di].cang;
			avg_rot_str->aver_str[conn_data[di].num_quad] +=
			                      conn_data[di].cstr;
			avg_rot_str->aver_num_ang[conn_data[di].num_quad]++;
		}
	}

	for ( i=0; i<quadrants->n_rigid_groups; i++ ) {
		if ( avg_rot_str->aver_num_ang[i] > 0 ) {
			avg_rot_str->aver_ang[i] /=
					   (double)avg_rot_str->aver_num_ang[i];
			avg_rot_str->aver_str[i] /=
					   (double)avg_rot_str->aver_num_ang[i];
		}
	}

	for ( di=0; di<connected->n_rigid_groups; di++ ) {

		if ( conn_data[di].n_peaks_in_conn <
		     gparams->min_num_pix_per_conn_group ) {
			if ( avg_rot_str->aver_num_ang[conn_data[di].num_quad]
			     > 0) {
				conn_data[di].cang =
				  avg_rot_str->aver_ang[conn_data[di].num_quad];
				conn_data[di].cstr =
				  avg_rot_str->aver_str[conn_data[di].num_quad];
				STATUS("Connected group %s has only %i useful "
				       "pixels. Using average angle: %0.4f "
				       "degrees\n", conn_data[di].name,
				       conn_data[di].n_peaks_in_conn,
				       conn_data[di].cang);
			} else {
				STATUS("Connected group %s does not have enough "
				       " peaks (%i). It will not be moved.\n",
				       conn_data[di].name,
				       conn_data[di].n_peaks_in_conn);
			}
		}
	}

	return 0;
}


static void correct_corner_coordinates(struct rg_collection *connected,
				       struct connected_data *conn_data)
{

	int di, ip;

	for ( di=0; di<connected->n_rigid_groups; di++ ) {
		for ( ip=0; ip<connected->rigid_groups[di]->n_panels; ip++ ) {

			struct panel *p;

			p = connected->rigid_groups[di]->panels[ip];

			if ( ip == 0 ) {

				p->cnx *= conn_data[di].cstr;
				p->cny *= conn_data[di].cstr;

			} else {

				struct panel *p0;
				double delta_x, delta_y;

				p0 = connected->rigid_groups[di]->panels[0];

				delta_x = p->cnx-p0->cnx/conn_data[di].cstr;
				delta_y = p->cny-p0->cny/conn_data[di].cstr;

				p->cnx = p0->cnx + delta_x;
				p->cny = p0->cny + delta_y;

			}
		}
	}
}


static void correct_rotation_and_stretch(struct rg_collection *connected,
                                         struct detector *det,
                                         struct connected_data *conn_data,
                                         double clen_to_use,
                                         double stretch_coeff,
                                         struct geoptimiser_params *gparams)
{

	int di, ip;

	STATUS("Applying rotation and stretch corrections.\n");
	for ( di=0; di<connected->n_rigid_groups; di++ ) {
		for ( ip=0; ip<connected->rigid_groups[di]->n_panels; ip++ ) {

			struct panel *p;
			double new_fsx, new_fsy, new_ssx, new_ssy;

			p = connected->rigid_groups[di]->panels[ip];

			new_fsx = p->fsx*cos(conn_data[di].cang)-
			          p->fsy*sin(conn_data[di].cang);
			new_fsy = p->fsx*sin(conn_data[di].cang)+
			          p->fsy*cos(conn_data[di].cang);
			new_ssx = p->ssx*cos(conn_data[di].cang)-
			          p->ssy*sin(conn_data[di].cang);
			new_ssy = p->ssx*sin(conn_data[di].cang)+
			          p->ssy*cos(conn_data[di].cang);
			p->fsx = new_fsx;
			p->fsy = new_fsy;
			p->ssx = new_ssx;
			p->ssy = new_ssy;
		}
	}


	if ( gparams->individual_coffset ) {

		STATUS("Using individual distances for rigid panels.\n");
		for ( di=0; di<connected->n_rigid_groups; di++ ) {
			for ( ip=0; ip<connected->rigid_groups[di]->n_panels;
			      ip++ ) {

				struct panel *p;

				p = connected->rigid_groups[di]->panels[ip];
				p->coffset -= (1.0-conn_data[di].cstr)*clen_to_use;
			}
		}

	} else {

		int pi;

		for ( di=0; di<connected->n_rigid_groups; di++ ) {
			conn_data[di].cstr = stretch_coeff;
		}

		for ( pi=0; pi<det->n_panels; pi++) {
			det->panels[pi].coffset -= clen_to_use*(1.0-stretch_coeff);
		}
		STATUS("Using a single offset distance for the whole "
		       "detector: %f m. Stretch ceofficient: %0.4f\n",
		       det->panels[0].coffset, stretch_coeff);
	}

	correct_corner_coordinates(connected, conn_data);
}


static void adjust_panel(struct connected_data *conn_data, int di, int ip,
                         struct rg_collection *connected,
                         struct pixel_maps *pixel_maps,
                         struct pixel_maps *recomputed_pixel_maps,
                         struct avg_displacements *avg_displ,
                         double stretch_coeff,
                         struct enhanced_det *edet, int *num_pix_displ)
{
	double c_stretch;
	struct panel *p;
	int ifs, iss;

	c_stretch = conn_data[di].cstr;

	//TODO

	if ( fabs(c_stretch)<FLT_EPSILON ) c_stretch = stretch_coeff;

	p = connected->rigid_groups[di]->panels[ip];

	for ( ifs=p->min_fs; ifs<p->max_fs+1; ifs++ ) {
		for ( iss=p->min_ss; iss<p->max_ss+1; iss++ ) {

			int pix_index = ifs+edet->width*iss;

			recomputed_pixel_maps->pix_to_x[pix_index] /= c_stretch;
			recomputed_pixel_maps->pix_to_y[pix_index] /= c_stretch;
			if ( num_pix_displ[pix_index] >=
			     conn_data[di].num_peaks_per_pixel) {

				avg_displ->displ_x[pix_index] -=
				               (pixel_maps->pix_to_x[pix_index]-
				    recomputed_pixel_maps->pix_to_x[pix_index]);
				avg_displ->displ_y[pix_index] -=
				               (pixel_maps->pix_to_y[pix_index]-
				    recomputed_pixel_maps->pix_to_y[pix_index]);
			}
		}
	}
}


static void adjust_displ_for_stretch(struct rg_collection *connected,
                                     struct pixel_maps *pixel_maps,
                                     struct pixel_maps *recomputed_pixel_maps,
                                     struct connected_data *conn_data,
                                     double stretch_coeff,
                                     struct enhanced_det *edet,
                                     struct avg_displacements *avg_displ,
                                     int *num_pix_displ)
{

	int di, ip;

	for ( di=0; di<connected->n_rigid_groups; di++ ) {
		for (ip=0; ip<connected->rigid_groups[di]->n_panels; ip++) {

			adjust_panel(conn_data, di, ip, connected,
			             pixel_maps, recomputed_pixel_maps,
			             avg_displ, stretch_coeff, edet,
			             num_pix_displ);
		}
	}
}


static void fill_av_conn(struct rg_collection *connected, int di,
                         struct connected_data *conn_data,
                         int *num_pix_displ, struct enhanced_det *edet,
                         double *list_displ_in_conn_fs,
                         double *list_displ_in_conn_ss,
                         struct avg_displacements *avg_displ)
{
	struct panel *p;
	int ifs, iss, ip;

	int counter = 0;

	for (ip=0; ip<connected->rigid_groups[di]->n_panels; ip++) {

		p = connected->rigid_groups[di]->panels[ip];

		for ( ifs=p->min_fs; ifs<p->max_fs+1; ifs++ ) {
			for ( iss=p->min_ss; iss<p->max_ss+1; iss++ ) {

				int pix_index = ifs+edet->width*iss;

				if ( num_pix_displ[pix_index] >=
				     conn_data[di].num_peaks_per_pixel ) {
					list_displ_in_conn_fs[counter] =
					          avg_displ->displ_x[pix_index];
					list_displ_in_conn_ss[counter] =
					          avg_displ->displ_y[pix_index];
					counter++;
				}
			}
		}
	}

	if ( counter != conn_data[di].n_peaks_in_conn ) {
		printf("counter: %i n_peaks_in_conn: %i\n",
		       counter, conn_data[di].n_peaks_in_conn);
		exit(0);
	}
}


static void fill_conn_data_sh(struct connected_data *conn_data,
                             double *av_in_panel_fs,
                             double *av_in_panel_ss, int di,
                             double max_peak_distance)
{
	conn_data[di].sh_x = comp_median(av_in_panel_fs,
	                                 conn_data[di].n_peaks_in_conn);
	conn_data[di].sh_y = comp_median(av_in_panel_ss,
	                                 conn_data[di].n_peaks_in_conn);
	STATUS("Panel %s, num pixels: %i, shifts (in pixels) X,Y: %0.8f, %0.8f\n",
	       conn_data[di].name, conn_data[di].n_peaks_in_conn,
	       conn_data[di].sh_x, conn_data[di].sh_y);
	if ( modulus2d(conn_data[di].sh_x, conn_data[di].sh_y ) >
	     0.8*max_peak_distance ) {
		STATUS("  WARNING: absolute shift is: %0.1f > 0.8*%0.1f pixels."
		       "  Increase the value of the max_peak_distance parameter!\n",
		       modulus2d(conn_data[di].sh_x, conn_data[di].sh_y),
		       max_peak_distance);
	}
}


static int compute_shift(struct rg_collection *connected,
                         struct connected_data *conn_data,
                         int *num_pix_displ,
                         struct enhanced_det *edet,
                         struct geoptimiser_params *gparams,
                         struct avg_displacements *avg_displ)
{
        STATUS("Computing shift corrections.\n");

	int di;

	for ( di=0; di<connected->n_rigid_groups; di++ ) {

		double *list_displ_in_conn_fs;
		double *list_displ_in_conn_ss;

		list_displ_in_conn_fs = malloc(conn_data[di].n_peaks_in_conn*
		                        sizeof(double));
		if  ( list_displ_in_conn_fs == NULL ) {
			ERROR("Failed to allocate memory for computing shifts\n");
			return 1;
		}
		list_displ_in_conn_ss = malloc(conn_data[di].n_peaks_in_conn*
		                        sizeof(double));
		if  ( list_displ_in_conn_ss == NULL ) {
			ERROR("Failed to allocate memory for computing shifts\n");
			free(list_displ_in_conn_fs);
			return 1;
		}

		fill_av_conn(connected, di, conn_data, num_pix_displ, edet,
		             list_displ_in_conn_fs, list_displ_in_conn_ss,
		             avg_displ);

		if ( conn_data[di].n_peaks_in_conn >=
		     gparams->min_num_pix_per_conn_group ) {

			fill_conn_data_sh(conn_data, list_displ_in_conn_fs,
			                  list_displ_in_conn_ss, di,
			                  gparams->max_peak_dist);
		} else {
			conn_data[di].sh_x = -10000.0;
			conn_data[di].sh_y = -10000.0;
		}

		free(list_displ_in_conn_fs);
		free(list_displ_in_conn_ss);
	}
	return 0;
}


struct avg_shift *initialize_avg_shift(int num_rigid_groups)
{
	struct avg_shift *avg_sh;
	int i;

	avg_sh = malloc(sizeof(struct avg_shift));
	if ( avg_sh == NULL ) {
		ERROR("Failed to allocate memory to compute shifts for "
		      "empty panels.\n");
		return NULL;
	}

	avg_sh->aver_x = malloc(num_rigid_groups*sizeof(double));
	if ( avg_sh->aver_x == NULL ) {
		ERROR("Failed to allocate memory to compute shifts for "
		      "empty panels.\n");
		free(avg_sh);
		return NULL;
	}
	avg_sh->aver_y = malloc(num_rigid_groups*sizeof(double));
	if ( avg_sh->aver_y == NULL ) {
		ERROR("Failed to allocate memory to compute shifts for "
		      "empty panels.\n");
		free(avg_sh->aver_x);
		free(avg_sh);
		return NULL;
	}
	avg_sh->aver_num_sh = malloc(num_rigid_groups*sizeof(int));
	if ( avg_sh->aver_num_sh == NULL ) {
		ERROR("Failed to allocate memory to compute shifts for "
		      "empty panels.\n");
		free(avg_sh->aver_x);
		free(avg_sh->aver_y);
		free(avg_sh);
		return NULL;
	}

	for ( i=0; i<num_rigid_groups; i++ ) {
		avg_sh->aver_x[i] = 0.0;
		avg_sh->aver_y[i] = 0.0;
		avg_sh->aver_num_sh[i] = 0;
	}

	return avg_sh;
}


void free_avg_shift(struct avg_shift *av_sh) {
	free(av_sh->aver_x);
	free(av_sh->aver_y);
	free(av_sh->aver_num_sh);
	free(av_sh);
}


static int compute_shift_for_empty_panels(struct rg_collection *quadrants,
                                            struct rg_collection *connected,
                                            struct connected_data *conn_data,
                                            struct geoptimiser_params* gparams)
{

	struct avg_shift *av_sh;
	int di, i;

	STATUS("Computing shift corrections for groups without the required "
	       "number of measurements.\n");
	av_sh = initialize_avg_shift(quadrants->n_rigid_groups);
	if ( av_sh == NULL ) return 1;

	for ( di=0; di<connected->n_rigid_groups; di++ ) {
		if ( conn_data[di].n_peaks_in_conn >=
		     gparams->min_num_pix_per_conn_group ) {
			av_sh->aver_x[conn_data[di].num_quad] +=
			              conn_data[di].sh_x;
			av_sh->aver_y[conn_data[di].num_quad] +=
			              conn_data[di].sh_y;
			av_sh->aver_num_sh[conn_data[di].num_quad]++;
		}
	}

	for ( i=0; i<quadrants->n_rigid_groups; i++ ) {
		if (av_sh->aver_num_sh[i]>0) {
			av_sh->aver_x[i] /= (double)av_sh->aver_num_sh[i];
			av_sh->aver_y[i] /= (double)av_sh->aver_num_sh[i];
		}
	}

	for (di=0; di<connected->n_rigid_groups; di++) {
		if ( conn_data[di].n_peaks_in_conn <
		     gparams->min_num_pix_per_conn_group ) {
			if ( av_sh->aver_num_sh[conn_data[di].num_quad] > 0 ) {
				conn_data[di].sh_x =
				          av_sh->aver_x[conn_data[di].num_quad];
				conn_data[di].sh_y =
				          av_sh->aver_y[conn_data[di].num_quad];
				STATUS("Panel %s doesn't not have eough (%i) "
				       "peaks. Using average shifts (in pixels) "
				       "X,Y: %0.2f,%0.2f\n", conn_data[di].name,
				       conn_data[di].n_peaks_in_conn,
				       conn_data[di].sh_x, conn_data[di].sh_y);
			} else {
				STATUS("Panel %s has not enough (%i) peaks. "
				       "It will not be moved.\n",
				       conn_data[di].name,
				       conn_data[di].n_peaks_in_conn);
			}
		}
	}

	free_avg_shift(av_sh);

	return 0;
}


static void correct_shift(struct rg_collection *connected,
                          struct connected_data *conn_data,
                          double clen_to_use)
{

	int di;
	int ip;

	STATUS("Applying shift corrections.\n");

	for ( di=0;di<connected->n_rigid_groups;di++ ) {
		for (ip=0; ip<connected->rigid_groups[di]->n_panels; ip++) {

			struct panel *p;

			p = connected->rigid_groups[di]->panels[ip];

			if ( conn_data[di].sh_x > -9999.0 &&
			     conn_data[di].sh_y > -9999.0 ) {

				p->cnx -= conn_data[di].sh_x;
				p->cny -= conn_data[di].sh_y;

			} else {
				STATUS("For some reason shifts for panel %s has "
				       "not been computed!\n", p->name);
			}
		}
	}
}


static struct connected_stretch_and_angles *initialize_connected_stretch_angles(
						 struct rg_collection *connected)
{

	struct connected_stretch_and_angles *csaa;

	csaa = malloc(sizeof(struct connected_stretch_and_angles));
	if ( csaa == NULL ) {
		ERROR("Failed to allocate memory to compute angles and "
		      "stretch.\n");
		return NULL;
	}
	csaa->stretch_coeff = malloc(connected->n_rigid_groups*sizeof(double));
	if ( csaa->stretch_coeff == NULL ) {
		ERROR("Failed to allocate memory to compute angles and "
		      "stretch.\n");
		free(csaa);
		return NULL;
	}
	csaa->num_angles = malloc(connected->n_rigid_groups*sizeof(long));
	if ( csaa->num_angles == NULL ) {
		ERROR("Failed to allocate memory to compute angles and "
		      "stretch.\n");
		free(csaa->stretch_coeff);
		free(csaa);
		return NULL;
	}
	csaa->num_coeff=0;
	return csaa;
}


static void scan_p1(int ip0, int ip1, struct rg_collection *connected,
                    struct avg_displacements *avg_displ,
                    struct connected_data *conn_data,
                    struct enhanced_det *edet,
                    struct pixel_maps *pixel_maps,
                    int* num_pix_displ,
                    int di, double min_dist,
                    long *num_ang, int ifs0, int iss0,
                    double c_x0, double c_y0, double cd_x0, double cd_y0,
                    int compute, double *angles, double *stretches)
{

	int iss1, ifs1;

	struct panel *p1 = connected->rigid_groups[di]->panels[ip1];

	int min_ss_p1, min_fs_p1;

	if ( ip0 == ip1 ) {
		min_fs_p1 = ifs0;
		min_ss_p1 = iss0;
	} else {
		min_fs_p1 = p1->min_fs;
		min_ss_p1 = p1->min_ss;
	}

	for ( iss1=min_ss_p1; iss1<p1->max_ss+1; iss1++ ) {

		for ( ifs1=min_fs_p1; ifs1<p1->max_fs+1; ifs1++ ) {

			int pix_index1 = ifs1+edet->width*iss1;

			if ( num_pix_displ[pix_index1]>=
			     conn_data[di].num_peaks_per_pixel ) {

				double dist;
				double c_x1, c_y1, cd_x1, cd_y1;
				double d_c_x, d_c_y, d_cd_x, d_cd_y;
				double len1, len2;

				c_x1 = pixel_maps->pix_to_x[pix_index1];
				c_y1 = pixel_maps->pix_to_y[pix_index1];
				cd_x1 = c_x1 - avg_displ->displ_x[pix_index1];
				cd_y1 = c_y1 - avg_displ->displ_y[pix_index1];
				d_c_x = c_x1-c_x0;
				d_c_y = c_y1-c_y0;
				d_cd_x = cd_x1-cd_x0;
				d_cd_y = cd_y1-cd_y0;

				dist = modulus2d(d_c_x,d_c_y);
				if ( dist < min_dist ) continue;

				len1 = modulus2d(d_c_x, d_c_y);
				len2 = modulus2d(d_cd_x, d_cd_y);
				if ( len1<FLT_EPSILON || len2<FLT_EPSILON ) {
					continue;
				}

				if (compute) {

					double scal_m;
					double multlen;

					scal_m = d_c_x * d_cd_x+
					         d_c_y * d_cd_y -
						 FLT_EPSILON;

					multlen = len1*len2;
					if ( fabs(scal_m)>=multlen ) {
						angles[*num_ang] = 0.0;
					} else {

						angles[*num_ang] = acos(scal_m
						                   /multlen);

						if (d_c_y * d_cd_x -
						    d_c_x * d_cd_y < 0) {
							angles[*num_ang] *= -1.;
						}

					}

					stretches[*num_ang] = len1/len2;
				}

				*num_ang = *num_ang+1;
			}
		}
	}
}


static void scan_p0(int ip0, struct rg_collection *connected,
                    struct avg_displacements *avg_displ,
                    struct connected_data *conn_data,
                    struct enhanced_det *edet,
                    struct pixel_maps *pixel_maps,
                    int* num_pix_displ,
                    int di, double min_dist,
                    long *num_ang, int compute,
                    double *angles, double *stretches)
{

	int iss0, ifs0, ip1;

	struct panel *p0 = connected->rigid_groups[di]->panels[ip0];

	for ( iss0=p0->min_ss; iss0<p0->max_ss+1; iss0++ ) {

		for ( ifs0=p0->min_fs; ifs0<p0->max_fs+1; ifs0++ ) {

			int pix_index0 = ifs0+edet->width*iss0;

			if ( num_pix_displ[pix_index0]>=
			     conn_data[di].num_peaks_per_pixel ) {

				double c_x0, c_y0, cd_x0, cd_y0;

				c_x0 = pixel_maps->pix_to_x[pix_index0];
				c_y0 = pixel_maps->pix_to_y[pix_index0];
				cd_x0 = c_x0 - avg_displ->displ_x[pix_index0];
				cd_y0 = c_y0 - avg_displ->displ_y[pix_index0];

				for ( ip1=ip0;
				      ip1<connected->rigid_groups[di]->n_panels;
				      ip1++ ) {
					scan_p1(ip0, ip1, connected, avg_displ,
					        conn_data, edet, pixel_maps,
					        num_pix_displ, di, min_dist,
					        num_ang, ifs0, iss0, c_x0,
					        c_y0, cd_x0, cd_y0, compute,
					        angles, stretches);
				}
			}
		}
	}
}


static double compute_rotation_and_stretch(struct rg_collection *connected,
                                           struct connected_data *conn_data,
                                           struct enhanced_det *edet,
                                           int *num_pix_displ,
                                           struct pixel_maps *pixel_maps,
                                           struct avg_displacements *avg_displ,
                                           double dist_coeff_for_rot_str,
                                           struct geoptimiser_params *gparams)
{
	int di;
	double stretch_cf;

	struct connected_stretch_and_angles *csaa;

	STATUS("Computing rotation and stretch corrections.\n");

	csaa = initialize_connected_stretch_angles(connected);
	if ( csaa == NULL ) {
		return -1.0;
	}

	for ( di=0; di<connected->n_rigid_groups; di++ ) {
		if ( conn_data[di].n_peaks_in_conn <
		     gparams->min_num_pix_per_conn_group ) {
			continue;
		}

		long max_num_ang = 0;

		double min_dist;
		double* angles;
		double* stretches;

		struct panel *first_p;
		long num_ang = 0;
		int ip0;
		int num_pix_first_p;

		first_p = connected->rigid_groups[di]->panels[0];

		num_pix_first_p = (first_p->max_fs+1-first_p->min_fs)*
				  (first_p->max_ss+1-first_p->min_ss);

		// TODO: MINRAD HERE IS NOT UNIVERSAL
		min_dist = dist_coeff_for_rot_str*sqrt(num_pix_first_p*
		           connected->rigid_groups[di]->n_panels);

		for ( ip0=0; ip0<connected->rigid_groups[di]->n_panels; ip0++ ) {

			scan_p0(ip0, connected, avg_displ, conn_data, edet,
			        pixel_maps, num_pix_displ, di, min_dist,
			        &num_ang, 0, NULL, NULL);
		}

		max_num_ang = num_ang+1;

		angles = malloc(max_num_ang*sizeof(double));
		if ( angles == NULL ) {
			ERROR("Error in allocating memory for angle "
			      "optimization\n");
			free(csaa->stretch_coeff);
			free(csaa->num_angles);
			free(csaa);
			return -1.0;
		}
		stretches = malloc(max_num_ang*sizeof(double));
		if ( stretches == NULL ) {
			ERROR("Error in allocating memory for stretch "
			      "optimization\n");
			free(angles);
			free(csaa->stretch_coeff);
			free(csaa->num_angles);
			free(csaa);
			return -1.0;
		}

		num_ang = 0;

		for ( ip0=0; ip0<connected->rigid_groups[di]->n_panels; ip0++ ) {

			scan_p0(ip0, connected, avg_displ, conn_data, edet,
			        pixel_maps, num_pix_displ, di, min_dist,
			        &num_ang, 1, angles, stretches);
		}

		if ( num_ang < 1 ) continue;
		conn_data[di].cang = -comp_median(angles, num_ang);
		conn_data[di].cstr = comp_median(stretches, num_ang);

		STATUS("Panel %s, num: %li, angle: %0.4f deg, stretch coeff: "
		       "%0.4f\n", conn_data[di].name, num_ang, conn_data[di].cang,
		       conn_data[di].cstr);

		csaa->stretch_coeff[csaa->num_coeff] = conn_data[di].cstr;
		csaa->num_angles[csaa->num_coeff] = num_ang;
		csaa->num_coeff++;

		free(angles);
		free(stretches);
	}

	stretch_cf = 1.0;

	printf("Computing overall stretch coefficient.\n");

	if ( csaa->num_coeff>0 ) {

		int peaks_per_p;

		peaks_per_p = gparams->min_num_peaks_per_pix;

		while ( peaks_per_p>=0 ) {

			double total_num;
			long di;

			stretch_cf = 0;
			total_num = 0;
			for ( di=0; di<csaa->num_coeff; di++ ) {
				if ( conn_data[di].num_peaks_per_pixel >=
				     peaks_per_p ) {
					stretch_cf += csaa->stretch_coeff[di]*
					         (double)csaa->num_angles[di];
					total_num += csaa->num_angles[di];
				}
			}

			if ( total_num > 0 ) {

				stretch_cf /= total_num;

				printf("(Using only connected groups for which "
				       "the minimum number of measurements per "
				       "pixel is %i).\n", peaks_per_p);
				break;
			}
			peaks_per_p--;
		}
	}

	if ( stretch_cf < FLT_EPSILON ) {
		stretch_cf = 1.0;
	}

	STATUS("The global stretch coefficient for the patterns is %0.4f\n",
	       stretch_cf);
	if ( gparams->nostretch ) {
		stretch_cf = 1.0;
		STATUS("However, distance offset will not be optimized, so the "
		       "stretching coefficient has been set to 1.0\n");
		for ( di=0; di<connected->n_rigid_groups; di++ ) {
			conn_data[di].cstr = stretch_cf;
		}

	}

	free(csaa->stretch_coeff);
	free(csaa->num_angles);
	free(csaa);

	return stretch_cf;
}


#ifdef HAVE_SAVE_TO_PNG

static void draw_panel(struct image *image, cairo_t *cr,
                       cairo_matrix_t *basic_m, GdkPixbuf **pixbufs, int i) {
	struct panel p = image->det->panels[i];
	int w = gdk_pixbuf_get_width(pixbufs[i]);
	int h = gdk_pixbuf_get_height(pixbufs[i]);
	cairo_matrix_t m;

	/* Start with the basic coordinate system */
	cairo_set_matrix(cr, basic_m);

	/* Move to the right location */
	cairo_translate(cr, p.cnx, p.cny);

	/* Twiddle directions according to matrix */
	cairo_matrix_init(&m, p.fsx, p.fsy, p.ssx, p.ssy, 0.0, 0.0);
	cairo_transform(cr, &m);

	gdk_cairo_set_source_pixbuf(cr, pixbufs[i], 0.0, 0.0);
	cairo_rectangle(cr, 0.0, 0.0, w, h);
}


struct rectangle
{
	int width, height;
	double min_x, min_y, max_x, max_y;
};


static int draw_detector(cairo_surface_t *surf, struct image *image,
                         struct rectangle rect)
{
	cairo_t *cr;
	cairo_matrix_t basic_m;
	cairo_matrix_t m;
	GdkPixbuf **pixbufs;
	int n_pixbufs;

	cr = cairo_create(surf);

	pixbufs = render_panels(image, 1, SCALE_GEOPTIMISER, 1, &n_pixbufs);

	/* Blank grey background */
	cairo_rectangle(cr, 0.0, 0.0, rect.width, rect.height);
	cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
	cairo_fill(cr);

	/* Set up basic coordinate system
	 *  - origin in the centre, y upwards. */
	cairo_identity_matrix(cr);
	cairo_matrix_init(&m, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0);
	cairo_translate(cr, -rect.min_x , rect.max_y);
	cairo_transform(cr, &m);
	cairo_get_matrix(cr, &basic_m);

	if (pixbufs != NULL) {

		int i;

		for (i = 0; i < image->det->n_panels; i++) {
			draw_panel(image, cr, &basic_m, pixbufs, i);
			cairo_fill(cr);
		}

	}

	/* Free old pixbufs */
	if (pixbufs != NULL) {
		int i;
		for (i = 0; i < n_pixbufs; i++) {
			g_object_unref(pixbufs[i]);
		}
		free(pixbufs);
	}

	return 0;

}


static int save_data_to_png(char *filename, struct enhanced_det *edet,
                            double *data)
{
	struct image im;
	int i;
	struct rectangle rect;
	GdkPixbuf *col_scale;
	cairo_t *cr;

	cairo_status_t r;
	cairo_surface_t *surf;

	im.det = edet->det;
	im.width = edet->width;
	im.height = edet->height;
	im.bad = NULL;
	im.dp = malloc(edet->det->n_panels*sizeof(float *));
	if ( im.dp == NULL ) {
		ERROR("Failed to allocate data\n");
		return 1;
	}
	for ( i=0; i<edet->det->n_panels; i++ ) {

		int fs, ss;
		struct panel *p;

		p = &edet->det->panels[i];

		im.dp[i] = calloc(p->w * p->h, sizeof(float));
		if ( im.dp[i] == NULL ) {
			ERROR("Failed to allocate data\n");
			return 1;
		}

		for ( ss=0; ss<p->h; ss++ ) {
		for ( fs=0; fs<p->w; fs++ ) {

			int idx;
			int cfs, css;
			float val;

			cfs = fs+p->min_fs;
			css = ss+p->min_ss;
			idx = cfs + css*edet->width;

			if ( data[idx] == -10000.0) {
				val = 0.0;
			} else if ( data[idx] > 1.0) {
				val = 1.0;
			} else {
				val = (float)data[idx];
			}
			val *= 10.0; /* render_panels sets this as max */

			im.dp[i][fs+p->w*ss] = val;

		}
		}
	}

	get_pixel_extents(im.det, &rect.min_x, &rect.min_y, &rect.max_x,
	                  &rect.max_y);

	if (rect.min_x > 0.0) rect.min_x = 0.0;
	if (rect.max_x < 0.0) rect.max_x = 0.0;
	if (rect.min_y > 0.0) rect.min_y = 0.0;
	if (rect.max_y < 0.0) rect.max_y = 0.0;

	rect.width = rect.max_x - rect.min_x;
	rect.height = rect.max_y - rect.min_y;

	/* Add a thin border */
	rect.width += 2.0;
	rect.height += 2.0;
	surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, rect.width + 20,
	                                  rect.height);

	draw_detector(surf, &im, rect);

	col_scale = render_get_colour_scale(20, rect.height, SCALE_GEOPTIMISER);

	cr = cairo_create(surf);
	cairo_identity_matrix(cr);
	cairo_translate(cr, rect.width, 0.0);
	cairo_rectangle(cr, 0.0, 0.0, 20.0, rect.height);
	gdk_cairo_set_source_pixbuf(cr, col_scale, 0.0, 0.0);
	cairo_fill(cr);
	cairo_destroy(cr);

	for ( i=0; i<edet->det->n_panels; i++ ) {
		free(im.dp[i]);
	}
	free(im.dp);

	r = cairo_surface_write_to_png(surf, filename);
	if ( r != CAIRO_STATUS_SUCCESS ) return 1;

	return 0;
}

#endif /* HAVE_SAVE_TO_PNG */


int check_and_enforce_cspad_dist(struct geoptimiser_params *gparams,
				 struct detector *det)
{
	int np = 0;
	int num_errors_found = 0;

	double dist_to_check = 197.0;
	double tol = 0.2;

	for ( np=0; np<det->n_panels; np = np+2 ) {

		double dist2;
                double cnx2, cny2;

		struct panel *ep = &det->panels[np];
		struct panel *op = &det->panels[np+1];

                cnx2 = ep->cnx + 197.0*ep->fsx;
                cny2 = ep->cny + 197.0*ep->fsy;

		dist2 = (( cnx2 - op->cnx )*( cnx2 - op->cnx ) +
		         ( cny2 - op->cny )*( cny2 - op->cny ));

		if ( dist2 > (tol*tol)) {

			num_errors_found += 1;

			STATUS("Warning: distance between panels %s and %s "
			       "is outside acceptable margins (Corners are "
                               "more than 0.2 pixels away: %3.2f).\n", ep->name,
			       op->name, sqrt(dist2));

			if ( gparams->enforce_cspad_layout ) {

				double new_op_cx, new_op_cy;

				new_op_cx = ep->cnx + ep->fsx*dist_to_check;
				new_op_cy = ep->cny + ep->fsy*dist_to_check;

				op->cnx = new_op_cx;
				op->cny = new_op_cy;

				STATUS("Enforcing distance....\n");
			}

		}

		if ( ep->fsx != op->fsx || ep->ssx != op->ssx ||
		     ep->fsy != op->fsy || ep->ssx != op->ssx ) {

			num_errors_found += 1;

			STATUS("Warning: relative orientation between panels "
			       "%s and %s is incorrect.\n", ep->name, op->name);

			if ( gparams->enforce_cspad_layout ) {

				STATUS("Enforcing relative orientation....\n");

				op->fsx = ep->fsx;
				op->ssx = ep->ssx;
				op->fsy = ep->fsy;
				op->ssy = ep->ssy;

				op->xfs = ep->xfs;
				op->xss = ep->xss;
				op->yfs = ep->yfs;
				op->yss = ep->yss;
			}

		}

	}
        return num_errors_found;
}


static struct connected_data *initialize_conn_data(struct rg_collection *quadrants,
						   struct rg_collection *connected)
{

	struct connected_data *conn_data;

	conn_data = malloc(connected->n_rigid_groups*
	                   sizeof(struct connected_data));

	int di;

	for (di=0; di<connected->n_rigid_groups; di++) {

		conn_data[di].num_quad = find_quad_for_connected(
		                         connected->rigid_groups[di],
		                         quadrants);
		conn_data[di].cang = 0.0;
		conn_data[di].cstr = 1.0;
		conn_data[di].sh_x = -10000.0;
		conn_data[di].sh_y = -10000.0;
		conn_data[di].num_peaks_per_pixel = 1;
		conn_data[di].name = connected->rigid_groups[di]->name;
		conn_data[di].n_peaks_in_conn = 0;
	}

	return conn_data;
}


static struct pixel_displ_list *initialize_pixel_displacement_list(
				  struct enhanced_det *edet)
{

	struct pixel_displ_list *pix_displ_list;
	int ipx;

	pix_displ_list = malloc(sizeof(struct pixel_displ_list));

	pix_displ_list->pix_displ_list = calloc(edet->num_pix,
	                                 sizeof(struct single_pixel_displ));
	if ( pix_displ_list->pix_displ_list == NULL ) {
		ERROR("Error allocating memory for pixel displacement data.\n");
		free(pix_displ_list);
		return NULL;
	}
	pix_displ_list->curr_pix_displ = calloc(edet->num_pix,
	                                 sizeof(struct single_pixel_displ*));
	if ( pix_displ_list->curr_pix_displ == NULL ) {
		ERROR("Error allocating memory for pixel displacement data.\n");
		free(pix_displ_list->pix_displ_list);
		free(pix_displ_list);
		return NULL;
	}
	pix_displ_list->num_pix_displ = calloc(edet->num_pix, sizeof(int));
	if ( pix_displ_list->num_pix_displ == NULL ) {
		ERROR("Error allocating memory for pixel displacement data.\n");
		free(pix_displ_list->curr_pix_displ);
		free(pix_displ_list->pix_displ_list);
		free(pix_displ_list);
		return NULL;
	}

	for ( ipx=0; ipx<edet->num_pix; ipx++ ) {
		pix_displ_list->pix_displ_list[ipx].dfs = -10000.0;
		pix_displ_list->pix_displ_list[ipx].dss = -10000.0;
		pix_displ_list->pix_displ_list[ipx].ne = NULL;
		pix_displ_list->curr_pix_displ[ipx] = &pix_displ_list->pix_displ_list[ipx];
		pix_displ_list->num_pix_displ[ipx] = 0;
	}

	return pix_displ_list;
}


static void free_pixel_displacement_list(
                struct pixel_displ_list *pix_displ_list,
                struct enhanced_det *edet)
{
	int i;
	struct single_pixel_displ *curr = NULL;
	struct single_pixel_displ *next = NULL;

	for ( i=0; i<edet->num_pix; i++ ) {

		curr = &pix_displ_list->pix_displ_list[i];

		if ( curr->ne != NULL ) {
			curr = curr->ne;
			while ( curr != NULL ) {
				next = curr->ne;
				free(curr);
				curr = next;
			}
		}
	}

	free(pix_displ_list->curr_pix_displ);
	free(pix_displ_list->pix_displ_list);
	free(pix_displ_list->num_pix_displ);
	free(pix_displ_list);
}


static struct avg_displacements *initialize_average_displacement(
				 struct enhanced_det *edet)

{
	static struct avg_displacements *avg_displ;

	avg_displ = malloc(sizeof(struct avg_displacements));
	avg_displ->displ_x = calloc(edet->num_pix, sizeof(double));
	if ( avg_displ->displ_x == NULL ) {
		ERROR("Error allocating memory for pixel properties.\n");
		return NULL;
	}
	avg_displ->displ_y = calloc(edet->num_pix, sizeof(double));
	if ( avg_displ->displ_y == NULL ) {
		ERROR("Error allocating memory for pixel properties.\n");
		free(avg_displ->displ_x);
		free(avg_displ);
		return NULL;
	}
	avg_displ->displ_abs = calloc(edet->num_pix, sizeof(double));
	if ( avg_displ->displ_abs == NULL ) {
		ERROR("Error allocating memory for pixel properties.\n");
		free(avg_displ->displ_x);
		free(avg_displ->displ_y);
		free(avg_displ);
		return NULL;
	}

	return avg_displ;
}


static void free_avg_displacement(struct avg_displacements * avg_displ)
{
	free(avg_displ->displ_x);
	free(avg_displ->displ_y);
	free(avg_displ->displ_abs);
	free(avg_displ);
}


static void free_pattern_list(struct pattern_list *pattern_list)
{

	int pti;

	for ( pti=0; pti<pattern_list->n_patterns; pti++ ) {
		int nuc;

		image_feature_list_free(pattern_list->patterns[pti]->im_list);
		reflist_free(pattern_list->patterns[pti]->ref_list);
		for ( nuc=0; nuc<pattern_list->patterns[pti]->n_unit_cells;
		      nuc++) {
			cell_free(pattern_list->patterns[pti]->unit_cells[nuc]);
		}
		free(pattern_list->patterns[pti]->filename);
		free(pattern_list->patterns[pti]);
	}
	free(pattern_list);
}


static int *extract_num_pix_free_displ_list(struct pixel_displ_list *pix_displ_list,
                                             struct enhanced_det *edet)
{
	int *num_pix;

	int i;
	struct single_pixel_displ *curr = NULL;
	struct single_pixel_displ *next = NULL;

	for ( i=0; i<edet->num_pix; i++ ) {

		curr = &pix_displ_list->pix_displ_list[i];

		if ( curr->ne != NULL ) {
			curr = curr->ne;
			while ( curr != NULL ) {
				next = curr->ne;
				free(curr);
				curr = next;
			}
		}
	}

	num_pix = pix_displ_list->num_pix_displ;

	free(pix_displ_list->curr_pix_displ);
	free(pix_displ_list->pix_displ_list);
	free(pix_displ_list);

	return num_pix;
}


static void recompute_panel_avg_displ(struct rg_collection *connected,
                                      struct connected_data *conn_data,
                                      int *num_pix_displ, int di, int ip,
                                      struct enhanced_det *edet,
                                      struct avg_displacements *avg_displ)
{
	struct panel *p;
	int ifs, iss;

	if (conn_data[di].sh_x < -9999.0) return;

	p = connected->rigid_groups[di]->panels[ip];

	for ( ifs=p->min_fs; ifs<p->max_fs+1; ifs++ ) {
		for ( iss=p->min_ss; iss<p->max_ss+1; iss++ ) {

			int pix_index = ifs+edet->width*iss;

			if ( num_pix_displ[pix_index]>=
			     conn_data[di].num_peaks_per_pixel ) {
				avg_displ->displ_x[pix_index] -= conn_data[di].sh_x;
				avg_displ->displ_y[pix_index] -= conn_data[di].sh_y;
				avg_displ->displ_abs[pix_index] = modulus2d(
				                      avg_displ->displ_x[pix_index],
				                      avg_displ->displ_y[pix_index]
				                                  );
			} else {
				avg_displ->displ_abs[pix_index] = -10000.0;
			}
		}
	}
}


void recompute_avg_displ(struct rg_collection *connected,
                         struct connected_data *conn_data,
                         int *num_pix_displ,
                         struct enhanced_det *edet,
                         struct avg_displacements *avg_displ)
{

	int di, ip;

	for ( di=0;di<connected->n_rigid_groups;di++ ) {
		for (ip=0; ip<connected->rigid_groups[di]->n_panels; ip++) {

			recompute_panel_avg_displ(connected, conn_data,
			                          num_pix_displ, di, ip, edet,
			                          avg_displ);

		}
	}
}


int optimize_geometry(struct geoptimiser_params *gparams,
                      struct detector *det,
                      struct rg_collection *quadrants,
                      struct rg_collection *connected)
{
	int max_fs = 0;
	int max_ss = 0;
	int pi;
	int ret;
	int write_ret;
	int maybe_cspad = 0;
	int *num_pix_displ;

	double res_sum;
	double avg_res;
	double clen_to_use;

	// for angles and stretch calculation use
	// only pixels which are distco*size_panel
	// away
	double dist_coeff_for_rot_str = 0.2;
	double total_error;

	struct pixel_maps *pixel_maps;
	struct pixel_maps *recomputed_pixel_maps;
	double stretch_coeff = 1.0;

	struct connected_data *conn_data = NULL;
	struct pattern_list *pattern_list;
	struct cell_params *cell_params;
	struct enhanced_det edet;
	struct avg_displacements *avg_displ;
	struct pixel_displ_list *pix_displ_list;

	STATUS("Maximum distance between peaks: %0.1f pixels.\n",
	       gparams->max_peak_dist);

	STATUS("Minimum number of measurements for a pixel to be included in the "
	       "refinement: %i\n",
	       gparams->min_num_peaks_per_pix);
	STATUS("Minimum number of measurements for connected group for accurate "
	       "estimation of position/orientation: %i\n",
	       gparams->min_num_pix_per_conn_group);

	if (det->n_panels == 64 ) {
		maybe_cspad = 1;
	}

	if ( maybe_cspad && !gparams->no_cspad ) {

		int num_errors = 0;

		STATUS("It looks like the detector is a CSPAD. "
		       "Checking relative distance and orientation of "
		       "connected ASICS.\n");
		STATUS("If the detector is not a CSPAD, please rerun the "
		       "program with the --no-cspad option.\n");
		if ( gparams->enforce_cspad_layout ) {
			STATUS("Enforcing CSPAD layout...\n");
		}

		num_errors = check_and_enforce_cspad_dist(gparams, det);

		if ( gparams->enforce_cspad_layout ) {

			int geom_wr;

			STATUS("Saving geometry with enforced CSPAD layout.\n"
			       "Please restart geometry optimization using the "
			       "optimized geometry from this run as input geometry "
			       "file.\n");
			geom_wr = write_detector_geometry_2(
			                        gparams->geometry_filename,
			                        gparams->outfile, det,
			                        gparams->command_line, 1);
			if ( geom_wr != 0 ) {
				ERROR("Error in writing output geometry file.\n");
				return 1;
			}
			STATUS("All done!\n");
			return 0;
		}

		if ( !gparams->enforce_cspad_layout && num_errors > 0 ) {
			ERROR("Relative distance and orientation of connected "
			      "ASICS do not respect the CSPAD layout.\n"
			      "Geometry optimization cannot continue.\n"
			      "Please rerun the program with the "
			      "--enforce-cspad-layout option.\n");
			return 1;
		}
	}

	pattern_list = read_patterns_from_steam_file(gparams->infile, det);
	if ( pattern_list->n_patterns < 1 ) {
		ERROR("Error reading stream file\n");
		return 1;
	}

	cell_params = compute_avg_cell_parameters(pattern_list);
	if ( cell_params == NULL ) {
		free(pattern_list);
		return 1;
	}

	res_sum = 0;
	for ( pi=0; pi<det->n_panels; pi++ ) {

		if ( det->panels[pi].max_fs > max_fs ) {
			max_fs = det->panels[pi].max_fs;
		}
		if ( det->panels[pi].max_ss > max_ss ) {
			max_ss = det->panels[pi].max_ss;
		}
		res_sum += det->panels[pi].res;
	}
	avg_res = res_sum/det->n_panels;

	edet.det = det;
	edet.width = max_fs+1;
	edet.height = max_ss+1;
	edet.num_pix = (max_fs+1)*(max_ss+1);

	clen_to_use = pick_clen_to_use(gparams, pattern_list, avg_res,
	                               cell_params);

	if ( clen_to_use == -1.0 ) return 1;

	avg_displ = initialize_average_displacement(&edet);
	if ( avg_displ == NULL ) {
		free(cell_params);
		free(pattern_list);
		return 1;
	}

	pixel_maps = initialize_pixel_maps(&edet);
	if ( pixel_maps == NULL ) {
		ERROR("Failed to allocate memory for pixel maps.\n");
		free_avg_displacement(avg_displ);
		free(cell_params);
		free(pattern_list);
		return 1;
	}

	pix_displ_list = initialize_pixel_displacement_list(&edet);
	if ( pix_displ_list == NULL ) {
		ERROR("Error allocating memory for connected structure data.\n");
		free_avg_displacement(avg_displ);
		free_pixel_maps(pixel_maps);
		free(pattern_list);
		return 1;
	}

	conn_data = initialize_conn_data(quadrants, connected);
	if ( conn_data == NULL ) {
		ERROR("Error allocating memory for connected structure data.\n");
		free_avg_displacement(avg_displ);
		free_pixel_maps(pixel_maps);
		free_pixel_displacement_list(pix_displ_list, &edet);
		free(pattern_list);
		return 1;
	}


	ret = compute_pixel_displacements(pattern_list, &edet, connected,
	                                  gparams, clen_to_use, conn_data,
	                                  avg_displ, pix_displ_list);
	if ( ret != 0 ) {
		free_avg_displacement(avg_displ);
		free_pixel_maps(pixel_maps);
		free_pixel_displacement_list(pix_displ_list, &edet);
		free(conn_data);
		free(pattern_list);
		return 1;
	}

	free_pattern_list(pattern_list);

	adjust_min_peaks_per_conn(&edet, connected, gparams, conn_data,
	                          pix_displ_list);

	ret = compute_avg_displacements(&edet, connected, pix_displ_list,
	                                conn_data, avg_displ);
	if ( ret != 0 ) {
		free_avg_displacement(avg_displ);
		free_pixel_maps(pixel_maps);
		free_pixel_displacement_list(pix_displ_list, &edet);
		free(conn_data);
		return 1;
	}

	num_pix_displ = extract_num_pix_free_displ_list(pix_displ_list, &edet);

	STATUS("Computing error before correction.\n");
	total_error = compute_error(connected, &edet, conn_data,
	                            num_pix_displ, avg_displ->displ_abs);

	STATUS("Detector-wide error before correction <delta^2> = %0.4f pixels.\n",
	       total_error);

	if ( gparams->error_maps ) {
		STATUS("Saving error map before correction.\n");

#ifdef HAVE_SAVE_TO_PNG

		ret = save_data_to_png("error_map_before.png", &edet,
		                        avg_displ->displ_abs);
		if ( ret != 0 ) {
			ERROR("Error while writing data to file.\n");

			free_avg_displacement(avg_displ);
			free_pixel_maps(pixel_maps);
			free(conn_data);
			free(num_pix_displ);
			return 1;
		}

#else /* HAVE_SAVE_TO_PNG */

		ERROR("WARNING: geoptimiser was compiled without GTK and cairo "
		       "support. Error maps will not be saved.\n");

#endif /* HAVE_SAVE_TO_PNG */

	}


	stretch_coeff = compute_rotation_and_stretch(connected, conn_data,
	                                             &edet, num_pix_displ,
	                                             pixel_maps, avg_displ,
                                                     dist_coeff_for_rot_str,
                                                     gparams);
	if ( stretch_coeff < 0.0 ) {
		free_avg_displacement(avg_displ);
		free_pixel_maps(pixel_maps);
		free(conn_data);
		free(num_pix_displ);
		return 1;
	}

	ret = compute_rot_stretch_for_empty_panels(quadrants, connected,
	                                           gparams, conn_data);
	if ( ret != 0 ) {
		free_avg_displacement(avg_displ);
		free_pixel_maps(pixel_maps);
		free(conn_data);
		free(num_pix_displ);
		return 1;
	}


	correct_rotation_and_stretch(connected, edet.det, conn_data,
	                             clen_to_use, stretch_coeff,
	                             gparams);

	recomputed_pixel_maps = initialize_pixel_maps(&edet);
	if ( recomputed_pixel_maps == NULL ) {
		free_avg_displacement(avg_displ);
		free_pixel_maps(pixel_maps);
		free(conn_data);
		free(num_pix_displ);
		return 1;
	}

	adjust_displ_for_stretch(connected, pixel_maps,
	                         recomputed_pixel_maps,
	                         conn_data,
	                         stretch_coeff, &edet, avg_displ,
	                         num_pix_displ);

	ret = compute_shift(connected, conn_data, num_pix_displ, &edet,
	                    gparams, avg_displ);
	if ( ret != 0 ) {
		free_avg_displacement(avg_displ);
		free_pixel_maps(pixel_maps);
		free(conn_data);
		free(num_pix_displ);
		free_pixel_maps(recomputed_pixel_maps);
		return 1;
	}

	compute_shift_for_empty_panels(quadrants, connected, conn_data,
	                               gparams);

	correct_shift(connected, conn_data, clen_to_use);

	recompute_avg_displ(connected, conn_data,
	                    num_pix_displ, &edet,
	                    avg_displ);

	if ( gparams->error_maps ) {

#ifdef HAVE_SAVE_TO_PNG

		STATUS("Saving error map after correction.\n");

		ret = save_data_to_png("error_map_after.png", &edet,
		                       avg_displ->displ_abs);
		if ( ret !=0 ) {
			ERROR("Error while writing data to file.\n");

			free_avg_displacement(avg_displ);
			free_pixel_maps(pixel_maps);
			free(conn_data);
			free(num_pix_displ);
			free_pixel_maps(recomputed_pixel_maps);
			return 1;
		}

#else /* HAVE_SAVE_TO_PNG */

		STATUS("ERROR: geoptimiser was compiled without GTK and cairo "
		       "support.\n Error maps will not be saved.\n");

#endif /* HAVE_SAVE_TO_PNG */

	}

	STATUS("Computing errors after correction.\n");
	total_error = compute_error(connected, &edet, conn_data,
	              num_pix_displ, avg_displ->displ_abs);

	STATUS("Detector-wide error after correction <delta^2> = %0.4f pixels.\n",
	       total_error);

	write_ret = write_detector_geometry_2(gparams->geometry_filename,
	                                      gparams->outfile, edet.det,
	                                      gparams->command_line, 1);
	if ( write_ret != 0 ) {
		ERROR("Error in writing output geometry file.\n");
		return 1;
	}
	STATUS("All done!\n");
	if ( gparams->error_maps ) {

#ifdef HAVE_SAVE_TO_PNG

		STATUS("Be sure to inspect error_map_before.png and "
		       "error_map_after.png !!\n");

#endif /* HAVE_SAVE_TO_PNG */
	}

	free_avg_displacement(avg_displ);
	free_pixel_maps(pixel_maps);
	free(conn_data);
	free(num_pix_displ);
	free_pixel_maps(recomputed_pixel_maps);
	return 0;
}


int main(int argc, char *argv[])
{
	int c, i;
	int ret_val;
	char buffer[256];
	char command_line[1024];

	char *quadrant_coll_name = NULL;
	char *connected_coll_name = NULL;

	struct geoptimiser_params *gparams;
	struct detector *det = NULL;
	struct rg_collection *quadrants;
	struct rg_collection *connected;
	struct beam_params beam;

	gparams = malloc(sizeof(struct geoptimiser_params));

	gparams->outfile = NULL;
	gparams->infile = NULL;
	gparams->geometry_filename = NULL;
	gparams->min_num_peaks_per_pix = 3;
	gparams->min_num_pix_per_conn_group = 100;
	gparams->only_best_distance = 0;
	gparams->enforce_cspad_layout = 0;
	gparams->nostretch = 0;
	gparams->individual_coffset = 0;
	gparams->no_cspad = 0;
	gparams->error_maps = 1;
	gparams->max_peak_dist = 4.0;

	const struct option longopts[] = {

	/* Options with long and short versions */
	{"help",                             0, NULL,               'h'},
	{"version",                          0, NULL,               10 },
	{"input",                            1, NULL,               'i'},
	{"output",                           1, NULL,               'o'},
	{"geometry",                         1, NULL,               'g'},
	{"quadrants",                        1, NULL,               'q'},
	{"connected",                        1, NULL,               'c'},
	{"min-num-peaks-per-pixel",          1, NULL,               'x'},
	{"min-num-pixels-per-conn-group",    1, NULL,               'p'},
	{"most-few-clen",                    0, NULL,               'l'},
	{"max-peak-dist",                    1, NULL,               'm'},
	{"individual-dist-offset",           0, NULL,               's'},

	/* Long-only options with no arguments */
	{"no-stretch",                       0, &gparams->nostretch,            1},
	{"no-error-maps",                    0, &gparams->error_maps,           0},
	{"enforce-cspad-layout",             0, &gparams->enforce_cspad_layout, 1},
	{"no-cspad",                         0, &gparams->no_cspad,             1},

	/* Long-only options with arguments */
	{"min-num-peaks-per-panel",          1, NULL,                          11},


	{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ho:i:g:q:c:o:x:p:lsm:",
	                       longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 10 :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			case 'o' :
			gparams->outfile = strdup(optarg);
			break;

			case 'i' :
			gparams->infile = strdup(optarg);
			break;

			case 'g' :
			gparams->geometry_filename = strdup(optarg);
			det = get_detector_geometry(gparams->geometry_filename,
			                            &beam);
			if ( det == NULL ) {
				ERROR("Failed to read detector geometry from "
				      "'%s'\n", optarg);
				return 1;
			}
			break;

			case 'q' :
			quadrant_coll_name = strdup(optarg);
			break;

			case 'c' :
			connected_coll_name = strdup(optarg);
			break;

			case 'x' :
			gparams->min_num_peaks_per_pix = atoi(optarg);
			break;

			case 'p' :
			gparams->min_num_pix_per_conn_group = atoi(optarg);
			break;

			case 'l' :
			gparams->only_best_distance = 1;
			break;

			case 'm' :
			gparams->max_peak_dist = strtof(optarg, NULL);
			break;

			case 's' :
			gparams->individual_coffset = 1;
			break;

			case 11:
			ERROR("WARNING: The --min-num-peaks-per-panel option has been "
			      "renamed to --min-num-pixels-per-conn-group. The "
			      "old option has been deprecated and will "
			      "soon be removed. It is currently remapped to "
			      "--min-num-pixels-per-conn-group\n");
			gparams->min_num_pix_per_conn_group = atoi(optarg);
			break;

		}
	}

	if ( gparams->geometry_filename == NULL ) {
		ERROR("You must provide a geometry to optimize.\n");
		return 1;
	}

	if ( gparams->infile == NULL ) {
		ERROR("You must provide an input stream file.\n");
		return 1;
	}

	if ( gparams->outfile == NULL ) {
		ERROR("You must provide an output filename.\n");
		return 1;
	}

	if ( quadrant_coll_name == NULL ) {
		ERROR("You must provide a rigid group collection for "
		      "quadrants.\n");
		return 1;
	}

	if ( connected_coll_name == NULL ) {
		ERROR("You must provide a rigid group collection for connected "
		      "panels.\n");
		return 1;
	}

	strcpy(command_line, "\0");

	quadrants = find_rigid_group_collection_by_name(det, quadrant_coll_name);
	if ( quadrants == NULL ) {
		ERROR("Cannot find rigid group collection for quadrants: %s\n",
		      quadrant_coll_name);
		return 1;
	}

	connected = find_rigid_group_collection_by_name(det,
	                                                connected_coll_name);
	if ( connected == NULL ) {
		ERROR("Cannot find rigid group collection for connected "
		      "asics: %s\n", connected_coll_name);
		return 1;
	}

	for ( i=0; i<argc; i++ ) {
		if ( i > 0 ) strcat(command_line, " ");
		strcpy(buffer, argv[i]);
		strcat(command_line, buffer);
	}

	g_type_init();
	ret_val = optimize_geometry(gparams, det, quadrants, connected);

	return ret_val;
}
