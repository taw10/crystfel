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
#include <cairo.h>
#include <gdk/gdk.h>

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
"  -h, --help                  Display this help message.\n"
"\n"
"      --version                        Print CrystFEL version number and\n"
"                                        exit.\n"
"  -i, --input=<filename>               Specify stream file to be used for \n"
"                                        geometry optimization.\n"
"  -g. --geometry=<file>                Get detector geometry from file.\n"
"  -o, --output=<filename>              Output stream.\n"
"  -q, --quadrants=<rg_coll>            Rigid group collection for quadrants.\n"
"  -c, --connected=<rg_coll>            Rigid group collection for connected\n"
"                                        ASICs.\n"
"      --no-error-maps                  Do not generate error map PNGs.\n"
"  -x, --min-num-peaks-per-pixel=<num>  Minimum number of peaks per pixel.\n"
"                                        Default: 3. \n"
"  -p, --min-num-peaks-per-panel=<num>  Minimum number of peaks per pixel.\n"
"                                        Default: 100.\n"
"  -l, --most-freq-clen                 Use only the most frequent camera\n"
"                                        length.\n"
"  -s, --individual-dist-offset         Use a distance offset for each panel.\n"
"                                        Default: whole-detector offset.\n"
"      --no-stretch                     Do not optimize distance offset.\n"
"                                        Default: distance offset is optimized\n"
"  -m  --max-peak-dist=<num>            Maximum distance between predicted and\n"
"                                        detected peaks (in pixels)\n"
"                                        Default: 4.0 pixels.\n"
);
}


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


struct single_pix_displ
{
	double dfs;
	double dss;
	struct single_pix_displ* ne;
};


struct connected_stretch_and_angles
{
	double *stretch_coeff;
	unsigned int *num_angles;
	int num_coeff;
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
              refl = next_refl(refl, iter) )
	{
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
			return atoi(stuff->fields[i]+
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


static void compute_avg_cell_parameters(struct pattern_list *pattern_list,
                                        double *avcp)
{
	int numavc;
	int j, i;
	double minc[6];
	double maxc[6];

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
				avcp[j] += cpar[j];
				if (cpar[j]<minc[j]) minc[j] = cpar[j];
				if (cpar[j]>maxc[j]) maxc[j] = cpar[j];
			}
			numavc++;

		}

	}

	if ( numavc>0 ) {
		for ( j=0; j<6; j++ ) avcp[j] /= numavc;
	}

	STATUS("Average cell coordinates:\n");
	STATUS("Average a, b, c (in nm): %6.3f, %6.3f, %6.3f\n",
	       avcp[0],avcp[1],avcp[2]);
	STATUS("Minimum -Maximum a, b, c:\n"
	       "\t%6.3f - %6.3f,\n"
	       "\t%6.3f - %6.3f,\n"
	       "\t%6.3f - %6.3f\n",
	       minc[0], maxc[0], minc[1], maxc[1], minc[2], maxc[2]);
	STATUS("Average alpha,beta,gamma in degrees: %6.3f, %6.3f, %6.3f\n",
	       avcp[3], avcp[4], avcp[5]);
	STATUS("Minimum - Maximum alpha,beta,gamma in degrees:\n"
	       "\t%5.2f - %5.2f,\n"
	       "\t%5.2f - %5.2f,\n"
	       "\t%5.2f - %5.2f\n",
	       minc[3], maxc[3], minc[4], maxc[4], minc[5], maxc[5]);

}


static double compute_clen_to_use(struct pattern_list *pattern_list,
                                    double istep, double *avcp,
                                    double max_peak_distance,
                                    int only_best_distance)
{
	int cp, i, u;
	int num_clens;
	int max_clens;
	int best_clen;
	int *clens_population;
	double *clens;
	double *lambdas;
	double irecistep;
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
			cqu = get_q_from_xyz(1/istep, 0, clens[i], lambdas[i]);

			irecistep = 1/cqu.u;

			min_braggp_dist = fmin(fmin(irecistep/avcp[0],
			                       irecistep/avcp[1]),
			                       irecistep/avcp[2]);
			STATUS("Camera length %0.4f m was found %i times.\n"
			       "Minimum inter-bragg peak distance (based on "
			       "average cell parameters): %0.1f pixels.\n",
			       clens[i], clens_population[i],
			       min_braggp_dist);
			if ( min_braggp_dist<1.2*max_peak_distance ) {
				STATUS("WARNING: The distance between Bragg "
				       "peaks is too small: "
				       "%0.1f < 1.2*%0.1f pixels.\n",
			               min_braggp_dist,
			               max_peak_distance);
			}
			if ( clens_population[i] > clens_population[best_clen] ) {
				best_clen = i;
				clen_to_use = clens[best_clen];
			}
		}
	}

	if ( only_best_distance ) {
		STATUS("Only %i patterns with  CLEN=%0.4f m will be used.\n",
		       clens_population[best_clen], clen_to_use);
	}

	free(clens);
	free(lambdas);
	free(clens_population);

	return clen_to_use;
}


static double comp_median(double *arr, unsigned int n)
{

	int low, high, median, middle, ll, hh;
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


static void free_all_curr_pix_displ(struct single_pix_displ *all_pix_displ,
                                    struct single_pix_displ **curr_pix_displ,
                                    int num_pix_in_slab)
{
    int i;
    struct single_pix_displ *curr = NULL;
    struct single_pix_displ *next = NULL;

    for ( i=0; i<num_pix_in_slab; i++ ) {

        curr = &all_pix_displ[i];


        if ( curr->ne != NULL ) {
            curr = curr->ne;
            while ( curr != NULL ) {
                next = curr->ne;
                free(curr);
                curr = next;
            }
        }
    }

    free(curr_pix_displ);
    free(all_pix_displ);
}


static int fill_pixel_statistics(int *num_pix_displ,
                                 struct single_pix_displ** curr_pix_displ,
				 struct single_pix_displ* all_pix_displ,
				 struct connected_data *conn_data,
				 int ifs, int iss, int di, int aw, int dfv,
				 double *displ_x,
                                 double *displ_y, double *displ_abs)
{
	double *cPxAfs;
	double *cPxAss;
	int cnu = 0;

	cPxAfs = calloc(num_pix_displ[ifs+aw*iss], sizeof(double));
	if ( cPxAfs == NULL ) {
		ERROR("Failed to allocate memory for pixel statistics.\n");
		return 1;
	}
	cPxAss = calloc(num_pix_displ[ifs+aw*iss], sizeof(double));
	if ( cPxAss == NULL ) {
		ERROR("Failed to allocate memory for pixel statistics.\n");
		free(cPxAfs);
		return 1;
	}

	curr_pix_displ[ifs+aw*iss] = &all_pix_displ[ifs+aw*iss];

	while (1) {
		if (curr_pix_displ[ifs+aw*iss]->dfs == dfv) break;
		cPxAfs[cnu] = curr_pix_displ[ifs+aw*iss]->dfs;
		cPxAss[cnu] = curr_pix_displ[ifs+aw*iss]->dss;
		cnu++;
		if ( curr_pix_displ[ifs+aw*iss]->ne == NULL ) break;
		curr_pix_displ[ifs+aw*iss] = curr_pix_displ[ifs+aw*iss]->ne;
	}

	if ( cnu < 1 ) return 0;

	displ_x[ifs+aw*iss] = comp_median(cPxAfs, cnu);
	displ_y[ifs+aw*iss] = comp_median(cPxAss, cnu);
	displ_abs[ifs+aw*iss] = modulus2d(displ_x[ifs+aw*iss],
	                                  displ_y[ifs+aw*iss]);
	conn_data[di].n_peaks_in_conn++;

	free(cPxAfs);
	free(cPxAss);

	return 0;
}



static int compute_panel_statistics(struct rg_collection *connected,
                                    int *num_pix_displ,
				    struct single_pix_displ** curr_pix_displ,
				    struct single_pix_displ* all_pix_displ,
				    struct connected_data *conn_data,
				    int di, int ip, int np,
				    double dfv,
				    int aw, double *displ_x,
                                    double *displ_y, double *displ_abs)
{
	struct panel *p;
	int ifs, iss;

	p = connected->rigid_groups[di]->panels[ip];

	for ( ifs=p->min_fs; ifs<p->max_fs+1; ifs++ ) {
		for ( iss=p->min_ss; iss<p->max_ss+1; iss++ ) {
			if ( num_pix_displ[ifs+aw*iss]>=np ) {

				int ret;

				ret = fill_pixel_statistics(num_pix_displ,
				                            curr_pix_displ,
				                            all_pix_displ,
				                            conn_data,
				                            ifs, iss, di, aw,
				                            dfv, displ_x,
				                            displ_y, displ_abs);

				if ( ret == -2 ) {
					break;
				} else if ( ret != 0 ) {
					return ret;
				}

			} else {
				displ_x[ifs+aw*iss] = dfv;
				displ_y[ifs+aw*iss] = dfv;
				displ_abs[ifs+aw*iss] = dfv;
			}

		}
	}
	return 0;
}



static int allocate_next_element(struct single_pix_displ** curr_pix_displ,
                                 int ipx)
{
	curr_pix_displ[ipx]->ne = malloc(sizeof(struct single_pix_displ));
	if ( curr_pix_displ[ipx]->ne == NULL ) {
		ERROR("Failed to allocate memory for pixel statistics.\n");
		return 1;
	}

	curr_pix_displ[ipx] = curr_pix_displ[ipx]->ne;

	return 0;
}


static int compute_pixel_statistics(struct pattern_list *pattern_list,
                                    struct detector *det,
                                    struct rg_collection *connected,
                                    struct rg_collection *quadrants,
                                    int num_pix_in_slab,
                                    int max_peak_distance, int aw,
                                    double dfv,
                                    double min_num_peaks_per_pixel,
                                    double min_num_peaks_per_panel,
                                    int only_best_distance,
                                    double clen_to_use,
                                    double *slab_to_x, double *slab_to_y,
                                    struct connected_data *conn_data,
                                    double *displ_x,
                                    double *displ_y, double *displ_abs,
                                    struct single_pix_displ* all_pix_displ,
                                    struct single_pix_displ** curr_pix_displ,
                                    int *num_pix_displ)
{
	int ipx, cp, ich, di, ip, np;

	for (di=0; di<connected->n_rigid_groups; di++) {

		conn_data[di].num_quad = find_quad_for_connected(
		                                    connected->rigid_groups[di],
		                                    quadrants);
		conn_data[di].cang = 0.0;
		conn_data[di].cstr = 1.0;
		conn_data[di].sh_x = dfv;
		conn_data[di].sh_y = dfv;
		conn_data[di].num_peaks_per_pixel = 1;
		conn_data[di].name = connected->rigid_groups[di]->name;
		conn_data[di].n_peaks_in_conn = 0;
	}


	for ( ipx=0; ipx<num_pix_in_slab; ipx++ ) {
		all_pix_displ[ipx].dfs = dfv;
		all_pix_displ[ipx].dss = dfv;
		all_pix_displ[ipx].ne = NULL;
		curr_pix_displ[ipx] = &all_pix_displ[ipx];
		num_pix_displ[ipx] = 0;
	}

	for ( cp=0; cp<pattern_list->n_patterns; cp++ ) {

		ImageFeatureList *flist = pattern_list->patterns[cp]->im_list;

		if ( only_best_distance ) {
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
			double rfs, rss;
			double crx, cry;
			Reflection *refl;

			RefList *rlist = pattern_list->patterns[cp]->ref_list;

			struct imagefeature *imfe = image_get_feature(flist, ich);
			compute_x_y(det, imfe->fs, imfe->ss, &fx, &fy);

			refl = find_closest_reflection(rlist, fx, fy, det,
			                               &min_dist);

			if ( refl == NULL ) continue;

			if ( min_dist < max_peak_distance ) {

				int ipx = ((int)rint(imfe->fs) + aw*
				           (int)rint(imfe->ss));

				if ( num_pix_displ[ipx]>0 ) {

					int ret;

					ret = allocate_next_element(curr_pix_displ,
					                            ipx);

					if ( ret != 0) return ret;

				}

				get_detector_pos(refl, &rfs, &rss);
				compute_x_y(det, rfs, rss, &crx, &cry);
				get_detector_pos(refl, &rfs, &rss);
				curr_pix_displ[ipx]->dfs = (fx-crx);
				curr_pix_displ[ipx]->dss = (fy-cry);
				curr_pix_displ[ipx]->ne = NULL;
				num_pix_displ[ipx]++;
			}
		}
	}

	for ( np=min_num_peaks_per_pixel; np>0; np-- ) {
		for ( di=0; di<connected->n_rigid_groups; di++ ) {
			if ( conn_data[di].num_peaks_per_pixel>np ) {
				continue;
			}

			for (ip=0; ip<connected->rigid_groups[di]->n_panels;
			     ip++) {

				int ret;

				ret = compute_panel_statistics(connected,
				                               num_pix_displ,
				                               curr_pix_displ,
				                               all_pix_displ,
				                               conn_data, di, ip,
				                               np,
				                               dfv,
				                               aw,
				                               displ_x,
				                               displ_y,
							       displ_abs);

				if ( ret !=0 ) return ret;
			}


			if ( conn_data[di].n_peaks_in_conn >=
			     min_num_peaks_per_panel ) {
				conn_data[di].num_peaks_per_pixel = np;
			}
		}
	}

    return 0;
}


static double compute_error(struct rg_collection *connected,
                            int aw,
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
					if ( num_pix_displ[ifs+aw*iss]>=
					     conn_data[di].num_peaks_per_pixel ) {

						double cer;

						cer = displ_abs[ifs+aw*iss]*
						      displ_abs[ifs+aw*iss];
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


static void fill_coordinate_matrices(struct detector *det, int aw,
                                     double *slab_to_x, double *slab_to_y)
{
	int pi;

	for ( pi=0; pi<det->n_panels; pi++ ) {

		struct panel *p;
		int iss, ifs;

		p = &det->panels[pi];

        for (iss=p->min_ss; iss < p->max_ss+1; iss++) {
            for (ifs=p->min_fs; ifs < p->max_fs+1; ifs++) {

				double xs, ys;

				compute_x_y(det, ifs, iss, &xs, &ys);

				slab_to_x[iss*aw+ifs] = xs;
				slab_to_y[iss*aw+ifs] = ys;

			}
		}
	}
}


static int correct_empty_panels(struct rg_collection *quadrants,
                                 struct rg_collection *connected,
                                 int min_num_peaks_per_panel,
                                 struct connected_data *conn_data)
{
	double *aver_ang;
	double *aver_str;
	int *aver_num_ang;

	int di,i;

	aver_ang = malloc(quadrants->n_rigid_groups*sizeof(double));
	if ( aver_ang == NULL ) {
		ERROR("Failed to allocate memory to correct empty panels.|n");
		return 1;
	}
	aver_str = malloc(quadrants->n_rigid_groups*sizeof(double));
	if ( aver_str == NULL ) {
		ERROR("Failed to allocate memory to correct empty panels.|n");
		free(aver_ang);
		return 1;
	}
	aver_num_ang = malloc(quadrants->n_rigid_groups*sizeof(int));
	if ( aver_num_ang == NULL ) {
		ERROR("Failed to allocate memory to correct empty panels.|n");
		free(aver_ang);
		free(aver_str);
		return 1;
	}

	for (i=0; i<quadrants->n_rigid_groups; i++) {
		aver_ang[i] = 0;
		aver_str[i] = 0;
		aver_num_ang[i] = 0;
	}

	for (di=0; di<connected->n_rigid_groups; di++) {
		if ( conn_data[di].n_peaks_in_conn>=min_num_peaks_per_panel ) {
			aver_ang[conn_data[di].num_quad] += conn_data[di].cang;
			aver_str[conn_data[di].num_quad] += conn_data[di].cstr;
			aver_num_ang[conn_data[di].num_quad]++;
		}
	}

	for ( i=0; i<quadrants->n_rigid_groups; i++ ) {
		if ( aver_num_ang[i]>0 ) {
			aver_ang[i] /= (double)aver_num_ang[i];
			aver_str[i] /= (double)aver_num_ang[i];
		}
	}

	for ( di=0; di<connected->n_rigid_groups; di++ ) {

		if ( conn_data[di].n_peaks_in_conn<min_num_peaks_per_panel ) {
			if (aver_num_ang[conn_data[di].num_quad]>0) {
				conn_data[di].cang =
				               aver_ang[conn_data[di].num_quad];
				conn_data[di].cstr =
				               aver_str[conn_data[di].num_quad];
				STATUS("Connected group %s has not enough peaks "
				       "(%i). Using average angle: %0.4f degrees\n",
				       conn_data[di].name,
				       conn_data[di].n_peaks_in_conn,
				       conn_data[di].cang);
			} else {
				STATUS("Connected group %s has not enough peaks "
				       "(%i). Left unchanged\n",
				       conn_data[di].name,
				       conn_data[di].n_peaks_in_conn);
			}
		}
	}

	free(aver_ang);
	free(aver_str);
	free(aver_num_ang);

    return 0;
}


static void correct_angle_and_stretch(struct rg_collection *connected,
                                      struct detector *det,
                                      struct connected_data *conn_data,
                                      double use_clen, double stretch_coeff,
                                      int individual_coffset)
{

	int di, ip;

	for ( di=0; di<connected->n_rigid_groups; di++ ) {
		for ( ip=0; ip<connected->rigid_groups[di]->n_panels; ip++ ) {

			struct panel *p;
			double newx, newy;

			p = connected->rigid_groups[di]->panels[ip];

			newx =
			    p->fsx*cos(conn_data[di].cang)-
			    p->fsy*sin(conn_data[di].cang);
			newy =
			    p->fsx*sin(conn_data[di].cang)+
			    p->fsy*cos(conn_data[di].cang);
			p->fsx = newx;
			p->fsy = newy;
			newx =
			    p->ssx*cos(conn_data[di].cang)-
			    p->ssy*sin(conn_data[di].cang);
			newy =
			    p->ssx*sin(conn_data[di].cang)+
			    p->ssy*cos(conn_data[di].cang);
			p->ssx = newx;
			p->ssy = newy;
		}
	}


	if ( individual_coffset == 0 ) {

		int pi;

		for (pi=0; pi<det->n_panels; pi++) {
			det->panels[pi].coffset -= use_clen*(1.0-stretch_coeff);
		}
		STATUS("Using a single offset distance for the whole detector: "
		       "%f m.\n", det->panels[0].coffset);

		for ( di=0; di<connected->n_rigid_groups; di++ ) {
			conn_data[di].cstr = stretch_coeff;
		}

	} else {

		STATUS("Using individual distances for rigid panels.\n");
		for ( di=0; di<connected->n_rigid_groups; di++ ) {
			for ( ip=0; ip<connected->rigid_groups[di]->n_panels;
			      ip++ ) {

				struct panel *p;

				p = connected->rigid_groups[di]->panels[ip];
				p->coffset -= (1.0-conn_data[di].cstr)*p->clen;

			}
		}
	}
}


static void shift_panels(struct rg_collection *connected,
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

				delta_x = (p->cnx-p0->cnx)/conn_data[di].cstr;
				delta_y = (p->cny-p0->cny)/conn_data[di].cstr;

				p->cnx = p0->cnx + delta_x;
				p->cny = p0->cny + delta_y;

			}
		}
	}
}


static void recompute_panel(struct connected_data *conn_data, int di, int ip,
                           struct rg_collection *connected,
                           double *slab_to_x,
                           double *slab_to_y,
			   double *recomputed_slab_to_x,
                           double *recomputed_slab_to_y,
                           double *displ_x, double *displ_y,
			   double stretch_coeff, int aw, int *num_pix_displ)
{

	double c_stretch;
	struct panel *p;
	int ifs, iss;

	c_stretch = conn_data[di].cstr;

	if ( fabs(c_stretch)<FLT_EPSILON ) c_stretch =
	                                          stretch_coeff;

	p = connected->rigid_groups[di]->panels[ip];

	for ( ifs=p->min_fs; ifs<p->max_fs+1; ifs++ ) {
		for ( iss=p->min_ss; iss<p->max_ss+1; iss++ ) {
			recomputed_slab_to_x[ifs+aw*iss] /= c_stretch;
			recomputed_slab_to_y[ifs+aw*iss] /= c_stretch;
			if ( num_pix_displ[ifs+aw*iss] >=
			     conn_data[di].num_peaks_per_pixel) {

			displ_x[ifs+aw*iss] -= (slab_to_x[ifs+aw*iss]-
				              recomputed_slab_to_x[ifs+aw*iss]);
			displ_y[ifs+aw*iss] -= (slab_to_y[ifs+aw*iss]-
				              recomputed_slab_to_y[ifs+aw*iss]);
			}
		}
	}
}


static void recompute_differences(struct rg_collection *connected,
                             double *slab_to_x, double *slab_to_y,
                             double *recomputed_slab_to_x,
                             double *recomputed_slab_to_y,
                             struct connected_data *conn_data,
                             double stretch_coeff, int aw,
                             double *displ_x, double *displ_y,
                             int *num_pix_displ)
{

	int di, ip;

	for ( di=0; di<connected->n_rigid_groups; di++ ) {
		for (ip=0; ip<connected->rigid_groups[di]->n_panels; ip++) {


			recompute_panel(conn_data, di, ip, connected,
			                slab_to_x, slab_to_y,
			                recomputed_slab_to_x,
			                recomputed_slab_to_y,
			                displ_x, displ_y,
			                stretch_coeff, aw, num_pix_displ);
		}
	}
}


static void fill_av_in_panel(struct rg_collection *connected, int di, int ip,
                             struct connected_data *conn_data,
                             int *num_pix_displ, int aw,
                             double *av_in_panel_fs,
                             double *av_in_panel_ss,
	                     double *displ_x, double *displ_y)
{
	struct panel *p;
	int ifs, iss;

	p = connected->rigid_groups[di]->panels[ip];

	for ( ifs=p->min_fs; ifs<p->max_fs+1; ifs++ ) {
		for ( iss=p->min_ss; iss<p->max_ss+1; iss++ ) {

			if (num_pix_displ[ifs+aw*iss]>=
			    conn_data[di].num_peaks_per_pixel) {
				av_in_panel_fs[conn_data[di].n_peaks_in_conn] =
					      displ_x[ifs+aw*iss];
				av_in_panel_ss[conn_data[di].n_peaks_in_conn] =
					      displ_y[ifs+aw*iss];
				conn_data[di].n_peaks_in_conn++;
			}
		}
	}
}


static void fill_con_data_sh(struct connected_data *conn_data,
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
	if ( modulus2d(conn_data[di].sh_x, conn_data[di].sh_y) >
			               0.8*max_peak_distance ) {
		STATUS("  WARNING: absolute shift is: %0.1f > 0.8*%0.1f pixels."
		       "  Increase the value of the max_peak_distance parameter!\n",
		       modulus2d(conn_data[di].sh_x, conn_data[di].sh_y),
		       max_peak_distance);
	}
}


static int compute_shifts(struct rg_collection *connected,
                          struct connected_data *conn_data,
                          int *num_pix_displ, int aw,
                          int min_num_peaks_per_panel,
                          double dfv, double max_peak_distance,
                          double *displ_x, double *displ_y)
{
	STATUS("Median for panels.\n");

	int di, ip;

	for ( di=0; di<connected->n_rigid_groups; di++ ) {

		int cmaxfs;
		int cmaxss;
		int num_all_pixels;
		double *av_in_panel_fs;
		double *av_in_panel_ss;

		cmaxfs = connected->rigid_groups[di]->panels[0]->max_fs+1-
		         connected->rigid_groups[di]->panels[0]->min_fs;
		cmaxss = connected->rigid_groups[di]->panels[0]->max_ss+1-
	                 connected->rigid_groups[di]->panels[0]->min_ss;

		num_all_pixels = cmaxfs*cmaxss*connected->rigid_groups[di]->n_panels;

		av_in_panel_fs = malloc(num_all_pixels*sizeof(double));
		if  ( av_in_panel_fs == NULL ) {
			ERROR("Failed to allocate memory for computing shifts\n");
			return 1;
		}
		av_in_panel_ss = malloc(num_all_pixels*sizeof(double));
		if  ( av_in_panel_ss == NULL ) {
			ERROR("Failed to allocate memory for computing shifts\n");
			free(av_in_panel_fs);
			return 1;
		}

		conn_data[di].n_peaks_in_conn = 0;

		for (ip=0; ip<connected->rigid_groups[di]->n_panels; ip++) {

			fill_av_in_panel(connected, di, ip, conn_data,
			                 num_pix_displ, aw, av_in_panel_fs,
			                 av_in_panel_ss, displ_x, displ_y);

		}

		if ( conn_data[di].n_peaks_in_conn>=min_num_peaks_per_panel ) {

			fill_con_data_sh(conn_data, av_in_panel_fs,
			                 av_in_panel_ss, di,
			                 max_peak_distance);

		} else {
			conn_data[di].sh_x = dfv;
			conn_data[di].sh_y = dfv;
		}
		free(av_in_panel_fs);
		free(av_in_panel_ss);
	}

	return 0;
}


static int compute_shifts_for_empty_panels(struct rg_collection *quadrants,
                                            struct rg_collection *connected,
                                            struct connected_data *conn_data,
                                            int min_num_peaks_per_panel)
{

	double *aver_x;
	double *aver_y;
	int *num_aver;
	int di, i;

	// shifts for empty
	aver_x = malloc(quadrants->n_rigid_groups*sizeof(double));
	if ( aver_x == NULL ) {
		ERROR("Failed to allocate memory to compute shifts for "
		      "empty panels.\n");
		return 1;
	}
	aver_y = malloc(quadrants->n_rigid_groups*sizeof(double));
	if ( aver_y == NULL ) {
		ERROR("Failed to allocate memory to compute shifts for "
		      "empty panels.\n");
		free(aver_x);
		return 1;
	}
	num_aver = malloc(quadrants->n_rigid_groups*sizeof(int));
	if ( num_aver == NULL ) {
		ERROR("Failed to allocate memory to compute shifts for "
		      "empty panels.\n");
		free(aver_x);
		free(aver_y);
		return 1;
	}

	for ( i=0; i<quadrants->n_rigid_groups; i++ ) {
		aver_x[i] = 0;
		aver_y[i] = 0;
		num_aver[i] = 0;
	}

	for ( di=0; di<connected->n_rigid_groups; di++ ) {
		if ( conn_data[di].n_peaks_in_conn>=min_num_peaks_per_panel ) {
			aver_x[conn_data[di].num_quad] += conn_data[di].sh_x;
			aver_y[conn_data[di].num_quad] += conn_data[di].sh_y;
			num_aver[conn_data[di].num_quad]++;
		}
	}

	for ( i=0; i<quadrants->n_rigid_groups; i++ ) {
		if (num_aver[i]>0) {
			aver_x[i] /= (double)num_aver[i];
			aver_y[i] /= (double)num_aver[i];
		}
	}

	for (di=0; di<connected->n_rigid_groups; di++) {
		if ( conn_data[di].n_peaks_in_conn<min_num_peaks_per_panel ) {
			if ( num_aver[conn_data[di].num_quad]>0 ) {
				conn_data[di].sh_x = aver_x[conn_data[di].num_quad];
				conn_data[di].sh_y = aver_y[conn_data[di].num_quad];
				STATUS("Panel %s has not enough (%i) peaks. "
				       "Using average shifts (in pixels) X,Y: "
				       "%0.2f,%0.2f\n", conn_data[di].name,
				       conn_data[di].n_peaks_in_conn,
				       conn_data[di].sh_x, conn_data[di].sh_y);
			} else {
				STATUS("Panel %s has not enough (%i) peaks. "
				       "Left unchanged.\n",
				       conn_data[di].name,
				       conn_data[di].n_peaks_in_conn);
			}
		}
	}

	free(aver_x);
	free(aver_y);
	free(num_aver);

	return 0;
}


static void correct_shifts(struct rg_collection *connected,
                           struct connected_data *conn_data,
                           double dfv, double clen_to_use)
{

	int di;
	int ip;

	for ( di=0;di<connected->n_rigid_groups;di++ ) {
		for (ip=0; ip<connected->rigid_groups[di]->n_panels; ip++) {

			struct panel *p;

			p = connected->rigid_groups[di]->panels[ip];

			if ( conn_data[di].sh_x>dfv+1.0 &&
			     conn_data[di].sh_y > dfv+1.0 ) {

				p->cnx -= conn_data[di].sh_x;
				p->cny -= conn_data[di].sh_y;

			} else {
				STATUS("For some reason pannel %s is empty!\n",
				       p->name);
			}
		}
	}
}


static void a_s_counting_loop(int *num_pix_displ, int ifs, int iss,
                              int di, struct connected_data *conn_data,
                              int aw, double *slab_to_x,
                              double *slab_to_y, struct panel *p0,
                              struct panel *p1, double *displ_x,
                              double *displ_y, double minrad,
                              int *num_ang)
{

	double coX, coY, cdX, cdY;

	if ( num_pix_displ[ifs+aw*iss]>=
	     conn_data[di].num_peaks_per_pixel ) {

		int ifs1, iss1;
		int max_fs1_tmp = p0->max_fs;
		int max_ss1_tmp = p0->max_ss;

		coX = slab_to_x[ifs+aw*iss];
		coY = slab_to_y[ifs+aw*iss];
		cdX = coX - displ_x[ifs+aw*iss];
		cdY = coY - displ_y[ifs+aw*iss];

		for (ifs1=ifs+1; ifs1<max_fs1_tmp+1; ifs1++) {

			if ( ifs1 == max_fs1_tmp ) {
				max_fs1_tmp = p1->max_fs;
			}

			for (iss1=iss+1; iss1<max_ss1_tmp+1; iss1++) {

				if ( iss1 == max_ss1_tmp ) {
					max_ss1_tmp = p1->max_ss;
				}

				if ( num_pix_displ[ifs1+aw*iss1]>=
				     conn_data[di].num_peaks_per_pixel ) {

					double dist;
					double coX1, coY1, cdX1, cdY1;
					double len1, len2;

					dist = modulus2d(ifs-ifs1,iss-iss1);
					if ( dist < minrad ) continue;
					coX1 = slab_to_x[ifs1+aw*iss1];
					coY1 = slab_to_y[ifs1+aw*iss1];
					cdX1 =
					  coX1 - displ_x[ifs1+aw*iss1];
					cdY1 =
					  coY1 - displ_y[ifs1+aw*iss1];

					len1 = modulus2d(coX1-coX, coY1-coY);
					len2 = modulus2d(cdX1-cdX, cdY1-cdY);
					if ( len1<FLT_EPSILON ||
					     len2<FLT_EPSILON ) {
						continue;
					}

					*num_ang = *num_ang+1;
				}
			}
		}
	}
}


static int a_s_processing_loop(int *num_pix_displ, int ifs, int iss,
                              int di, struct connected_data *conn_data,
                              int aw, double *slab_to_x,
                              double *slab_to_y, struct panel *p0,
                              struct panel *p1, double *displ_x,
                              double *displ_y, double minrad,
                              int max_num_ang, int *num_ang,
			      double *angles, double *stretches)
{
	double coX, coY, cdX, cdY;

	if ( num_pix_displ[ifs+aw*iss]>=
	     conn_data[di].num_peaks_per_pixel ) {

		int ifs1, iss1;
		int max_fs1_tmp = p0->max_fs;
		int max_ss1_tmp = p0->max_ss;

		if ( *num_ang>=max_num_ang ) return -2;

		coX = slab_to_x[ifs+aw*iss];
		coY = slab_to_y[ifs+aw*iss];
		cdX = coX - displ_x[ifs+aw*iss];
		cdY = coY - displ_y[ifs+aw*iss];

		for (ifs1=ifs+1; ifs1<max_fs1_tmp+1; ifs1++) {

			if ( ifs1 == max_fs1_tmp ) {
				max_fs1_tmp = p1->max_fs;
			}

			for (iss1=iss+1; iss1<max_ss1_tmp+1; iss1++) {

				if ( iss1 == max_ss1_tmp ) {
					max_ss1_tmp = p1->max_ss;
				}

				if ( num_pix_displ[ifs1+aw*iss1]>=
				     conn_data[di].num_peaks_per_pixel ) {

					double dist;
					double coX1, coY1, cdX1, cdY1;
					double len1, len2;
					double scalM;
					double multlen;

					if ( *num_ang>=max_num_ang ) return -2;
					dist = modulus2d(ifs-ifs1,iss-iss1);
					if (dist<minrad) return 0;
					coX1 = slab_to_x[ifs1+aw*iss1];
					coY1 = slab_to_y[ifs1+aw*iss1];
					cdX1 =
					  coX1 - displ_x[ifs1+aw*iss1];
					cdY1 =
					  coY1 - displ_y[ifs1+aw*iss1];

					len1 = modulus2d(coX1-coX, coY1-coY);
					len2 = modulus2d(cdX1-cdX, cdY1-cdY);
					scalM = (coX1-coX)*(cdX1-cdX)+
						(coY1-coY)*(cdY1-cdY)-
					       FLT_EPSILON;
					if ( len1<FLT_EPSILON ||
					     len2<FLT_EPSILON ) {
						return 0;
					}

					multlen = len1*len2;
					if ( fabs(scalM)>=multlen ) {
						angles[*num_ang] = 0.0;
					} else {

						angles[*num_ang] = 1.0;

						angles[*num_ang] =
						    acos(scalM/multlen);

						if ((coY1-coY)*(cdX1-cdX)-
						   (coX1-coX)*(cdY1-cdY) < 0) {
							angles[*num_ang] *= -1.;
						}

					}

					stretches[*num_ang] = len1/len2;

					*num_ang = *num_ang+1;
				}
			}
		}
	}
	return 0;
}




static int compute_angles_and_stretch(struct rg_collection *connected,
               struct connected_data *conn_data,
               int *num_pix_displ,
               double *slab_to_x,
               double *slab_to_y,
               double *displ_x,
               double *displ_y,
               int aw,
               int min_num_peaks_per_panel,
               double dist_coeff_ang_str,
               int num_peaks_per_pixel,
               double man_stretching_coeff,
               double *stretch_coeff)
{
	int di;
	int num_coeff;
	double stretch_cf;

	struct connected_stretch_and_angles *csaa;

	csaa = malloc(sizeof(struct connected_stretch_and_angles));
	if ( csaa == NULL ) {
		ERROR("Failed to allocate memory to compute angles and "
		      "stretch.\n");
		return 1;
	}
	csaa->stretch_coeff = malloc(connected->n_rigid_groups*sizeof(double));
	if ( csaa->stretch_coeff == NULL ) {
		ERROR("Failed to allocate memory to compute angles and "
		      "stretch.\n");
		free(csaa);
		return 1;
	}
	csaa->num_angles = malloc(connected->n_rigid_groups*sizeof(unsigned int));
	if ( csaa->num_angles == NULL ) {
		ERROR("Failed to allocate memory to compute angles and "
		      "stretch.\n");
		free(csaa->stretch_coeff);
		free(csaa);
		return 1;
	}

	csaa->num_coeff=0;

	for ( di=0; di<connected->n_rigid_groups; di++ ) {
		if ( conn_data[di].n_peaks_in_conn<min_num_peaks_per_panel ) {
			continue;
		}

		unsigned int max_num_ang = 0;

		double* angles;
		double* stretches;

		int cmaxfs;
		int cmaxss;

		struct panel *p;
		double minrad;
		int num_ang = 0;
		int ip0, ip1;

		p = connected->rigid_groups[di]->panels[0];

		cmaxfs = p->max_fs+1-p->min_fs;
		cmaxss = p->max_ss+1-p->min_ss;

		// TODO: MINRAD HERE IS NOT UNIVERSAL
		minrad = dist_coeff_ang_str*sqrt(cmaxfs*cmaxss*
		         connected->rigid_groups[di]->n_panels);

		for ( ip0=0; ip0<connected->rigid_groups[di]->n_panels; ip0++ ) {

			struct panel *p0 = connected->rigid_groups[di]->panels[ip0];

			for ( ip1=0; ip1<connected->rigid_groups[di]->n_panels;
			      ip1++ ) {

				struct panel *p1 =
				      connected->rigid_groups [di]->panels[ip1];

				int ifs, iss;
				int min_fs_tmp = p0->min_fs;
				int max_fs_tmp = p0->max_fs;
				int min_ss_tmp = p0->min_ss;
				int max_ss_tmp = p0->max_ss;

				for (ifs=min_fs_tmp; ifs<max_fs_tmp+1; ifs++) {

					if ( ifs == max_fs_tmp ) {
						min_fs_tmp = p1->min_fs;
						max_fs_tmp = p1->max_fs;
					}

					for (iss=min_ss_tmp; iss<max_ss_tmp+1;
					     iss++) {

						if ( iss == max_ss_tmp ) {
							min_ss_tmp = p1->min_ss;
							max_ss_tmp = p1->max_ss;
						}

						a_s_counting_loop(num_pix_displ,
						       ifs, iss, di, conn_data,
						       aw, slab_to_x, slab_to_y,
						       p0, p1, displ_x,
						       displ_y, minrad,
						       &num_ang);

					}
				}
			}
		}

		max_num_ang = conn_data[di].n_peaks_in_conn*
		              conn_data[di].n_peaks_in_conn;
		max_num_ang = num_ang+1;

		angles = malloc(max_num_ang*sizeof(double));
		if ( angles == NULL ) {
			ERROR("Error in allocating memory for angle "
			      "optimization\n");
			free(csaa->stretch_coeff);
			free(csaa->num_angles);
			free(csaa);
			return 1;
		}
		stretches = malloc(max_num_ang*sizeof(double));
		if ( stretches == NULL ) {
			ERROR("Error in allocating memory for stretch "
			      "optimization\n");
			free(angles);
			free(csaa->stretch_coeff);
			free(csaa->num_angles);
			free(csaa);
			return 1;
		}

		num_ang = 0;

		for ( ip0=0; ip0<connected->rigid_groups[di]->n_panels; ip0++ ) {

			struct panel *p0 = connected->rigid_groups[di]->panels[ip0];

			for ( ip1=0; ip1<connected->rigid_groups[di]->n_panels;
			      ip1++ ) {

				struct panel *p1 =
				      connected->rigid_groups [di]->panels[ip1];

				int ifs, iss;
				int min_fs_tmp = p0->min_fs;
				int max_fs_tmp = p0->max_fs;
				int min_ss_tmp = p0->min_ss;
				int max_ss_tmp = p0->max_ss;

				for (ifs=min_fs_tmp; ifs<max_fs_tmp+1; ifs++) {

					if ( ifs == max_fs_tmp ) {
						min_fs_tmp = p1->min_fs;
						max_fs_tmp = p1->max_fs;
					}

					for (iss=min_ss_tmp; iss<max_ss_tmp+1;
					     iss++) {

						int ret;

						if ( iss == max_ss_tmp ) {
							min_ss_tmp = p1->min_ss;
							max_ss_tmp = p1->max_ss;
						}

						ret = a_s_processing_loop(
						           num_pix_displ,
						           ifs, iss, di,
						           conn_data,
						           aw, slab_to_x,
						           slab_to_y,
						           p0, p1, displ_x,
						           displ_y, minrad,
						           max_num_ang,
						           &num_ang, angles,
						           stretches);

						if ( ret == -2 ) break;
					}
				}
			}
		}

		if ( num_ang<1 ) continue;
		conn_data[di].cang = -comp_median(angles,num_ang);
		conn_data[di].cstr = comp_median(stretches,num_ang);

		STATUS("Panel %s, num: %i, angle: %0.4f deg, stretch coeff: "
		       "%0.4f\n",
		conn_data[di].name, num_ang, conn_data[di].cang,
                conn_data[di].cstr);

		csaa->stretch_coeff[csaa->num_coeff] = conn_data[di].cstr;
		csaa->num_angles[csaa->num_coeff] = num_ang;
		csaa->num_coeff++;

		free(angles);
		free(stretches);
	}


	num_coeff = csaa->num_coeff;

	stretch_cf = 1;
	if (num_coeff>0) {

		int ipp;

		for ( ipp=num_peaks_per_pixel; ipp>=0; ipp-- ) {

			double total_num;
			int di;

			total_num = 0;
			for ( di=0; di<num_coeff; di++ ) {
				if ( conn_data[di].num_peaks_per_pixel>=ipp ) {
					total_num += csaa->num_angles[di];
				}
			}
			if ( total_num>1 ) {
				total_num = 1./total_num;
			} else {
				continue;
			}
			stretch_cf = 0;
			for ( di=0; di<num_coeff; di++ ) {
				if ( conn_data[di].num_peaks_per_pixel>=ipp ) {
					stretch_cf +=
					      total_num*csaa->stretch_coeff[di]*
					           (double)csaa->num_angles[di];
				}
			}
			break;
		}
	}

	if ( stretch_cf<FLT_EPSILON ) {
		stretch_cf = 1.0;
	}

	STATUS("The stretch coefficient for the patterns is %0.4f\n",
	       stretch_cf);
	if ( man_stretching_coeff>FLT_EPSILON ) {
		stretch_cf = man_stretching_coeff;
		STATUS("Using manually set stretch coefficient: %0.4f\n",
		       stretch_cf);

		for ( di=0; di<connected->n_rigid_groups; di++ ) {
			conn_data[di].cstr = man_stretching_coeff;
		}

	}

	free(csaa->stretch_coeff);
	free(csaa->num_angles);
	free(csaa);

	*stretch_coeff = stretch_cf;

	return 0;
}


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


static int unpack_slab(struct image *image)
{
	struct detector *det = image->det;
	int pi;

	image->dp = malloc(det->n_panels * sizeof(float *));
	image->bad = malloc(det->n_panels * sizeof(int *));
	if ( (image->dp == NULL) || (image->bad == NULL) ) {
		ERROR("Failed to allocate panels.\n");
		return 1;
	}

	for ( pi=0; pi<det->n_panels; pi++ ) {

		struct panel *p;
		int fs, ss;

		p = &det->panels[pi];
		image->dp[pi] = malloc(p->w*p->h*sizeof(float));
		image->bad[pi] = calloc(p->w*p->h, sizeof(int));
		if ( (image->dp[pi] == NULL) || (image->bad[pi] == NULL) ) {
			ERROR("Failed to allocate panel\n");
			return 1;
		}

		for ( ss=0; ss<p->h; ss++ ) {
		for ( fs=0; fs<p->w; fs++ ) {

			int idx;
			int cfs, css;

			cfs = fs+p->min_fs;
			css = ss+p->min_ss;
			idx = cfs + css*image->width;

			image->dp[pi][fs+p->w*ss] = image->data[idx];
			image->bad[pi][fs+p->w*ss] = 0;

		}
		}
	}

	return 0;
}


static int draw_detector(cairo_surface_t *surf, struct image *image,
                         struct rectangle rect)
{
	cairo_t *cr;
	cairo_matrix_t basic_m;
	cairo_matrix_t m;
	GdkPixbuf **pixbufs;
	int n_pixbufs;

	cr = cairo_create(surf);

	unpack_slab(image);
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


static int save_data_to_png(char *filename, struct detector *det,
                            int max_fs, int max_ss, double default_fill_value,
                            double *data)
{
	struct image im;
	int i;
	struct rectangle rect;
	GdkPixbuf *col_scale;
	cairo_t *cr;

	cairo_status_t r;
	cairo_surface_t *surf;

	im.data = malloc((max_fs+1)*(max_ss+1)*sizeof(float));
	if ( im.data == NULL ) {
		ERROR("Failed to allocate memory to save data.\n");
		return 1;
	}
	im.det = det;
	im.width = max_fs+1;
	im.height = max_ss+1;
	im.flags = NULL;

	for ( i=0; i<(max_fs+1)*(max_ss+1); i++) {
		if ( data[i] == default_fill_value ) {
			im.data[i] = 0.0;
		} else if ( data[i] > 1.0) {
			im.data[i] = 1.0;
		} else {
			im.data[i] = (float)data[i];
		}
		im.data[i] *= 10.0; /* render_panels sets this as max */
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

	r = cairo_surface_write_to_png(surf, filename);
	if (r != CAIRO_STATUS_SUCCESS) {
		free(im.data);
		return 1;
	}

	free(im.data);

	return 0;
}


static void calculate_panel_correction(int di, int ip, int aw,
                                       int *num_pix_displ,
                                       struct rg_collection *connected,
                                       struct connected_data *conn_data)
{
	struct panel *p;
	int ifs, iss;

	p = connected->rigid_groups[di]->panels[ip];
	for (ifs=p->min_fs; ifs<p->max_fs+1; ifs++) {
		for (iss=p->min_ss; iss<p->max_ss+1; iss++) {
			if ( num_pix_displ[ifs+aw*iss]>=
					conn_data[di].num_peaks_per_pixel ) {
					conn_data[di].n_peaks_in_conn++;
			}
		}
	}

}


static void compute_abs_displ(struct rg_collection *connected,
                              struct connected_data *conn_data,
			      int *num_pix_displ,
			      double dfv, int di, int ip, int aw,
		              double *displ_x,
		              double *displ_y,
			      double *displ_abs)
{
	struct panel *p;
	int ifs, iss;

	if (conn_data[di].sh_x < dfv+1) return;

	p = connected->rigid_groups[di]->panels[ip];

	for (ifs=p->min_fs; ifs<p->max_fs+1; ifs++) {
		for (iss=p->min_ss; iss<p->max_ss+1; iss++) {
			if ( num_pix_displ[ifs+aw*iss]>=
			     conn_data[di].num_peaks_per_pixel ) {
				displ_x[ifs+aw*iss] -= conn_data[di].sh_x;
				displ_y[ifs+aw*iss] -= conn_data[di].sh_y;
				displ_abs[ifs+aw*iss] = modulus2d(
							   displ_x[ifs+aw*iss],
							   displ_y[ifs+aw*iss]
							   );
			} else {
				displ_abs[ifs+aw*iss] = dfv;
			}
		}
	}
}


int check_and_enforce_cspad_dist(struct detector *det, int enforce)
{
	int np = 0;
	int num_errors_found = 0;

	double dist_to_check = 197.0;
	double tol = 0.2;

	for ( np=0; np<det->n_panels; np = np+2 ) {

		double dist2;

		struct panel *ep = &det->panels[np];
		struct panel *op = &det->panels[np+1];

		dist2 = (( ep->cnx - op->cnx )*( ep->cnx - op->cnx ) +
			 ( ep->cny - op->cny )*( ep->cny - op->cny ));

		if ( dist2 > (dist_to_check+tol)*(dist_to_check+tol) ||
		     dist2 < (dist_to_check-tol)*(dist_to_check-tol) ) {

			num_errors_found += 1;

			STATUS("Warning: distance between panels %s and %s "
			       "is outside acceptable margins.\n", ep->name,
			       op->name);

			if ( enforce ) {

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

			if ( enforce ) {

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



int optimize_geometry(char *infile, char *outfile, char *geometry_filename,
                      struct detector *det, struct rg_collection* quadrants,
                      struct rg_collection* connected,
                      int min_num_peaks_per_pixel, int min_num_peaks_per_panel,
                      int only_best_distance, int nostretch,
		      int individual_coffset, int error_maps,
		      int enforce_cspad_layout, int no_cspad,
		      double max_peak_dist, const char *command_line)
{
	int num_pix_in_slab;
	int max_fs = 0;
	int max_ss = 0;
	int aw = 0;
	int pi, di, ip, pti;
	int ret1, ret2;
	int ret;
	int write_ret;
	int maybe_cspad = 0;

	int *num_pix_displ;

	double res_sum;
	double istep;
	double clen_to_use;
	double man_stretching_coeff = 0.0;
	double avc[6] = {0.,0.,0.,0.,0.,0.};
	double dfv = -10000.0;
	// for angles and stretch calculation use
        // only pixels which are distco*size_panel
        // away
	double dist_coeff_ang_str = 0.2;
	double *displ_x;
	double *displ_y;
	double *displ_abs;
	double totalError;

	double* slab_to_x;
	double* slab_to_y;
	double* recomputed_slab_to_x;
	double* recomputed_slab_to_y;
	double stretch_coeff = 1;

	struct single_pix_displ *all_pix_displ;
	struct single_pix_displ **curr_pix_displ;
	struct connected_data *conn_data = NULL;
	struct pattern_list *pattern_list;

	if ( nostretch ) man_stretching_coeff = 1.0;

	STATUS("Maximum distance between peaks: %0.1f pixels.\n", max_peak_dist);

	STATUS("Minimum number of measurements for pixel to be included in the "
	       "refinement: %i\n",
	       min_num_peaks_per_pixel);
	STATUS("Minimum number of measurements for panel for accurate estimation "
	       "of position/orientation: %i\n", min_num_peaks_per_panel);

	if ( det->n_panels == 64 ) {
		maybe_cspad = 1;
	}

	if ( maybe_cspad && !no_cspad ) {

		int num_errors = 0;

		STATUS("It looks like the detector is a CSPAD. "
		       "Checking relative distance and orientation of "
		       "connected ASICS.\n");
		STATUS("If the detector is not a CSPAD, please rerun the "
		       "program with the --no-cspad option.\n");
		if ( enforce_cspad_layout ) {
			STATUS("Enforcing CSPAD layout...\n");
		}

		num_errors = check_and_enforce_cspad_dist(det,
							  enforce_cspad_layout);

		if ( enforce_cspad_layout ) {

			int geom_wr;

			STATUS("Saving geometry with enforced CSPAD layout.\n"
			       "Please restart geometry optimization using the "
			       "optimized geometry from this run as input geometry "
			       "file.\n");
			geom_wr = write_detector_geometry_2(geometry_filename,
							    outfile, det,
							    command_line, 1);
			if ( geom_wr != 0 ) {
				ERROR("Error in writing output geometry file.\n");
				return 1;
			}
			STATUS("All done!\n");
			return 0;
		}

		if ( !enforce_cspad_layout && num_errors > 0 ) {
			ERROR("Relative distance and orientation of connected "
			      "ASICS do not respect the CSPAD layout.\n"
			      "Geometry optimization cannot continue.\n"
			      "Please rerun the program with the "
			      "--enforce-cspad-layout option.\n");
			return 1;
		}
	}

	pattern_list = read_patterns_from_steam_file(infile, det);
	if ( pattern_list->n_patterns < 1 ) {
		ERROR("Error reading stream file\n");
		return 1;
	}

	compute_avg_cell_parameters(pattern_list, avc);

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

	istep = res_sum/det->n_panels;

	aw = max_fs+1;

	clen_to_use = compute_clen_to_use(pattern_list, istep, avc,
	                                  max_peak_dist,
	                                  only_best_distance);

	if ( clen_to_use == -1.0 ) return 1;

	num_pix_in_slab = (max_fs+1)*(max_ss+1);
	displ_x = calloc(num_pix_in_slab,sizeof(double));
	if ( displ_x == NULL ) {
		ERROR("Error allocating memory for pixel properties.\n");
		return 1;
	}
	displ_y = calloc(num_pix_in_slab,sizeof(double));
	if ( displ_y == NULL ) {
		ERROR("Error allocating memory for pixel properties.\n");
		free(displ_x);
		return 1;
	}
	displ_abs = calloc(num_pix_in_slab,sizeof(double));
	if ( displ_abs == NULL ) {
		ERROR("Error allocating memory for pixel properties.\n");
		free(displ_x);
		free(displ_y);
		return 1;
	}

	slab_to_x = malloc(num_pix_in_slab*sizeof(double));
	slab_to_y = malloc(num_pix_in_slab*sizeof(double));
	if ( slab_to_x == NULL ) {
		ERROR("Failed to allocate memory for pixel maps.\n");
		free(displ_x);
		free(displ_y);
		free(displ_abs);
		return 1;
	}
	slab_to_y = malloc(num_pix_in_slab*sizeof(double));
	if ( slab_to_y == NULL ) {
		ERROR("Failed to allocate memory for pixel maps.\n");
		free(displ_x);
		free(displ_y);
		free(displ_abs);
		free(slab_to_x);
		return 1;
	}

	fill_coordinate_matrices(det, aw, slab_to_x, slab_to_y);

	all_pix_displ = calloc(num_pix_in_slab,
	                       sizeof(struct single_pix_displ));
	if ( all_pix_displ == NULL ) {
		ERROR("Error allocating memory for connected structure data.\n");
		free(displ_x);
		free(displ_y);
		free(displ_abs);
		free(slab_to_x);
		free(slab_to_y);
		return 1;
	}
	curr_pix_displ = calloc(num_pix_in_slab,
	                        sizeof(struct single_pix_displ*));
	if ( curr_pix_displ == NULL ) {
		ERROR("Error allocating memory for connected structure data.\n");
		free(displ_x);
		free(displ_y);
		free(displ_abs);
		free(slab_to_x);
		free(slab_to_y);
		free(all_pix_displ);
		return 1;
	}
	num_pix_displ = calloc(num_pix_in_slab, sizeof(int));
	if ( num_pix_displ == NULL ) {
		ERROR("Error allocating memory for connected structure data.\n");
		free(displ_x);
		free(displ_y);
		free(displ_abs);
		free(slab_to_x);
		free(slab_to_y);
		free(all_pix_displ);
		free(curr_pix_displ);
		return 1;
	}

	conn_data = malloc(connected->n_rigid_groups*
	                   sizeof(struct connected_data));
	if ( conn_data == NULL ) {
		ERROR("Error allocating memory for connected structure data.\n");
		free(displ_x);
		free(displ_y);
		free(displ_abs);
		free(slab_to_x);
		free(slab_to_y);
		free(all_pix_displ);
		free(curr_pix_displ);
		free(num_pix_displ);
		return 1;
	}

	STATUS("Computing pixel statistics.\n");
	ret = compute_pixel_statistics(pattern_list, det, connected, quadrants,
	                               num_pix_in_slab, max_peak_dist,
	                               aw, dfv,
	                               min_num_peaks_per_pixel,
	                               min_num_peaks_per_panel,
	                               only_best_distance,
	                               clen_to_use, slab_to_x,
	                               slab_to_y, conn_data,
	                               displ_x, displ_y, displ_abs,
	                               all_pix_displ,
	                               curr_pix_displ,
	                               num_pix_displ);
	if ( ret != 0 ) {
		free(displ_x);
		free(displ_y);
		free(displ_abs);
		free(slab_to_x);
		free(slab_to_y);
		free_all_curr_pix_displ(all_pix_displ, curr_pix_displ,
		                        num_pix_in_slab);
		free(num_pix_displ);
		free(conn_data);
		return 1;
	}

	free_all_curr_pix_displ(all_pix_displ, curr_pix_displ, num_pix_in_slab);
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

	if ( error_maps ) {
		STATUS("Saving displacements before corrections.\n");
		ret1 = save_data_to_png("error_map_before.png", det, max_fs, max_ss,
		                        dfv, displ_abs);
		if ( ret1!=0 ) {
			ERROR("Error while writing data to file.\n");
			free(conn_data);
			free(displ_x);
			free(displ_y);
			free(displ_abs);
			free(num_pix_displ);
			free(slab_to_x);
			free(slab_to_y);
			return 1;
		}
	}

	STATUS("Computing initial error.\n");
	totalError = compute_error(connected, aw, conn_data,
	                           num_pix_displ, displ_abs);

	STATUS("The total initial error <delta^2> = %0.4f pixels.\n",
	       totalError);
	STATUS("Now calculating corrections.\n");

	for ( di=0;di<connected->n_rigid_groups;di++ ) {

		conn_data[di].n_peaks_in_conn = 0;

		for (ip=0; ip<connected->rigid_groups[di]->n_panels; ip++) {

			calculate_panel_correction(di, ip, aw, num_pix_displ,
			                           connected, conn_data);

		}
	}

	STATUS("Calculating angles and elongations.\n");

	ret = compute_angles_and_stretch(connected, conn_data,
	                             num_pix_displ,
	                             slab_to_x, slab_to_y,
	                             displ_x, displ_y,
	                             aw,
	                             min_num_peaks_per_panel,
	                             dist_coeff_ang_str,
	                             min_num_peaks_per_pixel,
	                             man_stretching_coeff, &stretch_coeff);

	if ( ret != 0 ) {
		free(conn_data);
		free(displ_x);
		free(displ_y);
		free(displ_abs);
		free(num_pix_displ);
		free(slab_to_x);
		free(slab_to_y);
		return 1;
	}
	ret = correct_empty_panels(quadrants, connected, min_num_peaks_per_panel,
	                           conn_data);

	if ( ret != 0 ) {
		free(conn_data);
		free(displ_x);
		free(displ_y);
		free(displ_abs);
		free(num_pix_displ);
		free(slab_to_x);
		free(slab_to_y);
		return 1;
	}

	correct_angle_and_stretch(connected, det, conn_data, clen_to_use,
	                          stretch_coeff, individual_coffset);

	shift_panels(connected, conn_data);

	recomputed_slab_to_x = malloc(num_pix_in_slab*sizeof(double));
	if ( recomputed_slab_to_x == NULL ) {
		ERROR("Failed to allocate memory for pixel maps.\n");
		free(displ_x);
		free(displ_y);
		free(displ_abs);
		free(slab_to_x);
		free(slab_to_y);
		free(num_pix_displ);
		free(conn_data);
		return 1;
	}
	recomputed_slab_to_y = malloc(num_pix_in_slab*sizeof(double));
	if ( recomputed_slab_to_y == NULL ) {
		ERROR("Failed to allocate memory for pixel maps.\n");
		free(displ_x);
		free(displ_y);
		free(displ_abs);
		free(slab_to_x);
		free(slab_to_y);
		free(num_pix_displ);
		free(conn_data);
		free(recomputed_slab_to_x);
		return 1;
	}

	fill_coordinate_matrices(det, aw, recomputed_slab_to_x,
	                         recomputed_slab_to_y);

	recompute_differences(connected, slab_to_x, slab_to_y,
	                      recomputed_slab_to_x,
	                      recomputed_slab_to_y, conn_data,
	                      stretch_coeff, aw, displ_x, displ_y,
	                      num_pix_displ);

	ret = compute_shifts(connected, conn_data, num_pix_displ, aw,
	                     min_num_peaks_per_panel, dfv,
	                     max_peak_dist, displ_x, displ_y );

	if ( ret != 0 ) return 1;

	compute_shifts_for_empty_panels(quadrants, connected, conn_data,
	                                min_num_peaks_per_panel);

	for ( di=0;di<connected->n_rigid_groups;di++ ) {
		for (ip=0; ip<connected->rigid_groups[di]->n_panels; ip++) {

			compute_abs_displ(connected, conn_data,
			                  num_pix_displ, dfv, di, ip, aw,
			                  displ_x, displ_y, displ_abs);

		}
	}

	correct_shifts(connected, conn_data, dfv, clen_to_use);

	if ( error_maps ) {
		STATUS("Saving displacements after corrections.\n");
		ret2 = save_data_to_png("error_map_after.png", det,  max_fs, max_ss,
		                        dfv, displ_x);
		if ( ret2 !=0 ) {
			ERROR("Error while writing data to file.\n");
			free(conn_data);
			free(displ_x);
			free(displ_y);
			free(displ_abs);
			free(num_pix_displ);
			free(slab_to_x);
			free(slab_to_y);
			free(recomputed_slab_to_x);
			free(recomputed_slab_to_y);
			return 1;
		}
	}

	STATUS("Computing final error.\n");
	totalError = compute_error(connected, aw, conn_data, num_pix_displ,
	                           displ_abs);

	STATUS("The total final error <delta^2> = %0.4f pixels.\n",totalError);

	write_ret = write_detector_geometry_2(geometry_filename, outfile, det,
	                                      command_line, 1);
	if ( write_ret != 0 ) {
		ERROR("Error in writing output geometry file.\n");
		return 1;
	}
	STATUS("All done!\n");
	if ( error_maps ) {
		STATUS("Be sure to inspect error_map_before.png and "
		       "error_map_after.png !!\n");
	}

	free(conn_data);
	free(displ_x);
	free(displ_y);
	free(displ_abs);
	free(num_pix_displ);
	free(slab_to_x);
	free(slab_to_y);
	free(recomputed_slab_to_x);
	free(recomputed_slab_to_y);
	return 0;
}

int main(int argc, char *argv[])
{
	int c, i;
	int ret_val;
	char buffer[256];
	char command_line[1024];
	char *outfile = NULL;
	char *infile = NULL;
	char *geometry_filename = NULL;
	char *quadrant_coll_name = NULL;
	char *connected_coll_name = NULL;
	int min_num_peaks_per_pixel = 3;
	int min_num_peaks_per_panel = 100;
	int only_best_distance = 0;
	int enforce_cspad_layout = 0;
	int nostretch = 0;
	int individual_coffset = 0;
	int no_cspad = 0;
	int error_maps = 1;
	double max_peak_dist = 4.0;

	struct detector *det = NULL;
	struct rg_collection *quadrants;
	struct rg_collection *connected;
	struct beam_params beam;

	const struct option longopts[] = {

        /* Options with long and short versions */
        {"help",                     0, NULL,               'h'},
        {"version",                  0, NULL,               10 },
        {"input",                    1, NULL,               'i'},
        {"output",                   1, NULL,               'o'},
        {"geometry",                 1, NULL,               'g'},
        {"quadrants",                1, NULL,               'q'},
        {"connected",                1, NULL,               'c'},
        {"min-num-peaks-per-pixel",  1, NULL,               'x'},
        {"min-num-peaks-per-panel",  1, NULL,               'p'},
        {"most-few-clen",            0, NULL,               'l'},
        {"max-peak-dist",            1, NULL,               'm'},
        {"individual-dist-offset",   0, NULL,               's'},

        /* Long-only options with no arguments */
        {"no-stretch",               0, &nostretch,                1},
        {"no-error-maps",            0, &error_maps,               0},
        {"enforce-cspad-layout",     0, &enforce_cspad_layout,     1},
        {"no-cspad",                 0, &no_cspad,                 1},

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
			outfile = strdup(optarg);
			break;

			case 'i' :
			infile = strdup(optarg);
			break;

			case 'g' :
			geometry_filename = strdup(optarg);
			det = get_detector_geometry(geometry_filename, &beam);
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
			min_num_peaks_per_pixel = atoi(optarg);
			break;

			case 'p' :
			min_num_peaks_per_panel = atoi(optarg);
			break;

			case 'l' :
			only_best_distance = 1;
			break;

			case 'm' :
			max_peak_dist = strtof(optarg, NULL);
			break;

			case 's' :
			individual_coffset = 1;
			break;

		}
	}

	if ( geometry_filename == NULL ) {
		ERROR("You must provide a geometry to optimize.\n");
		return 1;
	}

	if ( infile == NULL ) {
		ERROR("You must provide an input stream file.\n");
		return 1;
	}

	if ( outfile == NULL ) {
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
	ret_val = optimize_geometry(infile, outfile, geometry_filename, det,
	                        quadrants, connected, min_num_peaks_per_pixel,
	                        min_num_peaks_per_panel, only_best_distance,
				nostretch, individual_coffset, error_maps,
				enforce_cspad_layout, no_cspad,
	                        max_peak_dist, command_line);

	return ret_val;
}
