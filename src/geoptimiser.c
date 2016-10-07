/*
 * geoptimiser.c
 *
 * Refine detector geometry
 *
 * Copyright © 2014-2016 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2014-2015 Oleksandr Yefanov
 *   2014-2015 Valerio Mariani
 *   2014-2016 Thomas White <taw@physics.org>
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
#include <assert.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort.h>

#ifdef HAVE_CAIRO
#ifdef HAVE_GTK
#define HAVE_SAVE_TO_PNG 1
#include <cairo.h>
#include <gdk/gdk.h>
#endif /* HAVE_GTK */
#endif /* HAVE_CAIRO */

#include "detector.h"
#include "stream.h"
#include "version.h"
#include "crystal.h"
#include "image.h"
#include "utils.h"
#include "render.h"

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
"  -p, --min-num-pixels-per-conn-group=<num>    Minimum number of useful pixels per\n"
"                                                connected group.\n"
"                         f                       Default: 100.\n"
"  -l, --most-freq-clen                         Use only the most frequent camera\n"
"                                                length.\n"
"  -s, --individual-dist-offset                 Use a distance offset for each panel.\n"
"                                                Default: whole-detector offset.\n"
"      --no-stretch                             Do not optimize distance offset.\n"
"                                                Default: distance offset is optimized\n"
"  -m  --max-peak-dist=<num>                    Maximum distance between predicted and\n"
"                                                detected peaks (in pixels)\n"
"                                                 Default: half of minimal inter-Bragg distance\n"
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


struct single_pixel_displ
{
	double dx;
	double dy;
	struct single_pixel_displ *ne;
};


struct gpanel
{
	struct panel               *p;

	/* Individual pixel displacements */
	struct single_pixel_displ  *pix_displ_list;
	struct single_pixel_displ **curr_pix_displ;
	int                        *num_pix_displ;

	/* Average displacements for each pixel */
	double                     *avg_displ_x;
	double                     *avg_displ_y;
	double                     *avg_displ_abs;
};


static void compute_x_y(double fs, double ss, struct panel *p,
                        double *x, double *y)
{
	double xs, ys;

	xs = fs*p->fsx + ss*p->ssx;
	ys = fs*p->fsy + ss*p->ssy;

	*x = xs + p->cnx;
	*y = ys + p->cny;
}


static Reflection *find_closest_reflection(struct image *image,
                                           double fx, double fy,
                                           double *d)
{
	double dmin = HUGE_VAL;
	Reflection *closest = NULL;
	int i;

	for ( i=0; i<image->n_crystals; i++ ) {

		Reflection *refl;
		RefListIterator *iter;
		RefList *rlist = crystal_get_reflections(image->crystals[i]);

		for ( refl = first_refl(rlist, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			double ds;
			double rfs, rss;
			double rx, ry;

			get_detector_pos(refl, &rfs, &rss);

			compute_x_y(rfs, rss, get_panel(refl), &rx, &ry);

			ds = distance(rx, ry, fx, fy);

			if ( ds < dmin ) {
				dmin = ds;
				closest = refl;
			}

		}

	}

	if ( closest == NULL ) {
		*d = +INFINITY;
	} else {
		*d = dmin;
	}
	return closest;
}


static double get_average_clen(struct image *image)
{
	int i;
	struct stuff_from_stream *stuff = image->stuff_from_stream;

	if ( stuff == NULL ) {
		ERROR("No 'stuff' from stream!\n");
		return -1.0;
	}

	for ( i=0; i<stuff->n_fields; i++ ) {

		if ( strncmp(stuff->fields[i], "average_camera_length = ",
		             24) == 0 )
		{
			return atof(stuff->fields[i]+24);
		}
	}

	ERROR("Failed to recover average camera length from stream file\n");
	return -1.0;
}


static struct image *read_patterns_from_stream(const char *infile,
                                               struct detector *det, int *n)
{
	Stream *st;
	struct image *images;
	int n_chunks = 0;
	int max_images = 1024;
	int n_images = 0;

	images = malloc(max_images * sizeof(struct image));
	if ( images == NULL ) {
		ERROR("Failed to allocate memory for images.\n");
		return NULL;
	}

	st = open_stream_for_read(infile);
	if ( st == NULL ) {
		ERROR("Failed to open input stream '%s'\n", infile);
		free(images);
		return NULL;
	}

	do {

		images[n_images].det = det;

		if ( read_chunk_2(st, &images[n_images],
		                  STREAM_READ_REFLECTIONS
		                | STREAM_READ_PEAKS
		                | STREAM_READ_UNITCELL) != 0 ) break;

		n_chunks++; /* Number of chunks processed */

		/* Reject if there are no crystals (not indexed) */
		if ( images[n_images].n_crystals == 0 ) continue;

		images[n_images].avg_clen = get_average_clen(&images[n_images]);

		n_images++;  /* Number of images accepted */

		if ( n_images == max_images ) {

			struct image *images_new;

			images_new = realloc(images,
			      (max_images+1024)*sizeof(struct image));
			if ( images_new == NULL ) {
				ERROR("Failed to allocate memory for "
				      "patterns.\n");
				free(images);
				return NULL;
			}

			max_images += 1024;
			images = images_new;
		}

		if ( n_images % 1000 == 0 ) {
			STATUS("Loaded %i indexed patterns from %i total "
			       "patterns.\n", n_images, n_chunks);
		}


	} while ( 1 );

	close_stream(st);
	*n = n_images;

	STATUS("Found %i indexed patterns in file %s (from a total of %i).\n",
	       n_images, infile, n_chunks);

	return images;
}


static struct rvec get_q_from_xyz(double rx, double ry, double dist, double l)
{

	struct rvec q;
	double r = sqrt(rx*rx + ry*ry);
	double twotheta = atan2(r, dist);
	double az = atan2(ry, rx);

	q.u = 1.0/l * sin(twotheta)*cos(az);
	q.v = 1.0/l * sin(twotheta)*sin(az);
	q.w = 1.0/l * (cos(twotheta) - 1.0);

	return q;
}


static UnitCell *compute_avg_cell_parameters(struct image *images, int n)
{
	int numavc;
	int j, i;
	double minc[6];
	double maxc[6];
	double avg_cpar[6] = {0, 0, 0, 0, 0, 0};
	UnitCell *avg;

	for ( j=0; j<6; j++ ) {
		minc[j] = 1e100;
		maxc[j] = -1e100;
	}

	numavc = 0;
	for ( i=0; i<n; i++ ) {

		struct image *image;
		double cpar[6];
		int j, cri;

		image = &images[i];

		for ( cri=0; cri<images->n_crystals; cri++ ) {

			UnitCell *cell = crystal_get_cell(image->crystals[cri]);

			cell_get_parameters(cell,
			                    &cpar[0],  // a
			                    &cpar[1],  // b
			                    &cpar[2],  // c
			                    &cpar[3],  // alpha
			                    &cpar[4],  // beta
			                    &cpar[5]); // gamma

			for ( j=0; j<6; j++ ) {
				avg_cpar[j] += cpar[j];
				if ( cpar[j]<minc[j] ) minc[j] = cpar[j];
				if ( cpar[j]>maxc[j] ) maxc[j] = cpar[j];
			}
			numavc++;

		}

	}

	if ( numavc > 0 ) {
		for ( j=0; j<6; j++ ) avg_cpar[j] /= numavc;
	}

	avg = cell_new();

	cell_set_parameters(avg, avg_cpar[0], avg_cpar[1], avg_cpar[2],
	                         avg_cpar[3], avg_cpar[4], avg_cpar[5]);

	STATUS("Average cell parameters:\n");
	STATUS("Average a, b, c (in A): %6.3f, %6.3f, %6.3f\n",
	       avg_cpar[0]*1e10, avg_cpar[1]*1e10, avg_cpar[2]*1e10);
	STATUS("Minimum -Maximum a, b, c:\n"
	       "\t%6.3f - %6.3f A,\n"
	       "\t%6.3f - %6.3f A,\n"
	       "\t%6.3f - %6.3f A\n",
	       minc[0]*1e10, maxc[0]*1e10, minc[1]*1e10,
	       maxc[1]*1e10, minc[2]*1e10, maxc[2]*1e10);
	STATUS("Average alpha,beta,gamma in degrees: %6.3f, %6.3f, %6.3f\n",
	       rad2deg(avg_cpar[3]), rad2deg(avg_cpar[4]), rad2deg(avg_cpar[5]));
	STATUS("Minimum - Maximum alpha,beta,gamma in degrees:\n"
	       "\t%5.2f - %5.2f,\n"
	       "\t%5.2f - %5.2f,\n"
	       "\t%5.2f - %5.2f\n",
	       rad2deg(minc[3]), rad2deg(maxc[3]),
	       rad2deg(minc[4]), rad2deg(maxc[4]),
	       rad2deg(minc[5]), rad2deg(maxc[5]));

	return avg;
}


static double pick_clen_to_use(struct geoptimiser_params *gparams,
                               struct image *images, int n,
                               double avg_res, UnitCell *avg)
{
	int cp, i, u;
	int num_clens;
	int best_clen;
	int *clens_population;
	double *clens;
	double *lambdas;
	double min_braggp_dist;
	double clen_to_use;
	struct rvec cqu;
	double a, b, c, al, be, ga;

	/* These need to be big enough for the number of different camera
	 * lengths in the data set.  There are probably only a few, but assume
	 * the worst case here - a unique camera length for each frame */
	clens = calloc(n, sizeof(double));
	clens_population = calloc(n, sizeof(int));
	lambdas = calloc(n, sizeof(double));
	if ((lambdas == NULL) || (clens == NULL) || (clens_population == NULL))
	{
		ERROR("Failed to allocate memory for clen calculation.\n");
		free(lambdas);
		free(clens);
		free(clens_population);
		return -1.0;
	}

	num_clens = 0;

	for ( cp=0; cp<n; cp++ ) {

		int i;
		int found = 0;

		for ( i=0; i<num_clens; i++ ) {
			if ( fabs(images[cp].avg_clen - clens[i]) <0.0001 ) {
				clens_population[i]++;
				lambdas[i] += images[cp].lambda;
				found = 1;
				break;
			}
		}

		if ( found ) continue;

		clens[num_clens] = images[cp].avg_clen;
		clens_population[num_clens] = 1;
		lambdas[num_clens] = images[cp].lambda;
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
	cell_get_parameters(avg, &a, &b, &c, &al, &be, &ga);
	for ( i=0; i<num_clens; i++ ) {

		assert(clens_population[i] > 0);

		cqu = get_q_from_xyz(1.0/avg_res, 0, clens[i], lambdas[i]);

		min_braggp_dist = fmin(fmin((1.0/cqu.u)/a,
		                            (1.0/cqu.u)/b),
		                            (1.0/cqu.u)/c);

		STATUS("Camera length %0.4f m was found %i times.\n"
		       "Minimum inter-bragg peak distance (based on "
		       "average cell parameters): %0.1f pixels.\n",
		       clens[i], clens_population[i], min_braggp_dist);

		if ( min_braggp_dist<1.2*gparams->max_peak_dist ) {
			STATUS("WARNING: The distance between Bragg peaks is "
			       "too small: %0.1f < 1.2*%0.1f pixels.\n",
		               min_braggp_dist, gparams->max_peak_dist);
		}

		if ( gparams->max_peak_dist==0.0 ) {
                        gparams->max_peak_dist = 0.5*min_braggp_dist;
			STATUS("WARNING: Maximum distance between peaks is "
			       "set to: %0.1f pixels.\n",
			       gparams->max_peak_dist);
		}

		if ( clens_population[i] > clens_population[best_clen] ) {
			best_clen = i;
			clen_to_use = clens[best_clen];
		}

	}

	if ( gparams->only_best_distance ) {
		STATUS("Only %i patterns with camera length %0.4f m will be "
		       "used.\n", clens_population[best_clen], clen_to_use);
	}

	free(clens);
	free(lambdas);
	free(clens_population);

	return clen_to_use;
}


static double comp_median(double *arr, long n)
{
	gsl_sort(arr, 1, n);
	return gsl_stats_median_from_sorted_data(arr, 1, n);
}


static int find_quad_for_connected(struct rigid_group *rg,
                                   struct rg_collection *quadrants)
{
	struct panel *p;
	int qi;

	/* The quadrant for a group of connected panels is the quadrant to which
	 * the first panel in the connected set belongs */
	p = rg->panels[0];

	for ( qi=0; qi<quadrants->n_rigid_groups; qi++ ) {
		if ( panel_is_in_rigid_group(quadrants->rigid_groups[qi], p) ) {
			return qi;
		}
	}

	/* Hopefully never reached */
	ERROR("Couldn't find quadrant for connected group!\n");
	abort();
}


/* Take all the (valid) displacements for pixel "i" in panel "gp", calculate
 * the median displacements in each direction and the modulus */
static int fill_avg_pixel_displ(struct gpanel *gp, int i)
{
	double *list_dx;
	double *list_dy;
	int count = 0;
	int ei;

	list_dx = calloc(gp->num_pix_displ[i], sizeof(double));
	list_dy = calloc(gp->num_pix_displ[i], sizeof(double));
	if ( (list_dx == NULL) || (list_dy == NULL) ) {
		ERROR("Failed to allocate memory for pixel statistics.\n");
		free(list_dx);
		free(list_dy);
		return 1;
	}

	gp->curr_pix_displ[i] = &gp->pix_displ_list[i];

	for ( ei=0; ei<gp->num_pix_displ[i]; ei++ ) {

		struct single_pixel_displ *pix;

		pix = gp->curr_pix_displ[i];

		if ( pix->dx == -10000.0 ) break;
		list_dx[count] = pix->dx;
		list_dy[count] = pix->dy;
		count++;
		if ( pix->ne == NULL ) {
			break;
		} else {
			gp->curr_pix_displ[i] = gp->curr_pix_displ[i]->ne;
		}
	}

	if ( count < 1 ) {
		free(list_dx);
		free(list_dy);
		return 0;
	}

	gp->avg_displ_x[i] = comp_median(list_dx, count);
	gp->avg_displ_y[i] = comp_median(list_dy, count);
	gp->avg_displ_abs[i] = modulus2d(gp->avg_displ_x[i],
	                                 gp->avg_displ_y[i]);

	free(list_dx);
	free(list_dy);
	return 0;
}


static int allocate_next_element(struct single_pixel_displ **curr_pix_displ,
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


static int add_distance_to_list(struct gpanel *gp,
				struct imagefeature *imfe,
				Reflection *refl, double fx, double fy)
{
	int pix_index;
	double rfs, rss;
	double crx, cry;

	pix_index = ((int)rint(imfe->fs) + gp->p->w*(int)rint(imfe->ss));

	if ( gp->num_pix_displ[pix_index] > 0 ) {

		int ret;

		ret = allocate_next_element(gp->curr_pix_displ, pix_index);

		if ( ret != 0 ) return ret;

	}

	get_detector_pos(refl, &rfs, &rss);
	compute_x_y(rfs, rss, get_panel(refl), &crx, &cry);
	gp->curr_pix_displ[pix_index]->dx = fx - crx;
	gp->curr_pix_displ[pix_index]->dy = fy - cry;
	gp->curr_pix_displ[pix_index]->ne = NULL;
	gp->num_pix_displ[pix_index]++;

	return 0;
}


static int count_pixels_with_min_peaks(struct gpanel *gp, int min_num_peaks)
{
	int pixel_count = 0;
	int ifs, iss;

	for ( iss=0; iss<gp->p->h; iss++ ) {
	for ( ifs=0; ifs<gp->p->w; ifs++ ) {

		int idx = ifs+gp->p->w*iss;
		if ( gp->num_pix_displ[idx] >= min_num_peaks ) {
			pixel_count += 1;
		}

	}
	}
	return pixel_count;
}


static void adjust_min_peaks_per_conn(struct rg_collection *connected,
                                      struct gpanel *gpanels,
                                      struct detector *det,
                                      struct geoptimiser_params *gparams,
                                      struct connected_data *conn_data)
{
	int min_num_peaks, di, ip;

	STATUS("Adjusting the minimum number of measurements per pixel in "
	       "order to have enough measurements for each connected group.\n");

	for ( di=0; di<connected->n_rigid_groups; di++ ) {

		for ( min_num_peaks=gparams->min_num_peaks_per_pix;
		      min_num_peaks>0; min_num_peaks-- )
		{
			int di_count = 0;

			for ( ip=0; ip<connected->rigid_groups[di]->n_panels;
			      ip++ )
			{
				int pix_count;
				struct panel *p;
				struct gpanel *gp;

				p = connected->rigid_groups[di]->panels[ip];
				gp = &gpanels[panel_number(det, p)];

				pix_count = count_pixels_with_min_peaks(gp,
				                                 min_num_peaks);

				di_count += pix_count;

			}

			conn_data[di].n_peaks_in_conn = di_count;
			if ( di_count >= gparams->min_num_pix_per_conn_group ) {
				conn_data[di].num_peaks_per_pixel = min_num_peaks;
				break;
			}
		}

		STATUS("Minimum number of measurements per pixel for connected "
		       "group %s has been set to %i\n", conn_data[di].name,
		       conn_data[di].num_peaks_per_pixel);
	}
}


static int compute_pixel_displacements(struct image *images, int n_images,
                                       struct gpanel *gpanels,
                                       struct detector *det,
                                       struct rg_collection *connected,
                                       struct geoptimiser_params *gparams,
                                       double clen_to_use,
                                       struct connected_data *conn_data)
{
	int cp;

	STATUS("Computing pixel displacements.\n");

	for ( cp=0; cp<n_images; cp++ ) {

		int fi;
		ImageFeatureList *flist = images[cp].features;

		if ( gparams->only_best_distance ) {
			if ( fabs(images[cp].avg_clen - clen_to_use) > 0.0001 ) {
				continue;
			}
		}

		for ( fi=0; fi<image_feature_count(images[cp].features); fi++ ) {

			double min_dist;
			double fx, fy;
			Reflection *refl;
			struct imagefeature *imfe;

			imfe = image_get_feature(flist, fi);
			if ( imfe == NULL ) continue;

			compute_x_y(imfe->fs, imfe->ss, imfe->p, &fx, &fy);

			/* Find the closest reflection (from all crystals) */
			refl = find_closest_reflection(&images[cp], fx, fy,
			                               &min_dist);
			if ( refl == NULL ) continue;

			if ( min_dist < gparams->max_peak_dist ) {

				struct gpanel *gp;
				int r;
				gp = &gpanels[panel_number(det, imfe->p)];

				r = add_distance_to_list(gp, imfe, refl, fx, fy);
				if ( r ) return r;

			}
		}
	}

	return 0;
}


static int compute_avg_pix_displ(struct gpanel *gp, int idx,
                                 int num_peaks_per_pixel)
{
	int ret;

	if ( gp->num_pix_displ[idx] >= num_peaks_per_pixel ) {

		ret = fill_avg_pixel_displ(gp, idx);
		if ( ret != 0 ) return ret;

	} else {

		gp->avg_displ_x[idx] = -10000.0;
		gp->avg_displ_y[idx] = -10000.0;
		gp->avg_displ_abs[idx] = -10000.0;

	}

	return 0;

}


static int compute_avg_displacements(struct detector *det,
                                     struct rg_collection *connected,
                                     struct connected_data *conn_data,
                                     struct gpanel *gpanels)
{
	int di, ip, ifs, iss;
	int pix_index, ret;
	struct panel *p;

	for ( di=0; di<connected->n_rigid_groups; di++ ) {

		for ( ip=0; ip<connected->rigid_groups[di]->n_panels; ip++ ) {

			int pp;
			struct gpanel *gp;

			p = connected->rigid_groups[di]->panels[ip];
			pp = panel_number(det, p);
			gp = &gpanels[pp];

			for ( iss=0; iss<p->h; iss++ ) {
			for ( ifs=0; ifs<p->w; ifs++ ) {

				pix_index = ifs+p->w*iss;

				ret = compute_avg_pix_displ(gp, pix_index,
				             conn_data[di].num_peaks_per_pixel);

				if ( ret != 0 ) return ret;

			}
			}
		}

	}
	return 0;
}


static double compute_error(struct rg_collection *connected,
                            struct detector *det,
                            struct connected_data *conn_data,
                            struct gpanel *gpanels)
{
	double total_error = 0;
	int num_total_error = 0;
	int di, ip;

	for ( di=0;di<connected->n_rigid_groups;di++ ) {

		struct panel *p;
		double connected_error = 0;
		int num_connected_error = 0;
		int ifs, iss;

		for ( ip=0; ip<connected->rigid_groups[di]->n_panels; ip++ ) {

			p = connected->rigid_groups[di]->panels[ip];

			for ( iss=0; iss<p->h; iss++ ) {
			for ( ifs=0; ifs<p->w; ifs++ ) {

				int pix_index;
				struct gpanel *gp;
				int pp = panel_number(det, p);

				gp = &gpanels[pp];
				pix_index = ifs+p->w*iss;

				if ( gp->num_pix_displ[pix_index]
				       >= conn_data[di].num_peaks_per_pixel )
				{
					double cer;

					cer = gp->avg_displ_abs[pix_index]
					    * gp->avg_displ_abs[pix_index];
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
			       "more than %d peaks: RMSD = %0.4f pixels.\n",
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


static int compute_rot_stretch_for_empty_panels(struct rg_collection *quads,
                                                struct rg_collection *conn,
                                                int min_pix,
                                                struct connected_data *conn_data)
{
	int di,i;
	double *aver_ang;
	double *aver_str;
	int *n;

	STATUS("Computing rotation and elongation corrections for groups "
	       "without the required number of measurements.\n");

	aver_ang = malloc(quads->n_rigid_groups*sizeof(double));
	aver_str = malloc(quads->n_rigid_groups*sizeof(double));
	n = malloc(quads->n_rigid_groups*sizeof(double));

	for ( i=0; i<quads->n_rigid_groups; i++ ) {
		aver_ang[i] = 0.0;
		aver_str[i] = 0.0;
		n[i] = 0;
	}

	/* Calculate the mean values for the groups which DO have
	 * enough measurements */
	for ( di=0; di<conn->n_rigid_groups; di++ ) {
		if ( conn_data[di].n_peaks_in_conn >= min_pix ) {
			aver_ang[conn_data[di].num_quad] += conn_data[di].cang;
			aver_str[conn_data[di].num_quad] += conn_data[di].cstr;
			n[conn_data[di].num_quad]++;
		}
	}

	/* Divide totals to get means */
	for ( i=0; i<quads->n_rigid_groups; i++ ) {
		aver_ang[i] /= (double)n[i];
		aver_str[i] /= (double)n[i];
	}

	for ( di=0; di<conn->n_rigid_groups; di++ ) {

		int qn = conn_data[di].num_quad;
		assert(qn < quads->n_rigid_groups);

		if ( conn_data[di].n_peaks_in_conn >= min_pix ) continue;

		if ( n[qn] > 0 ) {

			conn_data[di].cang = aver_ang[qn];
			conn_data[di].cstr = aver_str[qn];
			STATUS("Connected group %s has only %i useful "
			       "pixels. Using average angle: %0.4f "
			       "degrees\n", conn_data[di].name,
			       conn_data[di].n_peaks_in_conn,
			       rad2deg(conn_data[di].cang));

		} else {

			STATUS("Connected group %s does not have enough "
			       " peaks (%i). It will not be moved.\n",
			       conn_data[di].name,
			       conn_data[di].n_peaks_in_conn);
		}
	}

	free(aver_ang);
	free(aver_str);
	free(n);

	return 0;
}


static void correct_rotation_and_stretch(struct rg_collection *connected,
                                         struct detector *det,
                                         struct gpanel *gpanels,
                                         struct connected_data *conn_data,
                                         double clen_to_use,
                                         double stretch_coeff,
                                         struct geoptimiser_params *gparams)
{

	int di, ip;

	STATUS("Applying rotation and stretch corrections.\n");

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

	for ( di=0; di<connected->n_rigid_groups; di++ ) {

		int npp = conn_data[di].num_peaks_per_pixel;
		double cs = conn_data[di].cstr;

		if ( fabs(cs)<FLT_EPSILON ) cs = stretch_coeff;

		for ( ip=0; ip<connected->rigid_groups[di]->n_panels; ip++ ) {

			struct panel *p;
			double new_fsx, new_fsy, new_ssx, new_ssy;
			double new_cnx, new_cny;
			int fs, ss;
			struct gpanel *gp;

			p = connected->rigid_groups[di]->panels[ip];
			gp = &gpanels[panel_number(det, p)];

			new_fsx = p->fsx*cos(conn_data[di].cang)-
			          p->fsy*sin(conn_data[di].cang);
			new_fsy = p->fsx*sin(conn_data[di].cang)+
			          p->fsy*cos(conn_data[di].cang);
			new_ssx = p->ssx*cos(conn_data[di].cang)-
			          p->ssy*sin(conn_data[di].cang);
			new_ssy = p->ssx*sin(conn_data[di].cang)+
			          p->ssy*cos(conn_data[di].cang);

			/* Calculate corner adjustment
			 * NB All panels follow the first one */
			if ( ip == 0 ) {

				new_cnx = p->cnx * cs;
				new_cny = p->cny * cs;

			} else {

				struct panel *p0;
				double delta_x, delta_y;

				p0 = connected->rigid_groups[di]->panels[0];

				delta_x = p->cnx - p0->cnx / cs;
				delta_y = p->cny - p0->cny / cs;

				new_cnx = p0->cnx + delta_x;
				new_cny = p0->cny + delta_y;

			}

			/* The average displacements now need to be updated */
			for ( ss=0; ss<p->h; ss++ ) {
			for ( fs=0; fs<p->w; fs++ ) {

				int i = fs + p->w*ss;

				if ( gp->num_pix_displ[i] < npp ) continue;

				gp->avg_displ_x[i] += fs*(new_fsx/cs - p->fsx)
				                    + ss*(new_ssx/cs - p->ssx)
				                    + new_cnx/cs - p->cnx;

				gp->avg_displ_y[i] += fs*(new_fsy/cs - p->fsy)
				                    + ss*(new_ssy/cs - p->ssy)
				                    + new_cny/cs - p->cny;
			}
			}

			p->fsx = new_fsx;
			p->fsy = new_fsy;
			p->ssx = new_ssx;
			p->ssy = new_ssy;

			p->cnx = new_cnx;
			p->cny = new_cny;

		}
	}

}


/* Collect together all the offsets for each group in "connected"
 * Only offsets which have enough peaks per pixel will be used. */
static int collate_offsets_for_rg(struct rigid_group *group,
                                  struct detector *det,
                                  struct gpanel *gpanels,
                                  int num_peaks_per_pixel,
                                  double *list_dx, double *list_dy, int list_sz)
{
	int ip;
	int counter = 0;

	for ( ip=0; ip<group->n_panels; ip++ ) {

		int ifs, iss;
		struct panel *p = group->panels[ip];
		struct gpanel *gp = &gpanels[panel_number(det, p)];

		for ( iss=0; iss<p->h; iss++ ) {
		for ( ifs=0; ifs<p->w; ifs++ ) {

			int idx = ifs+p->w*iss;

			if ( gp->num_pix_displ[idx] < num_peaks_per_pixel ) {
				continue;
			}

			if ( list_sz == counter ) {
				ERROR("Too many offsets in group %s!\n",
				      group->name);
				return 0;
			}

			list_dx[counter] = gp->avg_displ_x[idx];
			list_dy[counter] = gp->avg_displ_y[idx];
			counter++;

		}
		}
	}

	return counter;
}


static void fill_conn_data_sh(struct connected_data *conn_data,
                              double *list_dx, double *list_dy,
                              int di, double mpd)
{
	conn_data[di].sh_x = comp_median(list_dx,
	                                 conn_data[di].n_peaks_in_conn);
	conn_data[di].sh_y = comp_median(list_dy,
	                                 conn_data[di].n_peaks_in_conn);

	STATUS("Group %s, num pixels: %i, shifts x,y: %0.8f, %0.8f px\n",
	       conn_data[di].name, conn_data[di].n_peaks_in_conn,
	       conn_data[di].sh_x, conn_data[di].sh_y);

	if ( modulus2d(conn_data[di].sh_x, conn_data[di].sh_y ) > 0.8*mpd ) {
		STATUS("WARNING: absolute shift is: %0.1f > 0.8*%0.1f px "
		       "Increase the value of max_peak_distance!\n",
		       modulus2d(conn_data[di].sh_x, conn_data[di].sh_y), mpd);
	}
}


static int compute_shift(struct rg_collection *connected,
                         struct connected_data *conn_data,
                         struct detector *det,
                         struct geoptimiser_params *gparams,
                         struct gpanel *gpanels)
{
	int di;

        STATUS("Computing shift corrections.\n");

	for ( di=0; di<connected->n_rigid_groups; di++ ) {

		double *list_dx;
		double *list_dy;
		int ct;

		list_dx = malloc(conn_data[di].n_peaks_in_conn*sizeof(double));
		list_dy = malloc(conn_data[di].n_peaks_in_conn*sizeof(double));
		if  ( (list_dx == NULL) || (list_dy == NULL) ) {
			ERROR("Failed to allocate memory for computing shifts\n");
			free(list_dx);
			free(list_dy);
			return 1;
		}

		ct = collate_offsets_for_rg(connected->rigid_groups[di], det,
		                            gpanels,
		                            conn_data[di].num_peaks_per_pixel,
		                            list_dx, list_dy,
		                            conn_data[di].n_peaks_in_conn);

		if ( ct != conn_data[di].n_peaks_in_conn ) {
			ERROR("Wrong number of peaks for group %s!\n",
			      connected->rigid_groups[di]->name);
			ERROR("Counter: %i n_peaks_in_conn: %i\n",
			       ct, conn_data[di].n_peaks_in_conn);
			abort();
		}

		if ( conn_data[di].n_peaks_in_conn
		          >= gparams->min_num_pix_per_conn_group )
		{
			fill_conn_data_sh(conn_data, list_dx, list_dy, di,
			                  gparams->max_peak_dist);

		} else {
			conn_data[di].sh_x = -10000.0;
			conn_data[di].sh_y = -10000.0;
		}

		free(list_dx);
		free(list_dy);
	}

	return 0;
}


static int compute_shift_for_empty_panels(struct rg_collection *quadrants,
                                            struct rg_collection *connected,
                                            struct connected_data *conn_data,
                                            int min_pix)
{
	int di, i;
	double *aver_x;
	double *aver_y;
	int *n;

	aver_x = malloc(quadrants->n_rigid_groups * sizeof(double));
	aver_y = malloc(quadrants->n_rigid_groups * sizeof(double));
	n = malloc(quadrants->n_rigid_groups * sizeof(int));

	for ( i=0; i<quadrants->n_rigid_groups; i++ ) {
		aver_x[i] = 0.0;
		aver_y[i] = 0.0;
		n[i] = 0;
	}

	STATUS("Computing shift corrections for groups without the required "
	       "number of measurements.\n");

	for ( di=0; di<connected->n_rigid_groups; di++ ) {
		if ( conn_data[di].n_peaks_in_conn >= min_pix ) {
			aver_x[conn_data[di].num_quad] += conn_data[di].sh_x;
			aver_y[conn_data[di].num_quad] += conn_data[di].sh_y;
			n[conn_data[di].num_quad]++;
		}
	}

	for ( i=0; i<quadrants->n_rigid_groups; i++ ) {
		aver_x[i] /= (double)n[i];
		aver_y[i] /= (double)n[i];
	}

	for ( di=0; di<connected->n_rigid_groups; di++ ) {
		if ( conn_data[di].n_peaks_in_conn >= min_pix ) continue;

		int qn = conn_data[di].num_quad;

		if ( n[qn] > 0 ) {

			conn_data[di].sh_x = aver_x[qn];
			conn_data[di].sh_y = aver_y[qn];
			STATUS("Panel %s doesn't not have enough (%i) "
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


	free(aver_x);
	free(aver_y);
	free(n);

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


static void scan_p1(int ip0, int ip1, struct rg_collection *connected,
                    struct connected_data *conn_data,
                    struct detector *det, struct gpanel *gpanels,
                    int di, double min_dist,
                    long *num_ang, int ifs0, int iss0,
                    double c_x0, double c_y0, double cd_x0, double cd_y0,
                    int compute, double *angles, double *stretches)
{

	int iss1, ifs1;
	struct panel *p1;
	struct gpanel *gp1;

	p1 = connected->rigid_groups[di]->panels[ip1];
	gp1 = &gpanels[panel_number(det, p1)];

	int min_ss_p1, min_fs_p1;

	if ( ip0 == ip1 ) {
		min_fs_p1 = ifs0;
		min_ss_p1 = iss0;
	} else {
		min_fs_p1 = 0;
		min_ss_p1 = 0;
	}

	for ( iss1=min_ss_p1; iss1<p1->h; iss1++ ) {
	for ( ifs1=min_fs_p1; ifs1<p1->w; ifs1++ ) {

		double dist;
		double c_x1, c_y1, cd_x1, cd_y1;
		double d_c_x, d_c_y, d_cd_x, d_cd_y;
		double len1, len2;
		int pix_index1 = ifs1+p1->w*iss1;

		if ( gp1->num_pix_displ[pix_index1]
		           < conn_data[di].num_peaks_per_pixel ) continue;

		compute_x_y(ifs1, iss1, p1, &c_x1, &c_y1);
		cd_x1 = c_x1 - gp1->avg_displ_x[pix_index1];
		cd_y1 = c_y1 - gp1->avg_displ_y[pix_index1];
		d_c_x = c_x1-c_x0;
		d_c_y = c_y1-c_y0;
		d_cd_x = cd_x1-cd_x0;
		d_cd_y = cd_y1-cd_y0;

		dist = modulus2d(d_c_x,d_c_y);
		if ( dist < min_dist ) continue;

		len1 = modulus2d(d_c_x, d_c_y);
		len2 = modulus2d(d_cd_x, d_cd_y);
		if ( len1<FLT_EPSILON || len2<FLT_EPSILON ) continue;

		if ( compute ) {

			double scal_m;
			double multlen;

			scal_m = d_c_x * d_cd_x+ d_c_y * d_cd_y - FLT_EPSILON;

			multlen = len1*len2;
			if ( fabs(scal_m)>=multlen ) {
				angles[*num_ang] = 0.0;
			} else {

				angles[*num_ang] = acos(scal_m/multlen);

				if (d_c_y * d_cd_x - d_c_x * d_cd_y < 0) {
					angles[*num_ang] *= -1.;
				}

			}

			stretches[*num_ang] = len1/len2;
		}

		*num_ang = *num_ang+1;
	}
	}
}


/* Executed for each panel in the connected group */
static void scan_p0(int ip0, struct rg_collection *connected,
                    struct connected_data *conn_data,
                    struct detector *det, struct gpanel *gpanels,
                    int di, double min_dist,
                    long *num_ang, int compute,
                    double *angles, double *stretches)
{
	int iss0, ifs0, ip1;
	struct gpanel *gp;
	struct panel *p0;

	p0 = connected->rigid_groups[di]->panels[ip0];
	gp = &gpanels[panel_number(det, p0)];

	for ( iss0=0; iss0<p0->h; iss0++ ) {
	for ( ifs0=0; ifs0<p0->w; ifs0++ ) {

		double c_x0, c_y0, cd_x0, cd_y0;
		int pix_index0 = ifs0+p0->w*iss0;

		if ( gp->num_pix_displ[pix_index0]
		                < conn_data[di].num_peaks_per_pixel ) continue;

		compute_x_y(ifs0, iss0, p0, &c_x0, &c_y0);
		cd_x0 = c_x0 - gp->avg_displ_x[pix_index0];
		cd_y0 = c_y0 - gp->avg_displ_y[pix_index0];

		for ( ip1=ip0; ip1<connected->rigid_groups[di]->n_panels;
		      ip1++ )
		{
			scan_p1(ip0, ip1, connected,
			        conn_data, det, gpanels, di, min_dist,
			        num_ang, ifs0, iss0, c_x0,
			        c_y0, cd_x0, cd_y0, compute,
			        angles, stretches);
		}

	}
	}
}


static double compute_rotation_and_stretch(struct rg_collection *connected,
                                           struct connected_data *conn_data,
                                           struct detector *det,
                                           struct gpanel *gpanels,
                                           double dist_coeff_for_rot_str,
                                           struct geoptimiser_params *gparams)
{
	int di;
	double stretch_cf;
	long num_coeff;
	long *num_angles;
	double *stretch_coeff;

	STATUS("Computing rotation and stretch corrections.\n");

	stretch_coeff = malloc(connected->n_rigid_groups*sizeof(double));
	num_angles = malloc(connected->n_rigid_groups*sizeof(long));
	num_coeff = 0;

	for ( di=0; di<connected->n_rigid_groups; di++ ) {

		long max_num_ang = 0;

		double min_dist;
		double *angles;
		double *stretches;

		struct panel *first_p;
		long num_ang = 0;
		int ip0;
		int num_pix_first_p;

		/* Enough peaks in this connected group? */
		if ( conn_data[di].n_peaks_in_conn
		     < gparams->min_num_pix_per_conn_group ) continue;

		first_p = connected->rigid_groups[di]->panels[0];

		num_pix_first_p = first_p->w * first_p->h;

		/* FIXME: minrad here is not universal */
		min_dist = dist_coeff_for_rot_str * sqrt(num_pix_first_p
		                       * connected->rigid_groups[di]->n_panels);

		for ( ip0=0; ip0<connected->rigid_groups[di]->n_panels; ip0++ ) {
			scan_p0(ip0, connected, conn_data, det, gpanels,
			        di, min_dist, &num_ang, 0, NULL, NULL);
		}

		max_num_ang = num_ang+1;

		angles = malloc(max_num_ang*sizeof(double));
		if ( angles == NULL ) {
			ERROR("Error allocating memory for angle "
			      "optimization\n");
			return -1.0;
		}
		stretches = malloc(max_num_ang*sizeof(double));
		if ( stretches == NULL ) {
			ERROR("Error allocating memory for stretch "
			      "optimization\n");
			return -1.0;
		}

		num_ang = 0;

		for ( ip0=0; ip0<connected->rigid_groups[di]->n_panels; ip0++ ) {
			scan_p0(ip0, connected, conn_data, det, gpanels,
			        di, min_dist, &num_ang, 1, angles, stretches);
		}

		if ( num_ang < 1 ) continue;

		conn_data[di].cang = -comp_median(angles, num_ang);
		conn_data[di].cstr = comp_median(stretches, num_ang);

		STATUS("Panel %s, num: %li, angle: %0.4f deg, stretch coeff: "
		       "%0.4f\n", conn_data[di].name, num_ang,
		                  rad2deg(conn_data[di].cang),
		                  conn_data[di].cstr);

		stretch_coeff[num_coeff] = conn_data[di].cstr;
		num_angles[num_coeff] = num_ang;
		num_coeff++;

		free(angles);
		free(stretches);
	}

	stretch_cf = 1.0;

	printf("Computing overall stretch coefficient.\n");

	if ( num_coeff>0 ) {

		int peaks_per_p;

		peaks_per_p = gparams->min_num_peaks_per_pix;

		while ( peaks_per_p>=0 ) {

			double total_num;
			long di;

			stretch_cf = 0;
			total_num = 0;
			for ( di=0; di<num_coeff; di++ ) {
				if ( conn_data[di].num_peaks_per_pixel >=
				     peaks_per_p ) {
					stretch_cf += stretch_coeff[di]*
					         (double)num_angles[di];
					total_num += num_angles[di];
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

	free(stretch_coeff);
	free(num_angles);

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


static int save_data_to_png(char *filename, struct detector *det,
                            struct gpanel *gpanels)
{
	struct image im;
	int i;
	struct rectangle rect;
	GdkPixbuf *col_scale;
	cairo_t *cr;

	cairo_status_t r;
	cairo_surface_t *surf;

	im.det = det;
	im.bad = NULL;
	im.dp = malloc(det->n_panels*sizeof(float *));
	if ( im.dp == NULL ) {
		ERROR("Failed to allocate data\n");
		return 1;
	}
	for ( i=0; i<det->n_panels; i++ ) {

		int fs, ss;
		struct panel *p;

		p = &det->panels[i];

		im.dp[i] = calloc(p->w * p->h, sizeof(float));
		if ( im.dp[i] == NULL ) {
			ERROR("Failed to allocate data\n");
			return 1;
		}

		for ( ss=0; ss<p->h; ss++ ) {
		for ( fs=0; fs<p->w; fs++ ) {

			int idx;
			float val;

			idx = fs + ss*p->w;

			if ( gpanels[i].avg_displ_abs[idx] == -10000.0) {
				val = 0.0;
			} else if ( gpanels[i].avg_displ_abs[idx] > 1.0) {
				val = 1.0;
			} else {
				val = (float)gpanels[i].avg_displ_abs[idx];
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

	for ( i=0; i<det->n_panels; i++ ) {
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


static int initialize_pixel_displacement_list(struct gpanel *gp)
{
	int ipx;

	gp->pix_displ_list = calloc(gp->p->w*gp->p->h,
	                            sizeof(struct single_pixel_displ));
	if ( gp->pix_displ_list == NULL ) {
		ERROR("Error allocating memory for pixel displacement data.\n");
		return 1;
	}

	gp->curr_pix_displ = calloc(gp->p->w*gp->p->h,
	                            sizeof(struct single_pixel_displ *));
	if ( gp->curr_pix_displ == NULL ) {
		ERROR("Error allocating memory for pixel displacement data.\n");
		free(gp->pix_displ_list);
		return 1;
	}
	gp->num_pix_displ = calloc(gp->p->w*gp->p->h, sizeof(int));
	if ( gp->num_pix_displ == NULL ) {
		ERROR("Error allocating memory for pixel displacement data.\n");
		free(gp->pix_displ_list);
		free(gp->curr_pix_displ);
		return 1;
	}

	for ( ipx=0; ipx<gp->p->w*gp->p->h; ipx++ ) {
		gp->pix_displ_list[ipx].dx = -10000.0;
		gp->pix_displ_list[ipx].dy = -10000.0;
		gp->pix_displ_list[ipx].ne = NULL;
		gp->curr_pix_displ[ipx] = &gp->pix_displ_list[ipx];
		gp->num_pix_displ[ipx] = 0;
	}

	return 0;
}


static void free_displ_lists(struct gpanel *gpanels, int n)
{
	int j;
	struct single_pixel_displ *curr = NULL;
	struct single_pixel_displ *next = NULL;

	for ( j=0; j<n; j++ ) {

		int i;
		struct gpanel *gp = &gpanels[j];

		for ( i=0; i<gp->p->w*gp->p->h; i++ ) {

			curr = &gp->pix_displ_list[i];

			if ( curr->ne != NULL ) {
				curr = curr->ne;
				while ( curr != NULL ) {
					next = curr->ne;
					free(curr);
					curr = next;
				}
			}
		}

		free(gp->curr_pix_displ);
		free(gp->pix_displ_list);

	}
}


static void recompute_panel_avg_displ(struct rg_collection *connected,
                                      struct connected_data *conn_data,
                                      struct panel *p, struct gpanel *gp,
                                      int num_peaks_per_pixel,
                                      double sh_x, double sh_y)
{
	int ifs, iss;

	for ( iss=0; iss<p->h; iss++ ) {
	for ( ifs=0; ifs<p->w; ifs++ ) {

		int pix_index = ifs+p->w*iss;

		if ( gp->num_pix_displ[pix_index] >= num_peaks_per_pixel ) {

			gp->avg_displ_x[pix_index] -= sh_x;
			gp->avg_displ_y[pix_index] -= sh_y;
			gp->avg_displ_abs[pix_index] = modulus2d(
			                            gp->avg_displ_x[pix_index],
			                            gp->avg_displ_y[pix_index]);

		} else {
			gp->avg_displ_abs[pix_index] = -10000.0;
		}

	}
	}
}


void recompute_avg_displ(struct rg_collection *connected,
                         struct connected_data *conn_data,
                         struct detector *det, struct gpanel *gpanels)
{

	int di, ip;

	for ( di=0;di<connected->n_rigid_groups;di++ ) {

		if (conn_data[di].sh_x < -9999.0) continue;

		for (ip=0; ip<connected->rigid_groups[di]->n_panels; ip++) {

			struct panel *p;
			struct gpanel *gp;

			p = connected->rigid_groups[di]->panels[ip];
			gp = &gpanels[panel_number(det, p)];
			recompute_panel_avg_displ(connected, conn_data, p, gp,
			                     conn_data[di].num_peaks_per_pixel,
			                     conn_data[di].sh_x,
			                     conn_data[di].sh_y);

		}
	}
}


int optimize_geometry(struct geoptimiser_params *gparams,
                      struct detector *det,
                      struct rg_collection *quadrants,
                      struct rg_collection *connected)
{
	int pi;
	int ret;
	int write_ret;

	double res_sum;
	double avg_res;
	double clen_to_use;

	// for angles and stretch calculation use
	// only pixels which are distco*size_panel
	// away
	double dist_coeff_for_rot_str = 0.2;
	double total_error;

	double stretch_coeff = 1.0;

	struct connected_data *conn_data = NULL;
	struct image *images;
	int n_images = 0;
	UnitCell *avg_cell;
	struct gpanel *gpanels;

	STATUS("Maximum distance between peaks: %0.1f pixels.\n",
	       gparams->max_peak_dist);

	STATUS("Minimum number of measurements for a pixel to be included in "
	       "the refinement: %i\n", gparams->min_num_peaks_per_pix);
	STATUS("Minimum number of measurements for connected group for "
	       "accurate estimation of position/orientation: %i\n",
	       gparams->min_num_pix_per_conn_group);

	if ( (det->n_panels == 64) && !gparams->no_cspad ) {

		int num_errors = 0;

		STATUS("It looks like the detector is a CSPAD. "
		       "Checking relative distance and orientation of "
		       "connected ASICS.\n");
		STATUS("If the detector is not a CSPAD, please rerun the "
		       "program with the --no-cspad option.\n");

		STATUS("Enforcing CSPAD layout...\n");
		num_errors = check_and_enforce_cspad_dist(gparams, det);

		if ( gparams->enforce_cspad_layout ) {

			int geom_wr;

			STATUS("Saving geometry with enforced CSPAD layout.\n"
			       "Please restart geometry optimization using the "
			       "optimized geometry from this run as input "
			       "geometry file.\n");
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

	images = read_patterns_from_stream(gparams->infile, det, &n_images);
	if ( (n_images < 1) || (images == NULL) ) {
		ERROR("Error reading stream file\n");
		return 1;
	}

	avg_cell = compute_avg_cell_parameters(images, n_images);
	if ( avg_cell == NULL ) {
		free(images);
		return 1;
	}

	res_sum = 0;
	for ( pi=0; pi<det->n_panels; pi++ ) {
		res_sum += det->panels[pi].res;
	}
	avg_res = res_sum/det->n_panels;

	clen_to_use = pick_clen_to_use(gparams, images, n_images, avg_res,
	                               avg_cell);
	if ( clen_to_use < 0.0 ) return 1;

	gpanels = calloc(det->n_panels, sizeof(struct gpanel));
	if ( gpanels == NULL ) {
		ERROR("Failed to allocate panels.\n");
		return 1;
	}
	for ( pi=0; pi<det->n_panels; pi++ ) {

		struct gpanel *gp = &gpanels[pi];
		int npx = det->panels[pi].w * det->panels[pi].h;

		gp->p = &det->panels[pi];

		gp->avg_displ_x = malloc(npx*sizeof(double));
		gp->avg_displ_y = malloc(npx*sizeof(double));
		gp->avg_displ_abs = malloc(npx*sizeof(double));

		if ( (gp->avg_displ_x == NULL) || (gp->avg_displ_y == NULL)
		  || (gp->avg_displ_abs == NULL) )
		{
			ERROR("Failed to allocate displacements.\n");
			return 1;
		}

		if ( initialize_pixel_displacement_list(gp) ) {
			ERROR("Failed to allocate lists.\n");
			return 1;
		}

	}

	conn_data = initialize_conn_data(quadrants, connected);
	if ( conn_data == NULL ) {
		ERROR("Error allocating memory for connected structure data.\n");
		return 1;
	}

	if ( compute_pixel_displacements(images, n_images, gpanels, det,
	                                 connected, gparams, clen_to_use,
	                                 conn_data) ) return 1;

	adjust_min_peaks_per_conn(connected, gpanels, det, gparams, conn_data);

	if ( compute_avg_displacements(det, connected, conn_data, gpanels) ) {
		free(conn_data);
		free(images);
		return 1;
	}

	free_displ_lists(gpanels, det->n_panels);
	/* gpanels[].num_pix_displ is still there */

	STATUS("Computing error before correction.\n");
	total_error = compute_error(connected, det, conn_data, gpanels);

	STATUS("Detector-wide error before correction: RMSD = %0.4f pixels.\n",
	       total_error);

	if ( gparams->error_maps ) {

		STATUS("Saving error map before correction.\n");

#ifdef HAVE_SAVE_TO_PNG

		if ( save_data_to_png("error_map_before.png", det, gpanels) ) {
			ERROR("Error while writing data to file.\n");
			free(conn_data);
			free(images);
			return 1;
		}

#else /* HAVE_SAVE_TO_PNG */

		ERROR("WARNING: geoptimiser was compiled without GTK and cairo "
		       "support. Error maps will not be saved.\n");

#endif /* HAVE_SAVE_TO_PNG */

	}

	stretch_coeff = compute_rotation_and_stretch(connected, conn_data,
	                                             det, gpanels,
                                                     dist_coeff_for_rot_str,
                                                     gparams);
	if ( stretch_coeff < 0.0 ) {
		free(conn_data);
		return 1;
	}

	ret = compute_rot_stretch_for_empty_panels(quadrants, connected,
	                        gparams->min_num_pix_per_conn_group, conn_data);
	if ( ret ) {
		free(conn_data);
		return 1;
	}

	correct_rotation_and_stretch(connected, det, gpanels, conn_data,
	                             clen_to_use, stretch_coeff,
	                             gparams);

	ret = compute_shift(connected, conn_data, det, gparams, gpanels);
	if ( ret != 0 ) {
		free(conn_data);
		return 1;
	}

	compute_shift_for_empty_panels(quadrants, connected, conn_data,
	                               gparams->min_num_pix_per_conn_group);

	correct_shift(connected, conn_data, clen_to_use);

	recompute_avg_displ(connected, conn_data, det, gpanels);

	if ( gparams->error_maps ) {

#ifdef HAVE_SAVE_TO_PNG

		STATUS("Saving error map after correction.\n");

		ret = save_data_to_png("error_map_after.png", det, gpanels);
		if ( ret ) {
			ERROR("Error while writing data to file.\n");
			free(conn_data);
			return 1;
		}

#else /* HAVE_SAVE_TO_PNG */

		STATUS("ERROR: geoptimiser was compiled without GTK and cairo "
		       "support.\n Error maps will not be saved.\n");

#endif /* HAVE_SAVE_TO_PNG */

	}

	STATUS("Computing errors after correction.\n");
	total_error = compute_error(connected, det, conn_data, gpanels);

	STATUS("Detector-wide error after correction: RMSD = %0.4f pixels.\n",
	       total_error);

	write_ret = write_detector_geometry_2(gparams->geometry_filename,
	                                      gparams->outfile, det,
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

	free(conn_data);
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
	gparams->max_peak_dist = 0.0;

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

#ifdef HAVE_SAVE_TO_PNG
	g_type_init();
#endif

	ret_val = optimize_geometry(gparams, det, quadrants, connected);

	return ret_val;
}
