/*
 * peaks.c
 *
 * Peak search and other image analysis
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_statistics_int.h>
#include <pthread.h>

#include "image.h"
#include "utils.h"
#include "index.h"
#include "peaks.h"
#include "detector.h"
#include "filters.h"
#include "diffraction.h"


#define MAX_HITS (1024)

/* How close a peak must be to an indexed position to be considered "close"
 * for the purposes of double hit detection and sanity checking. */
#define PEAK_CLOSE (30.0)

/* How close a peak must be to an indexed position to be considered "close"
 * for the purposes of integration. */
#define PEAK_REALLY_CLOSE (10.0)

/* Degree of polarisation of X-ray beam */
#define POL (1.0)


#define PEAK_WINDOW_SIZE (10)
#define MAX_PEAKS (2048)
#define INTEGRATION_RADIUS (10)


static int is_hot_pixel(struct image *image, int x, int y)
{
	int dx, dy;
	int w, v;

	w = image->width;
	v = (1*image->data[x+w*y])/2;

	if ( x+1 >= image->width ) return 0;
	if ( x-1 < 0 ) return 0;
	if ( y+1 >= image->height ) return 0;
	if ( y-1 < 0 ) return 0;

	/* Must be at least one adjacent bright pixel */
	for ( dx=-1; dx<=+1; dx++ ) {
	for ( dy=-1; dy<=+1; dy++ ) {
		if ( (dx==0) && (dy==0) ) continue;
		if ( image->data[(x+dx)+w*(y+dy)] >= v ) return 0;
	}
	}

	return 1;
}


static int cull_peaks_in_panel(struct image *image, struct panel *p)
{
	int i, n;
	int nelim = 0;

	n = image_feature_count(image->features);

	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		int j, ncol;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		if ( f->x < p->min_x ) continue;
		if ( f->x > p->max_x ) continue;
		if ( f->y < p->min_y ) continue;
		if ( f->y > p->max_y ) continue;

		/* How many peaks are in the same column? */
		ncol = 0;
		for ( j=0; j<n; j++ ) {

			struct imagefeature *g;

			if ( i==j ) continue;

			g = image_get_feature(image->features, j);
			if ( g == NULL ) continue;

			if ( p->badrow == 'x' ) {
				if ( fabs(f->y - g->y) < 2.0 ) ncol++;
			} else if ( p->badrow == 'y' ) {
				if ( fabs(f->x - g->x) < 2.0 ) ncol++;
			} else {
				ERROR("Invalid badrow direction.\n");
				abort();
			}

		}

		/* More than three? */
		if ( ncol <= 3 ) continue;

		/* Yes?  Delete them all... */
		nelim = 0;
		for ( j=0; j<n; j++ ) {
			struct imagefeature *g;
			g = image_get_feature(image->features, j);
			if ( g == NULL ) continue;
			if ( p->badrow == 'x' ) {
				if ( fabs(f->y - g->y) < 2.0 ) {
					image_remove_feature(image->features,
					                     j);
					nelim++;
				}
			} else if ( p->badrow == 'y' ) {
				if ( fabs(f->x - g->x) < 2.0 ) {
					image_remove_feature(image->features,
					                     j);
					nelim++;
				}
			} else {
				ERROR("Invalid badrow direction.\n");
				abort();
			}

		}

	}

	return nelim;
}


/* Post-processing of the peak list to remove noise */
static int cull_peaks(struct image *image)
{
	int nelim = 0;
	struct panel *p;
	int i;

	for ( i=0; i<image->det->n_panels; i++ ) {
		p = &image->det->panels[i];
		nelim += cull_peaks_in_panel(image, p);
	}

	return nelim;
}


/* Returns non-zero if peak has been vetoed.
 * i.e. don't use result if return value is not zero. */
static int integrate_peak(struct image *image, int xp, int yp,
                           float *xc, float *yc, float *intensity,
                           int do_polar, int do_sa)
{
	signed int x, y;
	const int lim = INTEGRATION_RADIUS * INTEGRATION_RADIUS;
	double total = 0;
	int xct = 0;
	int yct = 0;

	for ( x=-INTEGRATION_RADIUS; x<+INTEGRATION_RADIUS; x++ ) {
	for ( y=-INTEGRATION_RADIUS; y<+INTEGRATION_RADIUS; y++ ) {

		struct panel *p = NULL;
		double val, sa, pix_area, Lsq, dsq, proj_area;
		float tt = 0.0;
		double phi, pa, pb, pol;
		uint16_t flags;

		/* Circular mask */
		if ( x*x + y*y > lim ) continue;

		if ( ((x+xp)>=image->width) || ((x+xp)<0) ) continue;
		if ( ((y+yp)>=image->height) || ((y+yp)<0) ) continue;

		/* Veto this peak if we tried to integrate in a bad region */
		if ( image->flags != NULL ) {
			flags = image->flags[(x+xp)+image->width*(y+yp)];
			if ( !((flags & 0x01) && (flags & 0x04)) ) return 1;
		}

		val = image->data[(x+xp)+image->width*(y+yp)];

		if ( do_sa || do_polar ) {

			p = find_panel(image->det, x+xp, y+yp);
			if ( p == NULL ) return 1;

		}

		if ( do_sa ) {

			/* Area of one pixel */
			pix_area = pow(1.0/p->res, 2.0);
			Lsq = pow(p->clen, 2.0);

			/* Area of pixel as seen from crystal (approximate) */
			tt = get_tt(image, x+xp, y+yp);
			proj_area = pix_area * cos(tt);

			/* Calculate distance from crystal to pixel */
			dsq = pow(((double)(x+xp) - p->cx) / p->res, 2.0);
			dsq += pow(((double)(y+yp) - p->cy) / p->res, 2.0);

			/* Projected area of pixel divided by distance squared */
			sa = 1.0e7 * proj_area / (dsq + Lsq);

			val /= sa;

		}

		if ( do_polar ) {

			if ( !do_sa ) tt = get_tt(image, x+xp, y+yp);

			phi = atan2(y+yp, x+xp);
			pa = pow(sin(phi)*sin(tt), 2.0);
			pb = pow(cos(tt), 2.0);
			pol = 1.0 - 2.0*POL*(1-pa) + POL*(1.0+pb);

			val /= pol;

		}

		total += val;

		xct += val*(xp+x);
		yct += val*(yp+y);

	}
	}

	/* The centroid is excitingly undefined if there is no intensity */
	if ( total != 0 ) {
		*xc = (float)xct / total;
		*yc = (float)yct / total;
		*intensity = total;
	} else {
		*xc = (float)xp;
		*yc = (float)yp;
		*intensity = 0;
	}

	return 0;
}


void search_peaks(struct image *image)
{
	int x, y, width, height;
	float *data;
	double d;
	int idx;
	float fx = 0.0;
	float fy = 0.0;
	float intensity = 0.0;
	int nrej_dis = 0;
	int nrej_hot = 0;
	int nrej_pro = 0;
	int nrej_fra = 0;
	int nrej_bad = 0;
	int nacc = 0;
	int ncull;

	data = image->data;
	width = image->width;
	height = image->height;

	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();

	for ( x=1; x<image->width-1; x++ ) {
	for ( y=1; y<image->height-1; y++ ) {

		double dx1, dx2, dy1, dy2;
		double dxs, dys;
		double grad;
		int mask_x, mask_y;
		int sx, sy;
		double max;
		unsigned int did_something;
		int r;

		/* Overall threshold */
		if ( data[x+width*y] < 50 ) continue;

		/* Get gradients */
		dx1 = data[x+width*y] - data[(x+1)+width*y];
		dx2 = data[(x-1)+width*y] - data[x+width*y];
		dy1 = data[x+width*y] - data[(x+1)+width*(y+1)];
		dy2 = data[x+width*(y-1)] - data[x+width*y];

		/* Average gradient measurements from both sides */
		dxs = ((dx1*dx1) + (dx2*dx2)) / 2;
		dys = ((dy1*dy1) + (dy2*dy2)) / 2;

		/* Calculate overall gradient */
		grad = dxs + dys;

		if ( grad < 100000 ) continue;

		mask_x = x;
		mask_y = y;

		do {

			max = data[mask_x+width*mask_y];
			did_something = 0;

			for ( sy=biggest(mask_y-PEAK_WINDOW_SIZE/2, 0);
			      sy<smallest(mask_y+PEAK_WINDOW_SIZE/2, height-1);
			      sy++ ) {
			for ( sx=biggest(mask_x-PEAK_WINDOW_SIZE/2, 0);
			      sx<smallest(mask_x+PEAK_WINDOW_SIZE/2, width-1);
			      sx++ ) {

				if ( data[sx+width*sy] > max ) {
					max = data[sx+width*sy];
					mask_x = sx;
					mask_y = sy;
					did_something = 1;
				}

			}
			}

			/* Abort if drifted too far from the foot point */
			if ( distance(mask_x, mask_y, x, y) > 50.0 ) break;

		} while ( did_something );

		/* Too far from foot point? */
		if ( distance(mask_x, mask_y, x, y) > 50.0 ) {
			nrej_dis++;
			continue;
		}

		/* Should be enforced by bounds used above.  Muppet check. */
		assert(mask_x < image->width);
		assert(mask_y < image->height);
		assert(mask_x >= 0);
		assert(mask_y >= 0);

		/* Isolated hot pixel? */
		if ( is_hot_pixel(image, mask_x, mask_y) ) {
			nrej_hot++;
			continue;
		}

		/* Centroid peak and get better coordinates.
		 * Don't bother doing polarisation/SA correction, because the
		 * intensity of this peak is only an estimate at this stage. */
		r = integrate_peak(image, mask_x, mask_y,
		                   &fx, &fy, &intensity, 0, 0);
		if ( r ) {
			/* Bad region - don't detect peak */
			nrej_bad++;
			continue;
		}

		/* It is possible for the centroid to fall outside the image */
		if ( (fx < 0.0) || (fx > image->width)
		  || (fy < 0.0) || (fy > image->height) ) {
			nrej_fra++;
			continue;
		}

		/* Check for a nearby feature */
		image_feature_closest(image->features, fx, fy, &d, &idx);
		if ( d < 15.0 ) {
			nrej_pro++;
			continue;
		}

		/* Add using "better" coordinates */
		image_add_feature(image->features, fx, fy, image, intensity,
		                  NULL);
		nacc++;

	}
	}

	if ( image->det != NULL ) {
		ncull = cull_peaks(image);
		nacc -= ncull;
	} else {
		STATUS("Not culling peaks because I don't have a "
		       "detector geometry file.\n");
		ncull = 0;
	}

	STATUS("%i accepted, %i box, %i hot, %i proximity, %i outside frame, "
	       "%i in bad regions, %i badrow culled.\n",
	       nacc, nrej_dis, nrej_hot, nrej_pro, nrej_fra, nrej_bad, ncull);
}


void dump_peaks(struct image *image, pthread_mutex_t *mutex)
{
	int i;

	/* Get exclusive access to the output stream if necessary */
	if ( mutex != NULL ) pthread_mutex_lock(mutex);

	printf("Peaks from peak search in %s\n", image->filename);
	printf("  x/px     y/px   (1/d)/nm^-1    Intensity\n");

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		double q, rx, ry, rz;
		struct imagefeature *f;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		map_position(image, f->x, f->y, &rx, &ry, &rz);
		q = modulus(rx, ry, rz);

		printf("%8.3f %8.3f %8.3f    %12.3f\n",
		       f->x, f->y, q/1.0e9, f->intensity);

	}

	printf("\n");

	if ( mutex != NULL ) pthread_mutex_unlock(mutex);
}


int find_projected_peaks(struct image *image, UnitCell *cell)
{
	int x, y;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	struct reflhit *hits;
	int n_hits = 0;

	hits = malloc(sizeof(struct reflhit)*MAX_HITS);
	if ( hits == NULL ) return 0;

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double hd, kd, ld;  /* Indices with decimal places */
		double dh, dk, dl;  /* Distances in h,k,l directions */
		signed int h, k, l;
		struct rvec q;
		double dist;
		int found = 0;
		int j;

		q = get_q(image, x, y, 1, NULL, 1.0/image->lambda);

		hd = q.u * ax + q.v * ay + q.w * az;
		kd = q.u * bx + q.v * by + q.w * bz;
		ld = q.u * cx + q.v * cy + q.w * cz;

		h = (signed int)rint(hd);
		k = (signed int)rint(kd);
		l = (signed int)rint(ld);

		dh = hd - h;
		dk = kd - k;
		dl = ld - l;
		dist = sqrt(pow(dh, 2.0) + pow(dk, 2.0) + pow(dl, 2.0));
		if ( dist > 0.1 ) continue;

		for ( j=0; j<n_hits; j++ ) {
			if ( (hits[j].h == h) && (hits[j].k == k)
			                      && (hits[j].l == l) ) {

				if ( dist < hits[j].min_distance ) {
					hits[j].min_distance = dist;
					hits[j].x = x;
					hits[j].y = y;
				}

				found = 1;

			}
		}

		if ( !found ) {
			hits[n_hits].min_distance = dist;
			hits[n_hits].x = x;
			hits[n_hits].y = y;
			hits[n_hits].h = h;
			hits[n_hits].k = k;
			hits[n_hits].l = l;
			n_hits++;
			assert(n_hits < MAX_HITS);
		}

	}
	}

	STATUS("Found %i reflections\n", n_hits);
	image->hits = hits;
	image->n_hits = n_hits;

	return n_hits;
}


int peak_sanity_check(struct image *image, UnitCell *cell)
{
	int i;
	int n_sane = 0;

	find_projected_peaks(image, cell);
	if ( image->n_hits == 0 ) return 0;  /* Failed sanity check: no peaks */

	for ( i=0; i<image->n_hits; i++ ) {

		double d;
		int idx;
		struct imagefeature *f;

		f = image_feature_closest(image->features,
                                          image->hits[i].x, image->hits[i].y,
		                          &d, &idx);
		if ( (f != NULL) && (d < PEAK_CLOSE) ) {
			n_sane++;
		}

	}

	STATUS("Sanity factor: %f / %f = %f\n", (float)n_sane,
	                                  (float)image->n_hits,
                                          (float)n_sane / (float)image->n_hits);
	if ( (float)n_sane / (float)image->n_hits < 0.8 ) return 0;

	return 1;
}


void output_intensities(struct image *image, UnitCell *cell,
                        pthread_mutex_t *mutex, int polar, int sa)
{
	int i;
	int n_found;
	int n_indclose = 0;
	int n_foundclose = 0;
	int n_veto = 0;
	int n_veto_second = 0;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double a, b, c, al, be, ga;

	if ( image->n_hits == 0 ) find_projected_peaks(image, cell);
	if ( image->n_hits == 0 ) return;

	/* Get exclusive access to the output stream if necessary */
	if ( mutex != NULL ) pthread_mutex_lock(mutex);

	/* Explicit printf() used here (not normally allowed) because
	 * we really want to output to stdout */
	printf("Reflections from indexing in %s\n", image->filename);
	printf("Orientation (wxyz): %7.5f %7.5f %7.5f %7.5f\n",
	       image->orientation.w, image->orientation.x,
	       image->orientation.y, image->orientation.z);
	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
	printf("Cell parameters %7.5f %7.5f %7.5f nm, %7.5f %7.5f %7.5f deg\n",
	       a*1.0e9, b*1.0e9, c*1.0e9,
	       rad2deg(al), rad2deg(be), rad2deg(ga));
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	printf("astar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	       asx/1e9, asy/1e9, asz/1e9);
	printf("bstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	       bsx/1e9, bsy/1e9, bsz/1e9);
	printf("cstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	       csx/1e9, csy/1e9, csz/1e9);

	if ( image->f0_available ) {
		printf("f0 = %7.5f (arbitrary gas detector units)\n",
		       image->f0);
	} else {
		printf("f0 = invalid\n");
	}

	for ( i=0; i<image->n_hits; i++ ) {

		float x, y, intensity;
		double d;
		int idx;
		struct imagefeature *f;

		/* Wait.. is there a really close feature which was detected? */
		if ( image->features != NULL ) {
			f = image_feature_closest(image->features,
			                          image->hits[i].x,
			                          image->hits[i].y,
			                          &d, &idx);
		} else {
			f = NULL;
		}
		if ( (f != NULL) && (d < PEAK_REALLY_CLOSE) ) {

			int r;

			/* f->intensity was measured on the filtered pattern,
			 * so instead re-integrate using old coordinates.
			 * This will produce further revised coordinates. */
			r = integrate_peak(image, f->x, f->y, &x, &y,
			                   &intensity, polar, sa);
			if ( r ) {
				/* The original peak (which also went through
				 * integrate_peak(), but with the mangled
				 * image data) would have been rejected if it
				 * was in a bad region.  Integration of the same
				 * peak included a bad region this time. */
				n_veto_second++;
				continue;
			}
			intensity = f->intensity;

		} else {

			int r;

			r = integrate_peak(image,
			                   image->hits[i].x,image->hits[i].y,
			                   &x, &y, &intensity, polar, sa);
			if ( r ) {
				/* Plain old ordinary peak veto */
				n_veto++;
				continue;
			}

		}

		if ( (f != NULL) && (d < PEAK_CLOSE) ) {
			n_indclose++;
		}

		/* Write h,k,l, integrated intensity and centroid coordinates */
		printf("%3i %3i %3i %6f (at %5.2f,%5.2f)\n",
		       image->hits[i].h, image->hits[i].k, image->hits[i].l,
		       intensity, x, y);

		image->hits[i].x = x;
		image->hits[i].y = y;

	}

	n_found = image_feature_count(image->features);
	for ( i=0; i<n_found; i++ ) {

		struct imagefeature *f;
		int j;

		f = image_get_feature(image->features, i);

		if ( f == NULL ) continue;

		for ( j=0; j<image->n_hits; j++ ) {

			double d;

			d = pow(image->hits[j].x-f->x, 2.0)
			  + pow(image->hits[j].y-f->y, 2.0);

			if ( d < PEAK_CLOSE ) n_foundclose++;

		}

	}

	printf("Peak statistics: %i peaks found by the peak search out of "
	       "%i were close to indexed positions. "
	       "%i indexed positions out of %i were close to detected peaks.\n",
	       n_foundclose, n_found, n_indclose, image->n_hits);
	printf("%i integrations using indexed locations were aborted because "
	       "they hit one or more bad pixels.\n", n_veto);
	printf("%i integrations using peak search locations were aborted "
	       "because they hit one or more bad pixels.\n", n_veto_second);
	/* Blank line at end */
	printf("\n");

	if ( mutex != NULL ) pthread_mutex_unlock(mutex);
}
