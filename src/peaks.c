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

#include "image.h"
#include "utils.h"
#include "index.h"
#include "peaks.h"
#include "detector.h"
#include "filters.h"
#include "diffraction.h"


#define MAX_HITS (1024)

struct reflhit {
	signed int h;
	signed int k;
	signed int l;
	double min_distance;
	int x;
	int y;
};


#define PEAK_WINDOW_SIZE (10)
#define MAX_PEAKS (2048)
#define INTEGRATION_RADIUS (10)

static int in_streak(int x, int y)
{
	if ( (y>512) && (y<600) && (abs(x-489)<15) ) return 1;
	if ( (y>600) && (abs(x-480)<25) ) return 1;
	return 0;
}


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


/* Post-processing of the peak list to remove noise */
static void cull_peaks(struct image *image)
{
	int i, n;
	int nelim = 0;

	n = image_feature_count(image->features);

	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		int j, ncol;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		/* How many peaks are in exactly the same column? */
		ncol = 0;
		for ( j=0; j<n; j++ ) {

			struct imagefeature *g;

			if ( i==j ) continue;

			g = image_get_feature(image->features, j);
			if ( g == NULL ) continue;
			if ( f->x == g->x ) ncol++;

		}

		/* More than three? */
		if ( ncol <= 3 ) continue;

		/* Yes?  Delete them all... */
		nelim = 0;
		for ( j=0; j<n; j++ ) {
			struct imagefeature *g;
			g = image_get_feature(image->features, j);
			if ( g == NULL ) continue;
			if ( f->x == g->x ) {
				image_remove_feature(image->features, j);
				nelim++;
			}
		}

	}

	//STATUS("%i peaks eliminated\n", nelim);
}


static void integrate_peak(struct image *image, int xp, int yp,
                           float *xc, float *yc, float *intensity)
{
	signed int x, y;
	const int lim = INTEGRATION_RADIUS * INTEGRATION_RADIUS;
	int total = 0;
	int xct = 0;
	int yct = 0;

	for ( x=-INTEGRATION_RADIUS; x<+INTEGRATION_RADIUS; x++ ) {
	for ( y=-INTEGRATION_RADIUS; y<+INTEGRATION_RADIUS; y++ ) {

		int val;

		/* Circular mask */
		if ( x*x + y*y > lim ) continue;

		if ( ((x+xp)>=image->width) || ((x+xp)<0) ) continue;
		if ( ((y+yp)>=image->height) || ((y+yp)<0) ) continue;

		val = image->data[(x+xp)+image->width*(y+yp)];

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
	int nacc = 0;

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

		/* Overall threshold */
		if ( data[x+width*y] < 800 ) continue;

		/* Ignore streak */
		if ( in_streak(x, y) ) continue;

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

		/* Centroid peak and get better coordinates */
		integrate_peak(image, mask_x, mask_y, &fx, &fy, &intensity);

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

	STATUS("%i accepted, %i box, %i hot, %i proximity, %i outside frame\n",
	       nacc, nrej_dis, nrej_hot, nrej_pro, nrej_fra);

	cull_peaks(image);
}


void dump_peaks(struct image *image)
{
	int i;

	printf("x/px\ty/px\t(1/d)/nm^-1\n");

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		double q, rx, ry, rz;
		struct imagefeature *f;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		map_position(image, f->x, f->y, &rx, &ry, &rz);
		q = modulus(rx, ry, rz);

		printf("%7.3f\t%7.3f\t%7.3f\t%7.3f\n", f->x, f->y, q/1.0e9, 1.0);

	}
}


void output_intensities(struct image *image, UnitCell *cell)
{
	int x, y;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	struct reflhit hits[MAX_HITS];
	int n_hits = 0;
	int i;

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

	/* Explicit printf() used here (not normally allowed) because
	 * we really want to output to stdout */
	printf("New pattern: %s %7.5f %7.5f %7.5f %7.5f\n", image->filename,
	       image->orientation.w, image->orientation.x,
	       image->orientation.y, image->orientation.z);
	for ( i=0; i<n_hits; i++ ) {

		float x, y, intensity;
		double d;
		int idx;
		struct imagefeature *f;

		integrate_peak(image, hits[i].x, hits[i].y, &x, &y, &intensity);

		/* Wait.. is there a closer feature which was detected? */
		f = image_feature_closest(image->features, x, y, &d, &idx);
		if ( (d < 10.0) && (f != NULL) ) {

			/* f->intensity was measured on the filtered pattern,
			 * so instead re-integrate using old coordinates.
			 * This will produce further revised coordinates. */
			integrate_peak(image, f->x, f->y, &x, &y, &intensity);
			intensity = f->intensity;

		}

		/* Write h,k,l, integrated intensity and centroid coordinates */
		printf("%3i %3i %3i %6f (at %5.2f,%5.2f)\n",
		       hits[i].h, hits[i].k, hits[i].l, intensity, x, y);

	}

	/* Blank line at end */
	printf("\n");
}
