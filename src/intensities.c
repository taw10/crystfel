/*
 * intensities.c
 *
 * Extract intensities from patterns
 *
 * (c) 2007-2009 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "image.h"
#include "intensities.h"
#include "cell.h"
#include "sfac.h"


#define MAX_HITS (1024)


struct reflhit {
	signed int h;
	signed int k;
	signed int l;
	double min_distance;
	int x;
	int y;
};


static int sum_nearby_points(int16_t *data, int width, int x, int y)
{
	int dx, dy;
	int intensity = 0;

	for ( dx=-3; dx<=3; dx++ ) {
	for ( dy=-3; dy<=3; dy++ ) {
		intensity += data[(x+dx) + width*(y+dy)];
	}
	}

	return intensity;
}


void output_intensities(struct image *image)
{
	int x, y;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	struct reflhit hits[MAX_HITS];
	int n_hits = 0;
	int i;

	cell_get_cartesian(image->molecule->cell, &ax, &ay, &az,
		                                  &bx, &by, &bz,
		                                  &cx, &cy, &cz);

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double hd, kd, ld;  /* Indices with decimal places */
		double dh, dk, dl;  /* Distances in h,k,l directions */
		signed int h, k, l;
		struct rvec q;
		double dist;
		int found = 0;
		int j;

		q = image->qvecs[x + image->width*y];

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
	printf("New pattern: %7.5f %7.5f %7.5f %7.5f\n",
	       image->orientation.w, image->orientation.x,
	       image->orientation.y, image->orientation.z);
	for ( i=0; i<n_hits; i++ ) {

		int intensity;

		/* Bounds check */
		if ( hits[i].x + 3 >= image->width ) continue;
		if ( hits[i].x - 3 < 0 ) continue;
		if ( hits[i].y + 3 >= image->height ) continue;
		if ( hits[i].y - 3 < 0 ) continue;

		intensity = sum_nearby_points(image->data, image->width,
		                              hits[i].x, hits[i].y);

		printf("%3i %3i %3i %6i (at %i,%i)\n",
		       hits[i].h, hits[i].k, hits[i].l, intensity,
		       hits[i].x, hits[i].y);

	}

	/* Blank line at end */
	printf("\n");
}
