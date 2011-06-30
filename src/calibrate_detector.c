/*
 * calibrate_detector.c
 *
 * Attempt to refine detector geometry
 *
 * (c) 2011 Rick Kirian <rkirian@asu.edu>
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
#include <fenv.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

#include "utils.h"
#include "image.h"
#include "detector.h"
#include "index.h"
#include "hdf5-file.h"
#include "stream.h"
#include "peaks.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options] -i <file.h5>\n\n", s);
	printf(
"Stream-based optimisation of detector geometry.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"  -g. --geometry=<file>      Get detector geometry from file.\n"
"  -i, --input=<file>         Input filename.\n"
"  -m, --method=<method>      The calibration method.  Choose from:\n"
"               xy            Determine panel shifts in plane of detector\n"
"  -o, --output=<file>        Name of output geometry file.\n"
"  -n, --npeaks=<number>      Don't refine unless this many peaks are found\n"
"                              in the whole stream.\n"
"\n");
}


static int write_detector_geometry(const char *filename, struct detector *det)
{
	struct panel *p;
	int pi;
	FILE *fh;

	if ( filename == NULL ) return 2;
	if ( det->n_panels < 1 ) return 3;

	fh = fopen(filename, "w");
	if ( fh == NULL ) return 1;

	for ( pi=0; pi<det->n_panels; pi++) {

		p = &(det->panels[pi]);

		if ( p == NULL ) return 4;

		fprintf(fh, "%s/min_fs = %d\n", p->name, p->min_fs);
		fprintf(fh, "%s/min_ss = %d\n", p->name, p->min_ss);
		fprintf(fh, "%s/max_fs = %d\n", p->name, p->max_fs);
		fprintf(fh, "%s/max_ss = %d\n", p->name, p->max_ss);
		fprintf(fh, "%s/badrow_direction = %C\n", p->name, p->badrow);
		fprintf(fh, "%s/res = %g\n", p->name, p->res);
		fprintf(fh, "%s/peak_sep = %g\n", p->name, p->peak_sep);
		fprintf(fh, "%s/clen = %s\n", p->name, p->clen_from);
		fprintf(fh, "%s/fs = %+fx %+fy\n", p->name, p->fsx, p->fsy);
		fprintf(fh, "%s/ss = %+fx %+fy\n", p->name, p->ssx, p->ssy);
		fprintf(fh, "%s/corner_x = %g\n", p->name, p->cnx);
		fprintf(fh, "%s/corner_y = %g\n", p->name, p->cny);
		fprintf(fh, "%s/no_index = %d\n", p->name, p->no_index);

	}
	fclose(fh);

	return 0;
}


static int calculate_projected_peak(struct panel *panel, struct rvec q,
                                    double kk, double *fs, double *ss)
{
	if ( panel == NULL ) return 1;

	double xd, yd, cl;
	double plx, ply;
	const double den = kk + q.w;

	/* Camera length for this panel */
	cl = panel->clen;

	/* Coordinates of peak relative to central beam, in m */
	xd = cl * q.u / den;
	yd = cl * q.v / den;

	/* Convert to pixels */
	xd *= panel->res;
	yd *= panel->res;

	/* Convert to relative to the panel corner */
	plx = xd - panel->cnx;
	ply = yd - panel->cny;

	*fs = panel->xfs*plx + panel->yfs*ply;
	*ss = panel->xss*plx + panel->yss*ply;

	*fs += panel->min_fs;
	*ss += panel->min_ss;

	/* Now, is this on this panel? */
	if ( *fs < panel->min_fs ) return 3;
	if ( *fs > panel->max_fs ) return 3;
	if ( *ss < panel->min_ss ) return 3;
	if ( *ss > panel->max_ss ) return 3;

	return 0;
}


static struct rvec nearest_bragg(struct image *image, struct rvec q)
{
	struct rvec g;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double hd, kd, ld;
	double h, k, l;
	int s;
	double U[9];
	double hvec[3];

	/* Miller indices of nearest Bragg reflection */
	cell_get_cartesian(image->indexed_cell, &ax, &ay, &az,
	                                        &bx, &by, &bz,
	                                        &cx, &cy, &cz);

	hd = q.u * ax + q.v * ay + q.w * az;
	kd = q.u * bx + q.v * by + q.w * bz;
	ld = q.u * cx + q.v * cy + q.w * cz;

	h = lrint(hd);
	k = lrint(kd);
	l = lrint(ld);

	/* Now get scattering vector for reflection hkl
	 * by solving the equation U*q = h */
	U[0] = ax;  U[1] = ay;  U[2] = az;
	U[3] = bx;  U[4] = by;  U[5] = bz;
	U[6] = cx;  U[7] = cy;  U[8] = cz;

	hvec[0] = h;  hvec[1] = k;  hvec[2] = l;

	gsl_matrix_view m = gsl_matrix_view_array(U, 3, 3);
	gsl_vector_view b = gsl_vector_view_array(hvec, 3);
	gsl_vector *x = gsl_vector_alloc(3);

	gsl_permutation *perm = gsl_permutation_alloc(3);
	gsl_linalg_LU_decomp(&m.matrix, perm, &s);
	gsl_linalg_LU_solve(&m.matrix, perm, &b.vector, x);

	/* Outgoing wavevector for hkl */
	g.u = x->data[0];
	g.v = x->data[1];
	g.w = x->data[2];

	return g;
}


static void refine_xy(FILE *fh, struct image *image, int minpeaks,
                       const char *outfilename)
{
	double *weightedSumFS;
	double *weightedSumSS;
	double *summedWeights;
	double *meanShiftFS;
	double *meanShiftSS;
	int *peaksFound;
	double cnx, cny;
	double xsh, ysh;
	int pi;
	int nChunks;

	weightedSumFS = calloc(image->det->n_panels, sizeof(double));
	weightedSumSS = calloc(image->det->n_panels, sizeof(double));
	summedWeights = calloc(image->det->n_panels, sizeof(double));
	peaksFound = calloc(image->det->n_panels, sizeof(double));
	meanShiftFS = calloc(image->det->n_panels, sizeof(double));
	meanShiftSS = calloc(image->det->n_panels, sizeof(double));

	/* Initialize arrays  */
	for (pi=0; pi<image->det->n_panels; pi++) {
		weightedSumFS[pi] = 0;
		weightedSumSS[pi] = 0;
		summedWeights[pi] = 0;
		peaksFound[pi] = 0;
		meanShiftFS[pi] = 0;
		meanShiftSS[pi] = 0;
	}

	fesetround(1);  /* Round towards nearest */

	nChunks = 0;
	while ( 1 )  {

		int fail;
		int nFeatures;
		int i;

		/* Get next chunk */
		fail = read_chunk(fh, image);
		if ( fail ) break;
		nChunks += 1;

		/* Skip if no peaks found */
		if ( image->features == NULL ) {
			continue;
		}

		/* Loop through peaks to determine mean panel shift */
		nFeatures = image_feature_count(image->features);

		for ( i=0; i<nFeatures; i++ ) {

			struct panel *p;
			struct imagefeature *thisFeature;
			struct rvec q;
			double twotheta;
			struct rvec g;
			double fs, ss;   /* Observed peaks */
			double pfs, pss; /* Predicted peaks */
			double dfs, dss; /* Observed - predicted */
			int pi;
			double thisWeight;

			/* If we find a feature, determine peak
			 * position */
			thisFeature = image_get_feature(image->features, i);
			if ( thisFeature == NULL ) {
				continue;
			}

			fs = thisFeature->fs;
			ss = thisFeature->ss;

			p = find_panel(image->det, fs, ss);
			if ( p == NULL ) {
				continue;
			}
			if ( p->no_index ) continue;

			/* Now determine the predicted peak position.
			 * Scattering vector of this peak */
			q = get_q(image, fs, ss, &twotheta, 1.0/image->lambda);

			/* Scattering vector of nearest bragg peak */
			g = nearest_bragg(image, q);

			/* Coordinate of this predicted peak */
			fail = calculate_projected_peak(p, g, 1.0/image->lambda,
			                                &pfs, &pss);
			if ( fail != 0 ) continue;

			/* Finally, we have the shift for this peak */
			dfs = pfs - fs;
			dss = pss - ss;

			/* Add this shift to the weighted sum over shifts */
			pi = find_panel_number(image->det,fs,ss);
			thisWeight = 1.0;
			weightedSumFS[pi] += thisWeight*dfs;
			weightedSumSS[pi] += thisWeight*dss;
			summedWeights[pi] += thisWeight;
			peaksFound[pi] += 1;

		}

	}

	/* Calculate weighted average shift in peak positions */
	for ( pi=0; pi<image->det->n_panels; pi++ ) {
		meanShiftFS[pi] = weightedSumFS[pi]/summedWeights[pi];
		meanShiftSS[pi] = weightedSumSS[pi]/summedWeights[pi];
	}

	/* Populate the image structure with new geometry info */
	for ( pi=0; pi<image->det->n_panels; pi++ ) {

		struct panel *p;

		p = &image->det->panels[pi];

		xsh = 0;
		ysh = 0;

		if ( peaksFound[pi] >= minpeaks ) {

			/* Convert shifts from raw coords to lab frame */
			xsh = meanShiftFS[pi]*p->fsx + meanShiftSS[pi]*p->ssx;
			ysh = meanShiftFS[pi]*p->fsy + meanShiftSS[pi]*p->ssy;

			/* Add shifts to original panel corner locations */
			cnx = p->cnx + xsh;
			cny = p->cny + ysh;

		} else {

			/* Not refined? use original coordinates */
			cnx = p->cnx;
			cny = p->cny;

		}

		image->det->panels[pi].cnx = cnx;
		image->det->panels[pi].cny = cny;
		if ( peaksFound[pi] < minpeaks) {
			image->det->panels[pi].no_index = 1;
		}

		STATUS("Panel %s, # peaks: %10d, mean shifts: %f %f\n",
		       p->name, peaksFound[pi], xsh, ysh);

	}

	/* Write the new geometry file */
	write_detector_geometry(outfilename, image->det);
}


int main(int argc, char *argv[])
{
	char c;
	struct image image;
	char *filename = NULL;
	char *geometry = NULL;
	char *method = NULL;
	char *outputfile = NULL;
	FILE *fh = NULL;
	FILE *outfh = NULL;
	int minpeaks = 0;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"geometry",           1, NULL,               'g'},
		{"method",             1, NULL,               'm'},
		{"output",             1, NULL,               'o'},
		{"npeaks",             1, NULL,               'n'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:g:m:o:n:",
	                        longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;
		case 'i' :
			filename = strdup(optarg);
			break;
		case 'g' :
			geometry = strdup(optarg);
			break;
		case 'm' :
			method = strdup(optarg);
			break;
		case 'o' :
			outputfile = strdup(optarg);
			break;
		case 'n' :
			minpeaks = atoi(optarg);
			break;
		case 0 :
			break;
		default :
			return 1;
		}

	}

	if ( filename == NULL ) {
		ERROR("You must specify the input filename with -i\n");
		return 1;
	}

	fh = fopen(filename,"r");
	if ( fh == NULL ) {
		ERROR("Couldn't open file '%s'\n", filename);
		return 1;
	}

	if ( geometry == NULL ) {
		ERROR("You need to specify a geometry file with --geometry\n");
		return 1;
	}

	image.det = get_detector_geometry(geometry);
	if ( image.det == NULL ) {
		ERROR("Failed to read detector geometry from %s\n", geometry);
		return 1;
	}
	free(geometry);

	image.width = image.det->max_fs;
	image.height = image.det->max_ss;

	if ( outputfile != NULL ) {
		STATUS("Writing result to '%s'\n", outputfile);
		outfh = fopen(outputfile, "w");
	} else {
		ERROR("You need to specify an output file.\n");
		return 1;
	}

	if ( minpeaks == 0 ) {
		minpeaks = 1;
		STATUS("You did not specify minimum number of peaks.\n");
		STATUS("Using default value of %d\n", minpeaks);
	}

	if ( method == NULL ) {
		STATUS("You did not specify a refinement method "
		       "- using default.\n");
		method = strdup("xy");
	}

	if ( strcmp(method,"xy") == 0 ) {

		STATUS("Using refinement method %s\n", method);

		refine_xy(fh, &image, minpeaks, outputfile);

	} else {

		printf("Refinement method %s not recognized\n",method);
   		return 1;

	}

	fclose(outfh);

	return 0;
}
