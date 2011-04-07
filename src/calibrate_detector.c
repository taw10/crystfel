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
"  -m, --method=<method>      The calibration method.\n"
"               xy            Determine panel shifts in plane of detector\n"
"  -o, --output=<file>        Output results here"
"  -n, --npeaks=<number>      Don't refine unless this many peaks observed\n"
"                             (in the whole stream, not a single shot)\n"
"\n");
}

//FIXME: should this function be in detector.h?
int calculate_projected_peak(struct panel *panel, struct rvec q,
                             double kk, double *fs, double *ss)
{

	if (panel == NULL) return 1;

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

//FIXME: should this function be in detector.h?
struct rvec nearest_bragg(struct image image, struct rvec q)
{

	struct rvec g;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double hd, kd, ld;
	double h, k, l;

	/* miller indices of nearest Bragg reflection */
	cell_get_cartesian(image.indexed_cell, &ax, &ay, &az,
	                                       &bx, &by, &bz,
	                                       &cx, &cy, &cz);

	hd = q.u * ax + q.v * ay + q.w * az;
	kd = q.u * bx + q.v * by + q.w * bz;
	ld = q.u * cx + q.v * cy + q.w * cz;

	h = lrint(hd);
	k = lrint(kd);
	l = lrint(ld);

	/* now get scattering vector for reflectin [hkl]
	 * this means solving the equation U*q = h */
	double U[] = {ax, ay, az,
	              bx, by, bz,
	              cx, cy, cz};

	double hvec[] = {h,k,l};

	gsl_matrix_view m
		= gsl_matrix_view_array (U, 3, 3);

	gsl_vector_view b
		= gsl_vector_view_array (hvec, 3);

	gsl_vector *x = gsl_vector_alloc (3);

	int s;

	gsl_permutation * perm = gsl_permutation_alloc (3);
	gsl_linalg_LU_decomp (&m.matrix, perm, &s);
	gsl_linalg_LU_solve (&m.matrix, perm, &b.vector, x);

	/* outgoing wavevector for [hkl] */
	g.u = x->data[0];
	g.v = x->data[1];
	g.w = x->data[2];

	return g;

}


int main(int argc, char *argv[])
{


	char c;
	struct image image;
	struct panel p;
	char *filename = NULL;
	char *geometry = NULL;
	char *method = NULL;
	char *outputfile = NULL;
	FILE *fh = NULL;
	FILE *outfh = NULL;
	int nFeatures;
	int i;
	int fail;
	int nChunks;
	int minpeaks=0;

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
	while ((c = getopt_long(argc, argv, "hi:g:m:o:n:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;
		case 0 :
			break;
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
		ERROR("Problem opening file\n");
		return 1;
	}
	printf("Read stream file: %s\n",filename);

	if ( geometry == NULL ) {
		ERROR("You need to specify a geometry file with --geometry\n");
		return 1;
	}

	image.det = get_detector_geometry(geometry);
	if ( image.det == NULL ) {
		ERROR("Failed to read detector geometry from %s\n", geometry);
		return 1;
	}
	printf("Read geometry file: %s\n",geometry);
	image.width = image.det->max_fs;
	image.height = image.det->max_ss;
	free(geometry);

	if ( !(outputfile == NULL) ) {
		printf("Writing result to file: %s\n",outputfile);
      outfh = fopen(outputfile,"w");
   } else {
		ERROR("You did not specify an output file.\n");
		return 1;
	}
	if ( minpeaks == 0 ) {
		minpeaks = 1;
		printf("You did not specify minimum number of peaks."
             " Setting default value of %d\n",minpeaks);
	}

	if ( method == NULL ) {
		printf("You did not specify a refinement method-- setting default.\n");
		method = strdup("xy");
	}


	if ( !strcmp(method,"xy" ) ) {

		printf("Using refinement method %s\n",method);

		double * weightedSumFS, * weightedSumSS;
		double * summedWeights;
		double * meanShiftFS, * meanShiftSS;
		int * peaksFound;
		double cnx, cny;
		double xsh, ysh;

		weightedSumFS = (double *)calloc(sizeof(double), image.det->n_panels);
		weightedSumSS = (double *)calloc(sizeof(double), image.det->n_panels);
		summedWeights = (double *)calloc(sizeof(double), image.det->n_panels);
		peaksFound = (int *)calloc(sizeof(int), image.det->n_panels);
		meanShiftFS = (double *)calloc(sizeof(double), image.det->n_panels);
		meanShiftSS = (double *)calloc(sizeof(double), image.det->n_panels);

		/* initialize arrays (is there a standard function to do this?)  */
		int pi;
		for (pi=0; pi<image.det->n_panels; pi++) {
			weightedSumFS[pi] = 0;
	      weightedSumSS[pi] = 0;
   	   summedWeights[pi] = 0;
     		peaksFound[pi] = 0;
      	meanShiftFS[pi] = 0;
      	meanShiftSS[pi] = 0;
		}


		fesetround(1);  /* Round towards nearest */

		nChunks = 0;
		while (1) {

		   // check for peaks in next file
			fail = read_chunk(fh, &image);
		   if ( fail == 1 ) {
		      // FIXME:should check if this is EOF, or broken file handle
		      break;
		   }
			nChunks += 1;


			// move on to the next chunk if no peaks found
			if ( image.features == NULL ) {
				continue;
			}

			// now loop through peaks to determine mean panel shift
			nFeatures = image_feature_count(image.features);

			for (i=0; i<nFeatures; i++) {

				struct panel * p;
				struct imagefeature * thisFeature;
				struct rvec q;
				double twotheta;
				struct rvec g;
				double fs, ss;   /* observed peaks */
				double pfs, pss; /* predicted peaks */
				double dfs, dss; /* observed - predicted */
				int pi;
				double thisWeight;
				//int fail;

				/* if we find a feature, determine peak position */
				thisFeature = image_get_feature(image.features,i);
				if ( thisFeature == NULL ) {
					continue;
				}

				fs = thisFeature->fs;
				ss = thisFeature->ss;

				p = find_panel(image.det, fs, ss);
				if ( p == NULL ) {
					continue;
				}
				if ( p->no_index ) continue;

				/* now determine the predicted peak position */

				/* scattering vector of this peak */
				q = get_q(&image, fs, ss, &twotheta, 1.0/image.lambda);

				/* scattering vector of nearest bragg peak */
				g = nearest_bragg(image, q);

				/* coordinate of this predicted peak */
				fail = calculate_projected_peak(p,g,1/image.lambda,&pfs,&pss);

				/* check for error, e.g. out of panel */
				if ( fail != 0 ) continue;

				/* Finally, we have the shift in position of this peak */
				dfs = pfs - fs;
				dss = pss - ss;

				/* Add this shift to the weighted sum over shifts */
				pi = find_panel_number(image.det,fs,ss);
				thisWeight = 1; // FIXME: use real weighting some day
				weightedSumFS[pi] += thisWeight*dfs;
				weightedSumSS[pi] += thisWeight*dss;
				summedWeights[pi] += thisWeight;
				peaksFound[pi] += 1;

			} /* end loop through image features */

		} /* end loop through stream chunks */

		/* calculate weighted average shift in peak positions */
		for (pi=0; pi < image.det->n_panels; pi++) {
			meanShiftFS[pi] = weightedSumFS[pi]/summedWeights[pi];
			meanShiftSS[pi] = weightedSumSS[pi]/summedWeights[pi];
		}

		/* first populate the image structure with new geometry info */
		for (pi=0; pi < image.det->n_panels; pi++) {

			p = image.det->panels[pi];

			xsh = 0;
			ysh = 0;

			if ( peaksFound[pi] >= minpeaks ) {

				/* convert shifts from raw coords to lab frame */
				xsh = meanShiftFS[pi]*p.fsx + meanShiftSS[pi]*p.ssx;
				ysh = meanShiftFS[pi]*p.fsy + meanShiftSS[pi]*p.ssy;

				/* add shifts to original panel corner locations */
				cnx = p.cnx + xsh;
				cny = p.cny + ysh;

			} else {
				/* not refined? use original coordinates */
				cnx = p.cnx;
				cny = p.cny;
			}

			image.det->panels[pi].cnx = cnx;
			image.det->panels[pi].cny = cny;
			if ( peaksFound[pi] < minpeaks) image.det->panels[pi].no_index = 1;

			printf("panel %s, # peaks: %10d, mean shifts: %f %f\n",p.name, peaksFound[pi],xsh,ysh);

		}

		/* write the new geometry file */
		print_detector_geometry(outfh,&image);

	} else {

		printf("Refinement method %s not recognized\n",method);
   	return 1;

	}


	fclose(outfh);

	printf("Done!\n");

	return 0;
}
