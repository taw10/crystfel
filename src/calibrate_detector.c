/*
 * calibrate_detector.c
 *
 * Attempt to refine detector geometry
 *
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
"  -m, --method=<method>      The calibration method.  Possiblities are\n"
"               xy            Determine panel shifts in plane of detector\n"
"                             (i.e. ignore camera length calibration\n"
"  -o, --output=<file>        Output results here"
"  -n, --npeaks=<number>      Don't refine unless this many peaks observed\n"
"                             (in the whole stream, not a single shot)\n"
"\n");
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

				/* update reflection list since geometry may have changed 
            from what is in the data stream (e.g. if iterating this
            calibration procedure) */
 	 		//image.reflections = 
         //	find_projected_peaks(&image,image.indexed_cell,0, 0.1);

			//printf("chunk %d\n",nChunks);
			//cell_print(image.indexed_cell);

			// now loop through peaks to determine mean panel shift
			nFeatures = image_feature_count(image.features);

			for (i=0; i<nFeatures; i++) {

				struct panel * p;
				struct imagefeature * thisFeature;			
				double ax, ay, az;
				double bx, by, bz;
				double cx, cy, cz;
				double hd, kd, ld;  /* Indices with decimal places */
				signed int h, k, l;
				struct rvec q;
				double twotheta;
				double fs, ss;   /* observed peaks */
				double pfs, pss; /* predicted peaks */
				double dfs, dss; /* observed - predicted */
				int pi;
				double thisWeight;

	
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
				q = get_q(&image, fs, ss, &twotheta, 1.0/image.lambda);

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
     
				gsl_vector *xx = gsl_vector_alloc (3);
       
				int s;
     
				gsl_permutation * perm = gsl_permutation_alloc (3);
				gsl_linalg_LU_decomp (&m.matrix, perm, &s);
				gsl_linalg_LU_solve (&m.matrix, perm, &b.vector, xx);
   

				// outgoing wavevector 
				double x = xx->data[0];
				double y = xx->data[1];
				double z = xx->data[2];

				double kk;
				double xd, yd, cl;
				double plx, ply;

				kk = 1/image.lambda;
				const double den = kk + z;

				/* Camera length for this panel */
				cl = p->clen;

				/* Coordinates of peak relative to central beam, in m */
				xd = cl * x / den;
				yd = cl * y / den;
				
				/* Convert to pixels */
				xd *= p->res;
				yd *= p->res;
				
				/* Convert to relative to the panel corner */
				plx = xd - p->cnx;
				ply = yd - p->cny;
				
				pfs = p->xfs*plx + p->yfs*ply;
				pss = p->xss*plx + p->yss*ply;
				
				pfs += p->min_fs;
				pss += p->min_ss;

				/* Now, is this on this panel? */
				if ( fs < p->min_fs ) continue;
				if ( fs > p->max_fs ) continue;
				if ( ss < p->min_ss ) continue;
				if ( ss > p->max_ss ) continue;

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
		
		/* now generate a new geometry file */
		for (pi=0; pi < image.det->n_panels; pi++) {
	
			p = image.det->panels[pi];

			xsh = 0;
			ysh = 0;
	
			if ( peaksFound[pi] >= minpeaks ) {

				//printf("meanShift: %f %f\n",meanShiftFS[pi],meanShiftSS[pi]);

				/* convert shifts from raw coords to lab frame */
				xsh = meanShiftFS[pi]*p.fsx + meanShiftSS[pi]*p.ssx;
				ysh = meanShiftFS[pi]*p.fsy + meanShiftSS[pi]*p.ssy;
		
				/* add shifts to original panel corner locations */
				cnx = p.cnx + xsh;
				cny = p.cny + ysh;
	
				//printf("new panel shift: %f %f\n",cnx,cny);

	
			} else {

				/* not refined: use original coordinates */
				cnx = p.cnx;
				cny = p.cny;

			}

			printf("panel %s, # peaks: %10d, mean shifts: %f %f\n",p.name, peaksFound[pi],xsh,ysh);

			//FIXME: there should be a function in geometry.c to write 
			//       these values to text file, since it will be useful on 
			//       other places as well (e.g. writing the geometry to 
			//       the data stream)
			fprintf(outfh,"%s/min_fs = %d\n",p.name,p.min_fs);
			fprintf(outfh,"%s/min_ss = %d\n",p.name,p.min_ss);
			fprintf(outfh,"%s/max_fs = %d\n",p.name,p.max_fs);
	      fprintf(outfh,"%s/max_ss = %d\n",p.name,p.max_ss);
			fprintf(outfh,"%s/badrow_direction = %C\n",p.name,p.badrow);
			fprintf(outfh,"%s/res = %g\n",p.name,p.res);
			fprintf(outfh,"%s/peak_sep = %g\n",p.name,p.peak_sep);
			fprintf(outfh,"%s/clen = %s\n",p.name,p.clen_from);
			//FIXME: the following is sketchy, but it will work for now.  we need
			//       to generalise the parser in detector.c
			char coord;
			char sign;
			if (p.fsx != 0){
				if (p.fsx>0){sign='+';}else{sign='-';}
				coord = 'x';
			} else {
				if (p.fsy>0){sign='+';}else{sign='-';}
				coord = 'y';
			}
			fprintf(outfh,"%s/fs = %C%C\n",p.name, sign, coord);
			if (p.ssx != 0){
            if (p.ssx>0){sign='+';}else{sign='-';}
            coord = 'x';
         } else {
            if (p.ssy>0){sign='+';}else{sign='-';}
            coord = 'y';
         }
			fprintf(outfh,"%s/ss = %C%C\n",p.name, sign, coord);
			fprintf(outfh,"%s/corner_x = %g\n",p.name,cnx);
	      fprintf(outfh,"%s/corner_y = %g\n",p.name,cny);
			if ( peaksFound[pi] < minpeaks ) {
				fprintf(outfh,"%s/no_index = %d\n",p.name,1);
			} else {
				fprintf(outfh,"%s/no_index = %d\n",p.name,p.no_index);
			}
			fprintf(outfh,"\n\n");
	
		}

	} else {

		printf("Refinement method %s not recognized\n",method);
   	return 1;

	}
	
	if ( !(outputfile == NULL) ) {
      fclose(outfh);
   }

	

	printf("Done!\n");

	return 0;
}
