/*
 * mosflm.c
 *
 * Invoke the DPS auto-indexing algorithm through MOSFLM
 *
 * (c) 2010 Richard Kirian <rkirian@asu.edu>
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
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "image.h"
#include "mosflm.h"
#include "utils.h"
#include "sfac.h"
#include "peaks.h"


/*

todo

** properly read the newmat file (don't use fscanf-- spaces between numers
   are not guaranteed)

** "success" is indicated by existence of NEWMAT file written by mosflm.
   Better to interact with mosflm directly in order to somehow verify success.

** investigate how these keywords affect mosflms behavior:

   MOSAICITY
   DISPERSION
   DIVERGENCE
   POLARISATION
   POSTREF BEAM
   POSTREF USEBEAM OFF
   PREREFINE ON
   EXTRA ON
   POSTREF ON

   These did not seem to affect the results by my (Rick's) experience, probably
   because they are only used conjunction with image intensity data, but it's
   worth another look at the documentation.

*/



#define MOSFLM_VERBOSE 0


static int read_newmat(const char * filename, struct image *image)
{
	FILE * fh;
	float asx, asy, asz;
	float bsx, bsy, bsz;
	float csx, csy, csz;
	int n;
	double c;

	fh = fopen(filename,"r");
	if (fh == NULL){
		STATUS("found newmat.\n");
		return 1;
	}
	n  = fscanf(fh,"%f %f %f\n",&asx,&bsx,&csx);
	n += fscanf(fh,"%f %f %f\n",&asy,&bsy,&csy);
	n += fscanf(fh,"%f %f %f\n",&asz,&bsz,&csz);
	if (n != 9) {
		STATUS("<9 parameters.\n");
		return 1;
	}
	fclose(fh);

	/* mosflm A matrix is multiplied by lambda, so fix this */
	c = 1/image->lambda;

	image->candidate_cells[0] = cell_new();

	cell_set_reciprocal(image->candidate_cells[0],
                            asz*c, asy*c, asx*c,
                            bsz*c, bsy*c, bsx*c,
                            csz*c, csy*c, csx*c);

        image->ncells = 1;

        return 0;
}


/* write .spt file for mosflm */
/* need to sort mosflm peaks by intensity... */
struct sptline {
	double x; /* x coordinate of peak */
	double y; /* y coordinate of peak */
	double h; /* height of peak */
	double s; /* sigma of peak */
};


static int compare_vals(const void *ap, const void *bp)
{
	const struct sptline a = *(struct sptline *)ap;
	const struct sptline b = *(struct sptline *)bp;

	if ( a.h < b.h ) return 1;
	if ( a.h > b.h ) return -1;
	return 0;
}


static void write_spt(struct image *image)
{
	FILE *fh;
	int i;
	char filename[1024];
	double fclen=67.8;  /* fake camera length in mm */
	double fpix=0.075;  /* fake pixel size in mm */
	double pix;
	double height=100;
	double sigma=1;
	int nPeaks = image_feature_count(image->features);

	snprintf(filename, 1023, "xfel-%i.spt", image->id);

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file xfel.spt\n");
		return;
	}

	fprintf(fh, "%10d %10d %10.8f %10.6f %10.6f\n", 1, 1, fpix, 1.0, 0.0);
	fprintf(fh, "%10d %10d\n", 1, 1);
	fprintf(fh, "%10.5f %10.5f\n", 0.0, 0.0);

	struct sptline *sptlines;
	sptlines = malloc(sizeof(struct sptline)*nPeaks);

	for ( i=0; i<nPeaks; i++ ) {

		struct imagefeature *f;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		struct panel *pan;
		pan = find_panel(image->det,f->x,f->y);
		if ( pan == NULL ) continue;

		pix = 1000/pan->res; /* pixel size in mm */
		height = f->intensity;

		sptlines[i].x = (f->y - pan->cy)*pix*fclen/pan->clen/1000;
		sptlines[i].y = -(f->x - pan->cx)*pix*fclen/pan->clen/1000;
		sptlines[i].h = height;
		sptlines[i].s = sigma;

	}

	qsort(sptlines, nPeaks, sizeof(struct sptline), compare_vals);

	for ( i=0; i<nPeaks; i++ ) {

		fprintf(fh, "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
		        sptlines[i].x, sptlines[i].y,
		        0.0, 0.0,
		        sptlines[i].h, sptlines[i].s);

	}

	fprintf(fh,"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
	           -999.0,-999.0,-999.0,-999.0,-999.0,-999.0);
	fclose(fh);
}


/* Write a dummy 1x1 pixel image file for mosflm.  Without post refinement,
 * mosflm will ignore this, but it must be present. */
static void write_img(struct image *image)
{
	FILE *fh;
	char filename[1024];
	unsigned short int * intimage;

	intimage = malloc(sizeof(unsigned short int));
	intimage[0] = 1;

	snprintf(filename, 1023, "xfel-%i_001.img", image->id);

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		return;
	}

	fprintf(fh,"{\nHEADER_BYTES=512;\n");
	fprintf(fh,"BYTE_ORDER=little_endian;\n");
	fprintf(fh,"TYPE=unsigned_short;\n");
	fprintf(fh,"DIM=2;\n");
	fprintf(fh,"SIZE1=1;\n");
	fprintf(fh,"SIZE2=1;\n");
	fprintf(fh,"}\n");
	while ( ftell(fh) < 512 ) { fprintf(fh," "); };
	fwrite(intimage,sizeof(unsigned short int),1,fh);
	fclose(fh);
}


void run_mosflm(struct image *image, UnitCell *cell)
{
	int i,j;
	char mos_cmd[1024];
	char symm[64];
	const char *sg;
	double a,b,c,alpha,beta,gamma;
	double wavelength; /* angstrom */
	char newmatfile[128];
	int fail;
	pid_t pid;
	int r;

	write_spt(image);
	write_img(image); /* dummy image */

	wavelength = image->lambda*1e10;
	cell_get_parameters(cell, &a, &b, &c, &alpha, &beta, &gamma);
	sg = cell_get_spacegroup(cell);
	sprintf(newmatfile,"xfel-%i.newmat",image->id);

	/* need to remove white space from spacegroup... */
	j = 0;
	for(i = 0; i < strlen(sg);i++)
	{
		if (sg[i] != ' ') {
			symm[j] = sg[i];
			j++;
		}
	}
	symm[j] = '\0';


	/* build a script to run mosflm */
	sprintf(mos_cmd,"%s","ipmosflm << eof-mosflm > /dev/null\n");
	sprintf(mos_cmd,"%s%s",mos_cmd,
	                     "DETECTOR ROTATION HORIZONTAL ANTICLOCKWISE"
	                     " ORIGIN LL FAST HORIZONTAL RECTANGULAR\n");
	sprintf(mos_cmd,"%sCELL %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
	                   mos_cmd,
	                   a*1e10,b*1e10,c*1e10,
	                   rad2deg(alpha),rad2deg(beta),rad2deg(gamma));
	sprintf(mos_cmd,"%sSYMM %s\n",mos_cmd,symm);
	sprintf(mos_cmd,"%sDISTANCE %8.4f\n",mos_cmd,67.8);
	sprintf(mos_cmd,"%sBEAM %8.4f %8.4f\n",mos_cmd,0.0,0.0);
	sprintf(mos_cmd,"%sWAVELENGTH %10.5f\n",mos_cmd,wavelength);
	sprintf(mos_cmd,"%sNEWMAT %s\n",mos_cmd,newmatfile);
	sprintf(mos_cmd,"%sIMAGE xfel-%i_001.img phi 0 0\n",mos_cmd,image->id);
	sprintf(mos_cmd,"%sAUTOINDEX DPS FILE xfel-%i.spt IMAGE 1\n",
	                      mos_cmd,image->id);
	sprintf(mos_cmd,"%sGO\n",mos_cmd);
	sprintf(mos_cmd,"%s%s",mos_cmd,"eof-mosflm\n");

	/* remove the previous NEWMAT file prior to running mosflm */
	remove(newmatfile);

	/* Run the mosflm script */
	pid = fork();
	if ( !( (pid != 0) && (pid != -1) ) ) {
		if ( pid == -1 ) {
			ERROR("fork() failed.\n");
		} else {

			/* Forked successfully, child process */
			if ( system(mos_cmd) ) {
				ERROR("MOSFLM execution failed.\n");
				exit(0);
				return;
			}

			exit(0);

		}
	}
	waitpid(pid, &r, 0);

	/* Read the mosflm NEWMAT file and set cell candidate */
	/* Existence of this file means possible success. Pretty shady. */
	fail = read_newmat(newmatfile,image);
	if (fail) {
		printf("Failed to read mosflm NEWMAT file.\n");
		return;
	}
}
