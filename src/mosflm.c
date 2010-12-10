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
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <assert.h>
#include <sys/ioctl.h>

#if HAVE_FORKPTY_LINUX
#include <pty.h>
#elif HAVE_FORKPTY_BSD
#include <util.h>
#endif


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

	/* Run the mosflm script */
	fail = system(mos_cmd);
	if (fail) { 
		ERROR("mosflm execution failed.\n");
		return; 
	}

	/* Read the mosflm NEWMAT file and set cell candidate */ 
	/* Existence of this file means possible success. Pretty shady. */
	fail = read_newmat(newmatfile,image);
	if (fail) {
		printf("Failed to read mosflm NEWMAT file.\n"); 
		return;
	}

	/* remove the mosflm NEWMAT file */
	//remove(newmatfile);
	
	return;
}

