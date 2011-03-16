/*
 * povray.c
 *
 * Invoke POV-ray
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "cell.h"
#include "utils.h"
#include "symmetry.h"
#include "render_hkl.h"
#include "povray.h"


#define MAX_PROC (256)


int povray_render_animation(UnitCell *cell, RefList *list, unsigned int nproc,
                            const char *sym, int wght, double boost)
{
	FILE *fh;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	pid_t pids[MAX_PROC];
	float max;
	int i;
	Reflection *refl;
	RefListIterator *iter;

	if ( (nproc > MAX_PROC) || (nproc < 1) ) {
		ERROR("Number of processes must be a number between 1 and %i\n",
		      MAX_PROC);
		return 1;
	}

	fh = fopen("render.pov", "w");
	fprintf(fh, "/* POV-Ray scene written by CrystFEL */\n\n");
	fprintf(fh, "#include \"colors.inc\"\n");
	fprintf(fh, "#include \"textures.inc\"\n\n");
	fprintf(fh, "global_settings {\n");
	fprintf(fh, "	assumed_gamma 1.0\n");
	fprintf(fh, "	ambient_light 5.0\n");
	fprintf(fh, "}\n\n");

	/* First quarter */
	fprintf(fh, "#if ( (clock >= 0) & (clock <= 124) )\n");
	fprintf(fh, "camera { location <0.0, -3.0, 0.0>"
	            " sky z direction 1.1*y\n"
	            " right -x*(image_width/image_height)\n"
	            " look_at <0.0, 0.0, 0.0> }\n\n");
	fprintf(fh, "#end\n");

	/* Second quarter */
	fprintf(fh, "#if ( (clock >= 125) & (clock <= 249) )\n");
	fprintf(fh, "camera { location <0.0,"
	            " -(2.0+cos(radians((clock-125)*(180/125)))), 0.0>"
	            " sky z direction 1.1*y\n"
	            " right -x*(image_width/image_height)\n"
	            " look_at <0.0, 0.0, 0.0> }\n\n");
	fprintf(fh, "#end\n");

	/* Third quarter */
	fprintf(fh, "#if ( (clock >= 250) & (clock <= 374) )\n");
	fprintf(fh, "camera { location <0.0, -1.0, 0.0>"
	            " sky z direction 1.1*y\n"
	            " right -x*(image_width/image_height)\n"
	            " look_at <0.0, 0.0, 0.0> }\n\n");
	fprintf(fh, "#end\n");

	/* Fourth quarter */
	fprintf(fh, "#if ( (clock >= 375) & (clock <= 500) )\n");
	fprintf(fh, "camera { location <0.0,"
	            " -(2.0+cos(radians((clock-375)*(180/125)+180))), 0.0>"
	            " sky z direction 1.1*y\n"
	            " right -x*(image_width/image_height)\n"
	            " look_at <0.0, 0.0, 0.0> }\n\n");
	fprintf(fh, "#end\n");

	fprintf(fh, "light_source { <-3.0 -3.0 3.0> White }\n");
	fprintf(fh, "light_source { <+3.0 -3.0 3.0> White }\n");
	fprintf(fh, "light_source { <0.0, -3.0, 0.0> 2*White }\n");
	fprintf(fh, "plane {z,-2.0 pigment { rgb <0.0, 0.0, 0.1> } }\n");
	fprintf(fh, "plane {-z,-2.0 pigment { rgb <0.0, 0.0, 0.05> } }\n\n");

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	fprintf(fh, "#declare WCA = (720/19);\n");
	fprintf(fh, "#declare WCL = (360/8.5);\n");
	fprintf(fh, "#declare TA = (4.875);\n");
	fprintf(fh, "#declare TB = (1.125);\n");

	fprintf(fh, "#declare TRANS = \n");
	fprintf(fh, "transform {\n");

	/* First half */

	/* Acceleration */
	fprintf(fh, "#if ( clock <= 24 )\n"
	"rotate <0, 0, 0.5*WCA*(clock/25)*(clock/25)>\n"
	"#end\n"

	/* Cruise */
	"#if ( (clock >= 25) & (clock <= 224) )\n"
	"rotate <0, 0, (WCA/2)+WCA*((clock-25)/25)>\n"
	"#end\n"

	/* Overlap */

	/* Deceleration */
	"#if ( (clock >= 225) & (clock <= 274) )\n"
	"rotate <0, 0, 360-WCA + WCA*((clock-225)/25) "
	" - 0.5*(WCA/2)*((clock-225)/25)*((clock-225)/25) >\n"
	"#end\n"

	/* Acceleration */
	"#if ( (clock >= 225) & (clock <= 274) )\n"
	"rotate <0.5*(WCL/2)*((clock-225)/25)*((clock-225)/25), 0, 0>\n"
	"#end\n"

	/* Second half */

	/* Cruise */
	"#if ( (clock >= 275) & (clock <= 396) )\n"
	"rotate <WCL + WCL*((clock-275)/25), 0, 0>\n"
	"#end\n"

	/* Deceleration to pause */
	"#if ( (clock >= 397) & (clock <= 421) )\n"
	"rotate <(1+TA)*WCL+ WCL*((clock-397)/25) "
	" - 0.5*WCL*((clock-397)/25)*((clock-397)/25), 0, 0 >\n"
	"#end\n"

	/* Acceleration after pause */
	"#if ( (clock >= 422) & (clock <= 446) )\n"
	"rotate <(1.5+TA)*WCL"
	" + 0.5*WCL*((clock-422)/25)*((clock-422)/25), 0, 0>\n"
	"#end\n"

	/* Final Cruise */
	"#if ( (clock >= 447) & (clock <= 474) )\n"
	"rotate <(2+TA)*WCL + WCL*((clock-447)/25), 0, 0>\n"
	"#end\n"

	/* Final Deceleration */
	"#if ( (clock >= 475) & (clock <= 499) )\n"
	"rotate <(2+TA+TB)*WCL + WCL*((clock-475)/25) "
	" - 0.5*WCL*((clock-475)/25)*((clock-475)/25), 0, 0 >\n"
	"#end\n");

	fprintf(fh, "}\n");

	max = 0.0;
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		float val;
		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		switch ( wght ) {
		case WGHT_I :
			val = get_intensity(refl);
			break;
		case WGHT_SQRTI :
			val = get_intensity(refl);
			val = (val>0.0) ? sqrt(val) : 0.0;
			break;
		case WGHT_COUNTS :
			val = get_redundancy(refl);
			val /= (float)num_equivs(h, k, l, sym);
			break;
		case WGHT_RAWCOUNTS :
			val = get_redundancy(refl);
			break;
		default :
			ERROR("Invalid weighting.\n");
			abort();
		}

		if ( val > max ) max = val;

	}
	max /= boost;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;float radius;
		int s;
		float val, p, r, g, b, trans;
		int j;

		get_indices(refl, &h, &k, &l);

		switch ( wght ) {
		case WGHT_I :
			val = get_intensity(refl);
			break;
		case WGHT_SQRTI :
			val = get_intensity(refl);
			val = (val>0.0) ? sqrt(val) : 0.0;
			break;
		case WGHT_COUNTS :
			val = get_redundancy(refl);
			val /= (float)num_equivs(h, k, l, sym);
			break;
		case WGHT_RAWCOUNTS :
			val = get_redundancy(refl);
			break;
		default :
			ERROR("Invalid weighting.\n");
			abort();
		}

		s = val / (max/6);
		p = fmod(val, max/6);
		p /= (max/6);

		r = 0;	g = 0;	b = 0;

		if ( (val < 0.0) ) {
			s = 0;
			p = 1.0;
		}
		if ( (val > max) ) {
			s = 6;
		}
		switch ( s ) {
		case 0 :   /* Black to blue */
			r = 0.0;  g = 0.0;  b = p;
			break;
		case 1 :   /* Blue to green */
			r = 0.0;  g = p;  b = 1.0-p;
			break;
		case 2 :   /* Green to red */
			r =p;  g = 1.0-p;  b = 0.0;
			break;
		case 3 :   /* Red to Orange */
			r = 1.0;  g = 0.5*p;  b = 0.0;
			break;
		case 4 :   /* Orange to Yellow */
			r = 1.0;  g = 0.5 + 0.5*p;  b = 0.0;
			break;
		case 5 :   /* Yellow to White */
			r = 1.0;  g = 1.0;  b = 1.0*p;
			break;
		case 6 :   /* Pixel has hit the maximum value */
			r = 1.0;  g = 1.0;  b = 1.0;
			break;
		}

		if ( val <= 0.0 ) continue;
		radius = 0.01 * pow(val, 0.25)/pow(max, 0.25);
		trans = 1.0-(val/max);

		/* For each equivalent */
		for ( j=0; j<num_equivs(h, k, l, sym); j++ ) {

			signed int he, ke, le;
			float x, y, z;

			get_equiv(h, k, l, &he, &ke, &le, sym, j);

			x = asx*he + bsx*ke + csx*le;
			y = asy*he + bsy*ke + csy*le;
			z = asz*he + bsz*ke + csz*le;

			fprintf(fh, "sphere { <%.5f, %.5f, %.5f>, %.5f "
			            "texture{pigment{color rgb <%f, %f, %f>"
			            " transmit %f} "
			            "finish { reflection %f } } \n"
			            "transform { TRANS }\n"
			            "}\n",
			            x/1e9, y/1e9, z/1e9, radius,
			            r, g, b, trans, 0.1*(1.0-trans));

		}

	}

	fprintf(fh, "\n");
	fclose(fh);

	for ( i=0; i<nproc; i++ ) {

		pids[i] = fork();
		if ( !( (pids[i] != 0) && (pids[i] != -1) ) ) {
			if ( pids[i] == -1 ) {
				ERROR("fork() failed.\n");
			} else {

				char minf[256];
				char maxf[256];
				float nf, xf, nsec;

				nsec = 500.0 / (float)nproc;
				nf = nsec * (float)i;
				xf = (nsec * (float)i + nsec) - 1.0;

				snprintf(minf, 255, "+SF%i", (int)nf);
				snprintf(maxf, 255, "+EF%i", (int)xf);

				/* Forked successfully, child process */
				execlp("povray", "", "+W1024", "+H768",
				       "+Irender.pov", "+Orender.png",
				       "+KFI0", "+KFF499", "+KI0", "+KF499",
				       minf, maxf, "-D", NULL);

			}
		} /* else start the next one */
	}

	for ( i=0; i<nproc; i++ ) {
		int r;
		waitpid(pids[i], &r, 0);
	}

	return 0;
}
