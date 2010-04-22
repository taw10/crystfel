/*
 * render_hkl.c
 *
 * Draw pretty renderings of reflection lists
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
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
#include <sys/types.h>
#include <sys/wait.h>

#include "utils.h"
#include "reflections.h"


#define MAX_PROC (256)


static void show_help(const char *s)
{
	printf("Syntax: %s [options] <file.hkl>\n\n", s);
	printf(
"Render intensity lists using POV-ray.\n"
"\n"
"  -h, --help      Display this help message.\n"
"  -j <n>          Run <n> instances of POV-ray in parallel.\n"
"\n");
}


int main(int argc, char *argv[])
{
	int c;
	UnitCell *cell;
	char *infile;
	double *ref1;
	signed int h, k, l;
	unsigned int *c1;
	FILE *fh;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	pid_t pids[MAX_PROC];
	float max;
	int nproc = 1;
	int i;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hj:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' : {
			show_help(argv[0]);
			return 0;
		}

		case 'j' : {
			nproc = atoi(optarg);
			break;
		}

		case 0 : {
			break;
		}

		default : {
			return 1;
		}
		}

	}

	if ( (nproc > MAX_PROC) || (nproc < 1) ) {
		ERROR("Number of processes is invalid.\n");
		return 1;
	}

	infile = argv[optind];

	cell = load_cell_from_pdb("molecule.pdb");
	c1 = new_list_count();
	ref1 = read_reflections(infile, c1);
	if ( ref1 == NULL ) {
		ERROR("Couldn't open file '%s'\n", infile);
		return 1;
	}

	fh = fopen("render.pov", "w");
	fprintf(fh, "/* POV-Ray scene written by Synth3D */\n\n");
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
	fprintf(fh, "#if ( (clock >= 125) & (clock <= 250) )\n");
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

	fprintf(fh, "#declare WC = (720/19);\n");


	fprintf(fh, "#declare TRANS = \n");
	fprintf(fh, "transform {\n");

	/* First half */

	/* Acceleration */
	fprintf(fh, "#if ( clock <= 24 )\n"
	"rotate <0, 0, 0.5*WC*(clock/25)*(clock/25)>\n"
	"#end\n"

	/* Cruise */
	"#if ( (clock >= 25) & (clock <= 224) )\n"
	"rotate <0, 0, (WC/2)+WC*((clock-25)/25)>\n"
	"#end\n"

	/* Overlap */

	/* Deceleration */
	"#if ( (clock >= 225) & (clock <= 274) )\n"
	"rotate <0, 0, 360-WC + WC*((clock-225)/25) "
	" - 0.5*(WC/2)*((clock-225)/25)*((clock-225)/25) >\n"
	"#end\n"

	/* Acceleration */
	"#if ( (clock >= 225) & (clock <= 274) )\n"
	"rotate <0.5*(WC/2)*((clock-225)/25)*((clock-225)/25), 0, 0>\n"
	"#end\n"

	/* Second half */

	/* Cruise */
	"#if ( (clock >= 275) & (clock <= 474) )\n"
	"rotate <WC+WC*((clock-275)/25), 0, 0>\n"
	"#end\n"

	/* Deceleration */
	"#if ( (clock >= 475) & (clock <= 499) )\n"
	"rotate <360-(WC/2)+ WC*((clock-475)/25) "
	" - 0.5*WC*((clock-475)/25)*((clock-475)/25), 0, 0 >\n"
	"#end\n");

	fprintf(fh, "}\n");

	max = 0.5e6;
	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {
	for ( l=-INDMAX; l<INDMAX; l++ ) {

		float radius, x, y, z;
		int s;
		float val, p, r, g, b, trans;

		if ( !lookup_count(c1, h, k, l) ) continue;

		val = lookup_intensity(ref1, h, k, l);

		val = max-val;

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
			case 0 : {  /* Black to blue */
				r = 0.0;  g = 0.0;  b = p;
				break;
			}
			case 1 : {  /* Blue to green */
				r = 0.0;  g = p;  b = 1.0-p;
				break;
			}
			case 2 : {  /* Green to red */
				r =p;  g = 1.0-p;  b = 0.0;
				break;
			}
			case 3 : {  /* Red to Orange */
				r = 1.0;  g = 0.5*p;  b = 0.0;
				break;
			}
			case 4 : {  /* Orange to Yellow */
				r = 1.0;  g = 0.5 + 0.5*p;  b = 0.0;
				break;
			}
			case 5 : {  /* Yellow to White */
				r = 1.0;  g = 1.0;  b = 1.0*p;
				break;
			}
			case 6 : {  /* Pixel has hit the maximum value */
				r = 1.0;  g = 1.0;  b = 1.0;
				break;
			}
		}

		val = max-val;

		if ( val <= 0.0 ) continue;
		radius = 0.1 * sqrt(sqrt(val))/1e2;
		radius -= 0.005;
		if ( radius > 0.03 ) radius = 0.03;
		if ( radius <= 0.0 ) continue;
		trans = (0.03-radius)/0.03;
		radius += 0.002;

		x = asx*h + bsx*k + csx*l;
		y = asy*h + bsy*k + csy*l;
		z = asz*h + bsz*k + csz*l;

		fprintf(fh, "sphere { <%.5f, %.5f, %.5f>, %.5f "
		            "texture{pigment{color rgb <%f, %f, %f>"
		            " transmit %f} "
		            "finish { reflection 0.1 } } \n"
			    "transform { TRANS }\n"
			    "}\n",
			x/1e9, y/1e9, z/1e9, radius, r, g, b, trans);

	}
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
				int nf, xf, nsec;

				nsec = 500 / nproc;
				nf = nsec * i;
				xf = (nsec * i + nsec)-1;

				snprintf(minf, 255, "+SF%i", nf);
				snprintf(maxf, 255, "+EF%i", xf);

				/* Forked successfully, child process */
				execlp("povray", "", "+W1024", "+H768",
				       "+Irender.pov", "+Orender.png",
				       "+KFI0", "+KFF499", "+KI0", "+KF499",
				       minf, maxf, NULL);

			}
		} /* else start the next one */
	}

	for ( i=0; i<nproc; i++ ) {
		int r;
		waitpid(pids[i], &r, 0);
	}

	return 0;
}
