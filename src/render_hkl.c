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
#include <cairo.h>
#include <cairo-pdf.h>

#include "utils.h"
#include "reflections.h"


#define MAX_PROC (256)


static void show_help(const char *s)
{
	printf("Syntax: %s [options] <file.hkl>\n\n", s);
	printf(
"Render intensity lists using POV-ray.\n"
"\n"
"  -h, --help       Display this help message.\n"
"      --povray     Render a 3D animation using POV-ray.\n"
"      --zone-axis  Render a 2D zone axis pattern.\n"
"  -j <n>           Run <n> instances of POV-ray in parallel.\n"
"  -p, --pdb=<file> PDB file from which to get the unit cell.\n"
);
}


static void render_za(UnitCell *cell, double *ref, unsigned int *c)
{
	cairo_surface_t *surface;
	cairo_t *dctx;
	double max_u, max_v, max_res, max_intensity, scale;
	double sep_u, sep_v, max_r;
	double as, bs, theta;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	signed int h, k;
	float wh, ht;

	wh = 1024;
	ht = 1024;

	surface = cairo_pdf_surface_create("za.pdf", wh, ht);

	if ( cairo_surface_status(surface) != CAIRO_STATUS_SUCCESS ) {
		fprintf(stderr, "Couldn't create Cairo surface\n");
		cairo_surface_destroy(surface);
		return;
	}

	dctx = cairo_create(surface);

	/* Black background */
	cairo_rectangle(dctx, 0.0, 0.0, wh, ht);
	cairo_set_source_rgb(dctx, 0.0, 0.0, 0.0);
	cairo_fill(dctx);

	max_u = 0.0;  max_v = 0.0;  max_intensity = 0.0;
	max_res = 0.0;

	/* Work out reciprocal lattice spacings and angles for this cut */
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	theta = angle_between(asx, asy, asz, bsx, bsy, bsz);
	as = modulus(asx, asy, asz) / 1e9;
	bs = modulus(bsx, bsy, bsz) / 1e9;
	STATUS("theta=%f\n", rad2deg(theta));

	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {

		double u, v, intensity, res;
		int ct;

		ct = lookup_count(c, h, k, 0);
		if ( ct < 1 ) continue;

		intensity = lookup_intensity(ref, h, k, 0) / (float)ct;

		res = resolution(cell, h, k, 0);
		if ( res > max_res ) max_res = res;

		if ( intensity != 0 ) {
			u =  (double)h*as*sin(theta);
			v =  (double)h*as*cos(theta) + k*bs;
			if ( fabs(u) > fabs(max_u) ) max_u = fabs(u);
			if ( fabs(v) > fabs(max_v) ) max_v = fabs(v);
			if ( fabs(intensity) > fabs(max_intensity) )
						max_intensity = fabs(intensity);
		}

	}
	}

	max_res /= 1e9;
	max_u /= 0.5;
	max_v /= 0.5;
	printf("Maximum resolution is %f nm^-1\n", max_res);

	if ( max_intensity <= 0.0 ) {
		max_r = 4.0;
		goto out;
	}

	/* Choose whichever scaling factor gives the smallest value */
	scale = ((double)wh-50.0) / (2*max_u);
	if ( ((double)ht-50.0) / (2*max_v) < scale ) {
		scale = ((double)ht-50.0) / (2*max_v);
	}

	sep_u = as * scale * cos(theta);
	sep_v = bs * scale;
	max_r = ((sep_u < sep_v)?sep_u:sep_v);

	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {

		double u, v, intensity, val;
		int ct;

		ct = lookup_count(c, h, k, 0);
		if ( ct < 1 ) continue;

		intensity = lookup_intensity(ref, h, k, 0) / (float)ct;
		val = 10.0*intensity/max_intensity;

		u = (double)h*as*sin(theta);
		v = (double)h*as*cos(theta) + k*bs;

		cairo_arc(dctx, ((double)wh/2)+u*scale*2,
				((double)ht/2)+v*scale*2, max_r, 0, 2*M_PI);

		cairo_set_source_rgb(dctx, val, val, val);
		cairo_fill(dctx);

	}
	}

out:
	/* Centre marker */
	cairo_arc(dctx, (double)wh/2,
			(double)ht/2, max_r, 0, 2*M_PI);
	cairo_set_source_rgb(dctx, 1.0, 0.0, 0.0);
	cairo_fill(dctx);

	cairo_surface_finish(surface);
	cairo_destroy(dctx);
}


static void povray_render_animation(UnitCell *cell, double *ref,
                                    unsigned int *c, unsigned int nproc)
{
	FILE *fh;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	pid_t pids[MAX_PROC];
	float max;
	int i;
	signed int h, k, l;

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

	max = 0.5e6;
	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {
	for ( l=-INDMAX; l<INDMAX; l++ ) {

		float radius, x, y, z;
		int s;
		float val, p, r, g, b, trans;

		if ( !lookup_count(c, h, k, l) ) continue;

		val = lookup_intensity(ref, h, k, l);

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
}


int main(int argc, char *argv[])
{
	int c;
	UnitCell *cell;
	char *infile;
	double *ref;
	unsigned int *cts;
	int config_povray = 0;
	int config_zoneaxis = 0;
	unsigned int nproc = 1;
	char *pdb = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"povray",             0, &config_povray,      1},
		{"zone-axis",          0, &config_zoneaxis,    1},
		{"pdb",                1, NULL,               'p'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hj:p:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'j' :
			nproc = atoi(optarg);
			break;

		case 'p' :
			pdb = strdup(optarg);
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	if ( pdb == NULL ) {
		pdb = strdup("molecule.pdb");
	}

	if ( (nproc > MAX_PROC) || (nproc < 1) ) {
		ERROR("Number of processes is invalid.\n");
		return 1;
	}

	infile = argv[optind];

	cell = load_cell_from_pdb(pdb);
	if ( cell == NULL ) {
		ERROR("Couldn't load unit cell from %s\n", pdb);
		return 1;
	}
	cts = new_list_count();
	ref = read_reflections(infile, cts, NULL);
	if ( ref == NULL ) {
		ERROR("Couldn't open file '%s'\n", infile);
		return 1;
	}

	if ( config_povray ) {
		povray_render_animation(cell, ref, cts, nproc);
	} else if ( config_zoneaxis ) {
		render_za(cell, ref, cts);
	} else {
		ERROR("Try again with either --povray or --zone-axis.\n");
	}

	free(pdb);

	return 0;
}
