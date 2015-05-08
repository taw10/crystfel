/*
 * render.c
 *
 * Render a high dynamic range buffer in some sensible way
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2012,2014 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <math.h>
#include <stdint.h>


#include "hdf5-file.h"
#include "render.h"
#include "peaks.h"
#include "filters.h"
#include "utils.h"


static void render_rgb(double val, double max,
                       double *rp, double *gp, double *bp)
{
	int s;
	double p;
	double r, g, b;

	s = val / (max/6);
	p = fmod(val, max/6.0);
	p /= (max/6.0);

	r = 0.0;  g = 0.0;  b = 0.0;

	if ( (val < 0.0) ) {
		p = fabs(val) / (max/6.0);
		*rp = 0.0;
		*gp = 0.5*p;
		*bp = 0.0;
		return;
	}
	if ( (val > max) ) {
		s = 6;
	}
	switch ( s ) {

		case 0 :	/* Black to blue */
		r = 0;  g = 0;  b = p;
		break;

		case 1 :	/* Blue to pink */
		r = p;  g = 0;  b = 1.0;
		break;

		case 2 :	/* Pink to red */
		r = 1.0;  g = 0;  b = (1.0-p)*1.0;
		break;

		case 3 :	/* Red to Orange */
		r = 1.0;  g = 0.5*p;  b = 0;
		break;

		case 4 :	/* Orange to Yellow */
		r = 1.0;  g = 0.5 + 0.5*p;  b = 0;
		break;

		case 5 :	/* Yellow to White */
		r = 1.0;  g = 1.0;  b = 1.0*p;
		break;

		case 6 :	/* Pixel has hit the maximum value */
		r = 1.0;  g = 1.0;  b = 1.0;
		break;

	}

	*rp = r;
	*gp = g;
	*bp = b;
}


static void render_geoptimiser(double val, double max,
                               double *rp, double *gp, double *bp)
{
	double r;
	double p;

	r = val/max;

	if ( val < 0.0 ) {
		*rp = 0.0;
		*gp = 0.0;
		*bp = 0.0;
		return;
	}

	if ( r >= 0.0 && r < 0.059 ) {
		p = (r-0.0)/(0.059-0.0);
		*rp = 0.0;
		*gp = 0.0;
		*bp = ((91.0/256.0)-0.0)*p;
		return;
	}

	if ( r >= 0.059 && r < 0.220 ) {
		p = (r-0.059)/(0.220-0.059);
		*rp = ((122.0/256.0)-0.0)*p;
		*gp = 0.0;
		*bp = ((227.0/256.0)-(91.0/256.0))*p+(91.0/256.0);
		return;
	}

	if ( r >= 0.220 && r < 0.376 ) {
		p = (r-0.220)/(0.376-0.220);
		*rp = ((195.0/256.0)-(122.0/256.0))*p+(122.0/256.0);
		*gp = 0.0;
		*bp = ((93.0/256.0)-(227.0/256.0))*p+(227.0/256.0);
		return;
	}

	if ( r >= 0.376 && r < 0.498 ) {
		p = (r-0.376)/(0.498-0.376);
		*rp = ((238.0/256.0)-(195.0/256.0))*p+(195.0/256.0);
		*gp = ((76.0/256.0)-0.0)*p;
		*bp = (0.0-(93.0/256.0))*p+(93.0/256.0);
		return;
	}

	if ( r >= 0.498 && r < 0.564 ) {
		p = (r-0.498)/(0.564-0.498);
		*rp = (1.0-(238.0/256.0))*p+(238.0/256.0);
		*gp = ((117.0/256.0)-(76.0/256.0))*p+(76.0/256.0);
		*bp = 0.0;
		return;
	}

	if ( r >= 0.564 && r < 0.815 ) {
		p = (r-0.564)/(0.815-0.564);
		*rp = 1.0;
		*gp = ((234.0/256.0)-(117.0/256.0))*p+(117.0/256.0);
		*bp = 0.0;
		return;
	}

	if ( r >= 0.815 && r < 1.0 ) {
		p = (r-0.815)/(1.0-0.815);
		*rp = 1.0;
		*gp = (1.0-(234.0/256.0))*p+(234.0/256.0);
		*bp = (1.0-0.0)*p;
		return;
	}

	if ( r >= 1.0 ) {
		*rp = 1.0; *gp = 1.0; *bp = 1.0;
		return;
	}
}


static void render_ratio(double val, double max,
                         double *rp, double *gp, double *bp)
{
	if ( val <= 1.0 ) {
		render_rgb(val, 2.0, rp, gp, bp);
	} else {
		/* Your homework is to simplify this expression */
		val = ((val-1.0)/(max-1.0)) * (max/2.0) + max/2.0;
		render_rgb(val, max, rp, gp, bp);
	}
}


static void render_mono(double val, double max,
                        double *rp, double *gp, double *bp)
{
	double p;
	p = val / max;
	if ( val < 0.0 ) p = 0.0;
	if ( val > max ) p = 1.0;
	*rp = p;
	*gp = p;
	*bp = p;
}


static void render_invmono(double val, double max,
                           double *rp, double *gp, double *bp)
{
	double p;
	p = val / max;
	p = 1.0 - p;
	if ( val < 0.0 ) p = 1.0;
	if ( val > max ) p = 0.0;
	*rp = p;
	*gp = p;
	*bp = p;
}


void render_scale(double val, double max, int scale,
                  double *rp, double *gp, double *bp)
{
	switch ( scale ) {

		case SCALE_COLOUR :
		render_rgb(val, max, rp, gp, bp);
		break;

		case SCALE_MONO :
		render_mono(val, max, rp, gp, bp);
		break;

		case SCALE_INVMONO :
		render_invmono(val, max, rp, gp, bp);
		break;

		case SCALE_RATIO :
		render_ratio(val, max, rp, gp, bp);
		break;

		case SCALE_GEOPTIMISER :
		render_geoptimiser(val, max, rp, gp, bp);
		break;
	}
}
