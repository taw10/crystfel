/*
 * post-refinement.h
 *
 * Post refinement
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2015 Thomas White <taw@physics.org>
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

#ifndef POST_REFINEMENT_H
#define POST_REFINEMENT_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>

#include "image.h"
#include "utils.h"
#include "crystal.h"
#include "geometry.h"


struct prdata
{
	int refined;
	int n_filtered;

	/* Before refinement */
	double initial_residual;
	double initial_free_residual;
	double initial_log_residual;
	double initial_free_log_residual;

	/* After refinement */
	double final_residual;
	double final_free_residual;
	double final_log_residual;
	double final_free_log_residual;
};

enum prflag
{
	PRFLAG_OK = 0,
	PRFLAG_FEWREFL = 16,
	PRFLAG_SOLVEFAIL = 17,
	PRFLAG_EARLY = 18,
	PRFLAG_CC = 19,
	PRFLAG_BIGB = 20,
};


extern const char *str_prflag(enum prflag flag);

extern void refine_all(Crystal **crystals, int n_crystals,
                       RefList *full, int nthreads, PartialityModel pmodel,
                       int no_scale, int no_pr, double max_B,
                       double *initial_residual,
                       double *initial_free_residual,
                       double *initial_log_residual,
                       double *initial_free_log_residual,
                       double *final_residual,
                       double *final_free_residual,
                       double *final_log_residual,
                       double *final_free_log_residual);

/* Exported so it can be poked by tests/pr_p_gradient_check */
extern double gradient(Crystal *cr, int k, Reflection *refl,
                       PartialityModel pmodel);

#endif	/* POST_REFINEMENT_H */
