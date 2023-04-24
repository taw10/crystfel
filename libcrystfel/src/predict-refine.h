/*
 * predict-refine.h
 *
 * Prediction refinement
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2019 Thomas White <taw@physics.org>
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

#ifndef PREDICT_REFINE_H
#define PREDICT_REFINE_H

#include "crystal.h"

struct image;

typedef void *Mille;

/**
 * \file predict-refine.h
 * Prediction refinement: refinement of indexing solutions before integration.
 */

extern Mille *crystfel_mille_new(const char *outFileName,
                                 int asBinary,
                                 int writeZero);
extern void crystfel_mille_free(Mille *m);

extern int refine_prediction(struct image *image, Crystal *cr, Mille *mille);
extern int refine_radius(Crystal *cr, struct image *image);


#endif	/* PREDICT_REFINE_H */
