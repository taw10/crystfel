/*
 * crystfel-mille.h
 *
 * Interface to Millepede geometry refinement
 *
 * Copyright © 2023 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2023 Thomas White <taw@physics.org>
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

#ifndef CRYSTFEL_MILLE_H
#define CRYSTFEL_MILLE_H

typedef void *Mille;

#include "cell.h"
#include "image.h"
#include "predict-refine.h"
#include "geometry.h"

/**
 * \file crystfel-mille.h
 * Detector geometry refinement using Millepede
 */

extern Mille *crystfel_mille_new(const char *outFileName,
                                 int asBinary,
                                 int writeZero);

extern void crystfel_mille_free(Mille *m);

extern int mille_label(int hierarchy_level, int member_index, enum gparam param);

extern void write_mille(Mille *mille, int n, UnitCell *cell,
                        struct reflpeak *rps, struct image *image);

extern void crystfel_mille_delete_last_record(Mille *m);

extern void crystfel_mille_write_record(Mille *m);

#endif	/* CRYSTFEL_MILLE_H */
