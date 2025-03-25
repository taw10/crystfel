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

#include <gsl/gsl_matrix.h>

typedef struct mille Mille;

#include "cell.h"
#include "image.h"
#include "predict-refine.h"

/**
 * \file crystfel-mille.h
 * Detector geometry refinement using Millepede
 */

extern Mille *crystfel_mille_new(const char *outFileName);
extern Mille *crystfel_mille_new_fd(int fd);

extern void crystfel_mille_free(Mille *m);

extern int mille_label(int group_serial, enum gparam param);
extern enum gparam mille_unlabel(int n);

extern void write_mille(Mille *mille, int n, UnitCell *cell,
                        enum gparam *rvl, int nl,
                        struct reflpeak *rps, struct image *image,
                        int max_level, gsl_matrix **Minvs);

extern void crystfel_mille_delete_last_record(Mille *m);

extern void crystfel_mille_write_record(Mille *m);

#endif	/* CRYSTFEL_MILLE_H */
