/*
 * error_model.h
 *
 * Error modelling
 *
 * Copyright © 2026 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2026 Thomas White <taw@physics.org>
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

#ifndef ERROR_MODEL_H
#define ERROR_MODEL_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


/**
 * An ErrorModel describes a statistical model for the errors
 * in the merged data.
 **/
typedef struct error_model ErrorModel;

typedef enum {

	EMODEL_EQUIVS,   /**< Observed intensity spread (the "old" method). */
	EMODEL_EV11,     /**< Evans (2011) with sdFac, sdB and sdAdd. */
	EMODEL_EV06,     /**< Evans (2006), like EMODEL_EV11, but without sdB. */
	EMODEL_KH23,     /**< Khouchen et al. 2023. */

} ErrorModelType;

extern ErrorModel *error_model_new(ErrorModelType t);
extern void error_model_free(ErrorModel *emodel);
extern void refine_error_model(RefList *full, ErrorModel *emodel);
extern void normal_probability_plot(RefList *full, ErrorModel *emodel);
extern ErrorModelType parse_error_model(const char *str, int *err);
extern void print_error_model(ErrorModel *emodel);

#endif	/* ERROR_MODEL_H */
