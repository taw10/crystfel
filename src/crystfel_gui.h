/*
 * crystfel_gui.h
 *
 * CrystFEL's main graphical user interface
 *
 * Copyright © 2020 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Thomas White <taw@physics.org>
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

#ifndef CRYSTFEL_GUI_H
#define CRYSTFEL_GUI_H

#include "gui_project.h"

struct crystfel_backend {
	const char *name;
	int (*run_unitcell)(struct crystfelproject *proj,
	                    const char *algo);
	void (*cancel)(struct crystfelproject *proj);
	void (*init)(struct crystfelproject *proj);
	void (*shutdown)(struct crystfelproject *proj);
};


extern void remove_infobar(struct crystfelproject *proj);

extern GtkWidget *create_infobar(struct crystfelproject *proj, const char *task,
                                 const char *extra_button,
                                 void (*cbfunc)(struct crystfelproject *proj));

#endif
