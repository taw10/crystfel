/*
 * gui_project.h
 *
 * GUI project persistence
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

#ifndef GUI_PROJECT_H
#define GUI_PROJECT_H

enum match_type_id
{
	 MATCH_EVERYTHING,
	 MATCH_CHEETAH_LCLS_H5,
	 MATCH_CHEETAH_CXI,
	 MATCH_CBF,
	 MATCH_CBFGZ,
};

struct crystfelproject;


#include "crystfel_gui.h"


extern enum match_type_id decode_matchtype(const char *type_id);

extern int load_project(struct crystfelproject *proj);

extern void add_file_to_project(struct crystfelproject *proj,
                                const char *filename,
                                const char *event);

extern void clear_project_files(struct crystfelproject *proj);

#endif
