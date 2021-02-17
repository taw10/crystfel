/*
 * gui_merge.h
 *
 * Merging via CrystFEL GUI
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

#ifndef GUI_MERGE_H
#define GUI_MERGE_H

#include <gtk/gtk.h>

#include "gui_project.h"

extern gint merge_sig(GtkWidget *widget,
                      struct crystfelproject *proj);

extern char **merging_command_line(const char *n_thread_str,
                                   struct gui_indexing_result *input,
                                   struct merging_params *params);

extern int write_merge_script(const char *filename,
                              struct gui_indexing_result *input,
                              const char *n_thread_str,
                              struct merging_params *params,
                              const char *out_hkl);

#endif
