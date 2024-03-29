/*
 * gui_index.h
 *
 * Indexing parts of GUI
 *
 * Copyright © 2020-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020-2021 Thomas White <taw@physics.org>
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

#ifndef GUI_INDEX_H
#define GUI_INDEX_H

#include <gtk/gtk.h>

#include "gui_project.h"

extern gint index_one_sig(GtkWidget *widget,
                          struct crystfelproject *proj);

extern gint index_all_sig(GtkWidget *widget,
                          struct crystfelproject *proj);

extern void cell_explorer_sig(struct crystfelproject *proj);

extern int read_number_processed(const char *filename);

extern int write_indexamajig_script(const char *script_filename,
                                    const char *geom_filename,
                                    const char *n_thread_str,
                                    const char *files_list,
                                    const char *stream_filename,
                                    const char *stdout_filename,
                                    const char *stderr_filename,
                                    const char *harvest_filename,
                                    const char *mille_filename,
                                    const char *serial_start,
                                    struct peak_params *peak_search_params,
                                    struct index_params *indexing_params,
                                    double wavelength_estimate,
                                    double clen_estimate,
                                    const char *prologue);

#endif
