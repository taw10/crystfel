/*
 * gui_ambi.h
 *
 * Resolve indexing ambiguities via CrystFEL GUI
 *
 * Copyright © 2021 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2021 Thomas White <taw@physics.org>
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

#ifndef GUI_AMBI_H
#define GUI_AMBI_H

#include <gtk/gtk.h>

#include "gui_project.h"

extern gint ambi_sig(GtkWidget *widget,
                     struct crystfelproject *proj);

extern int write_ambigator_script(const char *filename,
                                  struct gui_indexing_result *input,
                                  const char *n_thread_str,
                                  struct ambi_params *params,
                                  const char *out_stream,
                                  const char *stdout_filename,
                                  const char *stderr_filename,
                                  const char *fg_filename,
                                  const char *intermediate_rel_filename,
                                  const char *harvest_filename,
                                  const char *prologue);

extern double read_ambigator_progress(char *logfile_str, int niter);

#endif
