/*
 * crystfel_gui.h
 *
 * CrystFEL's main graphical user interface
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

#ifndef CRYSTFEL_GUI_H
#define CRYSTFEL_GUI_H

#include "gui_project.h"


struct gui_job_notes_page
{
	GtkWidget *textview;
};


enum gui_job_type
{
	GUI_JOB_INDEXING,
	GUI_JOB_PROCESS_HKL,
	GUI_JOB_PROCESS_HKL_SCALE,
	GUI_JOB_PARTIALATOR,
	GUI_JOB_AMBIGATOR,
};

extern void error_box(struct crystfelproject *proj, const char *message);

extern void add_running_task(struct crystfelproject *proj,
                             const char *task_desc,
                             struct crystfel_backend *backend,
                             void *job_priv);

extern void update_imageview(struct crystfelproject *proj);
extern void select_result(struct crystfelproject *proj,
                          const char *result_name);
extern const char *selected_result(struct crystfelproject *proj);

extern char *get_crystfel_path_str(void);

extern struct gui_job_notes_page *add_job_notes_page(GtkWidget *notebook);

extern GFile *make_job_folder(const char *job_title, const char *job_notes);

extern void force_peaks_on(struct crystfelproject *proj);
extern void force_refls_on(struct crystfelproject *proj);

extern char *make_new_job_title(const char *old_title);

extern char *relative_to_cwd(GFile *workdir, const char *filename);

#endif
