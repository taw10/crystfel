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

#include <gtk/gtk.h>

#include <peaks.h>
#include <stream.h>

#define MAX_RUNNING_TASKS (16)

enum match_type_id
{
	 MATCH_EVERYTHING,
	 MATCH_CHEETAH_LCLS_H5,
	 MATCH_CHEETAH_CXI,
	 MATCH_CBF,
	 MATCH_CBFGZ,
};

struct peak_params {
	enum peak_search_method method;
	float threshold;                /* zaef, pf8 */
	float min_sq_gradient;          /* zaef */
	float min_snr;                  /* zaef, pf8 */
	int min_pix_count;              /* pf8 */
	int max_pix_count;              /* pf8 */
	int local_bg_radius;            /* pf8 */
	int min_res;                    /* pf8 */
	int max_res;                    /* pf8 */
	float min_snr_biggest_pix;      /* pf9 */
	float min_snr_peak_pix;         /* pf9 */
	float min_sig;                  /* pf9 */
	float min_peak_over_neighbour;  /* pf9 */
	float pk_inn;
	float pk_mid;
	float pk_out;
	int half_pixel_shift;           /* cxi, hdf5 */
	int revalidate;
};

struct index_params {

	/* Indexing */
	char *cell_file;
	char *indexing_methods;
	int multi;
	int no_refine;
	int no_retry;
	int no_peak_check;
	int no_cell_check;
	float tols[6];
	int min_peaks;

	/* Integration */
	char *integration_method;
	int overpredict;
	float push_res;
};

struct crystfel_backend {

	const char *name;
	const char *friendly_name;

	/* Backend should provide a GTK widget to set options */
	GtkWidget *(*make_indexing_parameters_widget)(void *opts_priv);

	/* Called to ask the backend to start indexing frames.
	 * It should return a void pointer representing this job */
	void *(*run_indexing)(const char *job_title,
	                      const char *job_notes,
	                      char **filenames,
	                      char **events,
	                      int n_frames,
	                      char *geom_filename,
	                      struct peak_params *peak_search_params,
	                      struct index_params *indexing_params,
	                      void *opts_priv);

	/* Called to ask the backend to cancel the job */
	void (*cancel_task)(void *job_priv);

	/* Called to get the status of a task */
	int (*task_status)(void *job_priv,
	                   int *running,
	                   float *fraction_complete);

	/* Called to ask the backend to write its indexing options */
	void (*write_indexing_opts)(void *opts_priv, FILE *fh);

	/* Called when reading a project from file */
	void (*read_indexing_opt)(void *opts_priv,
	                          const char *key,
	                          const char *val);

	/* Backend should store options for indexing here */
	void *indexing_opts_priv;

};

struct gui_task
{
	GtkWidget *info_bar;
	GtkWidget *cancel_button;
	GtkWidget *progress_bar;
	int running;
	struct crystfel_backend *backend;
	void *job_priv;
};

struct gui_result
{
	char *name;
	char **streams;
	int n_streams;
};

struct crystfelproject {

	GtkWidget *window;
	GtkUIManager *ui;
	GtkActionGroup *action_group;

	GtkWidget *imageview;
	GtkWidget *icons;      /* Drawing area for task icons */
	GtkWidget *report;     /* Text view at the bottom for messages */
	GtkWidget *main_vbox;
	GtkWidget *image_info;
	GtkWidget *results_combo;

	int unsaved;

	int cur_frame;
	struct image *cur_image;

	char *geom_filename;
	char *stream_filename;
	char *data_top_folder;   /* For convenience only.  Filenames in
	                          * 'filenames' list should be complete */
	enum match_type_id data_search_pattern;

	DataTemplate *dtempl;
	Stream *stream;
	int n_frames;
	int max_frames;
	char **filenames;
	char **events;

	int show_peaks;
	struct peak_params peak_search_params;

	int show_refls;
	struct index_params indexing_params;
	int indexing_backend_selected;
	GtkWidget *indexing_opts;
	char *indexing_new_job_title;

	GtkWidget *type_combo;
	GtkWidget *peak_vbox;     /* Box for peak search parameter widgets */
	GtkWidget *peak_params;   /* Peak search parameter widgets */
	struct peak_params original_params;

	/* All the backends available in this project */
	struct crystfel_backend *backends;
	int n_backends;

	struct gui_task tasks[MAX_RUNNING_TASKS];
	int n_running_tasks;

	struct gui_result *results;
	int n_results;
};

extern enum match_type_id decode_matchtype(const char *type_id);

extern int match_filename(const char *fn, enum match_type_id mt);

extern int load_project(struct crystfelproject *proj);

extern void default_project(struct crystfelproject *proj);

extern int save_project(struct crystfelproject *proj);

extern void add_file_to_project(struct crystfelproject *proj,
                                const char *filename,
                                const char *event);

extern void clear_project_files(struct crystfelproject *proj);

extern int add_result(struct crystfelproject *proj,
                      char *name,
                      char **streams,
                      int n_streams);

#endif
