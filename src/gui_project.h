/*
 * gui_project.h
 *
 * GUI project persistence
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

#ifndef GUI_PROJECT_H
#define GUI_PROJECT_H

#include <gtk/gtk.h>

#include <peaks.h>
#include <stream.h>

#define MAX_RUNNING_TASKS (16)

enum match_type_id
{
	 MATCH_EVERYTHING,
	 MATCH_H5,
	 MATCH_CHEETAH_LCLS_H5,
	 MATCH_CHEETAH_CXI,
	 MATCH_CBF,
	 MATCH_CBFGZ,
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

	/* PinkIndexer-specific */
	int pinkindexer_cpeaks;
	int pinkindexer_use_max_res;
	double pinkindexer_max_res;
	int pinkindexer_angle_density;
	int pinkindexer_refinement_type;
	double pinkindexer_tolerance;
	int pinkindexer_use_refl_radius;
	double pinkindexer_refl_radius;
	double pinkindexer_max_imbalance;

	/* Integration */
	char *integration_method;
	int overpredict;
	float push_res;
	float ir_inn;
	float ir_mid;
	float ir_out;
	float fix_profile_radius;
	int use_fix_profile_radius;
	float fix_divergence;

	/* Stream output */
	int exclude_nonhits;
	int exclude_peaks;
	int exclude_refls;
	char **metadata_to_copy;
	int n_metadata;
	int millepede;
	int max_mille_level;
};

struct merging_params {

	char *model;   /* "process_hkl" in addition to xsphere/unity etc */
	char *symmetry;
	int scale;
	int bscale;
	int postref;
	int niter;
	char *polarisation;
	int deltacchalf;
	int min_measurements;
	float max_adu;
	char *custom_split;
	int pr_logs;
	char *twin_sym;
	float min_res;
	float push_res;
};

struct ambi_params {
	int use_res;
	double res_min;  /* Angstroms */
	double res_max;  /* Angstroms */
	int niter;
	int use_ncorr;
	int ncorr;
	char *sym;
	char *source_sym;
	int use_operator;
	char *operator;
};

struct gui_indexing_result
{
	char *name;

	int n_streams;
	char **streams;
	StreamIndex **indices;
	int need_rescan;
};

struct gui_merge_result
{
	char *name;
	char *indexing_result_name;  /* Indexing result this was derived from */
	char *hkl;    /* Complete merged data */
	char *hkl1;   /* First half-split */
	char *hkl2;   /* Second half-split */
};

struct crystfelproject;

struct crystfel_backend {

	const char *name;
	const char *friendly_name;

	/* Called to ask the backend to cancel the job */
	void (*cancel_task)(void *job_priv);

	/* Called to ask the backend to free any resources in the job record */
	void (*free_task)(void *job_priv);

	/* Called to get the status of a task */
	int (*task_status)(void *job_priv,
	                   int *running,
	                   float *fraction_complete);

	/* ....................... Indexing ........................ */

	/* Backend should provide a GTK widget to set options */
	GtkWidget *(*make_indexing_parameters_widget)(void *opts_priv);

	/* Called to ask the backend to start indexing frames.
	 * It should return a void pointer representing this job */
	void *(*run_indexing)(const char *job_title,
	                      const char *job_notes,
	                      struct crystfelproject *proj,
	                      void *opts_priv,
	                      double wavelength_estimate,
	                      double clen_estimate);

	/* Called to ask the backend to write its indexing options */
	void (*write_indexing_opts)(void *opts_priv, FILE *fh);

	/* Called when reading a project from file */
	void (*read_indexing_opt)(void *opts_priv,
	                          const char *key,
	                          const char *val);

	/* Backend should store options for indexing here */
	void *indexing_opts_priv;

	/* ....................... Merging ........................ */

	/* Backend should provide a GTK widget to set options */
	GtkWidget *(*make_merging_parameters_widget)(void *opts_priv);

	/* Called to ask the backend to start merging data.
	 * It should return a void pointer representing this job */
	void *(*run_merging)(const char *job_title,
	                     const char *job_notes,
	                     struct crystfelproject *proj,
	                     struct gui_indexing_result *input,
	                     void *opts_priv);

	/* Called to ask the backend to write its merging options */
	void (*write_merging_opts)(void *opts_priv, FILE *fh);

	/* Called when reading a project from file */
	void (*read_merging_opt)(void *opts_priv,
	                         const char *key,
	                         const char *val);

	/* Backend should store options for merging here */
	void *merging_opts_priv;

	/* .................. Indexing ambiguity .................. */

	/* Backend should provide a GTK widget to set options */
	GtkWidget *(*make_ambi_parameters_widget)(void *opts_priv);

	/* Called to ask the backend to start resolving indexing ambiguity.
	 * It should return a void pointer representing this job */
	void *(*run_ambi)(const char *job_title,
	                  const char *job_notes,
	                  struct crystfelproject *proj,
	                  struct gui_indexing_result *input,
	                  void *opts_priv);

	/* Called to ask the backend to write its ambigator options */
	void (*write_ambi_opts)(void *opts_priv, FILE *fh);

	/* Called when reading a project from file */
	void (*read_ambi_opt)(void *opts_priv,
	                      const char *key,
	                      const char *val);

	/* Backend should store options for ambigator here */
	void *ambi_opts_priv;
};

struct gui_task
{
	GtkWidget *info_bar;
	GtkWidget *cancel_button;
	GtkWidget *progress_bar;
	int running;
	struct crystfel_backend *backend;
	struct crystfelproject *proj;
	void *job_priv;
};


#define N_RANDOM_HISTORY (16)

struct crystfelproject {

	GtkWidget *window;
	GtkUIManager *ui;
	GtkActionGroup *action_group;

	GtkWidget *imageview;
	GtkWidget *colscale;
	GtkWidget *icons;      /* Drawing area for task icons */
	GtkWidget *report;     /* Text view at the bottom for messages */
	GtkWidget *main_vbox;
	GtkWidget *image_info;
	GtkWidget *results_combo;
	GtkWidget *next_button;
	GtkWidget *prev_button;
	GtkWidget *first_button;
	GtkWidget *last_button;
	int range_set;

	int unsaved;

	int cur_frame;
	struct image *cur_image;
	int random_history[N_RANDOM_HISTORY];
	int n_random_history;

	char *geom_filename;
	char *data_top_folder;   /* For convenience only.  Filenames in
	                          * 'filenames' list should be complete */
	enum match_type_id data_search_pattern;

	DataTemplate *dtempl;
	int n_frames;
	int max_frames;
	char **filenames;
	char **events;
	int show_centre;
	int resolution_rings;
	int rescan_on_change;

	char *unique_files[20];
	int n_unique_files;

	int show_peaks;
	struct peak_params peak_search_params;

	int show_refls;
	int label_refls;
	struct index_params indexing_params;
	int indexing_backend_selected;
	GtkWidget *indexing_opts;
	char *indexing_new_job_title;

	struct merging_params merging_params;
	int merging_backend_selected;
	GtkWidget *merging_opts;
	char *merging_new_job_title;

	GtkWidget *type_combo;
	GtkWidget *peak_vbox;     /* Box for peak search parameter widgets */
	GtkWidget *peak_params;   /* Peak search parameter widgets */
	struct peak_params original_params;

	/* All the backends available in this project */
	struct crystfel_backend *backends;
	int n_backends;

	GSList *tasks;

	struct gui_indexing_result *results;
	int n_results;

	struct gui_merge_result *merge_results;
	int n_merge_results;

	double fom_res_min;  /* Angstroms */
	double fom_res_max;  /* Angstroms */
	int fom_nbins;
	double fom_min_snr;
	int fom_min_meas;
	char *fom_cell_filename;

	double export_res_min;  /* Angstroms */
	double export_res_max;  /* Angstroms */

	char *ambi_new_job_title;
	int ambi_backend_selected;
	GtkWidget *ambi_opts;
	struct ambi_params ambi_params;
};

extern enum match_type_id decode_matchtype(const char *type_id);

extern int match_filename(const char *fn, enum match_type_id mt);

extern int load_project(struct crystfelproject *proj);

extern int default_project(struct crystfelproject *proj);

extern int save_project(struct crystfelproject *proj);

extern void add_file_to_project(struct crystfelproject *proj,
                                const char *filename,
                                const char *event);

extern void clear_project_files(struct crystfelproject *proj);
extern void clear_indexing_results(struct crystfelproject *proj);

extern int add_indexing_result(struct crystfelproject *proj,
                               const char *name,
                               char **streams,
                               int n_streams);

extern struct image *find_indexed_image(struct crystfelproject *proj,
                                        const char *result_name,
                                        const char *filename,
                                        const char *event,
                                        int permit_rescan);

extern void update_result_index(struct gui_indexing_result *result);

extern struct gui_indexing_result *find_indexing_result_by_name(struct crystfelproject *proj,
                                                                const char *name);

extern int add_merge_result(struct crystfelproject *proj, const char *name,
                            const char *indexing_result_name,
                            const char *hkl, const char *hkl1, const char *hkl2);

extern struct gui_merge_result *find_merge_result_by_name(struct crystfelproject *proj,
                                                          const char *name);

extern const char *selected_result(struct crystfelproject *proj);

#endif
