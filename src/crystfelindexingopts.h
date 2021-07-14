/*
 * crystfelindexingopts.h
 *
 * A GTK widget to set indexing options
 *
 * Copyright © 2020-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2020 Thomas White <taw@physics.org>
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

#ifndef CRYSTFELINDEXINGOPTS_H
#define CRYSTFELINDEXINGOPTS_H

#include <gtk/gtk.h>
#include <glib-object.h>

#define CRYSTFEL_TYPE_INDEXING_OPTS (crystfel_indexing_opts_get_type())

#define CRYSTFEL_INDEXING_OPTS(obj) (G_TYPE_CHECK_INSTANCE_CAST((obj), \
                                 CRYSTFEL_TYPE_INDEXING_OPTS, CrystFELIndexingOpts))

#define CRYSTFEL_IS_INDEXING_OPTS(obj) (G_TYPE_CHECK_INSTANCE_TYPE((obj), \
                                    CRYSTFEL_TYPE_INDEXING_OPTS))

#define CRYSTFEL_INDEXING_OPTS_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST((obj), \
                                         CRYSTFEL_TYPE_INDEXING_OPTS, CrystFELIndexingOpts))

#define CRYSTFEL_IS_INDEXING_OPTS_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((obj), \
                                            CRYSTFEL_TYPE_INDEXING_OPTS))

#define CRYSTFEL_INDEXING_OPTS_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS((obj), \
                                           CRYSTFEL_TYPE_INDEXING_OPTS, CrystFELIndexingOpts))

struct _crystfelindexingopts
{
	GtkNotebook parent_instance;

	/*< private >*/
	int show_stream_opts;
	GtkWidget *stream_params;  /* Stream output page */

	char *cell_file;
	GtkWidget *cell_chooser;
	GtkWidget *auto_indm;
	GtkWidget *indm_chooser;
	GtkListStore *indm_store;
	GtkWidget *multi;
	GtkWidget *refine;
	GtkWidget *retry;
	GtkWidget *check_peaks;
	GtkWidget *check_cell;
	GtkWidget *tols[6]; /* frac (not %) and radians */
	GtkWidget *enable_hitfind;
	GtkWidget *ignore_fewer_peaks;

	GtkWidget *integration_combo;
	GtkWidget *centering;
	GtkWidget *overpredict;
	GtkWidget *limit_res;
	GtkWidget *push_res;
	GtkWidget *ir_inn;
	GtkWidget *ir_mid;
	GtkWidget *ir_out;
	GtkWidget *fix_profile_radius_p;
	GtkWidget *fix_profile_radius;
	GtkWidget *fix_divergence;

	GtkWidget *pinkindexer_cpeaks;
	GtkWidget *pinkindexer_use_max_res;
	GtkWidget *pinkindexer_max_res;
	GtkWidget *pinkindexer_angle_density;
	GtkWidget *pinkindexer_refinement_type;
	GtkWidget *pinkindexer_tolerance;
	GtkWidget *pinkindexer_use_refl_radius;
	GtkWidget *pinkindexer_refl_radius;
	GtkWidget *pinkindexer_max_imbalance;

	GtkWidget *exclude_nonhits;
	GtkWidget *no_peaks_in_stream;
	GtkWidget *no_refls_in_stream;
	GtkListStore *copy_metadata_store;
};

struct _crystfelindexingoptsclass
{
	GtkNotebookClass parent_class;
};

typedef struct _crystfelindexingopts CrystFELIndexingOpts;
typedef struct _crystfelindexingoptsclass CrystFELIndexingOptsClass;

extern GType crystfel_indexing_opts_get_type(void);
extern GtkWidget *crystfel_indexing_opts_new(void);

extern char *crystfel_indexing_opts_get_cell_file(CrystFELIndexingOpts *opts);
extern char *crystfel_indexing_opts_get_indexing_method_string(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_multi_lattice(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_refine(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_retry(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_peak_check(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_cell_check(CrystFELIndexingOpts *opts);
extern void crystfel_indexing_opts_get_tolerances(CrystFELIndexingOpts *opts,
                                                  float *tols); /* frac (not %) and rad */
extern int crystfel_indexing_opts_get_min_peaks(CrystFELIndexingOpts *opts);

extern char *crystfel_indexing_opts_get_integration_method_string(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_overpredict(CrystFELIndexingOpts *opts);
extern float crystfel_indexing_opts_get_push_res(CrystFELIndexingOpts *opts);
extern void crystfel_indexing_opts_get_integration_radii(CrystFELIndexingOpts *opts,
                                                         float *ir_inn,
                                                         float *ir_mid,
                                                         float *ir_out);
extern int crystfel_indexing_opts_get_exclude_blanks(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_exclude_peaks(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_exclude_reflections(CrystFELIndexingOpts *opts);
extern char **crystfel_indexing_opts_get_metadata_to_copy(CrystFELIndexingOpts *opts,
                                                          int *n);
extern double crystfel_indexing_opts_get_fixed_profile_radius(CrystFELIndexingOpts *opts,
                                                              int *active);
extern double crystfel_indexing_opts_get_fixed_divergence(CrystFELIndexingOpts *opts);

extern int crystfel_indexing_opts_get_pinkindexer_cpeaks(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_pinkindexer_use_max_res(CrystFELIndexingOpts *opts);
extern double crystfel_indexing_opts_get_pinkindexer_max_res(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_pinkindexer_angle_density(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_pinkindexer_refinement_type(CrystFELIndexingOpts *opts);
extern double crystfel_indexing_opts_get_pinkindexer_tolerance(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_pinkindexer_use_refl_radius(CrystFELIndexingOpts *opts);
extern double crystfel_indexing_opts_get_pinkindexer_refl_radius(CrystFELIndexingOpts *opts);
extern double crystfel_indexing_opts_get_pinkindexer_max_imbalance(CrystFELIndexingOpts *opts);




extern void crystfel_indexing_opts_set_show_stream_opts(CrystFELIndexingOpts *opts,
                                                        int val);

extern void crystfel_indexing_opts_set_cell_file(CrystFELIndexingOpts *opts,
                                                 const char *cell_file);
extern void crystfel_indexing_opts_set_indexing_method_string(CrystFELIndexingOpts *opts,
                                                              const char *indm_str);
extern void crystfel_indexing_opts_set_multi_lattice(CrystFELIndexingOpts *opts,
                                                     int multi);
extern void crystfel_indexing_opts_set_refine(CrystFELIndexingOpts *opts,
                                              int refine);
extern void crystfel_indexing_opts_set_retry(CrystFELIndexingOpts *opts,
                                             int retry);
extern void crystfel_indexing_opts_set_peak_check(CrystFELIndexingOpts *opts,
                                                  int peak_check);
extern void crystfel_indexing_opts_set_cell_check(CrystFELIndexingOpts *opts,
                                                  int cell_check);
extern void crystfel_indexing_opts_set_tolerances(CrystFELIndexingOpts *opts,
                                                  float *tols); /* frac (not %) and rad */
extern void crystfel_indexing_opts_set_min_peaks(CrystFELIndexingOpts *opts,
                                                 int min_peaks);

extern void crystfel_indexing_opts_set_integration_method_string(CrystFELIndexingOpts *opts,
                                                                 const char *integr_str);
extern void crystfel_indexing_opts_set_overpredict(CrystFELIndexingOpts *opts,
                                                   int overpredict);
extern void crystfel_indexing_opts_set_push_res(CrystFELIndexingOpts *opts,
                                                float push_res);
extern void crystfel_indexing_opts_set_integration_radii(CrystFELIndexingOpts *opts,
                                                         float ir_inn,
                                                         float ir_mid,
                                                         float ir_out);
extern void crystfel_indexing_opts_set_metadata_to_copy(CrystFELIndexingOpts *opts,
                                                        char *const *headers,
                                                        int n_headers);
extern void crystfel_indexing_opts_set_exclude_blanks(CrystFELIndexingOpts *opts,
                                                      int flag);
extern void crystfel_indexing_opts_set_exclude_peaks(CrystFELIndexingOpts *opts,
                                                     int flag);
extern void crystfel_indexing_opts_set_exclude_reflections(CrystFELIndexingOpts *opts,
                                                           int flag);
extern void crystfel_indexing_opts_set_fixed_profile_radius(CrystFELIndexingOpts *opts,
                                                            int active,
                                                            double val);
extern void crystfel_indexing_opts_set_fixed_divergence(CrystFELIndexingOpts *opts,
                                                        double val);

extern void crystfel_indexing_opts_set_pinkindexer_cpeaks(CrystFELIndexingOpts *opts,
                                                          int val);
extern void crystfel_indexing_opts_set_pinkindexer_use_max_res(CrystFELIndexingOpts *opts,
                                                               int val);
extern void crystfel_indexing_opts_set_pinkindexer_max_res(CrystFELIndexingOpts *opts,
                                                           double val);
extern void crystfel_indexing_opts_set_pinkindexer_angle_density(CrystFELIndexingOpts *opts,
                                                                 int val);
extern void crystfel_indexing_opts_set_pinkindexer_refinement_type(CrystFELIndexingOpts *opts,
                                                                   int val);
extern void crystfel_indexing_opts_set_pinkindexer_tolerance(CrystFELIndexingOpts *opts,
                                                             double val);
extern void crystfel_indexing_opts_set_pinkindexer_use_refl_radius(CrystFELIndexingOpts *opts,
                                                                   int val);
extern void crystfel_indexing_opts_set_pinkindexer_refl_radius(CrystFELIndexingOpts *opts,
                                                               double val);
extern void crystfel_indexing_opts_set_pinkindexer_max_imbalance(CrystFELIndexingOpts *opts,
                                                                 double val);

#endif	/* CRYSTFELINDEXINGOPTS_H */
