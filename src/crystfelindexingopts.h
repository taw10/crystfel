/*
 * crystfelindexingopts.h
 *
 * A GTK widget to set indexing options
 *
 * Copyright © 2020 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
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
	GtkWidget *tols[6];
	GtkWidget *enable_hitfind;
	GtkWidget *ignore_fewer_peaks;

	GtkWidget *integration_combo;
	GtkWidget *centering;
	GtkWidget *overpredict;
	GtkWidget *limit_res;
	GtkWidget *push_res;
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
                                                  float *tols);
extern int crystfel_indexing_opts_get_min_peaks(CrystFELIndexingOpts *opts);

extern char *crystfel_indexing_opts_get_integration_method_string(CrystFELIndexingOpts *opts);
extern int crystfel_indexing_opts_get_overpredict(CrystFELIndexingOpts *opts);
extern float crystfel_indexing_opts_get_push_res(CrystFELIndexingOpts *opts);


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
                                                  float *tols);
extern void crystfel_indexing_opts_set_min_peaks(CrystFELIndexingOpts *opts,
                                                 int min_peaks);

extern void crystfel_indexing_opts_set_integration_method_string(CrystFELIndexingOpts *opts,
                                                                 const char *integr_str);
extern void crystfel_indexing_opts_set_overpredict(CrystFELIndexingOpts *opts,
                                                   int overpredict);
extern void crystfel_indexing_opts_set_push_res(CrystFELIndexingOpts *opts,
                                                float push_res);

#endif	/* CRYSTFELINDEXINGOPTS_H */
