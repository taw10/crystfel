/*
 * crystfelmergeopts.h
 *
 * A GTK widget to set merge options
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

#ifndef CRYSTFELMERGEOPTS_H
#define CRYSTFELMERGEOPTS_H

#include <gtk/gtk.h>
#include <glib-object.h>

#define CRYSTFEL_TYPE_MERGE_OPTS (crystfel_merge_opts_get_type())

#define CRYSTFEL_MERGE_OPTS(obj) (G_TYPE_CHECK_INSTANCE_CAST((obj), \
                                 CRYSTFEL_TYPE_MERGE_OPTS, CrystFELMergeOpts))

#define CRYSTFEL_IS_MERGE_OPTS(obj) (G_TYPE_CHECK_INSTANCE_TYPE((obj), \
                                    CRYSTFEL_TYPE_MERGE_OPTS))

#define CRYSTFEL_MERGE_OPTS_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST((obj), \
                                         CRYSTFEL_TYPE_MERGE_OPTS, CrystFELMergeOpts))

#define CRYSTFEL_IS_MERGE_OPTS_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((obj), \
                                            CRYSTFEL_TYPE_MERGE_OPTS))

#define CRYSTFEL_MERGE_OPTS_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS((obj), \
                                           CRYSTFEL_TYPE_MERGE_OPTS, CrystFELMergeOpts))

struct _crystfelmergeopts
{
	GtkNotebook parent_instance;

	/*< private >*/
	GtkWidget *model_combo;
	GtkWidget *symmetry;
	GtkWidget *scale;
	GtkWidget *bscale;
	GtkWidget *postref;
	GtkWidget *niter;
	GtkWidget *polarisation;
	GtkWidget *deltacchalf;
	GtkWidget *min_measurements;
	GtkWidget *use_max_adu;
	GtkWidget *max_adu;
	GtkWidget *custom_split;
	GtkWidget *custom_split_file;
	GtkWidget *detwin;
	GtkWidget *detwin_sym;
	GtkWidget *min_res;
	GtkWidget *min_res_val;
	GtkWidget *limit_res;  /* check box for push-res < infinity */
	GtkWidget *push_res;   /* entry for push-res value */
};

struct _crystfelmergeoptsclass
{
	GtkNotebookClass parent_class;
};

typedef struct _crystfelmergeopts CrystFELMergeOpts;
typedef struct _crystfelmergeoptsclass CrystFELMergeOptsClass;

extern GType crystfel_merge_opts_get_type(void);
extern GtkWidget *crystfel_merge_opts_new(void);

extern void crystfel_merge_opts_set_model(CrystFELMergeOpts *opts,
                                          const char *model);
extern void crystfel_merge_opts_set_symmetry(CrystFELMergeOpts *opts,
                                             const char *sym);
extern void crystfel_merge_opts_set_scale(CrystFELMergeOpts *opts,
                                          int scale);
extern void crystfel_merge_opts_set_bscale(CrystFELMergeOpts *opts,
                                           int bscale);
extern void crystfel_merge_opts_set_postref(CrystFELMergeOpts *opts,
                                            int postref);
extern void crystfel_merge_opts_set_niter(CrystFELMergeOpts *opts,
                                          int niter);
extern void crystfel_merge_opts_set_polarisation(CrystFELMergeOpts *opts,
                                                 const char *polar);
extern void crystfel_merge_opts_set_deltacchalf(CrystFELMergeOpts *opts,
                                                int deltacchalf);
extern void crystfel_merge_opts_set_min_measurements(CrystFELMergeOpts *opts,
                                                     int min_measurements);
extern void crystfel_merge_opts_set_max_adu(CrystFELMergeOpts *opts,
                                            float max_adu);
extern void crystfel_merge_opts_set_custom_split(CrystFELMergeOpts *opts,
                                                 const char *custom_split_file);
extern void crystfel_merge_opts_set_twin_sym(CrystFELMergeOpts *opts,
                                             const char *twin_sym);
extern void crystfel_merge_opts_set_min_res(CrystFELMergeOpts *opts,
                                            float min_res);
extern void crystfel_merge_opts_set_push_res(CrystFELMergeOpts *opts,
                                             float push_res);

extern const char *crystfel_merge_opts_get_model(CrystFELMergeOpts *opts);
extern const char *crystfel_merge_opts_get_symmetry(CrystFELMergeOpts *opts);
extern int crystfel_merge_opts_get_scale(CrystFELMergeOpts *opts);
extern int crystfel_merge_opts_get_bscale(CrystFELMergeOpts *opts);
extern int crystfel_merge_opts_get_postref(CrystFELMergeOpts *opts);
extern int crystfel_merge_opts_get_niter(CrystFELMergeOpts *opts);
extern const char *crystfel_merge_opts_get_polarisation(CrystFELMergeOpts *opts);
extern int crystfel_merge_opts_get_deltacchalf(CrystFELMergeOpts *opts);
extern int crystfel_merge_opts_get_min_measurements(CrystFELMergeOpts *opts);
extern float crystfel_merge_opts_get_max_adu(CrystFELMergeOpts *opts);
extern const char *crystfel_merge_opts_get_custom_split(CrystFELMergeOpts *opts);
extern const char *crystfel_merge_opts_get_twin_sym(CrystFELMergeOpts *opts);
extern float crystfel_merge_opts_get_min_res(CrystFELMergeOpts *opts);
extern float crystfel_merge_opts_get_push_res(CrystFELMergeOpts *opts);

#endif	/* CRYSTFELMERGEOPTS_H */
