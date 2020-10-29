/*
 * crystfelmergeopts.h
 *
 * A GTK widget to set merge options
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <gtk/gtk.h>
#include <glib-object.h>
#include <errno.h>
#include <math.h>

#include <integration.h>
#include <index.h>

#include "crystfelmergeopts.h"
#include "crystfelsymmetryselector.h"


G_DEFINE_TYPE(CrystFELMergeOpts,
              crystfel_merge_opts,
              GTK_TYPE_NOTEBOOK)


static void crystfel_merge_opts_class_init(CrystFELMergeOptsClass *klass)
{
}


static void crystfel_merge_opts_init(CrystFELMergeOpts *mo)
{
}


static GtkWidget *merge_parameters(CrystFELMergeOpts *mo)
{
	GtkWidget *box;
	GtkWidget *hbox;
	GtkWidget *label;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Model:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	mo->model_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->model_combo),
	                   FALSE, FALSE, 0);
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(mo->model_combo),
	                          "process_hkl",
	                          "Simple merging (process_hkl)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(mo->model_combo),
	                          "unity",
	                          "No partialities (unity)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(mo->model_combo),
	                          "xsphere",
	                          "Bandwidth integral (xsphere)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(mo->model_combo),
	                          "offset",
	                          "Monochromatic Ewald sphere offset (offset)");

	/* Symmetry */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Symmetry:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	mo->symmetry = crystfel_symmetry_selector_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->symmetry),
	                   FALSE, FALSE, 0);

	/* Scale, Bscale, post-ref on/off */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	mo->scale = gtk_check_button_new_with_label("Scale intensities");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->scale),
	                   FALSE, FALSE, 0);
	mo->bscale = gtk_check_button_new_with_label("Debye-Waller scaling");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->bscale),
	                   FALSE, FALSE, 0);
	mo->postref = gtk_check_button_new_with_label("Post-refinement");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->postref),
	                   FALSE, FALSE, 0);

	/* Number of iterations */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Number of scaling/post-refinement cycles:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	mo->niter = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(mo->niter), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->niter),
	                   FALSE, FALSE, 0);

	/* Polarisation horiz/vert/unpolarized beam/no correction */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Polarisation:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	mo->polarisation = gtk_combo_box_text_new();
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(mo->polarisation), "horiz",
	                          "Horizontal e-field (most synchrotrons and FELs)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(mo->polarisation), "vert",
	                          "Vertical e-field (LCLS-II hard X-ray undulator)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(mo->polarisation), "horiz50",
	                          "Unpolarized incident beam");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(mo->polarisation), "none",
	                          "No correction");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->polarisation),
	                   FALSE, FALSE, 0);

	/* deltaCChalf */
	mo->deltacchalf = gtk_check_button_new_with_label("Reject bad patterns according to ΔCC½");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(mo->deltacchalf),
	                   FALSE, FALSE, 0);

	/* Detector saturation value */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	mo->use_max_adu = gtk_check_button_new_with_label("Detector saturation cutoff:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->use_max_adu),
	                   FALSE, FALSE, 0);
	mo->max_adu = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(mo->max_adu), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->max_adu),
	                   FALSE, FALSE, 0);

	/* Minimum measurements */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Minimum number of measurements per merged reflection:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	mo->min_measurements = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(mo->min_measurements), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->min_measurements),
	                   FALSE, FALSE, 0);

	/* Custom split */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	mo->custom_split = gtk_check_button_new_with_label("Split datasets after merging");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->custom_split),
	                   FALSE, FALSE, 0);
	mo->custom_split_file = gtk_file_chooser_button_new("Dataset assignments",
	                                                    GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->custom_split_file),
	                   FALSE, FALSE, 0);

	/* Minimum pattern resolution */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	mo->min_res = gtk_check_button_new_with_label("Require minimum estimated pattern resolution");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->min_res),
	                   FALSE, FALSE, 0);
	mo->min_res_val = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(mo->min_res_val), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->min_res_val),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Å");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	/* push-res */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	mo->limit_res = gtk_check_button_new_with_label("Exclude measurements more than");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->limit_res),
	                   FALSE, FALSE, 0);
	mo->push_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(mo->push_res), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->push_res),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("nm^-1 above resolution limit");
	gtk_label_set_markup(GTK_LABEL(label), "nm<sup>-1</sup> above resolution limit");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	/* Detwin */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	mo->detwin = gtk_check_button_new_with_label("Refine indexing assignments");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->detwin),
	                   FALSE, FALSE, 0);
	mo->detwin_sym = crystfel_symmetry_selector_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(mo->detwin_sym),
	                   FALSE, FALSE, 0);

	return box;
}


GtkWidget *crystfel_merge_opts_new()
{
	CrystFELMergeOpts *mo;

	mo = g_object_new(CRYSTFEL_TYPE_MERGE_OPTS, NULL);

	gtk_notebook_append_page(GTK_NOTEBOOK(mo),
	                         merge_parameters(mo),
	                         gtk_label_new("Merging"));

	return GTK_WIDGET(mo);
}


static int set_if(GtkWidget *combo, const char *a, const char *new_id)
{
	if ( strcmp(a, new_id) == 0 ) {
		gtk_combo_box_set_active_id(GTK_COMBO_BOX(combo),
		                            new_id);
		return 1;
	}

	return 0;
}


static void set_int(GtkWidget *entry, int val)
{
	char tmp[64];
	snprintf(tmp, 63, "%i", val);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);
}


static void set_float(GtkWidget *entry, float val)
{
	char tmp[64];
	snprintf(tmp, 63, "%f", val);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);
}


void crystfel_merge_opts_set_model(CrystFELMergeOpts *opts,
                                   const char *model)
{
	int done = 0;
	done += set_if(opts->model_combo, model, "process_hkl");
	done += set_if(opts->model_combo, model, "unity");
	done += set_if(opts->model_combo, model, "xsphere");
	done += set_if(opts->model_combo, model, "offset");
	done += set_if(opts->model_combo, model, "ggpm");

	if ( done == 0 ) {
		ERROR("Unrecognised model '%s'\n", model);
	}
}


void crystfel_merge_opts_set_polarisation(CrystFELMergeOpts *opts,
                                          const char *polar)
{
	int done = 0;
	done += set_if(opts->polarisation, polar, "horiz");
	done += set_if(opts->polarisation, polar, "vert");
	done += set_if(opts->polarisation, polar, "none");
	done += set_if(opts->polarisation, polar, "horiz50");

	if ( done == 0 ) {
		ERROR("Unrecognised polarisation '%s'\n", polar);
	}
}


void crystfel_merge_opts_set_symmetry(CrystFELMergeOpts *opts,
                                      const char *sym)
{
	crystfel_symmetry_selector_set_group_symbol(CRYSTFEL_SYMMETRY_SELECTOR(opts->symmetry),
	                                            sym);
}


void crystfel_merge_opts_set_scale(CrystFELMergeOpts *opts,
                                   int scale)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->scale),
	                             scale);
}


void crystfel_merge_opts_set_bscale(CrystFELMergeOpts *opts,
                                    int bscale)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->bscale),
	                             bscale);
}


void crystfel_merge_opts_set_postref(CrystFELMergeOpts *opts,
                                     int postref)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->postref),
	                             postref);
}


void crystfel_merge_opts_set_niter(CrystFELMergeOpts *opts,
                                   int niter)
{
	set_int(opts->niter, niter);
}


void crystfel_merge_opts_set_deltacchalf(CrystFELMergeOpts *opts,
                                         int deltacchalf)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->deltacchalf),
	                             deltacchalf);
}


void crystfel_merge_opts_set_min_measurements(CrystFELMergeOpts *opts,
                                              int min_measurements)
{
	set_int(opts->min_measurements, min_measurements);
}


void crystfel_merge_opts_set_max_adu(CrystFELMergeOpts *opts,
                                     float max_adu)
{
	set_float(opts->max_adu, max_adu);
}


void crystfel_merge_opts_set_custom_split(CrystFELMergeOpts *opts,
                                          const char *custom_split_file)
{
	if ( custom_split_file != NULL ) {
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(opts->custom_split_file),
		                              custom_split_file);
	}
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->custom_split),
	                             (custom_split_file != NULL));
}


void crystfel_merge_opts_set_twin_sym(CrystFELMergeOpts *opts,
                                      const char *twin_sym)
{
	crystfel_symmetry_selector_set_group_symbol(CRYSTFEL_SYMMETRY_SELECTOR(opts->detwin_sym),
	                                            twin_sym);
}


void crystfel_merge_opts_set_min_res(CrystFELMergeOpts *opts,
                                     float min_res)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->min_res),
	                             !isinf(min_res));
	set_float(opts->min_res_val, min_res);
}


float crystfel_merge_opts_get_push_res(CrystFELMergeOpts *opts)
{
	if ( !gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->limit_res)) ) {
		return INFINITY;
	} else {
		const gchar *text;
		float push_res;
		char *rval;
		text = gtk_entry_get_text(GTK_ENTRY(opts->push_res));
		errno = 0;
		push_res = strtof(text, &rval);
		if ( *rval != '\0' ) {
			printf("Invalid value for push-res (%s)\n",
			      rval);
			return INFINITY;
		}
		return push_res;
	}
}


void crystfel_merge_opts_set_push_res(CrystFELMergeOpts *mo,
                                      float push_res)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(mo->limit_res),
	                             !isinf(push_res));
	set_float(mo->push_res, push_res);
}
