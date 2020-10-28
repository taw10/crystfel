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
	mo->symmetry = gtk_combo_box_text_new(); //crystfel_symmetry_selector_new();
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
	label = gtk_label_new("Detector saturation value:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
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
	mo->detwin_sym = gtk_combo_box_text_new(); //crystfel_symmetry_selector_new();
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

	gtk_widget_show_all(GTK_WIDGET(mo));
	return GTK_WIDGET(mo);
}


//int crystfel_merge_opts_get_multi_lattice(CrystFELMergeOpts *opts)
//{
//	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->multi));
//}


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
	char tmp[64];

	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(mo->limit_res),
	                             !isinf(push_res));

	snprintf(tmp, 63, "%f", push_res);
	gtk_entry_set_text(GTK_ENTRY(mo->push_res), tmp);
}
