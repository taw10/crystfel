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
	/* Scale on/off */
	/* B-factor scaling */
	/* Post-refinement on/off */
	/* Polarisation horiz/vert/unpolarized beam/no correction */
	/* Number of iterations */
	/* deltaCChalf */
	/* Detector saturation value */
	/* Minimum pattern resolution */
	/* push-res */
	/* Minimum measurements */
	/* Custom split file */
	/* Detwin */

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
