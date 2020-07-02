/*
 * gui_peaksearch.c
 *
 * Peak search parts of GUI
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdarg.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms-compat.h>
#include <assert.h>

#include <datatemplate.h>
#include <peaks.h>

#include "crystfel_gui.h"
#include "crystfelimageview.h"


void update_peaks(struct crystfelproject *proj)
{

	if ( proj->n_frames == 0 ) return;

	if ( proj->show_peaks ) {

		struct image *image;

		image = crystfel_image_view_get_image_struct(CRYSTFEL_IMAGE_VIEW(proj->imageview));
		if ( image == NULL ) return;

		image_feature_list_free(image->features);
		image->features = NULL;

		switch ( proj->peak_search_params.method ) {

		case PEAK_ZAEF:
			search_peaks(image,
			             proj->peak_search_params.threshold,
			             proj->peak_search_params.min_sq_gradient,
			             proj->peak_search_params.min_snr,
			             proj->peak_search_params.pk_inn,
			             proj->peak_search_params.pk_mid,
			             proj->peak_search_params.pk_out,
			             1);
			break;

		case PEAK_PEAKFINDER8:
			search_peaks_peakfinder8(image, 2048,
			                         proj->peak_search_params.threshold,
			                         proj->peak_search_params.min_snr,
			                         proj->peak_search_params.min_pix_count,
			                         proj->peak_search_params.max_pix_count,
			                         proj->peak_search_params.local_bg_radius,
			                         proj->peak_search_params.min_res,
			                         proj->peak_search_params.max_res,
			                         1);
			break;

		case PEAK_HDF5:
		case PEAK_CXI:
			image->features = image_read_peaks(proj->dtempl,
			                                   image->filename,
			                                   image->ev,
			                                   proj->peak_search_params.half_pixel_shift);
			if ( proj->peak_search_params.revalidate ) {
				validate_peaks(image,
				               proj->peak_search_params.min_snr,
				               proj->peak_search_params.pk_inn,
				               proj->peak_search_params.pk_mid,
				               proj->peak_search_params.pk_out,
				               1, 0);
			}
			break;

		default:
			ERROR("This peak detection method not implemented!\n");
			break;

		}

		crystfel_image_view_set_peaks(CRYSTFEL_IMAGE_VIEW(proj->imageview),
		                              image->features, 0);

	} else {
		crystfel_image_view_set_peaks(CRYSTFEL_IMAGE_VIEW(proj->imageview),
		                              NULL, 0);
	}
}


struct param_callback_vals
{
	float *pfval;
	int *pival;
	struct crystfelproject *proj;
};


static void check_param_callback(GtkWidget *checkbox,
                                 struct param_callback_vals *vals)
{
	int val = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbox));
	*(vals->pival) = val;
	update_peaks(vals->proj);
}


static void int_param_callback(GtkWidget *entry,
                               struct param_callback_vals *cbvals)
{
	const char *text;
	int val;

	text = gtk_entry_get_text(GTK_ENTRY(entry));
	if (sscanf(text, "%i", &val) != 1) {
		ERROR("Invalid value\n");
		return;
	}

	*(cbvals->pival) = val;
	cbvals->proj->unsaved = 1;
	update_peaks(cbvals->proj);
}


static void float_param_callback(GtkWidget *entry,
                                 struct param_callback_vals *cbvals)
{
	const char *text;
	float val;

	text = gtk_entry_get_text(GTK_ENTRY(entry));
	if (sscanf(text, "%f", &val) != 1) {
		ERROR("Invalid value\n");
		return;
	}

	*(cbvals->pfval) = val;
	cbvals->proj->unsaved = 1;
	update_peaks(cbvals->proj);
}


static void free_callback_params(gpointer cbvals,
                                 GClosure *closure)
{
	free(cbvals);
}


static void add_int_param(GtkWidget *params_box, const char *labeltext,
                          int *pval, struct crystfelproject *proj)
{
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *entry;
	struct param_callback_vals *cbvals;
	char tmp[64];

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(params_box),
	                   GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new(labeltext);
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 2.0);
	entry = gtk_entry_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   TRUE, TRUE, 2.0);
	snprintf(tmp, 63, "%i", *pval);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);

	cbvals = malloc(sizeof(struct param_callback_vals));
	if ( cbvals != NULL ) {
		cbvals->proj = proj;
		cbvals->pival = pval;
		g_signal_connect_data(G_OBJECT(entry), "activate",
		                      G_CALLBACK(int_param_callback),
		                      cbvals, free_callback_params, 0);
	} else {
		ERROR("Failed to connect parameter callback\n");
	}
}


static void add_float_param(GtkWidget *params_box, const char *labeltext,
                            float *pval, struct crystfelproject *proj)
{
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *entry;
	struct param_callback_vals *cbvals;
	char tmp[64];

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(params_box),
	                   GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new(labeltext);
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 2.0);
	entry = gtk_entry_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   TRUE, TRUE, 2.0);
	snprintf(tmp, 63, "%.2f", *pval);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);

	cbvals = malloc(sizeof(struct param_callback_vals));
	if ( cbvals != NULL ) {
		cbvals->proj = proj;
		cbvals->pfval = pval;
		g_signal_connect_data(G_OBJECT(entry), "activate",
		                      G_CALLBACK(float_param_callback),
		                      cbvals, free_callback_params, 0);
	} else {
		ERROR("Failed to connect parameter callback\n");
	}
}


static void add_check_param(GtkWidget *params_box, const char *labeltext,
                            int *pval, struct crystfelproject *proj)
{
	GtkWidget *checkbox;
	struct param_callback_vals *cbvals;

	checkbox = gtk_check_button_new_with_label(labeltext);
	gtk_box_pack_start(GTK_BOX(params_box),
	                   GTK_WIDGET(checkbox), FALSE, FALSE, 8.0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbox),
	                             *pval);

	cbvals = malloc(sizeof(struct param_callback_vals));
	if ( cbvals != NULL ) {
		cbvals->proj = proj;
		cbvals->pival = pval;
		g_signal_connect_data(G_OBJECT(checkbox), "toggled",
		                      G_CALLBACK(check_param_callback),
		                      cbvals, free_callback_params, 0);
	} else {
		ERROR("Failed to connect parameter callback\n");
	}
}


static void add_radii(GtkWidget *params_box,
                      struct crystfelproject *proj)
{
	add_float_param(params_box, "Peak radius (inner):",
	                &proj->peak_search_params.pk_inn, proj);
	add_float_param(params_box, "Peak radius (middle):",
	                &proj->peak_search_params.pk_mid, proj);
	add_float_param(params_box, "Peak radius (outer):",
	                &proj->peak_search_params.pk_out, proj);
}


static void peaksearch_algo_changed(GtkWidget *combo,
                                    struct crystfelproject *proj)
{
	const char *algo_id;

	if ( proj->peak_params != NULL ) {
		gtk_container_remove(GTK_CONTAINER(proj->peak_vbox),
		                     proj->peak_params);
		proj->peak_params = NULL;
	}

	proj->peak_params = gtk_vbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(proj->peak_vbox),
	                   GTK_WIDGET(proj->peak_params),
	                   FALSE, FALSE, 8.0);

	algo_id = gtk_combo_box_get_active_id(GTK_COMBO_BOX(combo));

	if ( strcmp(algo_id, "zaef") == 0 ) {

		proj->peak_search_params.method = PEAK_ZAEF;

		add_float_param(proj->peak_params, "Threshold:",
		                &proj->peak_search_params.threshold, proj);
		add_float_param(proj->peak_params, "Minimum squared gradient:",
		                &proj->peak_search_params.min_sq_gradient, proj);
		add_float_param(proj->peak_params, "Minimum signal/noise ratio:",
		                &proj->peak_search_params.min_snr, proj);
		add_radii(proj->peak_params, proj);

	} else if ( strcmp(algo_id, "peakfinder8") == 0 ) {

		proj->peak_search_params.method = PEAK_PEAKFINDER8;

		add_float_param(proj->peak_params, "Threshold:",
		                &proj->peak_search_params.threshold, proj);
		add_float_param(proj->peak_params, "Minimum signal/noise ratio:",
		                &proj->peak_search_params.min_snr, proj);
		add_int_param(proj->peak_params, "Minimum number of pixels:",
		              &proj->peak_search_params.min_pix_count, proj);
		add_int_param(proj->peak_params, "Maximum number of pixels:",
		              &proj->peak_search_params.max_pix_count, proj);
		add_int_param(proj->peak_params, "Local background radius:",
		              &proj->peak_search_params.local_bg_radius, proj);
		add_int_param(proj->peak_params, "Minimum resolution (pixels):",
		              &proj->peak_search_params.min_res, proj);
		add_int_param(proj->peak_params, "Maximum resolution (pixels):",
		              &proj->peak_search_params.max_res, proj);

	} else if ( strcmp(algo_id, "hdf5") == 0 ) {

		proj->peak_search_params.method = PEAK_HDF5;

		add_check_param(proj->peak_params, "Half pixel shift",
		                &proj->peak_search_params.half_pixel_shift,
		                proj);
		add_check_param(proj->peak_params, "Check peaks first",
		                &proj->peak_search_params.revalidate,
		                proj);
		add_radii(proj->peak_params, proj);

	} else if ( strcmp(algo_id, "cxi") == 0 ) {
		ERROR("algo_id should be hdf5, not cxi\n");

	} else {
		ERROR("Unrecognised peak search '%s'\n", algo_id);
	}

	/* FIXME: Radii */

	gtk_widget_show_all(proj->peak_vbox);
	proj->unsaved = 1;
	update_peaks(proj);
}


static void peaksearch_response_sig(GtkWidget *dialog, gint resp,
                                    struct crystfelproject *proj)
{

	if ( (resp == GTK_RESPONSE_DELETE_EVENT)
	  || (resp == GTK_RESPONSE_CANCEL) )
	{
		proj->peak_search_params = proj->original_params;
		update_peaks(proj);
	}

	gtk_widget_destroy(dialog);
	proj->peak_vbox = NULL;
	proj->peak_params = NULL;
}


gint peaksearch_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *combo;
	GtkWidget *w;

	if ( proj->peak_params != NULL ) return FALSE;

	proj->show_peaks = 1;

	w =  gtk_ui_manager_get_widget(proj->ui, "/ui/mainwindow/view/peaks");
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(w), 1);

	/* Take a copy of the original parameters, for reverting */
	proj->original_params = proj->peak_search_params;

	dialog = gtk_dialog_new_with_buttons("Peak search",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Discard changes", GTK_RESPONSE_CANCEL,
	                                     "Confirm", GTK_RESPONSE_ACCEPT,
	                                     NULL);

	g_signal_connect(G_OBJECT(dialog), "response",
	                 G_CALLBACK(peaksearch_response_sig), proj);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);
	proj->peak_vbox = vbox;

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Peak search algorithm");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(combo), TRUE, TRUE, 2.0);
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "zaef",
	                "Zaefferer gradient search (zaef)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "hdf5",
	                "Use the peak lists in the data files (hdf5/cxi)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "peakfinder8",
	                "Radial background estimation (peakfinder8)");
	#ifdef HAVE_FDIP
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "peakfinder9",
	                "Local background estimation (peakfinder9)");
	#endif
	g_signal_connect(G_OBJECT(combo), "changed",
	                 G_CALLBACK(peaksearch_algo_changed), proj);
	gtk_combo_box_set_active_id(GTK_COMBO_BOX(combo),
	                            str_peaksearch(proj->peak_search_params.method));
	proj->type_combo = combo;

	gtk_widget_show_all(dialog);

	return FALSE;
}
