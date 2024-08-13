/*
 * crystfel_gui.c
 *
 * CrystFEL's main graphical user interface
 *
 * Copyright © 2020-2023 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020-2023 Thomas White <taw@physics.org>
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
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <ctype.h>

#include <datatemplate.h>
#include <peaks.h>
#include <cell-utils.h>

#include "crystfelimageview.h"
#include "crystfelcolourscale.h"
#include "crystfel_gui.h"
#include "gui_import.h"
#include "gui_peaksearch.h"
#include "gui_index.h"
#include "gui_merge.h"
#include "gui_fom.h"
#include "gui_export.h"
#include "gui_align.h"
#include "gui_ambi.h"
#include "gui_project.h"
#include "version.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [data.stream]\n\n", s);
	printf(
"CrystFEL graphical user interface.\n"
"\n"
" -h, --help              Display this help message.\n"
"     --version           Print CrystFEL version number and exit.\n"

);
}

static int confirm_exit(struct crystfelproject *proj)
{
	int r;
	GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(proj->window),
	                                           0,
	                                           GTK_MESSAGE_QUESTION,
	                                           GTK_BUTTONS_NONE,
	                                           "Do you want to save the session?");
	gtk_dialog_add_buttons(GTK_DIALOG(dialog),
	                       "Save", GTK_RESPONSE_YES,
	                       "Don't save", GTK_RESPONSE_NO,
	                       "Cancel", GTK_RESPONSE_CANCEL,
	                       NULL);
	r = gtk_dialog_run(GTK_DIALOG(dialog));
	gtk_widget_destroy(dialog);
	if ( r == GTK_RESPONSE_YES ) {
		save_project(proj);
		return 1;
	}
	if ( r == GTK_RESPONSE_NO ) return 1;
	return 0;
}


void error_box(struct crystfelproject *proj, const char *message)
{
	GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(proj->window),
	                                           0,
	                                           GTK_MESSAGE_WARNING,
	                                           GTK_BUTTONS_NONE,
	                                           "%s",
	                                           message);
	gtk_dialog_add_buttons(GTK_DIALOG(dialog),
	                       "OK", GTK_RESPONSE_OK,
	                       NULL);
	gtk_dialog_run(GTK_DIALOG(dialog));
	gtk_widget_destroy(dialog);
}


/* Main window destroyed */
static gboolean delete_event_sig(GtkWidget *da, GdkEvent *event,
                                 struct crystfelproject *proj)
{
	if ( proj->unsaved ) {
		if ( !confirm_exit(proj) ) return TRUE;
	}
	gtk_main_quit();
	return FALSE;
}


static void add_ui_sig(GtkUIManager *ui, GtkWidget *widget,
                       GtkContainer *container)
{
	gtk_box_pack_start(GTK_BOX(container), widget, FALSE, FALSE, 0);
	if ( GTK_IS_TOOLBAR(widget) ) {
		gtk_toolbar_set_show_arrow(GTK_TOOLBAR(widget), TRUE);
	}
}


static void swap_data_arrays(struct image *a, struct image *b)
{
	float **swap;
	int **swap_bad;

	if ( (a==NULL) || (b==NULL) ) return;

	swap = a->dp;
	a->dp = b->dp;
	b->dp = swap;

	swap = a->sat;
	a->sat = b->sat;
	b->sat = swap;

	swap_bad = a->bad;
	a->bad = b->bad;
	b->bad = swap_bad;
}


void select_result(struct crystfelproject *proj,
                   const char *result_name)
{
	gtk_combo_box_set_active_id(GTK_COMBO_BOX(proj->results_combo),
	                            result_name);
}


const char *selected_result(struct crystfelproject *proj)
{
	return gtk_combo_box_get_active_id(GTK_COMBO_BOX(proj->results_combo));
}


/* Return non-zero if it makes sense to re-scan the streams for results.
 * e.g. if there are any jobs running, or one finished since the last scan.
 * Possible future improvement: exclude jobs which don't produce streams
 *  (i.e. merging) */
static int should_rescan_streams(struct crystfelproject *proj)
{
	GSList *item = proj->tasks;

	if ( !proj->rescan_on_change ) return 0;

	while ( item != NULL ) {
		struct gui_task *task = item->data;
		if ( task->running ) return 1;
		item = item->next;
	}

	return 0;
}


/* Bring the image view up to date after changing the selected image */
void update_imageview(struct crystfelproject *proj)
{
	char tmp[1024];
	char *ev_str;
	char *ev_sep;
	struct image *image;
	const gchar *results_name;

	if ( proj->image_info == NULL ) return;

	if ( (proj->dtempl == NULL)
	  || (proj->filenames == NULL)
	  || (proj->events == NULL)
	  || (proj->n_frames == 0) )
	{
		gtk_label_set_text(GTK_LABEL(proj->image_info),
		                   "Ready to load images");
		return;
	}

	if ( file_exists(proj->filenames[proj->cur_frame]) ) {
		image = image_read(proj->dtempl,
		                   proj->filenames[proj->cur_frame],
		                   proj->events[proj->cur_frame],
		                   0, 0, NULL);
	} else {
		STATUS("Image data file not present.\n");
		image = NULL;
	}

	if ( proj->events[proj->cur_frame] != NULL ) {
		ev_str = proj->events[proj->cur_frame];
		ev_sep = " ";
	} else {
		ev_str = "";
		ev_sep = "";
	}
	snprintf(tmp, 1023, "%s%s%s (frame %i of %i%s)",
	         proj->filenames[proj->cur_frame],
	         ev_sep,
	         ev_str,
	         proj->cur_frame+1,
	         proj->n_frames,
		 (image==NULL)?", load error":"");
	gtk_label_set_text(GTK_LABEL(proj->image_info), tmp);

	/* Give CrystFELImageView a chance to free resources */
	crystfel_image_view_set_image(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                              NULL);
	image_free(proj->cur_image);
	proj->cur_image = image;

	/* Look up results, if applicable */
	results_name = gtk_combo_box_get_active_id(GTK_COMBO_BOX(proj->results_combo));
	if ( strcmp(results_name, "crystfel-gui-internal") == 0 ) {
		update_peaks(proj);
	} else {
		struct image *res_im;

		res_im = find_indexed_image(proj,
		                            results_name,
		                            proj->filenames[proj->cur_frame],
		                            proj->events[proj->cur_frame],
		                            should_rescan_streams(proj));
		if ( res_im != NULL ) {
			swap_data_arrays(image, res_im);
			image_free(proj->cur_image);
			proj->cur_image = res_im;
		} else {
			ERROR("Failed to load chunk from stream.  "
			      "Just displaying the image.\n");
		}
	}

	crystfel_image_view_set_show_reflections(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                         proj->show_refls);
	crystfel_image_view_set_label_reflections(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                          proj->label_refls);
	crystfel_image_view_set_refl_box_size(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                      proj->indexing_params.ir_inn);
	crystfel_image_view_set_show_peaks(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                   proj->show_peaks);
	crystfel_image_view_set_show_centre(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                   proj->show_centre);
	crystfel_image_view_set_resolution_rings(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                         proj->resolution_rings);
	crystfel_image_view_set_image(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                              proj->cur_image);

	crystfel_colour_scale_scan_image(CRYSTFEL_COLOUR_SCALE(proj->colscale),
	                                 proj->cur_image);
	if ( !proj->range_set ) {
		crystfel_colour_scale_auto_range(CRYSTFEL_COLOUR_SCALE(proj->colscale));
	}

	gtk_widget_set_sensitive(proj->next_button,
	                         !(proj->cur_frame == proj->n_frames-1));
	gtk_widget_set_sensitive(proj->last_button,
	                         !(proj->cur_frame == proj->n_frames-1));
	gtk_widget_set_sensitive(proj->prev_button,
	                         !(proj->cur_frame == 0));
	gtk_widget_set_sensitive(proj->first_button,
	                         !(proj->cur_frame == 0));
}


static void push_random_frame(struct crystfelproject *proj, int fr)
{
	memmove(&proj->random_history[1],
	        &proj->random_history[0],
	        (N_RANDOM_HISTORY-1)*sizeof(int));
	proj->random_history[0] = fr;
	proj->n_random_history++;
	if ( proj->n_random_history > N_RANDOM_HISTORY ) {
		proj->n_random_history = N_RANDOM_HISTORY;
	}
}


static int pop_random_frame(struct crystfelproject *proj)
{
	int fr;
	assert(proj->n_random_history > 0);
	fr = proj->random_history[0];
	memmove(&proj->random_history[0],
	        &proj->random_history[1],
	        (N_RANDOM_HISTORY-1)*sizeof(int));
	proj->n_random_history--;
	return fr;
}


static int goto_frame(struct crystfelproject *proj,
                      const char *filename,
                      const char *ev)
{
	int i;
	for ( i=0; i<proj->n_frames; i++ ) {
		if ( strcmp(proj->events[i], ev) == 0 ) {
			if ( strcmp(proj->filenames[i], filename) == 0 ) {
				proj->cur_frame = i;
				return 0;
			}
		}
	}
	return 1;
}


struct goto_frame_stuff
{
	GtkWidget *framenumchk;
	GtkWidget *frameidchk;
	GtkWidget *framenum;
	GtkWidget *filename;
	GtkWidget *entry;  /* Frame ID */
	struct crystfelproject *proj;
};


static void free_goto_frame_stuff(gpointer stuff,
                                  GClosure *closure)
{
	free(stuff);
}


static void goto_frame_activate_sig(GtkEntry *entry, GtkDialog *dialog)
{
	gtk_dialog_response(dialog, GTK_RESPONSE_OK);
}


static void goto_frame_response_sig(GtkDialog *dialog, gint response_id,
                                    struct goto_frame_stuff *ctx)
{
	if ( response_id != GTK_RESPONSE_OK ) {
		gtk_widget_destroy(GTK_WIDGET(dialog));
		return;
	}

	if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(ctx->framenumchk)) ) {
		int id = atoi(gtk_entry_get_text(GTK_ENTRY(ctx->framenum)));
		if ( (id<=0) || (id>ctx->proj->n_frames) ) {
			error_box(ctx->proj, "Invalid frame number");
		} else {
			ctx->proj->cur_frame = id-1;
			update_imageview(ctx->proj);
			gtk_widget_destroy(GTK_WIDGET(dialog));
		}

	} else {
		int idx = gtk_combo_box_get_active(GTK_COMBO_BOX(ctx->filename));
		if ( (idx==-1) || goto_frame(ctx->proj, ctx->proj->unique_files[idx],
		                             gtk_entry_get_text(GTK_ENTRY(ctx->entry))) )
		{
			error_box(ctx->proj, "Frame ID not found");
		} else {
			gtk_widget_destroy(GTK_WIDGET(dialog));
			update_imageview(ctx->proj);
		}
	}
}


static void goto_toggle_sig(GtkWidget *radio, struct goto_frame_stuff *ctx)
{
	gtk_widget_set_sensitive(ctx->framenum, FALSE);
	gtk_widget_set_sensitive(ctx->filename, FALSE);
	gtk_widget_set_sensitive(ctx->entry, FALSE);

	if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(ctx->framenumchk)) ) {
		gtk_widget_set_sensitive(ctx->framenum, TRUE);
	} else {
		gtk_widget_set_sensitive(ctx->filename, TRUE);
		gtk_widget_set_sensitive(ctx->entry, TRUE);
	}
}


static gint goto_frame_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *window;
	GtkWidget *hbox;
	GtkWidget *vbox;
	GtkWidget *label;
	GtkWidget *content_area;
	GSList *grp;
	char tmp[64];
	int i;
	struct goto_frame_stuff *stuff;

	if ( proj->events == NULL ) return 0;

	stuff = malloc(sizeof(struct goto_frame_stuff));
	if ( stuff == NULL ) return 0;

	stuff->proj = proj;

	window = gtk_dialog_new_with_buttons("Jump to frame",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     GTK_STOCK_CANCEL, GTK_RESPONSE_CLOSE,
	                                     GTK_STOCK_OK, GTK_RESPONSE_OK, NULL);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(window));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	/* Jump to frame number (sequence) */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	stuff->framenumchk = gtk_radio_button_new_with_label(NULL, "Frame number:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(stuff->framenumchk),
	                   FALSE, FALSE, 4.0);
	stuff->framenum = gtk_entry_new();
	snprintf(tmp, 64, "%i", proj->cur_frame);
	gtk_entry_set_text(GTK_ENTRY(stuff->framenum), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(stuff->framenum),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(stuff->framenum), "activate",
			 G_CALLBACK(goto_frame_activate_sig), window);

	/* Jump to filename+event ID */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 2.0);
	gtk_widget_set_margin_top(hbox, 6.0);
	grp = gtk_radio_button_get_group(GTK_RADIO_BUTTON(stuff->framenumchk));
	stuff->frameidchk = gtk_radio_button_new_with_label(grp, "Filename:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(stuff->frameidchk),
	                   FALSE, FALSE, 4.0);
	stuff->filename = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(stuff->filename),
	                   FALSE, FALSE, 2.0);
	for ( i=0; i<proj->n_unique_files; i++ ) {
		char *bn = safe_basename(proj->unique_files[i]);
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(stuff->filename), NULL, bn);
		free(bn);
		if ( strcmp(proj->unique_files[i], proj->filenames[proj->cur_frame]) == 0 ) {
			gtk_combo_box_set_active(GTK_COMBO_BOX(stuff->filename), i);
		}
	}
	if ( proj->n_unique_files >= 20 ) {
		label = gtk_label_new("Only the first 20 filenames are listed here");
		gtk_label_set_markup(GTK_LABEL(label),
		                     "<b><i>Only the first 20 filenames are "
		                     "listed here</i></b>");
		gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(label), FALSE, FALSE, 8.0);
	}
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Frame ID:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 4.0);
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_widget_set_margin_start(label, 32);
	stuff->entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(stuff->entry),
	                   proj->events[proj->cur_frame]);
	gtk_box_pack_start(GTK_BOX(hbox), stuff->entry, TRUE, TRUE, 3);
	gtk_widget_set_margin_bottom(hbox, 6.0);

	g_signal_connect(stuff->framenumchk, "toggled", G_CALLBACK(goto_toggle_sig), stuff);
	g_signal_connect(stuff->frameidchk, "toggled", G_CALLBACK(goto_toggle_sig), stuff);
	goto_toggle_sig(stuff->framenumchk, stuff);
	g_signal_connect(G_OBJECT(stuff->entry), "activate",
			 G_CALLBACK(goto_frame_activate_sig), window);
	g_signal_connect_data(G_OBJECT(window), "response",
	                      G_CALLBACK(goto_frame_response_sig), stuff,
	                      free_goto_frame_stuff, 0);

	gtk_widget_show_all(window);
	gtk_widget_grab_focus(GTK_WIDGET(stuff->entry));

	return 0;
}


/* File->Quit */
static gint quit_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	if ( proj->unsaved ) {
		if ( !confirm_exit(proj) ) return TRUE;
	}
	gtk_main_quit();
	return FALSE;
}


static gint save_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	if ( save_project(proj) == 0 ) {
		STATUS("Saved project.\n");
	} else {
		ERROR("Could not save project.\n");
	}
	return FALSE;
}


static struct gui_indexing_result *current_result(struct crystfelproject *proj)
{
	const char *results_name;
	struct gui_indexing_result *res;

	results_name = gtk_combo_box_get_active_id(GTK_COMBO_BOX(proj->results_combo));
	if ( strcmp(results_name, "crystfel-gui-internal") == 0 ) return NULL;

	res = find_indexing_result_by_name(proj, results_name);
	if ( res == NULL ) {
		ERROR("Couldn't find result '%s'\n", results_name);
	}

	return res;
}


static gint rescan_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	struct gui_indexing_result *res = current_result(proj);
	if ( res != NULL ) {
		update_result_index(res);
	}
	return FALSE;
}


static gint peakogram_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	struct gui_indexing_result *res = current_result(proj);
	if ( res != NULL ) {
		GError *error = NULL;
		const gchar *args[128];
		GSubprocess *sp;
		int i;
		args[0] = "peakogram-stream";
		for ( i=0; i<MIN(res->n_streams, 60); i++ ) {
			args[1+2*i] = "-i";
			args[2+2*i] = res->streams[i];
		}
		args[1+2*res->n_streams] = NULL;

		sp = g_subprocess_newv(args, G_SUBPROCESS_FLAGS_NONE, &error);
		if ( sp == NULL ) {
			ERROR("Failed to invoke peakogram-stream: %s\n",
			      error->message);
			g_error_free(error);
		}
	} else {
		ERROR("Select indexing result first!\n");
	}
	return FALSE;
}


static gint reset_zoom_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	crystfel_image_view_reset_zoom(CRYSTFEL_IMAGE_VIEW(proj->imageview));
	return FALSE;
}


static gint reset_range_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	crystfel_colour_scale_auto_range(CRYSTFEL_COLOUR_SCALE(proj->colscale));
	return FALSE;
}


static gint first_frame_sig(GtkWidget *widget,
                            struct crystfelproject *proj)
{
	proj->cur_frame = 0;
	update_imageview(proj);
	return FALSE;
}


static gint prev_frame_sig(GtkWidget *widget,
                           struct crystfelproject *proj)
{
	if ( proj->cur_frame == 0 ) return FALSE;
	proj->cur_frame--;
	update_imageview(proj);
	return FALSE;
}


static gint random_frame_sig(GtkWidget *widget,
                             GdkEventButton *event,
                             struct crystfelproject *proj)
{
	if ( event->state & GDK_SHIFT_MASK ) {
		if ( proj->n_random_history > 0 ) {
			proj->cur_frame = pop_random_frame(proj);
			update_imageview(proj);
		}
	} else {
		push_random_frame(proj, proj->cur_frame);
		proj->cur_frame = random()*proj->n_frames / RAND_MAX;
		update_imageview(proj);
	}
	return FALSE;
}


static gint next_frame_sig(GtkWidget *widget,
                           struct crystfelproject *proj)
{
	if ( proj->cur_frame == proj->n_frames - 1 ) return FALSE;
	proj->cur_frame++;
	update_imageview(proj);
	return FALSE;
}


static gint last_frame_sig(GtkWidget *widget,
                           struct crystfelproject *proj)
{
	proj->cur_frame = proj->n_frames - 1;
	update_imageview(proj);
	return FALSE;
}


static gint about_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *window;

	const gchar *authors[] = {
		"Thomas White <taw@physics.org>",
		NULL
	};

	/* Alphabetical by family name */
	const gchar *contributors[] = {
		"Steve Aplin <steve.aplin@desy.de>",
		"Andrew Aquila <andrew.aquila@cfel.de>",
		"Anton Barty <anton.barty@desy.de>",
		"Kenneth Beyerlein <kenneth.beyerlein@desy.de>",
		"Wolfgang Brehm <wolfgang.brehm@gmail.com>",
		"Robert Bücker <robert.buecker@cssb-hamburg.de>",
		"Fedor Chervinskii <fedor.chervinskii@gmail.com>",
		"Nicholas Devenish <ndevenish@gmail.com>",
		"Lorenzo Galli <lorenzo.galli@desy.de>",
		"Cornelius Gati <cornelius.gati@cfel.de>",
		"Yaroslav Gevorkov <yaroslav.gevorkov@desy.de>",
		"Helen Ginn <helen@strubi.ox.ac.uk>",
		"Thomas Grant <tgrant@hwi.buffalo.edu>",
		"Pascal Hogan-Lamarre <pascal.hogan.lamarre@mail.utoronto.ca>",
		"Richard Kirian <rkirian@asu.edu>",
		"Valerio Mariani <<alerio.mariani@desy.de>",
		"Andrew Martin <andrew.martin@desy.de>",
		"Omri Mor <omor1@asu.edu>",
		"Takanori Nakane <nakane.t@gmail.com>",
		"Karol Nass <karol.nass@desy.de>",
		"Nicolas Riebesel <nicolas.riebesel@tuhh.de>",
		"Mamoru Suzuki <mamoru.suzuki@protein.osaka-u.ac.jp>",
		"Alexandra Tolstikova <alexandra.tolstikova@desy.de>",
		"Parker de Waal <Parker.deWaal@vai.org>",
		"Keitaro Yamashita <k.yamashita@spring8.or.jp>",
		"Oleksandr Yefanov <oleksandr.yefanov@desy.de>",
		"Chun Hong Yoon <chun.hong.yoon@desy.de>",
		"Nadia Zatsepin <nadia.zatsepin@asu.edu>",
		NULL
	};

	const gchar *uthash[] = {
		"Troy D. Hanson",
		NULL
	};

	const gchar *gpl = "CrystFEL is free software: you can redistribute it and/or modify\n"
	                   "it under the terms of the GNU General Public License as published by\n"
	                   "the Free Software Foundation, either version 3 of the License, or\n"
	                   "(at your option) any later version.\n\n"
	                   "CrystFEL is distributed in the hope that it will be useful,\n"
	                   "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
	                   "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
	                   "GNU General Public License for more details.\n\n"
	                   "You should have received a copy of the GNU General Public License\n"
	                   "along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.";

	window = gtk_about_dialog_new();
	gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(proj->window));

	gtk_about_dialog_set_logo_icon_name(GTK_ABOUT_DIALOG(window), "crystfel");
	gtk_about_dialog_set_program_name(GTK_ABOUT_DIALOG(window),
	        "CrystFEL graphical user interface");
	gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(window),
	                             crystfel_version_string());
	gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG(window),
		"© 2020-2023 Deutsches Elektronen-Synchrotron DESY, "
		"a research centre of the Helmholtz Association.");
	gtk_about_dialog_set_website(GTK_ABOUT_DIALOG(window),
		"https://www.desy.de/~twhite/crystfel");

	gtk_about_dialog_set_authors(GTK_ABOUT_DIALOG(window), authors);
	gtk_about_dialog_add_credit_section(GTK_ABOUT_DIALOG(window),
	                                    "Contributors", contributors);
	gtk_about_dialog_add_credit_section(GTK_ABOUT_DIALOG(window),
	                                    "Incorporates uthash by", uthash);

	gtk_about_dialog_set_license(GTK_ABOUT_DIALOG(window), gpl);

	g_signal_connect(window, "response", G_CALLBACK(gtk_widget_destroy),
			 NULL);

	gtk_widget_show_all(window);

	return 0;
}


static gint image_info_clicked_sig(GtkWidget *widget,
                                   struct crystfelproject *proj)
{
	GtkWidget *popover;
	GtkWidget *grid;
	GtkWidget *label;
	char tmp[64];

	if ( proj->cur_image == NULL ) return FALSE;

	popover = gtk_popover_new(widget);
	gtk_popover_set_position(GTK_POPOVER(popover),
	                         GTK_POS_BOTTOM);

	grid = gtk_grid_new();
	gtk_grid_set_row_spacing(GTK_GRID(grid), 4);
	gtk_grid_set_column_spacing(GTK_GRID(grid), 4);
	gtk_container_set_border_width(GTK_CONTAINER(grid), 6);

	label = gtk_label_new("Number of peaks:");
	gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
	snprintf(tmp, 63, "%i", image_feature_count(proj->cur_image->features));
	label = gtk_label_new(tmp);
	gtk_grid_attach(GTK_GRID(grid), label, 1, 0, 1, 1);

	label = gtk_label_new("Number of crystals:");
	gtk_grid_attach(GTK_GRID(grid), label, 0, 1, 1, 1);
	snprintf(tmp, 63, "%i", proj->cur_image->n_crystals);
	label = gtk_label_new(tmp);
	gtk_grid_attach(GTK_GRID(grid), label, 1, 1, 1, 1);

	gtk_container_add(GTK_CONTAINER(popover), grid);
	gtk_widget_show_all(grid);

#if GTK_CHECK_VERSION(3,22,0)
	gtk_popover_popup(GTK_POPOVER(popover));
#else
	gtk_widget_show_all(GTK_WIDGET(popover));
#endif

	return FALSE;
}


static gint results_combo_changed_sig(GtkComboBox *w,
                                      struct crystfelproject *proj)
{
	update_imageview(proj);
	return FALSE;
}


static gint rescan_on_change_sig(GtkWidget *w, struct crystfelproject *proj)
{
	proj->rescan_on_change = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(w));
	return FALSE;
}


static gint show_centre_sig(GtkWidget *w, struct crystfelproject *proj)
{
	proj->show_centre = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(w));
	crystfel_image_view_set_show_centre(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                    proj->show_centre);
	return FALSE;
}


static gint resolution_rings_sig(GtkWidget *w, struct crystfelproject *proj)
{
	int tr = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(w));

	if ( tr && proj->cur_image == NULL ) {
		error_box(proj, "Load some images first!\n");
		return FALSE;
	}

	if ( tr && isnan(detgeom_mean_camera_length(proj->cur_image->detgeom)) ) {
		error_box(proj, "Cannot display resolution rings accurately "
		          "because detector is not flat enough");
		tr = FALSE;
		gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(w), FALSE);
	}

	proj->resolution_rings = tr;
	crystfel_image_view_set_resolution_rings(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                         tr);

	return FALSE;
}


static gint show_peaks_sig(GtkWidget *w, struct crystfelproject *proj)
{
	proj->show_peaks = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(w));
	crystfel_image_view_set_show_peaks(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                   proj->show_peaks);
	return FALSE;
}


static gint show_refls_sig(GtkWidget *w, struct crystfelproject *proj)
{
	proj->show_refls = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(w));
	crystfel_image_view_set_show_reflections(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                         proj->show_refls);
	return FALSE;
}


static gint label_refls_sig(GtkWidget *w, struct crystfelproject *proj)
{
	proj->label_refls = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(w));
	crystfel_image_view_set_label_reflections(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                          proj->label_refls);
	return FALSE;
}


static void add_menu_bar(struct crystfelproject *proj, GtkWidget *vbox)
{
	GError *error = NULL;

	const char *ui = "<ui> <menubar name=\"mainwindow\">"
		"<menu name=\"file\" action=\"FileAction\">"
		"	<menuitem name=\"save\" action=\"SaveAction\" />"
		"	<menuitem name=\"quit\" action=\"QuitAction\" />"
		"</menu>"
		"<menu name=\"view\" action=\"ViewAction\" >"
		"	<menuitem name=\"peaks\" action=\"PeaksAction\" />"
		"	<menuitem name=\"refls\" action=\"ReflsAction\" />"
		"	<menuitem name=\"labelrefls\" action=\"LabelReflsAction\" />"
		"	<menuitem name=\"centre\" action=\"CentreAction\" />"
		"	<menuitem name=\"resrings\" action=\"ResolutionRingsAction\" />"
		"       <separator />"
		"	<menuitem name=\"resetzoom\" action=\"ResetZoomAction\" />"
		"	<menuitem name=\"resetrange\" action=\"ResetRangeAction\" />"
		"</menu>"
		"<menu name=\"tools\" action=\"ToolsAction\" >"
		"	<menuitem name=\"rescanonchange\" action=\"RescanOnChangeAction\" />"
		"	<menuitem name=\"rescan\" action=\"RescanAction\" />"
		"	<menuitem name=\"jumpframe\" action=\"JumpFrameAction\" />"
		"       <separator />"
		"	<menuitem name=\"peakogram\" action=\"PeakogramAction\" />"
		"</menu>"
		"<menu name=\"help\" action=\"HelpAction\">"
		"	<menuitem name=\"about\" action=\"AboutAction\" />"
		"</menu>"
		"</menubar></ui>";

	GtkActionEntry entries[] = {

		{ "FileAction", NULL, "_File", NULL, NULL, NULL },
		{ "SaveAction", GTK_STOCK_SAVE, "_Save", NULL, NULL,
			G_CALLBACK(save_sig) },
		{ "QuitAction", GTK_STOCK_QUIT, "_Quit", NULL, NULL,
			G_CALLBACK(quit_sig) },

		{ "ViewAction", NULL, "_View", NULL, NULL, NULL },
		{ "ResetZoomAction", NULL, "Reset zoom", NULL, NULL,
			G_CALLBACK(reset_zoom_sig) },
		{ "ResetRangeAction", NULL, "Reset colour scale", NULL, NULL,
			G_CALLBACK(reset_range_sig) },

		{ "ToolsAction", NULL, "_Tools", NULL, NULL, NULL },
		{ "RescanAction", NULL, "Rescan streams", NULL, NULL,
			G_CALLBACK(rescan_sig) },
		{ "JumpFrameAction", NULL, "Jump to frame", NULL, NULL,
			G_CALLBACK(goto_frame_sig) },
		{ "PeakogramAction", NULL, "Check detector saturation (peakogram)", NULL, NULL,
			G_CALLBACK(peakogram_sig) },

		{ "HelpAction", NULL, "_Help", NULL, NULL, NULL },
		{ "AboutAction", GTK_STOCK_ABOUT, "_About", NULL, NULL,
			G_CALLBACK(about_sig) },

	};

	GtkToggleActionEntry toggles[] = {
		{ "PeaksAction", NULL, "Peak detection results", NULL, NULL,
		  G_CALLBACK(show_peaks_sig), FALSE },
		{ "ReflsAction", NULL, "Calculated reflection positions", NULL, NULL,
		  G_CALLBACK(show_refls_sig), FALSE },
		{ "LabelReflsAction", NULL, "Show reflection indices", NULL, NULL,
		  G_CALLBACK(label_refls_sig), FALSE },
		{ "CentreAction", NULL, "Beam centre", NULL, NULL,
		  G_CALLBACK(show_centre_sig), FALSE },
		{ "ResolutionRingsAction", NULL, "Resolution rings", NULL, NULL,
		  G_CALLBACK(resolution_rings_sig), FALSE },
		{ "RescanOnChangeAction", NULL, "Rescan streams when changing frame", NULL, NULL,
		  G_CALLBACK(rescan_on_change_sig), FALSE },
	};

	proj->action_group = gtk_action_group_new("cellwindow");
	gtk_action_group_add_actions(proj->action_group, entries,
	                             G_N_ELEMENTS(entries), proj);
	gtk_action_group_add_toggle_actions(proj->action_group, toggles,
	                                    G_N_ELEMENTS(toggles), proj);

	proj->ui = gtk_ui_manager_new();
	gtk_ui_manager_insert_action_group(proj->ui, proj->action_group, 0);
	g_signal_connect(proj->ui, "add_widget", G_CALLBACK(add_ui_sig), vbox);
	if ( gtk_ui_manager_add_ui_from_string(proj->ui, ui, -1, &error) == 0 )
	{
		fprintf(stderr, "Error loading message window menu bar: %s\n",
			error->message);
		return;
	}

	gtk_window_add_accel_group(GTK_WINDOW(proj->window),
				   gtk_ui_manager_get_accel_group(proj->ui));
	gtk_ui_manager_ensure_update(proj->ui);
}


static void add_button(GtkWidget *vbox, const char *label, const char *imagen,
                       GCallback callback, struct crystfelproject *proj)
{
	GtkWidget *button;
	GtkWidget *image;

	button = gtk_button_new_with_label(label);
	g_object_set(G_OBJECT(button), "image-position", GTK_POS_TOP, NULL);
	image = gtk_image_new_from_resource(imagen);
	g_object_set(G_OBJECT(button), "image", image, NULL);
	g_object_set(G_OBJECT(button), "always-show-image", TRUE, NULL);
	gtk_button_set_relief(GTK_BUTTON(button), GTK_RELIEF_NONE);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(button), FALSE, FALSE, 4);

	if ( callback != NULL ) {
		g_signal_connect(G_OBJECT(button), "clicked",
		                 G_CALLBACK(callback), proj);
	}
}


static void add_task_buttons(GtkWidget *vbox, struct crystfelproject *proj)
{
	add_button(vbox, "Load data", "/de/desy/CrystFEL/crystfel-load-data.svg",
	           G_CALLBACK(import_sig), proj);
	add_button(vbox, "Peak detection", "/de/desy/CrystFEL/crystfel-peak-detection.svg",
	           G_CALLBACK(peaksearch_sig), proj);
	add_button(vbox, "Index this frame", "/de/desy/CrystFEL/crystfel-index-frame.svg",
	           G_CALLBACK(index_one_sig), proj);
	add_button(vbox, "Index all frames", "/de/desy/CrystFEL/crystfel-index-all.svg",
	           G_CALLBACK(index_all_sig), proj);
	add_button(vbox, "Determine unit cell", "/de/desy/CrystFEL/crystfel-unitcell.svg",
	           G_CALLBACK(cell_explorer_sig), proj);
	add_button(vbox, "Refine detector geometry", "/de/desy/CrystFEL/crystfel-geometry.svg",
	           G_CALLBACK(align_sig), proj);
	add_button(vbox, "Indexing ambiguity", "/de/desy/CrystFEL/crystfel-ambiguity.svg",
	           G_CALLBACK(ambi_sig), proj);
	add_button(vbox, "Merge", "/de/desy/CrystFEL/crystfel-merge.svg",
	           G_CALLBACK(merge_sig), proj);
	add_button(vbox, "Figures of merit", "/de/desy/CrystFEL/crystfel-fom.svg",
	           G_CALLBACK(fom_sig), proj);
	add_button(vbox, "Export data", "/de/desy/CrystFEL/crystfel-export.svg",
	           G_CALLBACK(export_sig), proj);
}


static void add_gui_message(enum log_msg_type type, const char *msg,
                            void *vp)
{
	GtkTextBuffer *buf;
	GtkTextIter iter;
	GtkTextMark *mark;
	struct crystfelproject *proj = vp;

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(proj->report));
	gtk_text_buffer_get_end_iter(buf, &iter);
	gtk_text_buffer_insert(buf, &iter, msg, -1);

	mark = gtk_text_mark_new(NULL, FALSE);
	gtk_text_buffer_get_end_iter(buf, &iter);
	gtk_text_buffer_add_mark(buf, mark, &iter);
	gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(proj->report),
	                             mark, 0.0, FALSE, 0.0, 0.0);
	gtk_text_buffer_delete_mark(buf, mark);
}


static void clear_log_sig(GtkMenuItem *widget,
                          struct crystfelproject *proj)
{
	GtkTextBuffer *buf;

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(proj->report));
	if ( buf != NULL ) {
		gtk_text_buffer_set_text(buf, "", -1);
	}
}


static void colscale_change_sig(CrystFELColourScale *colscale,
                                struct crystfelproject *proj)
{
	double lo, hi;
	crystfel_colour_scale_get_range(colscale, &lo, &hi);
	crystfel_image_view_set_colour_scale(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                     lo, hi);
	if ( proj->cur_image != NULL ) {
		proj->range_set = 1;
	}
}


static void add_log_menu_items(GtkTextView *textview,
                               GtkWidget *popup,
                               struct crystfelproject *proj)
{
	GtkWidget *item;

	if ( !GTK_IS_MENU(popup) ) return;

	item = gtk_separator_menu_item_new();
	gtk_menu_shell_append(GTK_MENU_SHELL(popup), item);
	gtk_widget_show(item);

	item = gtk_menu_item_new_with_label("Clear log");
	gtk_menu_shell_append(GTK_MENU_SHELL(popup), item);
	gtk_widget_show(item);
	g_signal_connect(item, "activate", G_CALLBACK(clear_log_sig), proj);
}


static void try_load_geom_from_result(struct crystfelproject *proj)
{
	struct gui_indexing_result *res;
	Stream *st;
	char *geom_str;
	const char *res_name = selected_result(proj);

	STATUS("No geometry file specified in project file.\n");
	if ( strcmp(res_name, "crystfel-gui-internal") == 0 ) {
		if ( proj->n_results == 0 ) {
			ERROR("No indexing results found.  "
			      "No geometry information!\n");
			return;
		}
		res_name = proj->results[0].name;
	}

	res = find_indexing_result_by_name(proj, res_name);
	if ( res == NULL ) {
		ERROR("Couldn't find results '%s'\n", res_name);
		return;
	}

	st = stream_open_for_read(res->streams[0]);
	if ( st == NULL ) {
		ERROR("Couldn't open stream '%s' from result '%s'\n",
		      res->streams[0], res_name);
		return;
	}

	geom_str = stream_geometry_file(st);
	if ( geom_str == NULL ) {
		ERROR("No geometry in stream '%s' from result '%s'\n",
		      res->streams[0], res_name);
		return;
	}

	proj->dtempl = data_template_new_from_string(geom_str);
	stream_close(st);

	STATUS("Using geometry from result '%s'\n", res_name);
}


int main(int argc, char *argv[])
{
	int c;
	struct crystfelproject proj;
	GtkWidget *vbox;
	GtkWidget *vpaned;
	GtkWidget *hpaned;
	GtkWidget *scroll;
	GtkWidget *frame;
	GtkWidget *main_vbox;
	GtkWidget *iv_hbox;
	GtkWidget *toolbar;
	GtkWidget *results_toolbar;
	GtkWidget *button;
	GtkWidget *label;
	int load_result;
	GtkAction *act;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,                1 },
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "h",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 1 :
			printf("CrystFEL: %s\n",
			       crystfel_version_string());
			printf("%s\n",
			       crystfel_licence_string());
			return 0;

			default :
			return 1;

		}

	}

	/* This isn't great, but necessary to make the command-line UI and file
	 * formats consistent with the other programs, which all use the C
	 * locale.  Better would be to have all the programs call
	 * setlocale(LC_ALL, "") and use the C locale temporarily when reading
	 * or writing a stream, reflection file, geometry file etc. */
	gtk_disable_setlocale();

	gtk_init(&argc, &argv);
	g_set_application_name("CrystFEL");

	gtk_icon_theme_add_resource_path(gtk_icon_theme_get_default(),
	                                 "/de/desy/CrystFEL");

	if ( default_project(&proj) ) {
		ERROR("Failed to set up default project\n");
		return 1;
	}

	proj.image_info = NULL;

	proj.window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(proj.window), "CrystFEL");
	gtk_window_set_default_icon_name("crystfel");
	g_signal_connect(G_OBJECT(proj.window), "delete-event",
	                 G_CALLBACK(delete_event_sig),
	                 &proj);

	vbox = gtk_vbox_new(FALSE, 0.0);
	gtk_container_add(GTK_CONTAINER(proj.window), vbox);
	add_menu_bar(&proj, vbox);

	vpaned = gtk_paned_new(GTK_ORIENTATION_VERTICAL);
	gtk_box_pack_end(GTK_BOX(vbox), vpaned, TRUE, TRUE, 0.0);

	hpaned = gtk_paned_new(GTK_ORIENTATION_HORIZONTAL);
	gtk_paned_pack1(GTK_PANED(vpaned), hpaned, TRUE, TRUE);

	toolbar = gtk_hbox_new(FALSE, 0.0);

	/* First */
	proj.first_button = gtk_button_new_from_icon_name("go-first", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), proj.first_button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(proj.first_button), "clicked",
	                 G_CALLBACK(first_frame_sig), &proj);
	gtk_widget_set_tooltip_text(proj.first_button, "First image");

	/* Prev */
	proj.prev_button = gtk_button_new_from_icon_name("go-previous", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), proj.prev_button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(proj.prev_button), "clicked",
	                 G_CALLBACK(prev_frame_sig), &proj);
	gtk_widget_set_tooltip_text(proj.prev_button, "Previous image");

	/* Random */
	button = gtk_button_new_from_icon_name("media-playlist-shuffle",
	                                       GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(button), "button-press-event",
	                 G_CALLBACK(random_frame_sig), &proj);
	gtk_widget_set_tooltip_text(button, "Jump to a random image.  "
	                                    "Shift-click to go back");

	/* Next */
	proj.next_button = gtk_button_new_from_icon_name("go-next", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), proj.next_button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(proj.next_button), "clicked",
	                 G_CALLBACK(next_frame_sig), &proj);
	gtk_widget_set_tooltip_text(proj.next_button, "Next image");

	/* Last */
	proj.last_button = gtk_button_new_from_icon_name("go-last", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), proj.last_button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(proj.last_button), "clicked",
	                 G_CALLBACK(last_frame_sig), &proj);
	gtk_widget_set_tooltip_text(proj.last_button, "Last image");

	/* Information about image */
	button = gtk_button_new_from_icon_name("document-properties",
	                                       GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_end(GTK_BOX(toolbar), button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(button), "clicked",
	                 G_CALLBACK(image_info_clicked_sig), &proj);

	results_toolbar = gtk_hbox_new(FALSE, 0.0);
	label = gtk_label_new("Show results from:");
	gtk_box_pack_start(GTK_BOX(results_toolbar), label,
	                   FALSE, FALSE, 4.0);
	proj.results_combo = gtk_combo_box_text_new();
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(proj.results_combo),
	                          "crystfel-gui-internal",
	                          "Calculations within GUI");
	gtk_combo_box_set_active(GTK_COMBO_BOX(proj.results_combo), 0);
	gtk_box_pack_start(GTK_BOX(results_toolbar), proj.results_combo,
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(proj.results_combo), "changed",
	                 G_CALLBACK(results_combo_changed_sig), &proj);

	/* Filename */
	proj.image_info = gtk_label_new("Ready to load images");
	gtk_label_set_selectable(GTK_LABEL(proj.image_info), TRUE);
	gtk_label_set_ellipsize(GTK_LABEL(proj.image_info),
	                        PANGO_ELLIPSIZE_START);
	gtk_box_pack_end(GTK_BOX(toolbar), proj.image_info, TRUE, TRUE, 0.0);

	main_vbox = gtk_vbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(main_vbox), toolbar,
	                   FALSE, FALSE, 2.0);
	gtk_box_pack_start(GTK_BOX(main_vbox), results_toolbar,
	                   FALSE, FALSE, 2.0);

	/* Main area stuff (toolbar and imageview) at right */
	proj.imageview = crystfel_image_view_new();
	proj.cur_frame = 0;
	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);

	iv_hbox = gtk_hbox_new(FALSE, 0.0);
	scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
	                               GTK_POLICY_ALWAYS, GTK_POLICY_ALWAYS);
	gtk_container_add(GTK_CONTAINER(scroll), GTK_WIDGET(proj.imageview));
	gtk_box_pack_start(GTK_BOX(main_vbox), iv_hbox, TRUE, TRUE, 0.0);
	gtk_box_pack_start(GTK_BOX(iv_hbox), scroll, TRUE, TRUE, 0.0);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(main_vbox));
	gtk_paned_pack2(GTK_PANED(hpaned), GTK_WIDGET(frame), TRUE, TRUE);
	proj.main_vbox = main_vbox;

	proj.colscale = crystfel_colour_scale_new();
	gtk_box_pack_start(GTK_BOX(iv_hbox), proj.colscale, FALSE, FALSE, 0.0);
	g_signal_connect(proj.colscale, "range-changed",
	                 G_CALLBACK(colscale_change_sig), &proj);

	/* Icon region at left */
	proj.icons = gtk_vbox_new(FALSE, 0.0);
	scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_container_set_border_width(GTK_CONTAINER(proj.icons), 16);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
	                               GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(scroll));
	gtk_container_add(GTK_CONTAINER(scroll), GTK_WIDGET(proj.icons));
	gtk_paned_pack1(GTK_PANED(hpaned), GTK_WIDGET(frame), FALSE, FALSE);
	add_task_buttons(proj.icons, &proj);

	/* Report (text) region at bottom */
	proj.report = gtk_text_view_new();
	gtk_text_view_set_editable(GTK_TEXT_VIEW(proj.report), FALSE);
	gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(proj.report), FALSE);
	scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
	                               GTK_POLICY_AUTOMATIC, GTK_POLICY_ALWAYS);
	g_signal_connect(proj.report, "populate-popup",
	                 G_CALLBACK(add_log_menu_items), &proj);
	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(scroll));
	gtk_container_add(GTK_CONTAINER(scroll), GTK_WIDGET(proj.report));
	gtk_paned_pack2(GTK_PANED(vpaned), GTK_WIDGET(frame), FALSE, FALSE);

	/* Send messages to report region */
	set_log_message_func(add_gui_message, &proj);

	if ( optind < argc ) {
		/* Create view of stream - probably temporary */
		load_result = load_stream(&proj, argv[optind++]);
	} else {
		/* Try to load state from disk */
		load_result = load_project(&proj);
	}

	if ( load_result == 0 ) {
		DataTemplate *dtempl;
		proj.cur_frame = 0;

		if ( proj.geom_filename != NULL ) {

			dtempl = data_template_new_from_file(proj.geom_filename);
			if ( dtempl != NULL ) {
				proj.dtempl = dtempl;
			}
		} else {
			try_load_geom_from_result(&proj);
		}

		update_imageview(&proj);
	}

	act = gtk_ui_manager_get_action(proj.ui, "/mainwindow/view/centre");
	gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(act),
	                             proj.show_centre);

	act = gtk_ui_manager_get_action(proj.ui, "/mainwindow/view/peaks");
	gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(act),
	                             proj.show_peaks);

	act = gtk_ui_manager_get_action(proj.ui, "/mainwindow/view/refls");
	gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(act),
	                             proj.show_refls);

	act = gtk_ui_manager_get_action(proj.ui, "/mainwindow/view/labelrefls");
	gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(act),
	                             proj.label_refls);

	act = gtk_ui_manager_get_action(proj.ui, "/mainwindow/view/resrings");
	gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(act),
	                             proj.resolution_rings);

	act = gtk_ui_manager_get_action(proj.ui, "/mainwindow/tools/rescanonchange");
	gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(act),
	                             proj.rescan_on_change);

	gtk_window_set_default_size(GTK_WINDOW(proj.window), 1024, 768);
	gtk_paned_set_position(GTK_PANED(hpaned), 172);
	gtk_paned_set_position(GTK_PANED(vpaned), 600);
	gtk_widget_show_all(proj.window);
	gtk_main();

	return 0;
}


struct infobar_data
{
	struct crystfelproject *proj;
	struct gui_task *task;
};


static void free_ib_callback_params(gpointer cbvals,
                                    GClosure *closure)
{
	free(cbvals);
}


static void remove_task(struct crystfelproject *proj,
                        struct gui_task *task)
{
	if ( task->running ) {
		ERROR("Attempt to remove a running task!\n");
		return;
	}

	if ( task->backend->free_task != NULL ) {
		task->backend->free_task(task->job_priv);
	}

	proj->tasks = g_slist_remove(proj->tasks, task);
	free(task);
}


static void infobar_response_sig(GtkInfoBar *infobar, gint resp,
                                 gpointer data)
{
	struct infobar_data *ibdata = data;

	if ( resp == GTK_RESPONSE_CANCEL ) {
		ibdata->task->backend->cancel_task(ibdata->task->job_priv);

	} else if ( resp == GTK_RESPONSE_CLOSE ) {

		gtk_widget_destroy(GTK_WIDGET(infobar));
		remove_task(ibdata->proj, ibdata->task);

	} else {
		ERROR("Unrecognised infobar response!\n");
	}
}


static gboolean update_info_bar(gpointer data)
{
	struct gui_task *task = data;
	int running;
	float frac_complete;

	if ( task->backend->task_status(task->job_priv,
	                                &running, &frac_complete) ) {
		ERROR("Error retrieving task status\n");
		return G_SOURCE_CONTINUE;
	}

	gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(task->progress_bar),
	                              frac_complete);

	if ( !running && task->running ) {

		int i;

		/* Task is no longer running */
		task->running = 0;
		gtk_widget_destroy(task->cancel_button);
		gtk_info_bar_set_show_close_button(GTK_INFO_BAR(task->info_bar),
		                                   TRUE);

		/* We don't have an easy way to get the result name from the
		 * task structure, so cheat by marking all of the results as
		 * "possibly out of date" */
		for ( i=0; i<task->proj->n_results; i++ ) {
			task->proj->results[i].need_rescan = 1;
		}

		return G_SOURCE_REMOVE;
	}

	return G_SOURCE_CONTINUE;
}


void add_running_task(struct crystfelproject *proj,
                      const char *task_desc,
                      struct crystfel_backend *backend,
                      void *job_priv)
{
	struct gui_task *task;
	struct infobar_data *ibdata;
	GtkWidget *bar_area;

	task = malloc(sizeof(struct gui_task));
	if ( task == NULL ) return;

	task->job_priv = job_priv;
	task->backend = backend;
	task->running = 1;
	task->proj = proj;
	proj->tasks = g_slist_append(proj->tasks, task);

	/* Progress info bar */
	task->info_bar = gtk_info_bar_new();
	gtk_info_bar_set_message_type(GTK_INFO_BAR(task->info_bar),
	                              GTK_MESSAGE_INFO);

	task->cancel_button = gtk_info_bar_add_button(GTK_INFO_BAR(task->info_bar),
	                                              GTK_STOCK_CANCEL,
	                                              GTK_RESPONSE_CANCEL);

	gtk_box_pack_end(GTK_BOX(proj->main_vbox), GTK_WIDGET(task->info_bar),
	                 FALSE, FALSE, 0.0);

	bar_area = gtk_info_bar_get_content_area(GTK_INFO_BAR(task->info_bar));

	/* Create progress bar */
	task->progress_bar = gtk_progress_bar_new();
	gtk_box_pack_start(GTK_BOX(bar_area),
	                   GTK_WIDGET(task->progress_bar),
	                   TRUE, TRUE, 0.0);
	gtk_progress_bar_set_text(GTK_PROGRESS_BAR(task->progress_bar),
	                          task_desc);
	gtk_progress_bar_set_show_text(GTK_PROGRESS_BAR(task->progress_bar),
	                               TRUE);

	ibdata = malloc(sizeof(struct infobar_data));
	if ( ibdata != NULL ) {
		ibdata->proj = proj;
		ibdata->task = task;
		g_signal_connect_data(G_OBJECT(task->info_bar), "response",
		                      G_CALLBACK(infobar_response_sig), ibdata,
		                      free_ib_callback_params, 0);
	}

	gtk_widget_show_all(task->info_bar);

#if GTK_CHECK_VERSION(3,22,29)
	gtk_info_bar_set_revealed(GTK_INFO_BAR(task->info_bar), TRUE);
#endif

	g_timeout_add(2000, update_info_bar, task);
}


static GFile *get_crystfel_path_gfile()
{
	GFile *self;
	GFileInfo *self_info;
	const char *self_target;
	GFile *tar;
	GFile *parent_dir;
	GError *error = NULL;

	self = g_file_new_for_path("/proc/self/exe");
	self_info = g_file_query_info(self, "standard",
	                              G_FILE_QUERY_INFO_NONE,
	                              NULL, &error);
	if ( self_info == NULL ) return NULL;

	self_target = g_file_info_get_symlink_target(self_info);
	if ( self_target == NULL ) return NULL;

	tar = g_file_new_for_path(self_target);
	if ( tar == NULL ) return NULL;

	parent_dir = g_file_get_parent(tar);
	if ( parent_dir == NULL ) return NULL;

	g_object_unref(self);
	g_object_unref(self_info);
	g_object_unref(tar);

	return parent_dir;
}


char *get_crystfel_path_str()
{
	char *path;
	GFile *crystfel_path = get_crystfel_path_gfile();
	if ( crystfel_path == NULL ) return NULL;
	path = g_file_get_path(crystfel_path);
	g_object_unref(crystfel_path);
	return path;
}


struct gui_job_notes_page *add_job_notes_page(GtkWidget *notebook)
{
	GtkWidget *box;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *scroll;
	struct gui_job_notes_page *notes;

	notes = malloc(sizeof(struct gui_job_notes_page));
	if ( notes == NULL ) return NULL;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	label = gtk_label_new("Whatever you enter here will be placed in "
	                      "the job's folder as 'notes.txt'");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   TRUE, TRUE, 0);
	notes->textview = gtk_text_view_new();
	gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(notes->textview),
	                            GTK_WRAP_WORD_CHAR);
	scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_container_add(GTK_CONTAINER(scroll), notes->textview);
	gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scroll),
	                                    GTK_SHADOW_ETCHED_IN);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(scroll),
	                   TRUE, TRUE, 2.0);

	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box,
	                         gtk_label_new("Notes"));
	return notes;
}


static int job_name_valid(const char *job_title)
{
	if ( strchr(job_title, '/') != NULL ) return 0;
	if ( strchr(job_title, '\\') != NULL ) return 0;
	return 1;
}


GFile *make_job_folder(const char *job_title, const char *job_notes)
{
	struct stat s;
	char *workdir;
	GFile *workdir_file;
	GFile *cwd_file;
	GFile *notes_file;
	char *notes_path;
	FILE *fh;

	if ( !job_name_valid(job_title) ) {
		ERROR("Invalid job name '%s'\n", job_title);
		return NULL;
	}

	workdir = strdup(job_title);
	if ( workdir == NULL ) return NULL;

	if ( stat(workdir, &s) != -1 ) {
		ERROR("Working directory already exists.  "
		      "Choose a different job name.\n");
		return NULL;
	}

	if ( mkdir(workdir, S_IRWXU) ) {
		ERROR("Failed to create working directory: %s\n",
		      strerror(errno));
		return NULL;
	}

	cwd_file = g_file_new_for_path(".");
	workdir_file = g_file_get_child(cwd_file, workdir);
	g_object_unref(cwd_file);

	/* Write the notes into notes.txt */
	notes_file = g_file_get_child(workdir_file, "notes.txt");
	notes_path = g_file_get_path(notes_file);
	fh = fopen(notes_path, "w");
	fputs(job_notes, fh);
	fclose(fh);
	g_free(notes_path);
	g_object_unref(notes_file);

	return workdir_file;
}


void force_peaks_on(struct crystfelproject *proj)
{
	GtkWidget *w;
	proj->show_peaks = 1;
	crystfel_image_view_set_show_peaks(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                   proj->show_peaks);

	w =  gtk_ui_manager_get_widget(proj->ui, "/ui/mainwindow/view/peaks");
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(w), 1);
}


void force_refls_on(struct crystfelproject *proj)
{
	GtkWidget *w;
	proj->show_refls = 1;
	crystfel_image_view_set_show_reflections(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                         proj->show_refls);

	w =  gtk_ui_manager_get_widget(proj->ui, "/ui/mainwindow/view/refls");
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(w), 1);
}


/* Given an "old" job title (possibly NULL), generate a new job title */
char *make_new_job_title(const char *orig_old_title)
{
	size_t len, i;
	char *old_title;

	if ( orig_old_title == NULL ) return NULL;

	old_title = strdup(orig_old_title);
	if ( old_title == NULL ) return NULL;

	len = strlen(old_title);

	for ( i=len-1; i>0; i-- ) {
		if ( !isdigit(old_title[i]) ) break;
	}
	if ( i == len-1 ) {
		/* No digits at end */
		char *new_title = malloc(len+3);
		if ( new_title == NULL ) return NULL;
		strcpy(new_title, old_title);
		strcat(new_title, "-2");
		return new_title;
	} else {
		/* Digits at end */
		int n = atoi(&old_title[i+1]);
		char *new_title = malloc(len+6);
		if ( new_title == NULL ) return NULL;
		old_title[i+1] = '\0';
		snprintf(new_title, len+6, "%s%i", old_title, n+1);
		return new_title;
	}
}


char *relative_to_cwd(GFile *workdir, const char *filename)
{
	GFile *current_dir;
	GFile *gfile;
	char *rel;

	current_dir = g_file_new_for_path(".");

	gfile = g_file_get_child(workdir, filename);
	rel = g_file_get_relative_path(current_dir, gfile);
	g_object_unref(gfile);
	g_object_unref(current_dir);

	return rel;
}
