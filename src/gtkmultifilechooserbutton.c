/*
 * gtkmultifilechooserbutton.c
 *
 * A GTK widget to select multiple files
 *
 * Copyright © 2021 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2021 Thomas White <taw@physics.org>
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

#include "gtkmultifilechooserbutton.h"


G_DEFINE_TYPE(GtkMultiFileChooserButton,
              gtk_multi_file_chooser_button,
              GTK_TYPE_BUTTON)


static void gtk_multi_file_chooser_button_finalize(GObject *obj)
{
	GtkMultiFileChooserButton *mfc = GTK_MULTI_FILE_CHOOSER_BUTTON(obj);
	free(mfc->text);
	G_OBJECT_CLASS(gtk_multi_file_chooser_button_parent_class)->finalize(obj);
}


static void gtk_multi_file_chooser_button_class_init(GtkMultiFileChooserButtonClass *klass)
{
	GObjectClass *object_class = G_OBJECT_CLASS(klass);
	object_class->finalize = gtk_multi_file_chooser_button_finalize;
}


static void gtk_multi_file_chooser_button_init(GtkMultiFileChooserButton *mfc)
{
	mfc->chooser = NULL;
	mfc->filenames = NULL;
}


static void chooser_destroy_sig(GtkWidget *widget, GtkMultiFileChooserButton *mfc)
{
	mfc->chooser = NULL;
}


static void delstring(gpointer data, gpointer user_data)
{
	g_free(data);
}


static void chooser_response_sig(GtkWidget *dialog, gint resp,
                                 GtkMultiFileChooserButton *mfc)
{
	if ( resp == GTK_RESPONSE_ACCEPT ) {

		int n;
		char tmp[64];

		if ( mfc->filenames != NULL ) {
			g_slist_foreach(mfc->filenames, delstring, NULL);
			g_slist_free(mfc->filenames);
		}

		mfc->filenames = gtk_file_chooser_get_filenames(GTK_FILE_CHOOSER(dialog));

		n = g_slist_length(mfc->filenames);
		snprintf(tmp, 63, "%i files selected", n);
		gtk_button_set_label(GTK_BUTTON(mfc), tmp);

	}

	gtk_widget_destroy(dialog);
}


static void click_sig(GtkMultiFileChooserButton *mfc, gpointer data)
{
	GtkWidget *parent;

	if ( mfc->chooser != NULL ) return;

	parent = gtk_widget_get_toplevel(GTK_WIDGET(mfc));
	if ( !GTK_IS_WINDOW(parent) ) {
		parent = NULL;
	}
	mfc->chooser = gtk_file_chooser_dialog_new(mfc->text,
	                                           GTK_WINDOW(parent),
	                                           GTK_FILE_CHOOSER_ACTION_OPEN,
	                                           GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
	                                           "Select", GTK_RESPONSE_ACCEPT,
	                                           NULL);
	gtk_file_chooser_set_select_multiple(GTK_FILE_CHOOSER(mfc->chooser), TRUE);
	g_signal_connect(mfc->chooser, "destroy",
	                 G_CALLBACK(chooser_destroy_sig), mfc);
	g_signal_connect(G_OBJECT(mfc->chooser), "response",
	                 G_CALLBACK(chooser_response_sig), mfc);
	gtk_widget_show(mfc->chooser);
}


GtkWidget *gtk_multi_file_chooser_button_new(const char *text)
{
	GtkMultiFileChooserButton *mfc;

	mfc = g_object_new(GTK_TYPE_MULTI_FILE_CHOOSER_BUTTON, NULL);
	mfc->text = strdup(text);
	gtk_button_set_label(GTK_BUTTON(mfc), mfc->text);
	g_signal_connect(mfc, "clicked", G_CALLBACK(click_sig), NULL);

	return GTK_WIDGET(mfc);
}


GSList *gtk_multi_file_chooser_button_get_filenames(GtkMultiFileChooserButton *mfc)
{
	return mfc->filenames;
}
