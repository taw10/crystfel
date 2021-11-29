/*
 * gtkmultifilechooserbutton.h
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

#ifndef GTKMULTIFILECHOOSERBUTTON_H
#define GTKMULTIFILECHOOSERBUTTON_H

#include <gtk/gtk.h>
#include <glib-object.h>

#define GTK_TYPE_MULTI_FILE_CHOOSER_BUTTON (gtk_multi_file_chooser_button_get_type())

#define GTK_MULTI_FILE_CHOOSER_BUTTON(obj) (G_TYPE_CHECK_INSTANCE_CAST((obj), \
                                 GTK_TYPE_MULTI_FILE_CHOOSER_BUTTON, GtkMultiFileChooserButton))

#define GTK_IS_MULTI_FILE_CHOOSER_BUTTON(obj) (G_TYPE_CHECK_INSTANCE_TYPE((obj), \
                                    GTK_TYPE_MULTI_FILE_CHOOSER_BUTTON))

#define GTK_MULTI_FILE_CHOOSER_BUTTON_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST((obj), \
                                         GTK_TYPE_MULTI_FILE_CHOOSER_BUTTON, GtkMultiFileChooser))

#define GTK_IS_MULTI_FILE_CHOOSER_BUTTON_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((obj), \
                                            GTK_TYPE_MULTI_FILE_CHOOSER_BUTTON))

#define GTK_MULTI_FILE_CHOOSER_BUTTON_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS((obj), \
                                           GTK_TYPE_MULTI_FILE_CHOOSER_BUTTON, GtkMultiFileChooser))

struct _gtkmultifilechooserbutton
{
	GtkButton parent_instance;

	/*< private >*/
	GtkWidget *chooser;
	char *text;
	GSList *filenames;
};

struct _gtkmultifilechooserbuttonclass
{
	GtkButtonClass parent_class;
	int dummy;
};

typedef struct _gtkmultifilechooserbutton GtkMultiFileChooserButton;
typedef struct _gtkmultifilechooserbuttonclass GtkMultiFileChooserButtonClass;

extern GType gtk_multi_file_chooser_button_get_type(void);
extern GtkWidget *gtk_multi_file_chooser_button_new(const char *label);
extern GSList *gtk_multi_file_chooser_button_get_filenames(GtkMultiFileChooserButton *mfc);

#endif	/* GTKMULTIFILECHOOSERBUTTON_H */
