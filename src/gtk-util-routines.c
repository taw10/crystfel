/*
 * gtk-util-routines.h
 *
 * GTK utilities
 *
 * Copyright © 2020-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020-2021 Thomas White <taw@physics.org>
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

#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <gtk/gtk.h>


char *get_all_text(GtkTextView *view)
{
	GtkTextBuffer *buf;
	GtkTextIter start, end;

	buf = gtk_text_view_get_buffer(view);

	gtk_text_buffer_get_start_iter(buf, &start);
	gtk_text_buffer_get_end_iter(buf, &end);

	return gtk_text_buffer_get_text(buf, &start, &end, FALSE);
}


float get_float(GtkWidget *entry)
{
	const gchar *text;
	char *rval;
	float val;
	text = gtk_entry_get_text(GTK_ENTRY(entry));
	errno = 0;
	val = strtof(text, &rval);
	if ( *rval != '\0' ) return NAN;
	return val;
}


unsigned int get_uint(GtkWidget *entry)
{
	const gchar *text;
	char *rval;
	unsigned long int val;
	text = gtk_entry_get_text(GTK_ENTRY(entry));
	errno = 0;
	val = strtoul(text, &rval, 10);
	if ( *rval != '\0' ) {
		printf("Invalid integer '%s'\n", text);
		return 0;
	}
	return val;
}


int get_bool(GtkWidget *widget)
{
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget));
}


int i_maybe_disable(GtkWidget *toggle, GtkWidget *widget)
{
	gtk_widget_set_sensitive(GTK_WIDGET(widget),
	                         gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggle)));
	return FALSE;
}


int i_maybe_disable_and_deselect(GtkWidget *toggle, GtkWidget *widget)
{
	int active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggle));
	gtk_widget_set_sensitive(GTK_WIDGET(widget), active);
	if ( !active ) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), FALSE);
	}
	return FALSE;
}


static int inv_maybe_disable(GtkWidget *toggle, GtkWidget *victim)
{
	if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggle)) ) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(victim), FALSE);
	}
	return FALSE;
}


void deselect_when_active(GtkWidget *toggle, GtkWidget *victim)
{
	g_signal_connect(G_OBJECT(toggle), "toggled",
	                 G_CALLBACK(inv_maybe_disable),
	                 victim);
	inv_maybe_disable(toggle, victim);
}


void set_active(GtkWidget *tb, int active)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tb), active);
}


void set_combo_id(GtkWidget *cb, const char *id)
{
	gtk_combo_box_set_active_id(GTK_COMBO_BOX(cb), id);
}


void redraw_widget(GtkWidget *wid)
{
	gint w, h;
	w = gtk_widget_get_allocated_width(GTK_WIDGET(wid));
	h = gtk_widget_get_allocated_height(GTK_WIDGET(wid));
	gtk_widget_queue_draw_area(GTK_WIDGET(wid), 0, 0, w, h);
}


const char *get_text_or_null(GtkEntry *entry)
{
	const char *text = gtk_entry_get_text(entry);
	if ( text[0] == '\0' ) return NULL;
	return text;
}
