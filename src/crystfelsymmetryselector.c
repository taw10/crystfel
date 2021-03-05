/*
 * crystfelsymmetryselector.c
 *
 * A GTK widget to choose a symmetry class
 *
 * Copyright © 2020-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2020-2021 Thomas White <taw@physics.org>
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

#include "crystfelsymmetryselector.h"


G_DEFINE_TYPE(CrystFELSymmetrySelector,
              crystfel_symmetry_selector,
              GTK_TYPE_BOX)


static void crystfel_symmetry_selector_class_init(CrystFELSymmetrySelectorClass *klass)
{
}


static void crystfel_symmetry_selector_init(CrystFELSymmetrySelector *mo)
{
}


GtkWidget *crystfel_symmetry_selector_new()
{
	CrystFELSymmetrySelector *sel;

	sel = g_object_new(CRYSTFEL_TYPE_SYMMETRY_SELECTOR, NULL);

	sel->entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(sel->entry), 10);
	gtk_box_pack_start(GTK_BOX(sel), GTK_WIDGET(sel->entry),
	                   FALSE, FALSE, 0.0);

	return GTK_WIDGET(sel);
}


char *crystfel_symmetry_selector_get_group_symbol(CrystFELSymmetrySelector *sel)
{
	const char *text = gtk_entry_get_text(GTK_ENTRY(sel->entry));
	if ( text == NULL ) return NULL;
	if ( text[0] == '\0' ) return NULL;
	return strdup(text);
}


int crystfel_symmetry_selector_set_group_symbol(CrystFELSymmetrySelector *sel,
                                                const char *pg_symbol)
{
	if ( pg_symbol != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(sel->entry), pg_symbol);
	}
	return 0;
}
