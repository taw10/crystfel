/*
 * gtk-util-routines.h
 *
 * GTK utilities
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

#ifndef GTK_UTIL_ROUTINES_H
#define GTK_UTIL_ROUTINES_H

#include <gtk/gtk.h>

extern char *get_all_text(GtkTextView *view);
extern int get_bool(GtkWidget *widget);
extern unsigned int get_uint(GtkWidget *entry);
extern float get_float(GtkWidget *entry);
extern int i_maybe_disable(GtkWidget *toggle, GtkWidget *widget);

#endif
