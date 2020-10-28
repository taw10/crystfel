/*
 * crystfelsymmetryselector.h
 *
 * A GTK widget to choose a symmetry class
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

#ifndef CRYSTFELSYMMETRYSELECTOR_H
#define CRYSTFELSYMMETRYSELECTOR_H

#include <gtk/gtk.h>
#include <glib-object.h>

#define CRYSTFEL_TYPE_SYMMETRY_SELECTOR (crystfel_symmetry_selector_get_type())

#define CRYSTFEL_SYMMETRY_SELECTOR(obj) (G_TYPE_CHECK_INSTANCE_CAST((obj), \
                                 CRYSTFEL_TYPE_SYMMETRY_SELECTOR, CrystFELSymmetrySelector))

#define CRYSTFEL_IS_SYMMETRY_SELECTOR(obj) (G_TYPE_CHECK_INSTANCE_TYPE((obj), \
                                    CRYSTFEL_TYPE_SYMMETRY_SELECTOR))

#define CRYSTFEL_SYMMETRY_SELECTOR_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST((obj), \
                                         CRYSTFEL_TYPE_SYMMETRY_SELECTOR, CrystFELSymmetrySelector))

#define CRYSTFEL_IS_SYMMETRY_SELECTOR_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((obj), \
                                            CRYSTFEL_TYPE_SYMMETRY_SELECTOR))

#define CRYSTFEL_SYMMETRY_SELECTOR_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS((obj), \
                                           CRYSTFEL_TYPE_SYMMETRY_SELECTOR, CrystFELSymmetrySelector))

struct _crystfelsymmetryselector
{
	GtkBox parent_instance;

	/*< private >*/
	GtkWidget *entry;
};

struct _crystfelsymmetryselectorclass
{
	GtkBoxClass parent_class;
};

typedef struct _crystfelsymmetryselector CrystFELSymmetrySelector;
typedef struct _crystfelsymmetryselectorclass CrystFELSymmetrySelectorClass;

extern GType crystfel_symmetry_selector_get_type(void);
extern GtkWidget *crystfel_symmetry_selector_new(void);

extern char *crystfel_symmetry_selector_get_group_symbol(CrystFELSymmetrySelector *sel);
extern int crystfel_symmetry_selector_set_group_symbol(CrystFELSymmetrySelector *sel,
                                                       const char *pg_symbol);

#endif	/* CRYSTFELSYMMETRYSELECTOR_H */
