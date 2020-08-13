/*
 * crystfelindexingopts.h
 *
 * A GTK widget to set indexing options
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

#ifndef CRYSTFELINDEXINGOPTS_H
#define CRYSTFELINDEXINGOPTS_H

#include <gtk/gtk.h>
#include <glib-object.h>

#define CRYSTFEL_TYPE_INDEXING_OPTS (crystfel_indexing_opts_get_type())

#define CRYSTFEL_INDEXING_OPTS(obj) (G_TYPE_CHECK_INSTANCE_CAST((obj), \
                                 CRYSTFEL_TYPE_INDEXING_OPTS, CrystFELIndexingOpts))

#define CRYSTFEL_IS_INDEXING_OPTS(obj) (G_TYPE_CHECK_INSTANCE_TYPE((obj), \
                                    CRYSTFEL_TYPE_INDEXING_OPTS))

#define CRYSTFEL_INDEXING_OPTS_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST((obj), \
                                         CRYSTFEL_TYPE_INDEXING_OPTS, CrystFELIndexingOpts))

#define CRYSTFEL_IS_INDEXING_OPTS_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((obj), \
                                            CRYSTFEL_TYPE_INDEXING_OPTS))

#define CRYSTFEL_INDEXING_OPTS_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS((obj), \
                                           CRYSTFEL_TYPE_INDEXING_OPTS, CrystFELIndexingOpts))

struct _crystfelindexingopts
{
	GtkNotebook       parent_instance;

	/*< private >*/
	int dummy;
};

struct _crystfelindexingoptsclass
{
	GtkNotebookClass parent_class;
	int dummy;
};

typedef struct _crystfelindexingopts CrystFELIndexingOpts;
typedef struct _crystfelindexingoptsclass CrystFELIndexingOptsClass;

extern GType crystfel_indexing_opts_get_type(void);
extern GtkWidget *crystfel_indexing_opts_new(void);

#endif	/* CRYSTFELINDEXINGOPTS_H */
