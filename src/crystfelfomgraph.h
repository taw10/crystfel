/*
 * crystfelfomgraph.h
 *
 * Figure of merit graph plot widget
 *
 * Copyright © 2020-2022 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2020-2022 Thomas White <taw@physics.org>
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

#ifndef CRYSTFELFOMGRAPH_H
#define CRYSTFELFOMGRAPH_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fom.h>

#define CRYSTFEL_TYPE_FOM_GRAPH (crystfel_fom_graph_get_type())

#define CRYSTFEL_FOM_GRAPH(obj) (G_TYPE_CHECK_INSTANCE_CAST((obj), \
                                 CRYSTFEL_TYPE_FOM_GRAPH, CrystFELFoMGraph))

#define CRYSTFEL_IS_FOM_GRAPH(obj) (G_TYPE_CHECK_INSTANCE_TYPE((obj), \
                                    CRYSTFEL_TYPE_FOM_GRAPH))

#define CRYSTFEL_FOM_GRAPH_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST((obj), \
                                         CRYSTFEL_TYPE_FOM_GRAPH, CrystFELFoMGraph))

#define CRYSTFEL_IS_FOM_GRAPH_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((obj), \
                                            CRYSTFEL_TYPE_FOM_GRAPH))

#define CRYSTFEL_FOM_GRAPH_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS((obj), \
                                           CRYSTFEL_TYPE_FOM_GRAPH, CrystFELFoMGraph))


#define COLSCALE_N_BINS (256)
#define COLSCALE_SAMPLE_SIZE (4096)

struct _crystfelfomgraph
{
	GtkDrawingArea       parent_instance;
	double               visible_width;
	double               visible_height;
};

struct _crystfelfomgraphclass
{
	GtkDrawingAreaClass parent_class;
};

typedef struct _crystfelfomgraph CrystFELFoMGraph;
typedef struct _crystfelfomgraphclass CrystFELFoMGraphClass;

extern GType crystfel_fom_graph_get_type(void);
extern GtkWidget *crystfel_fom_graph_new(void);

extern void crystfel_fom_graph_set_data(CrystFELFoMGraph *fg,
                                        double *shell_centers, int n_shells,
                                        enum fom_type *fom_types,
                                        double **fom_values, int n_foms);

#endif	/* CRYSTFELFOMGRAPH_H */
