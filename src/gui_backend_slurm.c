/*
 * gui_backend_slurm.c
 *
 * GUI backend for running jobs via SLURM
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

#include <pty.h>
#include <glib.h>
#include <sys/wait.h>
#include <gtk/gtk.h>

#include "crystfel_gui.h"


struct slurm_backend_priv
{
	int dummy;
};


static void init_backend(struct crystfelproject *proj)
{
}


static void shutdown_backend(struct crystfelproject *proj)
{
}


static void cancel(struct crystfelproject *proj)
{
}


static int run_unitcell(struct crystfelproject *proj,
                        const char *algo)
{
	return 0;
}


static GtkWidget *make_parameters(void)
{
	return gtk_label_new("SLURM params");
}


const struct crystfel_backend _backend_slurm =
	{
	 .name = "slurm",
	 .friendly_name = "SLURM",
	 .make_parameters = make_parameters,
	 .init = init_backend,
	 .shutdown = shutdown_backend,
	 .run_unitcell = run_unitcell,
	 .cancel = cancel,
	};

const struct crystfel_backend *backend_slurm = &_backend_slurm;
