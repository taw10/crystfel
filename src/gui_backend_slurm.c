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

#include <utils.h>

#include "gui_project.h"


struct slurm_indexing_opts
{
	char *partition;
	int block_size;
	char *email_address;
};


struct slurm_indexing_job
{
	double frac_complete;
	/* FIXME: List of SLURM job numbers to track */
};


static void cancel_indexing(void *job_priv)
{
	//struct slurm_indexing_job *job = job_priv;
}


static void *run_indexing(char **filenames,
                          char **events,
                          int n_frames,
                          char *geom_filename,
                          struct peak_params *peak_search_params,
                          struct index_params *indexing_params,
                          void *opts_priv)
{
	//struct slurm_indexing_opts *opts = opts_priv;
	return NULL;
}


static GtkWidget *make_indexing_parameters_widget(void *opts_priv)
{
	//struct slurm_indexing_opts *opts = opts_priv;

	return gtk_label_new("SLURM params");
}


static struct slurm_indexing_opts *make_default_slurm_opts()
{
	struct slurm_indexing_opts *opts = malloc(sizeof(struct slurm_indexing_opts));
	if ( opts == NULL ) return NULL;

	opts->partition = NULL;
	opts->block_size = 1000;
	opts->email_address = NULL;

	return opts;
}


static void write_indexing_opts(void *opts_priv, FILE *fh)
{
	struct slurm_indexing_opts *opts = opts_priv;

	fprintf(fh, "indexing.slurm.block_size %i\n",
	        opts->block_size);
	fprintf(fh, "indexing.slurm.partition %s\n",
	        opts->partition);
	fprintf(fh, "indexing.slurm.email_address %s\n",
	        opts->email_address);
}


static void read_indexing_opt(void *opts_priv,
                              const char *key,
                              const char *val)
{
	//struct slurm_indexing_opts *opts = opts_priv;

	STATUS("SLURM got %s = '%s'\n", key, val);
	/* FIXME: Parse and set */
}


int make_slurm_backend(struct crystfel_backend *be)
{
	be->name = "slurm";
	be->friendly_name = "SLURM";

	be->make_indexing_parameters_widget = make_indexing_parameters_widget;
	be->run_indexing = run_indexing;
	be->cancel_indexing = cancel_indexing;
	be->write_indexing_opts = write_indexing_opts;
	be->read_indexing_opt = read_indexing_opt;

	be->indexing_opts_priv = make_default_slurm_opts();
	if ( be->indexing_opts_priv == NULL ) return 1;

	return 0;
};
