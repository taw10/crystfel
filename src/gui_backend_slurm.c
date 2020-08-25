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


struct slurm_job
{
	double frac_complete;
	/* FIXME: List of SLURM job numbers to track */
};


static void cancel_task(void *job_priv)
{
	//struct slurm_job *job = job_priv;
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


static void block_size_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_indexing_opts *opts = data;
	convert_int(gtk_entry_get_text(entry), &opts->block_size);
}


static void partition_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_indexing_opts *opts = data;
	opts->partition = strdup(gtk_entry_get_text(entry));
}


static void email_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_indexing_opts *opts = data;
	opts->email_address = strdup(gtk_entry_get_text(entry));
}


static GtkWidget *make_indexing_parameters_widget(void *opts_priv)
{
	struct slurm_indexing_opts *opts = opts_priv;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *entry;
	GtkWidget *label;
	char tmp[64];

	vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Submit job to partition:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	entry = gtk_entry_new();
	if ( opts->partition != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(entry), opts->partition);
	}
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry), "maxwell");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(partition_activate_sig),
	                 opts);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Split job into blocks of");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	snprintf(tmp, 63, "%i", opts->block_size);
	entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(block_size_activate_sig),
	                 opts);
	label = gtk_label_new("frames");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Send notifications to:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	entry = gtk_entry_new();
	if ( opts->email_address != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(entry), opts->email_address);
	}
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry),
	                               "myself@example.org");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(email_activate_sig),
	                 opts);

	return vbox;
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

	if ( opts->partition != NULL) {
		fprintf(fh, "indexing.slurm.partition %s\n",
		        opts->partition);
	}

	if ( opts->email_address != NULL ) {
		fprintf(fh, "indexing.slurm.email_address %s\n",
		        opts->email_address);
	}
}


static void read_indexing_opt(void *opts_priv,
                              const char *key,
                              const char *val)
{
	struct slurm_indexing_opts *opts = opts_priv;

	if ( strcmp(key, "indexing.slurm.block_size") == 0 ) {
		if ( convert_int(val, &opts->block_size) ) {
			ERROR("Invalid block size: %s\n", val);
		}
	}

	if ( strcmp(key, "indexing.slurm.email_address") == 0 ) {
		opts->email_address = strdup(val);
	}

	if ( strcmp(key, "indexing.slurm.partition") == 0 ) {
		opts->partition = strdup(val);
	}
}


int make_slurm_backend(struct crystfel_backend *be)
{
	be->name = "slurm";
	be->friendly_name = "SLURM";

	be->make_indexing_parameters_widget = make_indexing_parameters_widget;
	be->run_indexing = run_indexing;
	be->write_indexing_opts = write_indexing_opts;
	be->read_indexing_opt = read_indexing_opt;
	be->cancel_task = cancel_task;

	be->indexing_opts_priv = make_default_slurm_opts();
	if ( be->indexing_opts_priv == NULL ) return 1;

	return 0;
};
