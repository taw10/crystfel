/*
 * gui_backend_local.c
 *
 * GUI backend for running jobs on the local machine
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

#include "crystfel_gui.h"


volatile sig_atomic_t at_zombies = 0;


struct local_backend_priv
{
	pid_t indexamajig_pid;
};


static void check_zombies(struct local_backend_priv *priv)
{
	pid_t r;
	int status;

	if ( !at_zombies ) return;

	r = waitpid(priv->indexamajig_pid, &status, WNOHANG);
	if ( r == -1 ) {
		ERROR("waitpid() failed\n");
	} else {
		STATUS("indexamajig process exited "
		       "with status %i\n", status);
	}
	priv->indexamajig_pid = 0;
	at_zombies = 0;
}


static gboolean index_readable(GIOChannel *source, GIOCondition cond,
                               void *vp)
{
	GIOStatus r;
	GError *err = NULL;
	struct crystfelproject *proj = vp;
	struct local_backend_priv *priv = proj->backend_private;
	gchar *line;
	double frac = 0.1;

	check_zombies(priv);

	r = g_io_channel_read_line(source, &line, NULL, NULL, &err);
	if ( r == G_IO_STATUS_EOF ) {
		STATUS("End of output.\n");
		return FALSE;
	}
	if ( r != G_IO_STATUS_NORMAL ) {
		STATUS("Read error?\n");
		return FALSE;
	}
	chomp(line);
	STATUS("Got line '%s'\n", line);

	gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(proj->progressbar),
	                              frac);

	g_free(line);

	return TRUE;
}


static int write_file_list(struct crystfelproject *proj)
{
	FILE *fh;
	int i;

	fh = fopen("files.lst", "w");
	if ( fh == NULL ) return 1;

	for ( i=0; i<proj->n_frames; i++ ) {
		fprintf(fh, "%s", proj->filenames[i]);
		if ( proj->events[i] != NULL ) {
			fprintf(fh, " %s\n", proj->events[i]);
		} else {
			fprintf(fh, "\n");
		}
	}

	fclose(fh);

	return 0;
}


static void add_arg(char **args, int pos, const char *label,
                    float val)
{
	char *str;

	str = malloc(64);
	if ( str == NULL ) return;

	snprintf(str, 63, "--%s=%f", label, val);
	args[pos] = str;
}


static int run_unitcell(struct crystfelproject *proj,
                        const char *algo)
{
	pid_t pid;
	int pty;
	GIOChannel *ioch;
	char *args[64];
	char index_str[64];
	char peaks_str[64];
	int n_args;
	struct local_backend_priv *priv = proj->backend_private;

	if ( at_zombies ) {
		STATUS("Try again in a moment!\n");
		return 1;
	}

	if ( write_file_list(proj) ) {
		STATUS("Failed to write list\n");
		return 1;
	}

	strcpy(index_str, "--indexing=");
	strncat(index_str, algo, 50);

	strcpy(peaks_str, "--peaks=");
	strncat(peaks_str,
	        str_peaksearch(proj->peak_search_params.method), 50);

	args[0] = "indexamajig";
	args[1] = "-i";
	args[2] = "files.lst";
	args[3] = "-g";
	args[4] = proj->geom_filename;
	args[5] = "-o";
	args[6] = "test.stream";
	args[7] = index_str;
	args[8] = "--no-check-cell";
	args[9] = "-j";
	args[10] = "1";
	args[11] = "--integration=none";
	args[12] = peaks_str;
	n_args = 13;

	if ( proj->peak_search_params.method == PEAK_ZAEF ) {
		add_arg(args, n_args++, "threshold",
		        proj->peak_search_params.threshold);
		add_arg(args, n_args++, "min-squared-gradient",
		        proj->peak_search_params.min_sq_gradient);
		add_arg(args, n_args++, "min-snr",
		        proj->peak_search_params.min_snr);
	} else if ( proj->peak_search_params.method == PEAK_PEAKFINDER8 ) {
		add_arg(args, n_args++, "threshold",
		        proj->peak_search_params.threshold);
		add_arg(args, n_args++, "min-snr",
		        proj->peak_search_params.min_snr);
		add_arg(args, n_args++, "min-pix-count",
		        proj->peak_search_params.min_pix_count);
		add_arg(args, n_args++, "max-pix-count",
		        proj->peak_search_params.max_pix_count);
		add_arg(args, n_args++, "local-bg-radius",
		        proj->peak_search_params.local_bg_radius);
		add_arg(args, n_args++, "min-res",
		        proj->peak_search_params.min_res);
		add_arg(args, n_args++, "max-res",
		        proj->peak_search_params.max_res);
	}

	args[n_args] = NULL;

	int i;
	for ( i=0; i<n_args; i++ ) {
		STATUS("%s ", args[i]);
	}
	STATUS("\n");

	pid = forkpty(&pty, NULL, NULL, NULL);
	if ( pid == -1 ) return 1;

	if ( pid == 0 ) {

		/* Child process */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		execvp("indexamajig", args);
		_exit(1);

	}

	priv->indexamajig_pid = pid;
	ioch = g_io_channel_unix_new(pty);
	g_io_add_watch(ioch, G_IO_IN | G_IO_ERR | G_IO_HUP,
	               index_readable, proj);

	/* FIXME: waitpid() on SIGCHLD */

	return 0;
}


static void cancel(struct crystfelproject *proj)
{
	struct local_backend_priv *priv = proj->backend_private;

	if ( priv->indexamajig_pid == 0 ) {
		ERROR("Indexamajig not running!\n");
		return;
	}

	ERROR("Cancelling!!!\n");
}


static void shutdown_backend(struct crystfelproject *proj)
{
	struct local_backend_priv *priv = proj->backend_private;
	free(priv);
}


static void sigchld_handler(int sig, siginfo_t *si, void *uc_v)
{
	at_zombies = 1;
}


static void init_backend(struct crystfelproject *proj)
{
	struct local_backend_priv *priv;
	struct sigaction sa;
	int r;

	priv = malloc(sizeof(struct local_backend_priv));
	if ( priv == NULL ) {
		ERROR("Failed to initialise backend\n");
		return;
	}

	priv->indexamajig_pid = 0;

	/* Signal handler */
	sa.sa_flags = SA_SIGINFO | SA_NOCLDSTOP | SA_RESTART;
	sigemptyset(&sa.sa_mask);
	sa.sa_sigaction = sigchld_handler;
	r = sigaction(SIGCHLD, &sa, NULL);
	if ( r == -1 ) {
	        ERROR("Failed to set signal handler!\n");
	        return;
	}

	/* Callback to check on signals */
	g_timeout_add_seconds(1, G_SOURCE_FUNC(check_zombies), priv);

	proj->backend_private = priv;
	STATUS("Local backend initialised.\n");
}


struct crystfel_backend _backend_local =
	{
	 .init = init_backend,
	 .shutdown = shutdown_backend,
	 .run_unitcell = run_unitcell,
	 .cancel = cancel,
	};

struct crystfel_backend *backend_local = &_backend_local;
