/*
 * gui_project.h
 *
 * GUI project persistence
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

#ifndef GUI_PROJECT_H
#define GUI_PROJECT_H

#include <gtk/gtk.h>

#include <peaks.h>

enum match_type_id
{
	 MATCH_EVERYTHING,
	 MATCH_CHEETAH_LCLS_H5,
	 MATCH_CHEETAH_CXI,
	 MATCH_CBF,
	 MATCH_CBFGZ,
};

struct peak_params {
	enum peak_search_method method;
	float threshold;                /* zaef, pf8 */
	float min_sq_gradient;          /* zaef */
	float min_snr;                  /* zaef, pf8 */
	int min_pix_count;              /* pf8 */
	int max_pix_count;              /* pf8 */
	int local_bg_radius;            /* pf8 */
	int min_res;                    /* pf8 */
	int max_res;                    /* pf8 */
	float min_snr_biggest_pix;      /* pf9 */
	float min_snr_peak_pix;         /* pf9 */
	float min_sig;                  /* pf9 */
	float min_peak_over_neighbour;  /* pf9 */
	float pk_inn;
	float pk_mid;
	float pk_out;
	int half_pixel_shift;           /* cxi, hdf5 */
	int revalidate;
};

struct crystfelproject {

	GtkWidget *window;
	GtkUIManager *ui;
	GtkActionGroup *action_group;

	GtkWidget *imageview;
	GtkWidget *icons;      /* Drawing area for task icons */
	GtkWidget *report;     /* Text view at the bottom for messages */
	GtkWidget *main_vbox;
	GtkWidget *image_info;

	int unsaved;

	int cur_frame;

	char *geom_filename;
	char *stream_filename;
	char *data_top_folder;   /* For convenience only.  Filenames in
	                          * 'filenames' list should be complete */
	enum match_type_id data_search_pattern;

	int n_frames;
	int max_frames;
	char **filenames;
	char **events;

	int show_peaks;
	struct peak_params peak_search_params;

	GtkWidget *type_combo;
	GtkWidget *peak_vbox;     /* Box for peak search parameter widgets */
	GtkWidget *peak_params;   /* Peak search parameter widgets */
	struct peak_params original_params;

	GtkWidget *unitcell_combo;

	GtkWidget *info_bar;
	void (*infobar_callback)(struct crystfelproject *proj);
	GtkWidget *progressbar;

	struct crystfel_backend *backend;
	void *backend_private;
};

extern enum match_type_id decode_matchtype(const char *type_id);

extern int match_filename(const char *fn, enum match_type_id mt);

extern int load_project(struct crystfelproject *proj);

extern int save_project(struct crystfelproject *proj);

extern void add_file_to_project(struct crystfelproject *proj,
                                const char *filename,
                                const char *event);

extern void clear_project_files(struct crystfelproject *proj);

#endif
