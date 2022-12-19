/*
 * gui_fom.c
 *
 * Figures of merit via CrystFEL GUI
 *
 * Copyright © 2020-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020-2021 Thomas White <taw@physics.org>
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

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <gtk/gtk.h>
#include <assert.h>

#include <utils.h>
#include <fom.h>
#include <reflist-utils.h>
#include <cell-utils.h>

#include "gui_project.h"
#include "crystfel_gui.h"
#include "gtk-util-routines.h"
#include "crystfelfomgraph.h"

#define MAX_DATASETS (64)


struct fom_window
{
	struct crystfelproject *proj;
	GtkWidget *min_res;
	GtkWidget *max_res;
	GtkWidget *num_bins;
	GtkWidget *min_snr;
	GtkWidget *min_meas;
	GtkWidget *cell_chooser;
	GtkWidget *graph;
	GtkWidget *input_combo;

	int n_foms;
	GtkWidget *fom_checkboxes[16];
	enum fom_type fom_types[16];
};


static int menu_selected(GtkWidget *w)
{
	return gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(w));
}


static int fom_selected(struct fom_window *f, int i)
{
	return menu_selected(f->fom_checkboxes[i]);
}


static int anomalous_foms_selected(struct fom_window *f)
{
	int i;
	for ( i=0; i<f->n_foms; i++ ) {
		if ( fom_selected(f, i) ) {
			if ( fom_is_anomalous(f->fom_types[i]) ) return 1;
		}
	}
	return 0;
}


static double *make_shell_centers(struct fom_shells *shells)
{
	int i;
	double *vals;

	vals = malloc(shells->nshells*sizeof(double));
	if ( vals == NULL ) return NULL;

	for ( i=0; i<shells->nshells; i++ ) {
		vals[i] = fom_shell_centre(shells, i);
	}

	return vals;
}


static double *make_fom_vals(struct fom_context *fctx,
                             struct fom_shells *shells)
{
	int i;
	double *vals;

	vals = malloc(shells->nshells*sizeof(double));
	if ( vals == NULL ) return NULL;

	for ( i=0; i<shells->nshells; i++ ) {
		vals[i] = fom_shell_value(fctx, i);
	}

	return vals;
}


static int load_dataset(struct gui_merge_result *result,
                        int need_ano, UnitCell *cell,
                        double min_res, double max_res,
                        int min_meas, double min_snr,
                        SymOpList **psym,
                        RefList **pall_refls,
                        RefList **pall_refls_anom,
                        RefList **ppart1,
                        RefList **ppart2,
                        RefList **ppart1_anom,
                        RefList **ppart2_anom)
{
	RefList *raw_refl;
	RefList *raw_part1;
	RefList *raw_part2;
	RefList *all_refls = NULL;
	RefList *all_refls_anom = NULL;
	RefList *part1 = NULL;
	RefList *part2 = NULL;
	RefList *part1_anom = NULL;
	RefList *part2_anom = NULL;
	SymOpList *sym;
	char *sym_str;
	char *sym_str_part1;
	char *sym_str_part2;

	raw_refl = read_reflections_2(result->hkl, &sym_str);
	if ( raw_refl == NULL ) {
		ERROR("Failed to load dataset %s (%s)\n",
		      result->name, result->hkl);
		return 1;
	}

	raw_part1 = read_reflections_2(result->hkl1, &sym_str_part1);
	if ( raw_part1 == NULL ) {
		ERROR("Failed to load part 1 dataset %s (%s)\n",
		      result->name, result->hkl1);
		return 1;
	}

	raw_part2 = read_reflections_2(result->hkl2, &sym_str_part2);
	if ( raw_part2 == NULL ) {
		ERROR("Failed to load part 2 dataset %s (%s)\n",
		      result->name, result->hkl2);
		return 1;
	}

	if ( (sym_str == NULL)
	  || (sym_str_part1 == NULL)
	  || (sym_str_part2 == NULL) )
	{
		ERROR("Reflection list has no point group\n");
		reflist_free(raw_refl);
		reflist_free(raw_part1);
		reflist_free(raw_part2);
		return 1;
	}

	if ( (strcmp(sym_str, sym_str_part1) != 0)
	  || (strcmp(sym_str, sym_str_part2) != 0) )
	{
		ERROR("Datasets do not have the same point group!\n");
		free(sym_str);
		free(sym_str_part1);
		free(sym_str_part2);
		reflist_free(raw_refl);
		reflist_free(raw_part1);
		reflist_free(raw_part2);
		return 1;
	}

	sym = get_pointgroup(sym_str);
	free(sym_str);
	free(sym_str_part1);
	free(sym_str_part2);

	fom_select_reflections(raw_refl, &all_refls,
	                       cell, sym,
	                       1e10/min_res, 1e10/max_res,
	                       min_snr, 0, 0, min_meas);
	if ( all_refls == NULL ) {
		ERROR("Failed to select reflections for dataset '%s'\n",
		      result->name);
		reflist_free(raw_refl);
		reflist_free(raw_part1);
		reflist_free(raw_part2);
		return 1;
	}

	fom_select_reflection_pairs(raw_part1, raw_part2,
	                            &part1, &part2,
	                            cell, sym, 0,
	                            1e10/min_res, 1e10/max_res,
	                            min_snr, 0, 0, min_meas);
	if ( (part1 == NULL) || (part2 == NULL) ) {
		ERROR("Failed to select reflection pairs for dataset '%s'\n",
		      result->name);
		reflist_free(all_refls);
		reflist_free(raw_refl);
		reflist_free(raw_part1);
		reflist_free(raw_part2);
		return 1;
	}

	if ( need_ano ) {

		fom_select_reflections(raw_refl, &all_refls_anom,
		                       cell, sym,
		                       1e10/min_res, 1e10/max_res,
		                       min_snr, 0, 0, min_meas);
		if ( all_refls_anom == NULL ) {
			ERROR("Failed to load dataset '%s'\n",
			      result->name);
			reflist_free(raw_refl);
			return 1;
		}

		fom_select_reflection_pairs(raw_part1, raw_part2,
		                            &part1_anom, &part2_anom,
		                            cell, sym, 1,
		                            1e10/min_res, 1e10/max_res,
		                            min_snr, 0, 0, min_meas);
		if ( (part1_anom == NULL) || (part2_anom == NULL) ) {
			ERROR("Failed to select anomalous reflection pairs "
			      "for dataset '%s'\n", result->name);
			reflist_free(part1);
			reflist_free(part2);
			reflist_free(all_refls);
			reflist_free(raw_refl);
			reflist_free(raw_part1);
			reflist_free(raw_part2);
			return 1;
		}
	}

	reflist_free(raw_refl);
	reflist_free(raw_part1);
	reflist_free(raw_part2);

	*pall_refls = all_refls;
	*pall_refls_anom = all_refls_anom;
	*ppart1 = part1;
	*ppart2 = part2;
	*ppart1_anom = part1_anom;
	*ppart2_anom = part2_anom;
	*psym = sym;
	return 0;
}


static struct fom_context *dispatch_fom(RefList *all_refls,
                                        RefList *all_refls_anom,
                                        RefList *part1,
                                        RefList *part2,
                                        RefList *part1_anom,
                                        RefList *part2_anom,
                                        UnitCell *cell,
                                        struct fom_shells *shells,
                                        const SymOpList *sym,
                                        enum fom_type fom)
{
	if ( fom_is_anomalous(fom) ) {
		if ( fom_is_comparison(fom) ) {
			if ( part1_anom == NULL ) return NULL;
			if ( part2_anom == NULL ) return NULL;
			return fom_calculate(part1_anom, part2_anom,
			                     cell, shells, fom, 1, sym);
		} else {
			if ( all_refls_anom == NULL ) return NULL;
			return fom_calculate(all_refls_anom, NULL,
			                     cell, shells, fom, 1, sym);
		}
	} else {
		if ( fom_is_comparison(fom) ) {
			if ( part1 == NULL ) return NULL;
			if ( part2 == NULL ) return NULL;
			return fom_calculate(part1, part2,
			                     cell, shells, fom, 1, sym);
		} else {
			if ( all_refls == NULL ) return NULL;
			return fom_calculate(all_refls, NULL,
			                     cell, shells, fom, 1, sym);
		}
	}
}


static void fom_response_sig(GtkWidget *dialog, gint resp,
                             struct fom_window *f)
{
	if ( resp != GTK_RESPONSE_APPLY ) {
		gtk_widget_destroy(dialog);
		return;
	}
}


static void update_fom(GtkWidget *widget, struct fom_window *f)
{
	int fom;
	int need_ano;
	UnitCell *cell;
	struct fom_shells *shells;

	f->proj->fom_res_min = get_float(f->min_res);
	f->proj->fom_res_max = get_float(f->max_res);
	f->proj->fom_min_snr = get_float(f->min_snr);
	f->proj->fom_min_meas = get_uint(f->min_meas);
	f->proj->fom_nbins = get_uint(f->num_bins);
	if ( isnan(f->proj->fom_res_min)
	  || isnan(f->proj->fom_res_max)
	  || isnan(f->proj->fom_min_snr)
	  || isinf(f->proj->fom_res_min)
	  || isinf(f->proj->fom_res_max) )
	{
		ERROR("Invalid parameters\n");
		return;
	}

	/* "Minimum resolution" should be the bigger number */
	if ( f->proj->fom_res_min < f->proj->fom_res_max ) {
		double tmp = f->proj->fom_res_min;
		f->proj->fom_res_min = f->proj->fom_res_max;
		f->proj->fom_res_max = tmp;
	}

	shells = fom_make_resolution_shells(1e10/f->proj->fom_res_min,
	                                    1e10/f->proj->fom_res_max,
	                                    f->proj->fom_nbins);
	if ( shells == NULL ) {
		ERROR("Failed to make resolution shells\n");
		return;
	}

	f->proj->fom_cell_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(f->cell_chooser));
	cell = load_cell_from_file(f->proj->fom_cell_filename);
	if ( cell == NULL ) {
		ERROR("Invalid cell file '%s'\n", f->proj->fom_cell_filename);
		return;
	}

	need_ano = anomalous_foms_selected(f);

	const char *name;
	struct gui_merge_result *result;
	SymOpList *sym = NULL;
	RefList *all_refls = NULL;
	RefList *all_refls_anom = NULL;
	RefList *part1 = NULL;
	RefList *part2 = NULL;
	RefList *part1_anom = NULL;
	RefList *part2_anom = NULL;

	/* Load dataset */
	name = gtk_combo_box_get_active_id(GTK_COMBO_BOX(f->input_combo));
	result = find_merge_result_by_name(f->proj, name);

	if ( load_dataset(result, need_ano, cell,
	                  f->proj->fom_res_min,
	                  f->proj->fom_res_max,
	                  f->proj->fom_min_meas,
	                  f->proj->fom_min_snr,
	                  &sym, &all_refls, &all_refls_anom,
	                  &part1, &part2, &part1_anom, &part2_anom) )
	{
		return;
	}

	double *shell_centers = malloc(shells->nshells*sizeof(double));
	double **fom_values = malloc(f->n_foms*sizeof(double *));
	enum fom_type *fom_types = malloc(f->n_foms*sizeof(enum fom_type));

	int fomi = 0;
	for ( fom=0; fom<f->n_foms; fom++ ) {

		struct fom_context *fctx;

		if ( !fom_selected(f, fom) ) continue;

		fctx = dispatch_fom(all_refls, all_refls_anom,
		                    part1, part2,
		                    part1_anom, part2_anom,
		                    cell, shells, sym,
		                    f->fom_types[fom]);
		if ( fctx == NULL ) {
			ERROR("Failed to calculate FoM %i for dataset %s\n",
			      f->fom_types[fom], name);
			continue;
		}

		fom_types[fomi] = f->fom_types[fom];
		fom_values[fomi] = make_fom_vals(fctx, shells);
		fomi++;

	}

	shell_centers = make_shell_centers(shells);
	crystfel_fom_graph_set_data(CRYSTFEL_FOM_GRAPH(f->graph),
	                            shell_centers, shells->nshells,
	                            fom_types, fom_values, fomi);

	reflist_free(all_refls);
	reflist_free(all_refls_anom);
	free_symoplist(sym);

}


static GtkWidget *add_item(GtkWidget *menu,
                           const char *text,
                           const char *markup,
                           struct fom_window *f)
{
	GtkWidget *label;
	GtkWidget *item;

	item = gtk_check_menu_item_new();

	label = gtk_label_new(text);
	if ( markup != NULL ) {
		gtk_label_set_markup(GTK_LABEL(label), markup);
	}
	gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);

	gtk_container_add(GTK_CONTAINER(item), label);
	gtk_widget_show_all(item);

	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	g_signal_connect(G_OBJECT(item), "toggled", G_CALLBACK(update_fom), f);

	return item;
}


static void add_separator(GtkWidget *menu)
{
	GtkWidget *item;
	item = gtk_separator_menu_item_new();
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	gtk_widget_show_all(item);
}


static GtkWidget *make_fom_menu(struct fom_window *fom)
{
	GtkWidget *menu;
	GtkWidget *item;

	menu = gtk_menu_new();

	item = gtk_tearoff_menu_item_new();
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	gtk_widget_show_all(item);

	/* Order of FoMs must match list below */
	fom->fom_checkboxes[0] = add_item(menu, "I/σ(I)", NULL, fom);
	fom->fom_checkboxes[1] = add_item(menu, "Completeness", NULL, fom);
	fom->fom_checkboxes[2] = add_item(menu, "Redundancy", NULL, fom);
	add_separator(menu);
	fom->fom_checkboxes[3] = add_item(menu, "Rsplit", "R<sub>split</sub>", fom);
	fom->fom_checkboxes[4] = add_item(menu, "CC", "CC<sub>½</sub>", fom);
	fom->fom_checkboxes[5] = add_item(menu, "CC*", "CC<sup>*</sup>", fom);
	add_separator(menu);
	fom->fom_checkboxes[6] = add_item(menu, "CCano", "CC<sub>ano</sub>", fom);
	fom->fom_checkboxes[7] = add_item(menu, "Rano", "R<sub>ano</sub>", fom);
	fom->fom_checkboxes[8] = add_item(menu, "Rano ÷ Rsplit",
	                                  "R<sub>ano</sub> ÷ R<sub>split</sub>", fom);
	fom->fom_checkboxes[9] = add_item(menu, "RMS anomalous correlation ratio", NULL, fom);
	add_separator(menu);
	fom->fom_checkboxes[10] = add_item(menu, "Fraction of differences within 1σ(I)", NULL, fom);
	fom->fom_checkboxes[11] = add_item(menu, "Fraction of differences within 2σ(I)", NULL, fom);

	/* Order must match the list above */
	fom->fom_types[0] = FOM_SNR;
	fom->fom_types[1] = FOM_COMPLETENESS;
	fom->fom_types[2] = FOM_REDUNDANCY;

	fom->fom_types[3] = FOM_RSPLIT;
	fom->fom_types[4] = FOM_CC;
	fom->fom_types[5] = FOM_CCSTAR;

	fom->fom_types[6] = FOM_CCANO;
	fom->fom_types[7] = FOM_RANO;
	fom->fom_types[8] = FOM_RANORSPLIT;
	fom->fom_types[9] = FOM_CRDANO;

	fom->fom_types[10] = FOM_D1SIG;
	fom->fom_types[11] = FOM_D2SIG;

	fom->n_foms = 12;

	return menu;
}


static int result_res_range(struct gui_merge_result *result,
                            UnitCell *cell,
                            double *lowres,
                            double *highres)
{
	RefList *rl;
	double rmin, rmax;

	rl = read_reflections(result->hkl);
	if ( rl == NULL ) return 1;

	resolution_limits(rl, cell, &rmin, &rmax);
	reflist_free(rl);

	*lowres = 1e10/rmin;
	*highres = 1e10/rmax;

	return 0;
}


static void res_range_to_data_sig(GtkButton *buton,
                                  struct fom_window *f)
{
	const char *name;
	struct gui_merge_result *result;
	gchar *cell_filename;
	UnitCell *cell;
	char tmp[64];
	double lowres;  /* Angstroms */
	double highres;

	cell_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(f->cell_chooser));
	if ( cell_filename == NULL ) {
		ERROR("You must choose the unit cell first.\n");
		return;
	}
	cell = load_cell_from_file(cell_filename);
	if ( cell == NULL ) {
		ERROR("Invalid cell file '%s'\n", cell_filename);
		g_free(cell_filename);
		return;
	}
	g_free(cell_filename);

	name = gtk_combo_box_get_active_id(GTK_COMBO_BOX(f->input_combo));
	result = find_merge_result_by_name(f->proj, name);

	if ( result_res_range(result, cell, &lowres, &highres) ) return;

	cell_free(cell);

	f->proj->fom_res_min = lowres;
	f->proj->fom_res_max = highres;
	snprintf(tmp, 64, "%.2f", f->proj->fom_res_min);
	gtk_entry_set_text(GTK_ENTRY(f->min_res), tmp);
	snprintf(tmp, 64, "%.2f", f->proj->fom_res_max);
	gtk_entry_set_text(GTK_ENTRY(f->max_res), tmp);

	update_fom(NULL, f);
}


static void cell_file_clear_sig(GtkButton *button,
                                struct fom_window *f)
{
	gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(f->cell_chooser),
	                              "(none)");
	update_fom(NULL, f);
}


gint fom_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *button;
	char tmp[64];
	struct fom_window *f;
	int i;

	f = malloc(sizeof(struct fom_window));
	if ( f == NULL ) return 0;

	f->proj = proj;
	f->n_foms = 0;

	dialog = gtk_dialog_new_with_buttons("Calculate figures of merit",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Close", GTK_RESPONSE_CLOSE,
	                                     NULL);

	g_signal_connect(G_OBJECT(dialog), "response",
	                 G_CALLBACK(fom_response_sig),
	                 f);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_box_pack_start(GTK_BOX(content_area), GTK_WIDGET(vbox), TRUE, TRUE, 0.0);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Results to show:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	f->input_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(f->input_combo),
	                   FALSE, FALSE, 4.0);
	for ( i=0; i<proj->n_merge_results; i++ ) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(f->input_combo),
		                          proj->merge_results[i].name,
		                          proj->merge_results[i].name);
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(f->input_combo),
	                         proj->n_merge_results-1);
	g_signal_connect(G_OBJECT(f->input_combo), "changed",
	                 G_CALLBACK(update_fom), f);

	label = gtk_label_new("Figures of merit to show:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	button = gtk_menu_button_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(button),
	                   FALSE, FALSE, 4.0);
	gtk_menu_button_set_popup(GTK_MENU_BUTTON(button),
	                          make_fom_menu(f));

	/* Unit cell */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Unit cell file:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	f->cell_chooser = gtk_file_chooser_button_new("Unit cell file",
	                                              GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_file_chooser_set_local_only(GTK_FILE_CHOOSER(f->cell_chooser),
	                                TRUE);
	if ( proj->fom_cell_filename != NULL ) {
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(f->cell_chooser),
		                              proj->fom_cell_filename);
	}
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(f->cell_chooser),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(f->cell_chooser), "file-set",
	                 G_CALLBACK(update_fom), f);
	button = gtk_button_new_from_icon_name("edit-clear",
	                                       GTK_ICON_SIZE_BUTTON);
	g_signal_connect(G_OBJECT(button), "clicked",
	                 G_CALLBACK(cell_file_clear_sig), f);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(button),
	                   FALSE, FALSE, 4.0);

	/* Resolution range */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Resolution range:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	f->min_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(f->min_res), 6);
	snprintf(tmp, 64, "%.2f", proj->fom_res_min);
	gtk_entry_set_text(GTK_ENTRY(f->min_res), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(f->min_res),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(f->min_res), "activate",
	                 G_CALLBACK(update_fom), f);
	label = gtk_label_new("to");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	f->max_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(f->max_res), 4);
	snprintf(tmp, 64, "%.2f", proj->fom_res_max);
	gtk_entry_set_text(GTK_ENTRY(f->max_res), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(f->max_res),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(f->max_res), "activate",
	                 G_CALLBACK(update_fom), f);
	label = gtk_label_new("Å");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	button = gtk_button_new_with_label("Reset to entire data");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(button),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(button), "clicked",
	                 G_CALLBACK(res_range_to_data_sig), f);

	/* Number of resolution bins */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Number of resolution bins:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	f->num_bins = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(f->num_bins), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(f->num_bins),
	                   FALSE, FALSE, 4.0);
	snprintf(tmp, 64, "%i", proj->fom_nbins);
	gtk_entry_set_text(GTK_ENTRY(f->num_bins), tmp);
	g_signal_connect(G_OBJECT(f->num_bins), "activate",
	                 G_CALLBACK(update_fom), f);

	/* Minimum I/sigI */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Minimum I/sigI:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	f->min_snr = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(f->min_snr), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(f->min_snr),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(f->min_snr), "activate",
	                 G_CALLBACK(update_fom), f);

	/* Min measurements per reflection */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	snprintf(tmp, 64, "%.2f", proj->fom_min_snr);
	gtk_entry_set_text(GTK_ENTRY(f->min_snr), tmp);
	label = gtk_label_new("Minimum number of measurements per reflection:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	f->min_meas = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(f->min_meas), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(f->min_meas),
	                   FALSE, FALSE, 4.0);
	snprintf(tmp, 64, "%i", proj->fom_min_meas);
	gtk_entry_set_text(GTK_ENTRY(f->min_meas), tmp);
	g_signal_connect(G_OBJECT(f->min_meas), "activate",
	                 G_CALLBACK(update_fom), f);

	f->graph = crystfel_fom_graph_new();
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(f->graph),
	                   TRUE, TRUE, 4.0);

	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(f->fom_checkboxes[3]), TRUE);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(f->fom_checkboxes[5]), TRUE);

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_CLOSE);
	gtk_widget_show_all(dialog);

	return FALSE;
}
