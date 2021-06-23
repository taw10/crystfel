/*
 * crystfelsymmetryselector.c
 *
 * A GTK widget to choose a symmetry class
 *
 * Copyright © 2020-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2020-2021 Thomas White <taw@physics.org>
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <gtk/gtk.h>
#include <glib-object.h>
#include <errno.h>

#include <cell.h>
#include <cell-utils.h>

#include "crystfelsymmetryselector.h"
#include "gtk-util-routines.h"


G_DEFINE_TYPE(CrystFELSymmetrySelector,
              crystfel_symmetry_selector,
              GTK_TYPE_BUTTON)


static void crystfel_symmetry_selector_class_init(CrystFELSymmetrySelectorClass *klass)
{
}


static void crystfel_symmetry_selector_init(CrystFELSymmetrySelector *mo)
{
}


static void selector_response_sig(GtkWidget *dialog, gint resp,
                                  CrystFELSymmetrySelector *sel)
{
	gtk_widget_destroy(sel->dialog);
	sel->dialog = NULL;
	free(sel->pointgroups);
	sel->n_pgs = 0;
}


static void add_lattice_types(GtkComboBoxText *combo)
{
	gtk_combo_box_text_append(combo, "triclinic", "Triclinic");
	gtk_combo_box_text_append(combo, "monoclinic", "Monoclinic");
	gtk_combo_box_text_append(combo, "orthorhombic", "Orthorhombic");
	gtk_combo_box_text_append(combo, "tetragonal", "Tetragonal");
	gtk_combo_box_text_append(combo, "rhombohedral", "Rhombohedral");
	gtk_combo_box_text_append(combo, "hexagonal", "Hexagonal");
	gtk_combo_box_text_append(combo, "cubic", "Cubic");
}


#define MAX_PG_WIDGETS 128

static int filter_pgs(GtkFlowBoxChild *child, gpointer data)
{
	int i;
	int idx;
	int centro, sohnke;
	const char *unique_axis;
	const char *s;
	LatticeType lattice_type;
	CrystFELSymmetrySelector *sel = data;

	centro = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(sel->centro_checkbox));
	sohnke = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(sel->sohnke_checkbox));
	unique_axis = gtk_combo_box_get_active_id(GTK_COMBO_BOX(sel->unique_axis_combo));

	s = gtk_combo_box_get_active_id(GTK_COMBO_BOX(sel->lattice_combo));

	lattice_type = lattice_from_str(s);
	if ( !has_unique_axis(lattice_type) ) {
		unique_axis = "*";
	}

	idx = MAX_PG_WIDGETS;
	for ( i=0; i<sel->n_pgs; i++ ) {
		if ( sel->pointgroups[i].w == gtk_bin_get_child(GTK_BIN(child)) ) {
			idx = i;
			break;
		}
	}

	if ( idx == MAX_PG_WIDGETS ) return FALSE;

	if ( (sel->pointgroups[idx].centro == centro)
	    && (sel->pointgroups[idx].sohnke == sohnke)
	    && (sel->pointgroups[idx].lattice_type == lattice_type)
	    && (sel->pointgroups[idx].unique_axis[0] == unique_axis[0]) )
	{
		return TRUE;
	} else {
		return FALSE;
	}
}


static void pg_press_sig(GtkWidget *button, CrystFELSymmetrySelector *sel)
{
	GtkWidget *label;
	const char *text;

	label = gtk_bin_get_child(GTK_BIN(button));
	text = gtk_label_get_text(GTK_LABEL(label));

	label = gtk_bin_get_child(GTK_BIN(sel));
	gtk_label_set_text(GTK_LABEL(label), text);
	sel->have_pg = 1;
	gtk_dialog_response(GTK_DIALOG(sel->dialog), GTK_RESPONSE_CANCEL);
}


static void add_pointgroup(CrystFELSymmetrySelector *sel,
                           const char *symbol,
                           int centro,
                           int sohnke,
                           char *unique_axis,
                           LatticeType latt)
{
	struct pointgroup_widget *pg;
	GtkWidget *label;

	pg = &sel->pointgroups[sel->n_pgs++];
	assert(sel->n_pgs < MAX_PG_WIDGETS);

	pg->w = gtk_button_new();
	label = gtk_label_new(symbol);
	gtk_container_add(GTK_CONTAINER(pg->w), label);
	gtk_flow_box_insert(GTK_FLOW_BOX(sel->flowbox), pg->w, -1);
	g_signal_connect(G_OBJECT(pg->w), "clicked",
	                 G_CALLBACK(pg_press_sig), sel);
	pg->symbol = symbol;
	pg->centro = centro;
	pg->sohnke = sohnke;
	pg->unique_axis = unique_axis;
	pg->lattice_type = latt;
}


static void update_opts_sig(GtkWidget *ignore, CrystFELSymmetrySelector *sel)
{
	int centro, sohnke, unique_axis;
	const char *s;
	LatticeType lattice_type;

	centro = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(sel->centro_checkbox));
	sohnke = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(sel->sohnke_checkbox));
	unique_axis = gtk_combo_box_get_active_id(GTK_COMBO_BOX(sel->unique_axis_combo))[0];

	s = gtk_combo_box_get_active_id(GTK_COMBO_BOX(sel->lattice_combo));
	lattice_type = lattice_from_str(s);
	if ( !has_unique_axis(lattice_type) ) {
		gtk_widget_hide(sel->unique_axis_combo);
		unique_axis = 'c';  /* Avoid showing 'weird' message */
	} else {
		gtk_widget_show(sel->unique_axis_combo);
	}

	if ( (lattice_type == L_HEXAGONAL) && (unique_axis != 'c') ) {
		gtk_widget_show(sel->weird);
	} else if ( unique_axis == 'a' ) {
		gtk_widget_show(sel->weird);
	} else if ( (lattice_type != L_MONOCLINIC) && (unique_axis == 'b') ) {
		gtk_widget_show(sel->weird);
	} else {
		gtk_widget_hide(sel->weird);
	}

	if ( lattice_type == L_RHOMBOHEDRAL ) {
		gtk_widget_show(sel->rhombo1);
		gtk_widget_show(sel->rhombo2);
	} else {
		gtk_widget_hide(sel->rhombo1);
		gtk_widget_hide(sel->rhombo2);
	}

	if ( !sohnke && !centro ) {
		gtk_widget_show(sel->nonbio);
	} else {
		gtk_widget_hide(sel->nonbio);
	}

	gtk_flow_box_invalidate_filter(GTK_FLOW_BOX(sel->flowbox));
}


static void open_selector(CrystFELSymmetrySelector *sel, gpointer data)
{
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *frame;
	int i;
	const char *pg;
	GtkWidget *parent;

	/* Don't open selector twice */
	if ( sel->dialog != NULL ) return;

	parent = gtk_widget_get_toplevel(GTK_WIDGET(sel));
	if ( !GTK_IS_WINDOW(parent) ) {
		parent = NULL;
	}
	sel->dialog = gtk_dialog_new_with_buttons("Choose point group",
	                                          GTK_WINDOW(parent),
	                                          GTK_DIALOG_DESTROY_WITH_PARENT,
	                                          "Cancel", GTK_RESPONSE_CANCEL,
	                                          NULL);

	vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);
	gtk_container_set_border_width(GTK_CONTAINER(vbox), 4);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(sel->dialog));
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	g_signal_connect(G_OBJECT(sel->dialog), "response",
	                 G_CALLBACK(selector_response_sig), sel);

	/* Lattice type chooser */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Lattice type:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	sel->lattice_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(sel->lattice_combo),
	                   FALSE, FALSE, 0);
	add_lattice_types(GTK_COMBO_BOX_TEXT(sel->lattice_combo));
	gtk_combo_box_set_active(GTK_COMBO_BOX(sel->lattice_combo), 0);

	/* Sohnke/centrosymmetric buttons */
	sel->sohnke_checkbox = gtk_check_button_new_with_label("Sohnke "
	              "(preserve anomalous signal for biological structures)");
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(sel->sohnke_checkbox),
	                   FALSE, FALSE, 4.0);
	sel->centro_checkbox = gtk_check_button_new_with_label("Centrosymmetric "
	              "(merge Friedel pairs for biological structures)");
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(sel->centro_checkbox),
	                   FALSE, FALSE, 4.0);
	deselect_when_active(sel->sohnke_checkbox, sel->centro_checkbox);
	deselect_when_active(sel->centro_checkbox, sel->sohnke_checkbox);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sel->sohnke_checkbox), TRUE);

	sel->nonbio = gtk_label_new("For biological structures, you "
	                            "should select one of the above");
	gtk_label_set_markup(GTK_LABEL(sel->nonbio),
	                     "<b><i>For biological structures, you "
	                     "should select one of the above</i></b>");
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(sel->nonbio),
	                   FALSE, FALSE, 4.0);

	/* Unique axis chooser */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Unique axis:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	sel->unique_axis_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(sel->unique_axis_combo),
	                   FALSE, FALSE, 0);
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(sel->unique_axis_combo),
	                          "a", "a (very uncommon)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(sel->unique_axis_combo),
	                          "b", "b (most common for monoclinic)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(sel->unique_axis_combo),
	                          "c", "c (most common for everything else)");
	gtk_combo_box_set_active_id(GTK_COMBO_BOX(sel->unique_axis_combo), "c");

	sel->weird = gtk_label_new("This unique axis choice is highly unconventional.");
	gtk_label_set_markup(GTK_LABEL(sel->weird),
	                     "<b><i>This unique axis choice is highly unconventional</i></b>");
	gtk_box_pack_start(GTK_BOX(vbox), sel->weird,
	                   FALSE, FALSE, 4.0);

	/* The point groups themselves */
	sel->pointgroups = malloc(MAX_PG_WIDGETS*sizeof(struct pointgroup_widget));
	if ( sel->pointgroups == NULL ) return;

	frame = gtk_frame_new("Possible point groups");
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), FALSE, FALSE, 4.0);

	sel->flowbox = gtk_flow_box_new();
	gtk_container_add(GTK_CONTAINER(frame), sel->flowbox);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);

	add_pointgroup(sel, "1",         0, 1, "*", L_TRICLINIC);
	add_pointgroup(sel, "-1",        1, 0, "*", L_TRICLINIC);
	add_pointgroup(sel, "2_uaa",     0, 1, "a", L_MONOCLINIC);
	add_pointgroup(sel, "2_uab",     0, 1, "b", L_MONOCLINIC);
	add_pointgroup(sel, "2",         0, 1, "c", L_MONOCLINIC);
	add_pointgroup(sel, "2/m_uaa",   1, 0, "a", L_MONOCLINIC);
	add_pointgroup(sel, "2/m_uab",   1, 0, "b", L_MONOCLINIC);
	add_pointgroup(sel, "2/m",       1, 0, "c", L_MONOCLINIC);
	add_pointgroup(sel, "m_uaa",     0, 0, "a", L_MONOCLINIC);
	add_pointgroup(sel, "m_uab",     0, 0, "b", L_MONOCLINIC);
	add_pointgroup(sel, "m",         0, 0, "c", L_MONOCLINIC);
	add_pointgroup(sel, "222",       0, 1, "*", L_ORTHORHOMBIC);
	add_pointgroup(sel, "mm2",       0, 0, "*", L_ORTHORHOMBIC);
	add_pointgroup(sel, "mmm",       1, 0, "*", L_ORTHORHOMBIC);

	add_pointgroup(sel, "4",         0, 1, "c", L_TETRAGONAL);
	add_pointgroup(sel, "4/m",       1, 0, "c", L_TETRAGONAL);
	add_pointgroup(sel, "422",       0, 1, "c", L_TETRAGONAL);
	add_pointgroup(sel, "4mm",       0, 0, "c", L_TETRAGONAL);
	add_pointgroup(sel, "-4",        0, 0, "c", L_TETRAGONAL);
	add_pointgroup(sel, "-42m",      0, 0, "c", L_TETRAGONAL);
	add_pointgroup(sel, "-4m2",      0, 0, "c", L_TETRAGONAL);
	add_pointgroup(sel, "4/mmm",     1, 0, "c", L_TETRAGONAL);

	add_pointgroup(sel, "4_uaa",     0, 1, "a", L_TETRAGONAL);
	add_pointgroup(sel, "4/m_uaa",   1, 0, "a", L_TETRAGONAL);
	add_pointgroup(sel, "422_uaa",   0, 1, "a", L_TETRAGONAL);
	add_pointgroup(sel, "4mm_uaa",   0, 0, "a", L_TETRAGONAL);
	add_pointgroup(sel, "-4_uaa",    0, 0, "a", L_TETRAGONAL);
	add_pointgroup(sel, "-42m_uaa",  0, 0, "a", L_TETRAGONAL);
	add_pointgroup(sel, "-4m2_uaa",  0, 0, "a", L_TETRAGONAL);
	add_pointgroup(sel, "4/mmm_uaa", 1, 0, "a", L_TETRAGONAL);

	add_pointgroup(sel, "4_uab",     0, 1, "b", L_TETRAGONAL);
	add_pointgroup(sel, "4/m_uab",   1, 0, "b", L_TETRAGONAL);
	add_pointgroup(sel, "422_uab",   0, 1, "b", L_TETRAGONAL);
	add_pointgroup(sel, "4mm_uab",   0, 0, "b", L_TETRAGONAL);
	add_pointgroup(sel, "-4_uab",    0, 0, "b", L_TETRAGONAL);
	add_pointgroup(sel, "-42m_uab",  0, 0, "b", L_TETRAGONAL);
	add_pointgroup(sel, "-4m2_uab",  0, 0, "b", L_TETRAGONAL);
	add_pointgroup(sel, "4/mmm_uab", 1, 0, "b", L_TETRAGONAL);

	add_pointgroup(sel, "3_R",       0, 1, "*", L_RHOMBOHEDRAL);
	add_pointgroup(sel, "32_R",      0, 1, "*", L_RHOMBOHEDRAL);
	add_pointgroup(sel, "3m_R",      0, 0, "*", L_RHOMBOHEDRAL);
	add_pointgroup(sel, "-3_R",      1, 0, "*", L_RHOMBOHEDRAL);
	add_pointgroup(sel, "-3m_R",     1, 0, "*", L_RHOMBOHEDRAL);

	add_pointgroup(sel, "3_H",       0, 1, "c", L_HEXAGONAL);
	add_pointgroup(sel, "-3_H",      1, 0, "c", L_HEXAGONAL);
	add_pointgroup(sel, "312_H",     0, 1, "c", L_HEXAGONAL);
	add_pointgroup(sel, "321_H",     0, 1, "c", L_HEXAGONAL);
	add_pointgroup(sel, "-31m_H",    1, 0, "c", L_HEXAGONAL);
	add_pointgroup(sel, "-3m1_H",    1, 0, "c", L_HEXAGONAL);
	add_pointgroup(sel, "3m1_H",     0, 0, "c", L_HEXAGONAL);
	add_pointgroup(sel, "31m_H",     0, 0, "c", L_HEXAGONAL);

	add_pointgroup(sel, "3_H_uaa",    0, 1, "a", L_HEXAGONAL);
	add_pointgroup(sel, "-3_H_uaa",   1, 0, "a", L_HEXAGONAL);
	add_pointgroup(sel, "312_H_uaa",  0, 1, "a", L_HEXAGONAL);
	add_pointgroup(sel, "321_H_uaa",  0, 1, "a", L_HEXAGONAL);
	add_pointgroup(sel, "-31m_H_uaa", 1, 0, "a", L_HEXAGONAL);
	add_pointgroup(sel, "-3m1_H_uaa", 1, 0, "a", L_HEXAGONAL);
	add_pointgroup(sel, "3m1_H_uaa",  0, 0, "a", L_HEXAGONAL);
	add_pointgroup(sel, "31m_H_uaa",  0, 0, "a", L_HEXAGONAL);

	add_pointgroup(sel, "3_H_uab",    0, 1, "b", L_HEXAGONAL);
	add_pointgroup(sel, "-3_H_uab",   1, 0, "b", L_HEXAGONAL);
	add_pointgroup(sel, "312_H_uab",  0, 1, "b", L_HEXAGONAL);
	add_pointgroup(sel, "321_H_uab",  0, 1, "b", L_HEXAGONAL);
	add_pointgroup(sel, "-31m_H_uab", 1, 0, "b", L_HEXAGONAL);
	add_pointgroup(sel, "-3m1_H_uab", 1, 0, "b", L_HEXAGONAL);
	add_pointgroup(sel, "3m1_H_uab",  0, 0, "b", L_HEXAGONAL);
	add_pointgroup(sel, "31m_H_uab",  0, 0, "b", L_HEXAGONAL);

	add_pointgroup(sel, "6",         0, 1, "c", L_HEXAGONAL);
	add_pointgroup(sel, "6/m",       1, 0, "c", L_HEXAGONAL);
	add_pointgroup(sel, "622",       0, 1, "c", L_HEXAGONAL);
	add_pointgroup(sel, "6/mmm",     1, 0, "c", L_HEXAGONAL);
	add_pointgroup(sel, "6mm",       0, 0, "c", L_HEXAGONAL);
	add_pointgroup(sel, "-6m2",      0, 0, "c", L_HEXAGONAL);
	add_pointgroup(sel, "-62m",      0, 0, "c", L_HEXAGONAL);

	add_pointgroup(sel, "6_uaa",     0, 1, "a", L_HEXAGONAL);
	add_pointgroup(sel, "6/m_uaa",   1, 0, "a", L_HEXAGONAL);
	add_pointgroup(sel, "622_uaa",   0, 1, "a", L_HEXAGONAL);
	add_pointgroup(sel, "6/mmm_uaa", 1, 0, "a", L_HEXAGONAL);
	add_pointgroup(sel, "6mm_uaa",   0, 0, "a", L_HEXAGONAL);
	add_pointgroup(sel, "-6m2_uaa",  0, 0, "a", L_HEXAGONAL);
	add_pointgroup(sel, "-62m_uaa",  0, 0, "a", L_HEXAGONAL);

	add_pointgroup(sel, "6_uab",     0, 1, "b", L_HEXAGONAL);
	add_pointgroup(sel, "6/m_uab",   1, 0, "b", L_HEXAGONAL);
	add_pointgroup(sel, "622_uab",   0, 1, "b", L_HEXAGONAL);
	add_pointgroup(sel, "6/mmm_uab", 1, 0, "b", L_HEXAGONAL);
	add_pointgroup(sel, "6mm_uab",   0, 0, "b", L_HEXAGONAL);
	add_pointgroup(sel, "-6m2_uab",  0, 0, "b", L_HEXAGONAL);
	add_pointgroup(sel, "-62m_uab",  0, 0, "b", L_HEXAGONAL);

	add_pointgroup(sel, "23",        0, 1, "*", L_CUBIC);
	add_pointgroup(sel, "m-3",       1, 0, "*", L_CUBIC);
	add_pointgroup(sel, "432",       0, 1, "*", L_CUBIC);
	add_pointgroup(sel, "-43m",      0, 0, "*", L_CUBIC);
	add_pointgroup(sel, "m-3m",      1, 0, "*", L_CUBIC);

	gtk_flow_box_set_filter_func(GTK_FLOW_BOX(sel->flowbox),
	                             filter_pgs, sel, NULL);
	g_object_set(G_OBJECT(sel->flowbox), "margin", 10, NULL);

	sel->rhombo1 = gtk_label_new("Looking for 'H3' / 'H32'?  Select "
	                            "hexagonal and use '3_H' / '321_H',");
	gtk_label_set_markup(GTK_LABEL(sel->rhombo1),
	                     "<i>Looking for 'H3' / 'H32'?  Select hexagonal "
	                     "and use '3_H' / '321_H'</i>");
	gtk_box_pack_start(GTK_BOX(vbox), sel->rhombo1, FALSE, FALSE, 2.0);
	sel->rhombo2 = gtk_label_new("respectively, or -3_H / -3m1_H to merge Friedel pairs");
	gtk_label_set_markup(GTK_LABEL(sel->rhombo2),
	                     "<i>respectively, or -3_H / -3m1_H to merge Friedel pairs</i>");
	gtk_box_pack_start(GTK_BOX(vbox), sel->rhombo2, FALSE, FALSE, 2.0);

	g_signal_connect(G_OBJECT(sel->unique_axis_combo), "changed",
	                 G_CALLBACK(update_opts_sig), sel);
	g_signal_connect(G_OBJECT(sel->sohnke_checkbox), "toggled",
	                 G_CALLBACK(update_opts_sig), sel);
	g_signal_connect(G_OBJECT(sel->centro_checkbox), "toggled",
	                 G_CALLBACK(update_opts_sig), sel);
	g_signal_connect(G_OBJECT(sel->lattice_combo), "changed",
	                 G_CALLBACK(update_opts_sig), sel);

	label = gtk_bin_get_child(GTK_BIN(sel));
	pg = gtk_label_get_text(GTK_LABEL(label));
	for ( i=0; i<sel->n_pgs; i++ ) {
		if ( strcmp(sel->pointgroups[i].symbol, pg) == 0 ) {
			set_active(sel->sohnke_checkbox,
			           sel->pointgroups[i].sohnke);
			set_active(sel->centro_checkbox,
			           sel->pointgroups[i].centro);
			set_combo_id(sel->unique_axis_combo,
			             sel->pointgroups[i].unique_axis);
			set_combo_id(sel->lattice_combo,
			             str_lattice(sel->pointgroups[i].lattice_type));
		}
	}

	gtk_widget_show_all(sel->dialog);
	update_opts_sig(NULL, sel);
}


GtkWidget *crystfel_symmetry_selector_new()
{
	CrystFELSymmetrySelector *sel;

	sel = g_object_new(CRYSTFEL_TYPE_SYMMETRY_SELECTOR, NULL);

	sel->dialog = NULL;
	sel->have_pg = 0;
	sel->n_pgs = 0;

	sel->label = gtk_label_new("Click to choose");
	gtk_container_add(GTK_CONTAINER(sel), sel->label);

	g_signal_connect(G_OBJECT(sel), "clicked",
	                 G_CALLBACK(open_selector), NULL);

	return GTK_WIDGET(sel);
}


char *crystfel_symmetry_selector_get_group_symbol(CrystFELSymmetrySelector *sel)
{
	if ( sel->have_pg ) {
		GtkWidget *label = gtk_bin_get_child(GTK_BIN(sel));
		return strdup(gtk_label_get_text(GTK_LABEL(label)));
	} else {
		return strdup("none");
	}
}


int crystfel_symmetry_selector_set_group_symbol(CrystFELSymmetrySelector *sel,
                                                const char *pg_symbol)
{
	if ( pg_symbol != NULL ) {
		GtkWidget *label = gtk_bin_get_child(GTK_BIN(sel));
		gtk_label_set_text(GTK_LABEL(label), pg_symbol);
		sel->have_pg = 1;
	}
	return 0;
}
