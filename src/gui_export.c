/*
 * gui_export.c
 *
 * Export data from CrystFEL GUI
 *
 * Copyright © 2021 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2021 Thomas White <taw@physics.org>
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

#ifdef HAVE_LIBCCP4
#include <ccp4/cmtzlib.h>
#include <ccp4/csymlib.h>
#include <ccp4/ccp4_parser.h>
#endif

#include <utils.h>
#include <reflist-utils.h>
#include <cell-utils.h>
#include <symmetry.h>

#include "version.h"
#include "gui_project.h"
#include "crystfel_gui.h"
#include "gtk-util-routines.h"


struct export_window
{
	struct crystfelproject *proj;
	GtkWidget *cell_chooser;
	GtkWidget *limit_res;
	GtkWidget *min_res;
	GtkWidget *max_res;
	GtkWidget *dataset;
	GtkWidget *format;
};


struct point_group_conversion
{
	char centering;
	const char *crystfel;
	int friedel;

	int xds_spgnum;

	const char *ccp4;
};


/* Table for converting CrystFEL's knowledge about centering, point group
 * and unique axis into something that can be recognised by external programs.
 * If xds_spgnum=0, ccp4=NULL, or something is missing form this table entirely,
 * it means that there is no way to represent the situation to that program
 * without re-indexing the dataset. */
struct point_group_conversion pg_conversions[] = {

	/* Triclinic  */
	{'P', "1",       0,      1,     "P 1"},
	{'P', "-1",      1,      1,     "P 1"},

	/* Monoclinic */
	{'P', "2_uaa",   0,      0,     "P211"},
	{'P', "m_uaa",   0,      0,     "Pm11"},
	{'P', "2/m_uaa", 1,      0,     "P211"},
	{'P', "2_uab",   0,      3,     "P121"},
	{'P', "m_uab",   0,      0,     "P1m1"},
	{'P', "2/m_uab", 1,      3,     "P121"},
	{'P', "2_uac",   0,      0,     "P112"},
	{'P', "m_uac",   0,      0,     "P11m"},
	{'P', "2/m_uac", 1,      0,     "P112"},
	{'P', "2",       0,      0,     "P121"}, /* unique axis c */
	{'P', "m",       0,      0,     "P11m"}, /* unique axis c */
	{'P', "2/m",     1,      0,     "P121"}, /* unique axis c */

	{'A', "2_uab",   0,      0,     "A121"},
	{'A', "m_uab",   0,      0,     "A1m1"},
	{'A', "2/m_uab", 1,      0,     "A121"},
	{'A', "2_uac",   0,      0,     "A112"},
	{'A', "m_uac",   0,      0,     "A11m"},
	{'A', "2/m_uac", 1,      0,     "A112"},
	{'A', "2",       0,      0,     "A121"}, /* unique axis c */
	{'A', "m",       0,      0,     "A11m"}, /* unique axis c */
	{'A', "2/m",     1,      0,     "A121"}, /* unique axis c */

	{'B', "2_uaa",   0,      0,     "B211"},
	{'B', "m_uaa",   0,      0,     "Bm11"},
	{'B', "2/m_uaa", 1,      0,     "B211"},
	{'B', "2_uac",   0,      0,     "B112"},
	{'B', "m_uac",   0,      0,     "B11m"},
	{'B', "2/m_uac", 1,      0,     "B112"},
	{'B', "2",       0,      0,     "B112"}, /* unique axis c */
	{'B', "m",       0,      0,     "B11m"}, /* unique axis c */
	{'B', "2/m",     1,      0,     "B112"}, /* unique axis c */

	{'C', "2_uaa",   0,      0,     "C211"},
	{'C', "m_uaa",   0,      0,     "Cm11"},
	{'C', "2/m_uaa", 1,      0,     "C211"},
	{'C', "2_uab",   0,      5,     "C121"},
	{'C', "m_uab",   0,      0,     "C1m1"},
	{'C', "2/m_uab", 1,      5,     "C121"},

	{'I', "2_uaa",   0,      0,     "I211"},
	{'I', "m_uaa",   0,      0,     "Im11"},
	{'I', "2/m_uaa", 1,      0,     "I211"},
	{'I', "2_uab",   0,      0,     "I121"},
	{'I', "m_uab",   0,      0,     "I1m1"},
	{'I', "2/m_uab", 1,      0,     "I121"},
	{'I', "2_uac",   0,      0,     "I112"},
	{'I', "m_uac",   0,      0,     "I11m"},
	{'I', "2/m_uac", 1,      0,     "I112"},
	{'I', "2",       0,      0,     "I121"}, /* unique axis c */
	{'I', "m",       0,      0,     "I11m"}, /* unique axis c */
	{'I', "2/m",     1,      0,     "I121"}, /* unique axis c */

	/* Orthorhombic */
	{'P', "222",       0,     16,     "P222"},
	{'P', "mmm",       1,     16,     "P222"},
	{'P', "mm2",       0,     25,     "Pmm2"},
	{'A', "222",       0,      0,     "A222"},
	{'A', "mmm",       1,      0,     "A222"},
	{'A', "mm2",       0,     38,     "Amm2"},
	{'B', "222",       0,      0,     "B222"},
	{'B', "mmm",       1,      0,     "B222"},
	{'B', "mm2",       0,      0,     "Bmm2"},
	{'C', "222",       0,     21,     "C222"},
	{'C', "mmm",       1,     21,     "C222"},
	{'C', "mm2",       0,     35,     "Cmm2"},
	{'F', "222",       0,     22,     "F222"},
	{'F', "mmm",       1,     22,     "F222"},
	{'F', "mm2",       0,     42,     "Fmm2"},
	{'I', "222",       0,     23,     "I222"},
	{'I', "mmm",       1,     23,     "I222"},
	{'I', "mm2",       0,     45,     "Imm2"},

	/* Tetragonal */
	{'P', "4",         0,     75,     "P4"},    /* unique axis c */
	{'P', "4/m",       1,     75,     "P4"},    /* unique axis c */
	{'P', "422",       0,     89,     "P422"},  /* unique axis c */
	{'P', "4/mmm",     1,     89,     "P422"},  /* unique axis c */
	{'P', "4mm",       0,     99,     "P4mm"},  /* unique axis c */
	{'P', "-4",        0,     81,     "P-4"},   /* unique axis c */
	{'P', "-42m",      0,    111,     "P-42m"}, /* unique axis c */
	{'P', "-4m2",      0,    115,     "P-4m2"}, /* unique axis c */
	{'P', "4_uac",     0,     75,     "P4"},
	{'P', "4/m_uac",   1,     75,     "P4"},
	{'P', "422_uac",   0,     89,     "P422"},
	{'P', "4/mmm_uac", 1,     89,     "P422"},
	{'P', "4mm_uac",   0,     99,     "P4mm"},
	{'P', "-4_uac",    0,     81,     "P-4"},
	{'P', "-42m_uac",  0,    111,     "P-42m"},
	{'P', "-4m2_uac",  0,    115,     "P-4m2"},
	{'I', "4",         0,     79,     "I4"},    /* unique axis c */
	{'I', "4/m",       1,     79,     "I4"},    /* unique axis c */
	{'I', "422",       0,     97,     "I422"},  /* unique axis c */
	{'I', "4/mmm",     1,     97,     "I422"},  /* unique axis c */
	{'I', "4mm",       0,    107,     "I4mm"},  /* unique axis c */
	{'I', "-4",        0,     82,     "I-4"},   /* unique axis c */
	{'I', "-42m",      0,    121,     "I-42m"}, /* unique axis c */
	{'I', "-4m2",      0,    119,     "I-4m2"}, /* unique axis c */
	{'I', "4_uac",     0,     79,     "I4"},
	{'I', "4/m_uac",   1,     79,     "I4"},
	{'I', "422_uac",   0,     97,     "I422"},
	{'I', "4/mmm_uac", 1,     97,     "I422"},
	{'I', "4mm_uac",   0,    107,     "I4mm"},
	{'I', "-4_uac",    0,     82,     "I-4"},
	{'I', "-42m_uac",  0,    121,     "I-42m"},
	{'I', "-4m2_uac",  0,    119,     "I-4m2"},

	/* Trigonal (rhombohedral) */
	{'R', "3_R",       0,      0,     "R3:R"},
	{'R', "-3_R",      1,      0,     "R3:R"},
	{'R', "32_R",      0,      0,     "R32:R"},
	{'R', "-3m_R",     1,      0,     "R32:R"},
	{'R', "3m_R",      0,      0,     "R3m:R"},

	/* Trigonal (rhombohedral on hexagonal axes) */
	{'H', "3_H",       0,    146,     "R3:H"},
	{'H', "-3_H",      1,    146,     "R3:H"},
	{'H', "32_H",      0,    155,     "R3:H"},
	{'H', "-3m_H",     1,    155,     "R3:H"},
	{'H', "3m_H",      0,      0,     "R3m:H"},

	/* Trigonal (hexagonal) */
	{'P', "3_H",       0,    143,     "P3"},
	{'P', "-3_H",      1,    143,     "P3"},
	{'P', "312_H",     0,    149,     "P312"},
	{'P', "-31m_H",    1,    149,     "P312"},
	{'P', "321_H",     0,    150,     "P321"},
	{'P', "-3m1_H",    1,    150,     "P321"},
	{'P', "3m1_H",     0,    156,     "P3m1"},
	{'P', "31m_H",     0,    157,     "P31m"},

	/* Hexagonal */
	{'P', "6",         0,    168,     "P6"},
	{'P', "6/m",       1,    168,     "P6"},
	{'P', "622",       0,    177,     "P622"},
	{'P', "6/mmm",     1,    177,     "P622"},
	{'P', "6mm",       0,    177,     "P6mm"},
	{'P', "-6m2",      0,    187,     "P-6m2"},
	{'P', "-62m",      0,    189,     "P-62m"},

	/* Cubic */
	{'P', "23",        0,    195,     "P23"},
	{'P', "m-3",       1,    195,     "P23"},
	{'P', "432",       0,    207,     "P432"},
	{'P', "m-3m",      1,    207,     "P432"},
	{'P', "-43m",      0,    215,     "P -4 3 m"},
	{'I', "23",        0,    197,     "I23"},
	{'I', "m-3",       1,    197,     "I23"},
	{'I', "432",       0,    211,     "I432"},
	{'I', "m-3m",      1,    211,     "I432"},
	{'I', "-43m",      0,    217,     "I -4 3 m"},
	{'F', "23",        0,    196,     "F23"},
	{'F', "m-3",       1,    196,     "F23"},
	{'F', "432",       0,    209,     "F432"},
	{'F', "m-3m",      1,    209,     "F432"},
	{'F', "-43m",      0,    216,     "F -4 3 m"},

	{'*', NULL,  0, 0, NULL}
};


static int space_group_for_xds(const char *sym_str, char cen)
{
	int i = 0;
	do {
		if ( (pg_conversions[i].centering == cen)
		  && (strcmp(sym_str, pg_conversions[i].crystfel) == 0) )
		{
			return pg_conversions[i].xds_spgnum;
		}
		i++;
	} while (pg_conversions[i].centering != '*');

	ERROR("Couldn't derive XDS representation of symmetry.\n");
	return 0;
}


static const char *space_group_for_mtz(const char *sym_str, char cen)
{
	int i = 0;
	do {
		if ( (pg_conversions[i].centering == cen)
		  && (strcmp(sym_str, pg_conversions[i].crystfel) == 0) )
		{
			return pg_conversions[i].ccp4;
		}
		i++;
	} while (pg_conversions[i].centering != '*');

	ERROR("Couldn't derive CCP4 representation of symmetry.\n");
	return NULL;
}


static int export_to_xds(struct gui_merge_result *result,
                         const char *filename, UnitCell *cell,
                         double min_res, double max_res)
{
	FILE *fh;
	RefList *reflist;
	RefListIterator *iter;
	Reflection *refl;
	double a, b, c, al, be,ga;
	char *sym_str;
	SymOpList *sym;
	int spg;

	fh = fopen(filename, "w");
	if ( fh == NULL ) return 1;

	reflist = read_reflections_2(result->hkl, &sym_str);
	if ( reflist == NULL ) return 1;
	if ( sym_str == NULL ) return 1;

	sym = get_pointgroup(sym_str);
	if ( sym == NULL ) return 1;

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);

	spg = space_group_for_xds(sym_str, cell_get_centering(cell));
	if ( spg == 0 ) return 1;

	fprintf(fh, "!FORMAT=XDS_ASCII MERGE=TRUE FRIEDEL'S_LAW=%s\n",
	        is_centrosymmetric(sym) ? "TRUE" : "FALSE");
	fprintf(fh, "!SPACE_GROUP_NUMBER=%i\n", spg);
	fprintf(fh, "!UNIT_CELL_CONSTANT= %.2f %.2f %.2f %.2f %.2f %.2f\n",
	        a*1e10, b*1e10, c*1e10, rad2deg(al), rad2deg(be), rad2deg(ga));
	fprintf(fh, "!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=5\n");
	fprintf(fh, "!ITEM_H=1\n");
	fprintf(fh, "!ITEM_K=2\n");
	fprintf(fh, "!ITEM_L=3\n");
	fprintf(fh, "!ITEM_IOBS=4\n");
	fprintf(fh, "!ITEM_SIGMA(IOBS)=5\n");
	fprintf(fh, "!END_OF_HEADER\n");

	for ( refl = first_refl(reflist, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double one_over_d;

		get_indices(refl, &h, &k, &l);

		one_over_d = 2.0*resolution(cell, h, k, l);
		if ( (one_over_d > min_res) && (one_over_d < max_res) ) {

			fprintf(fh, "%6i %6i %6i %9.2f %9.2f\n",
			        h, k, l,
			        get_intensity(refl),
			        get_esd_intensity(refl));

		}
	}

	fprintf(fh, "!END_OF_DATA\n");
	free_symoplist(sym);
	free(sym_str);
	reflist_free(reflist);

	fclose(fh);
	return 0;
}


#ifdef HAVE_LIBCCP4
static CCP4SPG *add_mtz_symmetry_header(MTZ *mtz, const char *spg_name)
{
	CCP4SPG *spg;
	float rsymx[192][4][4];
	char ltypex[2];
	int i;

	spg = ccp4spg_load_by_spgname(spg_name);
	if ( spg == NULL ) {
		ERROR("Couldn't look up CCP4 space group '%s'\n", spg_name);
		return NULL;
	}

	for ( i=0; i<spg->nsymop; i++ ) {
		rotandtrn_to_mat4(rsymx[i], spg->symop[i]);
	}
	ltypex[0] = spg->symbol_old[0];
	ltypex[1] = '\0';

	ccp4_lwsymm(mtz, spg->nsymop, spg->nsymop_prim,
	            rsymx, ltypex, spg->spg_ccp4_num, spg->symbol_old,
	            spg->point_group);

	return spg;
}
#endif


static int export_to_mtz(struct gui_merge_result *result,
                         const char *filename, UnitCell *cell,
                         double min_res, double max_res)
{
#ifdef HAVE_LIBCCP4
	MTZ *mtz;
	MTZXTAL *cr;
	MTZSET *ds;
	MTZCOL *columns[7];
	double a, b, c, al, be, ga;
	int r;
	char tmp[128];
	float cellp[6];
	int refl_i;
	RefList *reflist;
	Reflection *refl;
	RefListIterator *iter;
	char *sym_str = NULL;
	CCP4SPG *spg;
	const char *spg_name;
	SymOpList *sym;

	reflist = read_reflections_2(result->hkl, &sym_str);
	if ( reflist == NULL ) return 1;
	if ( sym_str == NULL ) return 1;

	sym = get_pointgroup(sym_str);
	if ( sym == NULL ) return 1;

	spg_name = space_group_for_mtz(sym_str, cell_get_centering(cell));
	if ( spg_name == NULL ) {
		reflist_free(reflist);
		return 1;
	}

	mtz = MtzMalloc(0, 0);

	snprintf(tmp, 128, "Data exported via CrystFEL GUI, version %s",
	         crystfel_version_string());
	ccp4_lwtitl(mtz, tmp, 0);

	mtz->refs_in_memory = 0;
	mtz->fileout = MtzOpenForWrite(filename);

	spg = add_mtz_symmetry_header(mtz, spg_name);
	if ( spg == NULL ) {
		return 1;
	}

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
	cellp[0] = a*1e10;
	cellp[1] = b*1e10;
	cellp[2] = c*1e10;
	cellp[3] = rad2deg(al);
	cellp[4] = rad2deg(be);
	cellp[5] = rad2deg(ga);

	/* FIXME: Proposed labelling:
	 *  title = as above
	 *  project = basename of folder containing crystfel.project
	 *  crystal = name of indexing results run
	 *  dataset = name of merge results run */
	cr = MtzAddXtal(mtz, "Crystal_name", "Project_name", cellp);
	ds = MtzAddDataset(mtz, cr, result->name, 0.0);
	columns[0] = MtzAddColumn(mtz, ds, "H", "H");
	columns[1] = MtzAddColumn(mtz, ds, "K", "H");
	columns[2] = MtzAddColumn(mtz, ds, "L", "H");
	columns[3] = MtzAddColumn(mtz, ds, "I+", "J");
	columns[4] = MtzAddColumn(mtz, ds, "SIGI+", "Q");
	columns[5] = MtzAddColumn(mtz, ds, "I-", "J");
	columns[6] = MtzAddColumn(mtz, ds, "SIGI-", "Q");

	refl_i = 1;
	for ( refl = first_refl(reflist, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double one_over_d;
		int isym;

		get_indices(refl, &h, &k, &l);

		one_over_d = 2.0*resolution(cell, h, k, l);
		if ( (one_over_d > min_res) && (one_over_d < max_res) ) {

			float refldata[7];
			signed int nh, nk, nl;
			signed int fh, fk, fl;
			Reflection *friedel;
			Reflection *refl_plus;
			Reflection *refl_minus;

			/* Look for Friedel partner */
			if ( find_equiv_in_list(reflist, -h, -k, -l,
			                        sym, &fh, &fk, &fl) )
			{
				friedel = find_refl(reflist, fh, fk, fl);
			} else {
				friedel = NULL;
			}

			/* Move to CCP4's idea of the ASU */
			isym = ccp4spg_put_in_asu(spg, h, k, l, &nh, &nk, &nl);

			/* Ok, do we have an I+ or an I- ? */
			if ( is_odd(isym) ) {
				/* I+ */
				refl_plus = refl;
				refl_minus = friedel;
			} else {
				/* I- */
				refl_minus = refl;
				refl_plus = friedel;
			}

			/* If we are looking at an I-, only write it out now
			 * if the corresponding I+ if not in 'reflist'.
			 * If I+ is present, then this I- will get written when
			 * the Friedel pair is processed. */
			if ( !is_odd(isym) && (refl_plus != NULL) ) continue;

			refldata[0] = nh;
			refldata[1] = nk;
			refldata[2] = nl;
			if ( refl_plus != NULL ) {
				refldata[3] = get_intensity(refl_plus);
				refldata[4] = get_esd_intensity(refl_plus);
			} else {
				refldata[3] = NAN;
				refldata[4] = NAN;
			}
			if ( refl_minus != NULL ) {
				refldata[5] = get_intensity(refl_minus);
				refldata[6] = get_esd_intensity(refl_minus);
			} else {
				refldata[5] = NAN;
				refldata[6] = NAN;
			}

			ccp4_lwrefl(mtz, refldata, columns, 7, refl_i++);

		}
	}

	r = MtzPut(mtz, " ");
	ccp4spg_free(&spg);
	MtzFree(mtz);
	reflist_free(reflist);
	return 1-r; /* Yes, really.  MtzPut return values are backwards */
#else
	return 1;
#endif
}


static int export_data(struct export_window *win, char *filename)
{
	gchar *cell_filename;
	const char *dataset;
	const char *format;
	struct gui_merge_result *result;
	UnitCell *cell;
	int r = 0;
	double min_res = 0;
	double max_res = +INFINITY;

	dataset = gtk_combo_box_get_active_id(GTK_COMBO_BOX(win->dataset));
	if ( dataset == NULL ) {
		ERROR("Please select the dataset to export.\n");
		return 1;
	}

	format = gtk_combo_box_get_active_id(GTK_COMBO_BOX(win->format));
	if ( format == NULL ) {
		ERROR("Please select the data format to use.\n");
		return 1;
	}

	cell_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(win->cell_chooser));
	if ( cell_filename == NULL ) {
		ERROR("Please choose the unit cell file.\n");
		return 1;
	}

	cell = load_cell_from_file(cell_filename);
	if ( cell == NULL ) {
		ERROR("Failed to load unit cell file %s\n", cell_filename);
		return 1;
	}

	if ( get_bool(win->limit_res) ) {
		min_res = 1e10/get_float(win->min_res);
		max_res = 1e10/get_float(win->max_res);
	}

	result = find_merge_result_by_name(win->proj, dataset);
	if ( result == NULL ) {
		ERROR("Couldn't find merged dataset '%s'\n", dataset);
		return 1;
	}

	STATUS("Exporting dataset %s to %s, in format %s, using unit cell %s,"
	       "%f to %f m^-1\n", dataset, filename, format, cell_filename,
	       min_res, max_res);

	if ( strcmp(format, "mtz") == 0 ) {
		r = export_to_mtz(result, filename, cell, min_res, max_res);
	} else if ( strcmp(format, "xds") == 0 ) {
		r = export_to_xds(result, filename, cell, min_res, max_res);
	} else {
		ERROR("Unrecognised export format '%s'\n", format);
		return 1;
	}

	if ( r ) {
		ERROR("Export failed\n");
	}

	g_free(cell_filename);

	return 0;
}


static void export_response_sig(GtkWidget *dialog, gint resp,
                                struct export_window *win)
{
	int r = 0;

	if ( resp == GTK_RESPONSE_ACCEPT ) {
		gchar *filename;
		filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
		r = export_data(win, filename);
		g_free(filename);
	}

	if ( !r ) gtk_widget_destroy(dialog);
}


gint export_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	char tmp[64];
	struct export_window *win;
	int i;

	win = malloc(sizeof(struct export_window));
	if ( win == NULL ) return 0;

	win->proj = proj;

	dialog = gtk_file_chooser_dialog_new("Export data",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_FILE_CHOOSER_ACTION_SAVE,
	                                     GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
	                                     GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
	                                     NULL);
	gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog),
	                                               TRUE);

	g_signal_connect(G_OBJECT(dialog), "response",
	                 G_CALLBACK(export_response_sig),
	                 win);

	vbox = gtk_vbox_new(FALSE, 0.0);
	gtk_file_chooser_set_extra_widget(GTK_FILE_CHOOSER(dialog),
	                                  GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(vbox), 4);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Results to export:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->dataset = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->dataset),
	                   FALSE, FALSE, 4.0);
	for ( i=0; i<proj->n_merge_results; i++ ) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->dataset),
		                          proj->merge_results[i].name,
		                          proj->merge_results[i].name);
	}
	label = gtk_label_new("Format");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->format = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->format),
	                   FALSE, FALSE, 4.0);
#ifdef HAVE_LIBCCP4
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->format), "mtz",
	                          "MTZ");
#endif
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->format), "xds",
	                          "XDS ASCII");
	gtk_combo_box_set_active(GTK_COMBO_BOX(win->format), 0);

	label = gtk_label_new("Unit cell file:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->cell_chooser = gtk_file_chooser_button_new("Unit cell file",
	                                              GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_file_chooser_set_local_only(GTK_FILE_CHOOSER(win->cell_chooser),
	                                TRUE);
	/* Use the "FoM" cell file because there should only be one
	 * point of truth for the "final" cell parameters.  Eventually, I hope
	 * to determine this automatically. */
	gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(win->cell_chooser),
	                              proj->fom_cell_filename);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->cell_chooser),
	                   FALSE, FALSE, 4.0);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	win->limit_res = gtk_check_button_new_with_label("Restrict resolution range:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->limit_res),
	                   FALSE, FALSE, 4.0);
	win->min_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->min_res), 4);
	snprintf(tmp, 64, "%.2f", proj->export_res_min);
	gtk_entry_set_text(GTK_ENTRY(win->min_res), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->min_res),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("to");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->max_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->max_res), 4);
	snprintf(tmp, 64, "%.2f", proj->export_res_max);
	gtk_entry_set_text(GTK_ENTRY(win->max_res), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->max_res),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Å");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(win->limit_res), "toggled",
	                 G_CALLBACK(i_maybe_disable), win->min_res);
	g_signal_connect(G_OBJECT(win->limit_res), "toggled",
	                 G_CALLBACK(i_maybe_disable), win->max_res);
	gtk_widget_set_sensitive(win->min_res, FALSE);
	gtk_widget_set_sensitive(win->max_res, FALSE);

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_CLOSE);
	gtk_widget_show_all(dialog);

	return FALSE;
}
