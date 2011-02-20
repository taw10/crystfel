/*
 * dw-geomatic.c
 *
 * GUI geometry calibration
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gtk/gtk.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cairo.h>
#include <gdk-pixbuf/gdk-pixbuf.h>

#include "dw-geomatic.h"
#include "render.h"
#include "hdf5-file.h"
#include "utils.h"
#include "cell.h"
#include "geometry.h"
#include "peaks.h"


static void displaywindow_error(DWGeomatic *dw, const char *message)
{
	GtkWidget *window;

	window = gtk_message_dialog_new(GTK_WINDOW(dw->window),
					GTK_DIALOG_DESTROY_WITH_PARENT,
					GTK_MESSAGE_WARNING,
					GTK_BUTTONS_CLOSE, message);
	gtk_window_set_title(GTK_WINDOW(window), "Error");

	g_signal_connect_swapped(window, "response",
				 G_CALLBACK(gtk_widget_destroy), window);
	gtk_widget_show(window);
}


static void displaywindow_update(DWGeomatic *dw)
{
	gint width;

	if ( dw->image != NULL ) {
		dw->width = dw->image->width;
		dw->height = dw->image->height;
	} else {
		dw->width = 320;
		dw->height = 320;
	}

	width = dw->width;
	if ( dw->show_col_scale ) width += 20;

	if ( dw->pixbuf != NULL ) {
		gdk_pixbuf_unref(dw->pixbuf);
	}
	if ( dw->image != NULL ) {
		dw->pixbuf = render_get_image(dw->image, 1, dw->scale,
		                              dw->boostint);
	} else {
		dw->pixbuf = NULL;
	}

	if ( dw->col_scale != NULL ) {
		gdk_pixbuf_unref(dw->col_scale);
	}
	dw->col_scale = render_get_colour_scale(20, dw->height, dw->scale);

	gdk_window_invalidate_rect(dw->drawingarea->window, NULL, FALSE);
}


/* Window closed - clean up */
static gint displaywindow_closed(GtkWidget *window, DWGeomatic *dw)
{
	if ( dw->hdfile != NULL ) {
		hdfile_close(dw->hdfile);
	}

	exit(0);  /* This program only handles one image at a time */
}


static gboolean displaywindow_expose(GtkWidget *da, GdkEventExpose *event,
				     DWGeomatic *dw)
{
	cairo_t *cr;
	Reflection *refl;
	RefListIterator *iter;

	cr = gdk_cairo_create(da->window);

	/* Blank white background */
	cairo_rectangle(cr, 0.0, 0.0, da->allocation.width,
			da->allocation.height);
	cairo_set_source_rgb(cr, 0.0, 0.0, 1.0);
	cairo_fill(cr);

	if ( dw->pixbuf != NULL ) {
		gdk_draw_pixbuf(da->window,
				da->style->bg_gc[GTK_WIDGET_STATE(da)],
				dw->pixbuf,
				0, 0, 0, 0, dw->width, dw->height,
				GDK_RGB_DITHER_NONE, 0, 0);
	}

	if ( (dw->show_col_scale) && (dw->col_scale != NULL) ) {
		gdk_draw_pixbuf(da->window,
				da->style->bg_gc[GTK_WIDGET_STATE(da)],
				dw->col_scale,
				0, 0, dw->width, 0, 20, dw->height,
				GDK_RGB_DITHER_NONE, 0, 0);
	}

	if ( dw->image->det != NULL ) {

		RefList *peaks;
		UnitCell *rot;

		dw->image->bw = 0.01;
		dw->image->div = 0.0;
		dw->image->profile_radius = 0.0;

		rot = rotate_cell(dw->cell, dw->pos_x/100000.0, dw->pos_y/100000.0, 0.0);
		peaks = find_intersections(dw->image, rot, 0);
		cell_free(rot);

		for ( refl = first_refl(peaks, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) ) {

			double x, y;

			get_detector_pos(refl, &x, &y);
			cairo_new_path(cr);
			cairo_arc(cr, x, y, 3.0, 0.0, 2.0*M_PI);
			cairo_set_source_rgb(cr, 0.0, 0.0, 1.0);
			cairo_stroke(cr);

		}

		reflist_free(peaks);

	}

	return FALSE;
}


static gint displaywindow_close(GtkWidget *widget, DWGeomatic *dw)
{
	gtk_widget_destroy(dw->window);
	return 0;
}


static gint displaywindow_set_boostint_response(GtkWidget *widget,
						gint response,
						DWGeomatic *dw)
{
	int done = 1;

	if ( response == GTK_RESPONSE_OK ) {

		const char *sboostint;
		float boostint;
		int scanval;

		sboostint = gtk_entry_get_text(
		 			GTK_ENTRY(dw->boostint_dialog->entry));
		scanval = sscanf(sboostint, "%f", &boostint);
		if ( (scanval != 1) || (boostint <= 0) ) {
			displaywindow_error(dw, "Please enter a positive "
					"value for the intensity boost "
					"factor.");
			done = 0;
		} else {
			dw->boostint = boostint;
			displaywindow_update(dw);
		}
	}

	if ( done ) {
		gtk_widget_destroy(dw->boostint_dialog->window);
	}

	return 0;
}


static gint displaywindow_set_boostint_destroy(GtkWidget *widget,
					       DWGeomatic *dw)
{
	free(dw->boostint_dialog);
	dw->boostint_dialog = NULL;
	return 0;
}


static gint displaywindow_set_boostint_response_ac(GtkWidget *widget,
						   DWGeomatic *dw)
{
	return displaywindow_set_boostint_response(widget, GTK_RESPONSE_OK, dw);
}


/* Create a window to ask the user for a new intensity boost factor */
static gint displaywindow_set_boostint(GtkWidget *widget, DWGeomatic *dw)
{
	struct gmdialog *bd;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *table;
	GtkWidget *label;
	char tmp[64];

	if ( dw->boostint_dialog != NULL ) {
		return 0;
	}

	if ( dw->hdfile == NULL ) {
		return 0;
	}

	bd = malloc(sizeof(struct gmdialog));
	if ( bd == NULL ) return 0;
	dw->boostint_dialog = bd;

	bd->window = gtk_dialog_new_with_buttons("Intensity Boost",
					GTK_WINDOW(dw->window),
					GTK_DIALOG_DESTROY_WITH_PARENT,
					GTK_STOCK_CANCEL, GTK_RESPONSE_CLOSE,
					GTK_STOCK_OK, GTK_RESPONSE_OK, NULL);

	vbox = gtk_vbox_new(FALSE, 0);
	hbox = gtk_hbox_new(TRUE, 0);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(bd->window)->vbox),
			   GTK_WIDGET(hbox), FALSE, FALSE, 7);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(vbox), FALSE, FALSE, 5);

	table = gtk_table_new(3, 2, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(table), 5);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(table), FALSE, FALSE, 0);

	label = gtk_label_new("Boost Factor:");
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(table), GTK_WIDGET(label),
				  1, 2, 3, 4);

	bd->entry = gtk_entry_new();
	snprintf(tmp, 63, "%.2f", dw->boostint);
	gtk_entry_set_text(GTK_ENTRY(bd->entry), tmp);
	gtk_table_attach_defaults(GTK_TABLE(table), GTK_WIDGET(bd->entry),
				  2, 3, 3, 4);

	g_signal_connect(G_OBJECT(bd->entry), "activate",
			 G_CALLBACK(displaywindow_set_boostint_response_ac),
			 dw);
	g_signal_connect(G_OBJECT(bd->window), "response",
			 G_CALLBACK(displaywindow_set_boostint_response), dw);
	g_signal_connect(G_OBJECT(bd->window), "destroy",
			 G_CALLBACK(displaywindow_set_boostint_destroy), dw);
	gtk_window_set_resizable(GTK_WINDOW(bd->window), FALSE);
	gtk_widget_show_all(bd->window);
	gtk_widget_grab_focus(GTK_WIDGET(bd->entry));

	return 0;
}


static gint displaywindow_about(GtkWidget *widget, DWGeomatic *dw)
{
	GtkWidget *window;

	const gchar *authors[] = {
		"Thomas White <taw@physics.org>",
		NULL
	};

	window = gtk_about_dialog_new();
	gtk_window_set_transient_for(GTK_WINDOW(window),
				     GTK_WINDOW(dw->window));

	gtk_about_dialog_set_name(GTK_ABOUT_DIALOG(window), "geomatic");
	gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(window), PACKAGE_VERSION);
	gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG(window),
		"(c) 2006-2011 Thomas White <taw@physics.org> and others");
	gtk_about_dialog_set_comments(GTK_ABOUT_DIALOG(window),
		"Quick viewer for HDF files");
	gtk_about_dialog_set_license(GTK_ABOUT_DIALOG(window),
		"(c) 2006-2011 Thomas White <taw@physics.org>\n");
	gtk_about_dialog_set_website(GTK_ABOUT_DIALOG(window),
		"http://www.bitwiz.org.uk/");
	gtk_about_dialog_set_authors(GTK_ABOUT_DIALOG(window), authors);

	g_signal_connect(window, "response", G_CALLBACK(gtk_widget_destroy),
			 NULL);

	gtk_widget_show_all(window);

	return 0;
}


static gint displaywindow_set_colscale(GtkWidget *widget, DWGeomatic *dw)
{
	dw->show_col_scale = 1 - dw->show_col_scale;
	displaywindow_update(dw);
	return 0;
}


static void displaywindow_addui_callback(GtkUIManager *ui, GtkWidget *widget,
					 GtkContainer *container)
{
	gtk_box_pack_start(GTK_BOX(container), widget, FALSE, FALSE, 0);

	/* Enable overflow menu if this is a toolbar */
	if ( GTK_IS_TOOLBAR(widget) ) {
		gtk_toolbar_set_show_arrow(GTK_TOOLBAR(widget), TRUE);
	}
}


static gint displaywindow_setscale(GtkWidget *widget, GtkRadioAction *action,
                                   DWGeomatic *dw)
{
	switch ( gtk_radio_action_get_current_value(action) )
	{
		case 0 : dw->scale = SCALE_COLOUR; break;
		case 1 : dw->scale = SCALE_MONO; break;
		case 2 : dw->scale = SCALE_INVMONO; break;
	}
	displaywindow_update(dw);

	return 0;
}


static gint displaywindow_loadgeom_response(GtkWidget *d, gint response,
                                            DWGeomatic *dw)
{
	if ( response == GTK_RESPONSE_ACCEPT ) {

		char *file;
		struct detector *det;

		file = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(d));

		det = get_detector_geometry(file);
		g_free(file);

		if ( det == NULL ) {
			displaywindow_error(dw, "Invalid geometry file");
			return 1;
		}

		/* Validate geometry */
		if ( (1+det->max_fs != dw->image->width)
		  || (1+det->max_ss != dw->image->height) ) {

			displaywindow_error(dw,
			                  "Geometry does not match image size");
			return 1;

		} else {

			if ( dw->image->det != NULL ) {
				free_detector_geometry(dw->image->det);
			}
			dw->image->det = det;

			displaywindow_update(dw);
		}

	}

	gtk_widget_destroy(d);

	return 0;
}



static gint displaywindow_loadgeom(GtkWidget *widget, DWGeomatic *dw)
{
	GtkWidget *d;

	d = gtk_file_chooser_dialog_new("Load Geometry File",
	                                GTK_WINDOW(dw->window),
	                                GTK_FILE_CHOOSER_ACTION_SAVE,
	                                GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
	                                GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
	                                NULL);

	g_signal_connect(G_OBJECT(d), "response",
	                 G_CALLBACK(displaywindow_loadgeom_response), dw);

	gtk_widget_show_all(d);

	return 0;
}


static void displaywindow_addmenubar(DWGeomatic *dw, GtkWidget *vbox,
                                     int colscale)
{
	GError *error = NULL;
	GtkActionEntry entries[] = {

		{ "FileAction", NULL, "_File", NULL, NULL, NULL },
		{ "LoadGeomAction", GTK_STOCK_OPEN, "_Load Geometry", NULL,
			NULL, G_CALLBACK(displaywindow_loadgeom) },
		{ "CloseAction", GTK_STOCK_CLOSE, "_Close", NULL, NULL,
			G_CALLBACK(displaywindow_close) },

		{ "ViewAction", NULL, "_View", NULL, NULL, NULL },
		{ "ImagesAction", NULL, "Images", NULL, NULL, NULL },
		{ "BoostIntAction", NULL, "Boost Intensity...", "F5", NULL,
			G_CALLBACK(displaywindow_set_boostint) },

		{ "HelpAction", NULL, "_Help", NULL, NULL, NULL },
		{ "AboutAction", GTK_STOCK_ABOUT, "_About Geomatic...",
			NULL, NULL,
			G_CALLBACK(displaywindow_about) },

	};
	guint n_entries = G_N_ELEMENTS(entries);

	GtkToggleActionEntry toggles[] = {
		{ "ColScaleAction", NULL, "Show Colour Scale", NULL, NULL,
			G_CALLBACK(displaywindow_set_colscale), FALSE },
	};
	guint n_toggles = G_N_ELEMENTS(toggles);
	GtkRadioActionEntry radios[] = {
		{ "ColAction", NULL, "Colour", NULL, NULL,
			SCALE_COLOUR },
		{ "MonoAction", NULL, "Monochrome", NULL, NULL,
			SCALE_MONO },
		{ "InvMonoAction", NULL, "Inverse Monochrome", NULL, NULL,
			SCALE_INVMONO },
	};
	guint n_radios = G_N_ELEMENTS(radios);

	dw->action_group = gtk_action_group_new("geomatic");
	gtk_action_group_add_actions(dw->action_group, entries, n_entries, dw);
	gtk_action_group_add_toggle_actions(dw->action_group, toggles,
					    n_toggles, dw);
	gtk_action_group_add_radio_actions(dw->action_group, radios, n_radios,
	                                   colscale,
	                                   G_CALLBACK(displaywindow_setscale),
	                                   dw);

	dw->ui = gtk_ui_manager_new();
	gtk_ui_manager_insert_action_group(dw->ui, dw->action_group, 0);
	g_signal_connect(dw->ui, "add_widget",
			 G_CALLBACK(displaywindow_addui_callback), vbox);
	if ( gtk_ui_manager_add_ui_from_file(dw->ui,
	     DATADIR"/crystfel/geomatic.ui", &error) == 0 ) {
		fprintf(stderr, "Error loading message window menu bar: %s\n",
			error->message);
		return;
	}

	gtk_window_add_accel_group(GTK_WINDOW(dw->window),
				   gtk_ui_manager_get_accel_group(dw->ui));
	gtk_ui_manager_ensure_update(dw->ui);
}


struct newhdf {
	DWGeomatic *dw;
	char name[1024];
};

static gint displaywindow_newhdf(GtkMenuItem *item, struct newhdf *nh)
{
	hdfile_set_image(nh->dw->hdfile, nh->name);
	hdf5_read(nh->dw->hdfile, nh->dw->image, 0, 0.0);
	gtk_widget_set_size_request(GTK_WIDGET(nh->dw->drawingarea),
	                            nh->dw->image->width,
	                            nh->dw->image->height);
	displaywindow_update(nh->dw);
	return 0;
}


static GtkWidget *displaywindow_addhdfgroup(struct hdfile *hdfile,
                                            const char *group,
                                            DWGeomatic *dw, GSList **rgp)
{
	char **names;
	int *is_group;
	int *is_image;
	GtkWidget *ms;
	int n, i;

	if ( hdfile == NULL ) return NULL;

	names = hdfile_read_group(hdfile, &n, group, &is_group, &is_image);
	if ( n == 0 ) return NULL;

	ms = gtk_menu_new();

	for ( i=0; i<n; i++ ) {

		GtkWidget *item;
		GtkWidget *sub;

		if ( names[i] == NULL ) return NULL;

		if ( is_group[i] ) {

			item = gtk_menu_item_new_with_label(names[i]);

			sub = displaywindow_addhdfgroup(hdfile, names[i],
			                                dw, rgp);
			gtk_menu_item_set_submenu(GTK_MENU_ITEM(item), sub);

		} else if ( is_image[i] ) {

			struct newhdf *nh;

			item = gtk_radio_menu_item_new_with_label(*rgp,
			                                          names[i]);

			if ( *rgp == NULL ) {
				gtk_check_menu_item_set_active(
				               GTK_CHECK_MENU_ITEM(item), TRUE);
			} else {
				gtk_check_menu_item_set_active(
				              GTK_CHECK_MENU_ITEM(item), FALSE);
			}

			*rgp = gtk_radio_menu_item_get_group(
		                                     GTK_RADIO_MENU_ITEM(item));

			nh = malloc(sizeof(struct newhdf));
			if ( nh != NULL ) {
				strncpy(nh->name, names[i], 1023);
				nh->dw = dw;
				g_signal_connect(G_OBJECT(item), "activate",
			                  G_CALLBACK(displaywindow_newhdf), nh);
			}

		} else {

			char *tmp;

			item = gtk_menu_item_new_with_label(names[i]);

			tmp = hdfile_get_string_value(hdfile, names[i]);
			if ( tmp != NULL ) {

				GtkWidget *ss;
				GtkWidget *mss;

				mss = gtk_menu_new();
				ss = gtk_menu_item_new_with_label(tmp);
				gtk_widget_set_sensitive(ss, FALSE);
				gtk_menu_shell_append(GTK_MENU_SHELL(mss), ss);
				gtk_menu_item_set_submenu(GTK_MENU_ITEM(item),
				                          mss);

			}


		}

		gtk_menu_shell_append(GTK_MENU_SHELL(ms), item);

		free(names[i]);


	}

	free(is_group);
	free(is_image);

	return ms;
}


static GtkWidget *displaywindow_createhdfmenus(struct hdfile *hdfile,
                                               DWGeomatic *dw)
{
	GSList *rg = NULL;

	return displaywindow_addhdfgroup(hdfile, "/", dw, &rg);
}


static void displaywindow_update_menus(DWGeomatic *dw)
{
	GtkWidget *ms;
	GtkWidget *w;

	ms = displaywindow_createhdfmenus(dw->hdfile, dw);

	if ( ms == NULL ) {

		/* Too bad.  You'd better hope that /data/data exists... */
		ERROR("Couldn't get list of images in HDF file\n");
		w = gtk_ui_manager_get_widget(dw->ui,
					      "/ui/displaywindow/view/images");
		gtk_widget_set_sensitive(GTK_WIDGET(w), FALSE);

		/* Add a dummy menu so that the user knows what's going on */
		ms = gtk_menu_new();
		w = gtk_ui_manager_get_widget(dw->ui,
		                              "/ui/displaywindow/view/images");
		gtk_menu_item_set_submenu(GTK_MENU_ITEM(w), ms);

		return;

	}

	/* Make new menu be the submenu for File->Images */
	w = gtk_ui_manager_get_widget(dw->ui, "/ui/displaywindow/view/images");
	gtk_menu_item_set_submenu(GTK_MENU_ITEM(w), ms);

	gtk_widget_show_all(ms);
}


static gint displaywindow_release(GtkWidget *widget, GdkEventButton *event,
                                  DWGeomatic *dw)
{
	if ( (event->type == GDK_BUTTON_RELEASE) && (event->button == 1) ) {

		g_signal_handler_disconnect(GTK_OBJECT(dw->drawingarea),
		                            dw->motion_callback);
		dw->motion_callback = 0;

	}

	return 0;
}


static gint displaywindow_motion(GtkWidget *widget, GdkEventMotion *event,
                                 DWGeomatic *dw)
{
	double x, y;

	x = event->x - dw->motion_origx;
	y = event->y - dw->motion_origy;

	dw->pos_x += x;
	dw->pos_y += y;

	/* Schedule redraw */
	gtk_widget_queue_draw_area(dw->drawingarea, 0, 0,
	                           dw->width, dw->height);

	return 0;
}


static gint displaywindow_press(GtkWidget *widget, GdkEventButton *event,
                                DWGeomatic *dw)
{
	if ( dw->motion_callback != 0 ) {
		return 0;
	}

	if ( (event->type == GDK_BUTTON_PRESS) && (event->button == 1) ) {

		dw->motion_origx = event->x;
		dw->motion_origy = event->y;

		/* Connect motion callback */
		dw->motion_callback = g_signal_connect(
		                               GTK_OBJECT(dw->drawingarea),
		                               "motion-notify-event",
		                               G_CALLBACK(displaywindow_motion),
		                               dw);

	}

	return 0;

}


DWGeomatic *geomatic_open(const char *filename)
{
	DWGeomatic *dw;
	char *title;
	GtkWidget *vbox;
	GtkWidget *sw;
	int wr, hr;

	dw = calloc(1, sizeof(DWGeomatic));
	if ( dw == NULL ) return NULL;
	dw->pixbuf = NULL;
	dw->show_col_scale = 0;
	dw->col_scale = NULL;
	dw->boostint_dialog = NULL;
	dw->boostint = 1;
	dw->motion_callback = 0;
	dw->image = NULL;
	dw->scale = SCALE_COLOUR;
	dw->pos_x = 0.0;
	dw->pos_y = 0.0;
	dw->pos_z = 0.0;

	dw->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);

	if ( filename == NULL ) {
		title = strdup("No file - geomatic");
	} else {
		char *bn = safe_basename(filename);
		title = malloc(strlen(bn)+14);
		sprintf(title, "%s - geomatic", bn);
		free(bn);
	}
	gtk_window_set_title(GTK_WINDOW(dw->window), title);
	free(title);

	g_signal_connect(G_OBJECT(dw->window), "destroy",
			 G_CALLBACK(displaywindow_closed), dw);

	vbox = gtk_vbox_new(FALSE, 0);
	gtk_container_add(GTK_CONTAINER(dw->window), vbox);
	displaywindow_addmenubar(dw, vbox, dw->scale);

	/* Open the file, if any */
	if ( filename != NULL ) {

		dw->hdfile = hdfile_open(filename);
		if ( dw->hdfile == NULL ) {
			ERROR("Couldn't open file '%s'\n", filename);
			return NULL;
		} else if ( hdfile_set_first_image(dw->hdfile, "/") ) {
			ERROR("Couldn't select path\n");
			return NULL;
		} else {
			dw->image = calloc(1, sizeof(struct image));
			hdf5_read(dw->hdfile, dw->image, 0, 0.0);
		}

	} else {
		return NULL;
	}

	dw->drawingarea = gtk_drawing_area_new();
	gtk_widget_set_size_request(GTK_WIDGET(dw->drawingarea),
	                            dw->image->width,
	                            dw->image->height);
	sw = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
	                               GTK_POLICY_AUTOMATIC,
	                               GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(sw),
	                                      dw->drawingarea);
	gtk_box_pack_start(GTK_BOX(vbox), sw, TRUE, TRUE, 0);

	wr = dw->image->width;
	hr = dw->image->height;
	if ( wr > 640 ) wr = 640;
	if ( hr > 640 ) hr = 640;
	gtk_widget_set_size_request(GTK_WIDGET(dw->window), wr, hr);

	g_signal_connect(GTK_OBJECT(dw->drawingarea), "expose-event",
			 G_CALLBACK(displaywindow_expose), dw);

	gtk_window_set_resizable(GTK_WINDOW(dw->window), TRUE);
	gtk_widget_show_all(dw->window);

	dw->cell = load_cell_from_pdb("1JB0.pdb");

	dw->boostint = 1.0;
	displaywindow_update(dw);

	gtk_widget_add_events(GTK_WIDGET(dw->drawingarea),
	                      GDK_BUTTON_PRESS_MASK
	                      | GDK_BUTTON_RELEASE_MASK
	                      | GDK_BUTTON1_MOTION_MASK);
	g_object_set(G_OBJECT(dw->drawingarea), "can-focus", TRUE, NULL);

	g_signal_connect(GTK_OBJECT(dw->drawingarea), "button-press-event",
	                 G_CALLBACK(displaywindow_press), dw);
	g_signal_connect(GTK_OBJECT(dw->drawingarea), "button-release-event",
	                 G_CALLBACK(displaywindow_release), dw);

	if ( dw->hdfile != NULL ) displaywindow_update_menus(dw);

	return dw;
}
