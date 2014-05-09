/*
 * hdf5-file.h
 *
 * Read/write HDF5 data files
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2012 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani

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

#ifndef HDF5_H
#define HDF5_H

struct event_list;

#include <stdint.h>
#include <hdf5.h>
#include "image.h"
#include "events.h"

struct hdfile;
struct copy_hdf5_field;

#include "image.h"
#include "events.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int hdf5_write(const char *filename, const void *data,
                      int width, int height, int type);

extern int hdf5_write_image(const char *filename, struct image *image,
                            char *element);

extern int hdf5_read(struct hdfile *f, struct image *image,
                     const char* element, int satcorr);

extern int hdf5_read2(struct hdfile *f, struct image *image,
			   struct event *ev, int satcorr);

extern int check_path_existence(hid_t fh, const char *path);

extern struct hdfile *hdfile_open(const char *filename);
int hdfile_set_image(struct hdfile *f, const char *path,
                     struct panel *p);
extern int16_t *hdfile_get_image_binned(struct hdfile *hdfile,
                                         int binning, int16_t *maxp);
extern char **hdfile_read_group(struct hdfile *f, int *n, const char *parent,
                                int **p_is_group, int **p_is_image);
extern int hdfile_set_first_image(struct hdfile *f, const char *group);
extern void hdfile_close(struct hdfile *f);

extern int hdfile_is_scalar(struct hdfile *f, const char *name, int verbose);
char *hdfile_get_string_value(struct hdfile *f, const char *name,
                              struct event* ev);
extern int get_peaks(struct image *image, struct hdfile *f, const char *p);
extern double get_value(struct hdfile *f, const char *name);

extern double get_ev_based_value(struct hdfile *f, const char *name,
                                 struct event *ev);

extern struct copy_hdf5_field *new_copy_hdf5_field_list(void);
extern void free_copy_hdf5_field_list(struct copy_hdf5_field *f);

extern void copy_hdf5_fields(struct hdfile *f,
                             const struct copy_hdf5_field *copyme,
                             FILE *fh, struct event *ev);
extern void add_copy_hdf5_field(struct copy_hdf5_field *copyme,
                                const char *name);
extern struct event_list *fill_event_list(struct hdfile* hdfile,
                                          struct detector* det);

#ifdef __cplusplus
}
#endif

#endif	/* HDF5_H */
