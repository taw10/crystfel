/*
 * hdf5-file.h
 *
 * Read/write HDF5 data files
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2017 Thomas White <taw@physics.org>
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

struct copy_hdf5_field;

#include "image.h"
#include "events.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \file hdf5-file.h
 * HDF5 utility functions
 */

extern int hdf5_write_image(const char *filename, const struct image *image,
                            char *element);

extern int check_path_existence(hid_t fh, const char *path);

extern struct copy_hdf5_field *new_copy_hdf5_field_list(void);
extern void free_copy_hdf5_field_list(struct copy_hdf5_field *f);
extern void copy_hdf5_fields(struct hdfile *f,
                             const struct copy_hdf5_field *copyme,
                             FILE *fh, struct event *ev);
extern void add_copy_hdf5_field(struct copy_hdf5_field *copyme,
                                const char *name);
extern int hdfile_get_value(struct hdfile *f, const char *name,
                            struct event *ev, void *val, hid_t memtype);
extern int hdfile_is_scalar(struct hdfile *f, const char *name, int verbose);
extern char *hdfile_get_string_value(struct hdfile *f, const char *name,
                                     struct event *ev);

#ifdef __cplusplus
}
#endif

#endif	/* HDF5_H */
