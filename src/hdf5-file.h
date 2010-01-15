/*
 * hdf5.h
 *
 * Read/write HDF5 data files
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef HDF5_H
#define HDF5_H

#include <stdint.h>

#include "image.h"

struct hdfile;


extern int hdf5_write(const char *filename, const int16_t *data,
                      int width, int height);

extern int hdf5_read(struct hdfile *f, struct image *image);

extern struct hdfile *hdfile_open(const char *filename);
extern int hdfile_set_image(struct hdfile *f, const char *path);
extern int hdfile_get_width(struct hdfile *f);
extern int hdfile_get_height(struct hdfile *f);
extern int16_t *hdfile_get_image_binned(struct hdfile *hdfile,
                                         int binning, int16_t *maxp);
extern int hdfile_get_unbinned_value(struct hdfile *f, int x, int y,
                                     int16_t *val);
extern char **hdfile_read_group(struct hdfile *f, int *n, const char *parent,
                                int **p_is_group, int **p_is_image);
extern int hdfile_set_first_image(struct hdfile *f, const char *group);
extern void hdfile_close(struct hdfile *f);

extern char *hdfile_get_string_value(struct hdfile *f, const char *name);
extern struct image *hdfile_get_image(struct hdfile *f);

#endif	/* HDF5_H */
