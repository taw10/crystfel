/*
 * hdf5-file.h
 *
 * Read/write HDF5 data files
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
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
#include <hdf5.h>

#include "image.h"

struct hdfile;

struct copy_hdf5_field;


extern int hdf5_write(const char *filename, const void *data,
                      int width, int height, int type);

extern int hdf5_read(struct hdfile *f, struct image *image, int satcorr);

extern struct hdfile *hdfile_open(const char *filename);
extern int hdfile_set_image(struct hdfile *f, const char *path);
extern int16_t *hdfile_get_image_binned(struct hdfile *hdfile,
                                         int binning, int16_t *maxp);
extern char **hdfile_read_group(struct hdfile *f, int *n, const char *parent,
                                int **p_is_group, int **p_is_image);
extern int hdfile_set_first_image(struct hdfile *f, const char *group);
extern void hdfile_close(struct hdfile *f);

extern char *hdfile_get_string_value(struct hdfile *f, const char *name);
extern int get_peaks(struct image *image, struct hdfile *f, const char *p);
extern double get_value(struct hdfile *f, const char *name);

extern struct copy_hdf5_field *new_copy_hdf5_field_list(void);
extern void free_copy_hdf5_field_list(struct copy_hdf5_field *f);
extern void copy_hdf5_fields(struct hdfile *f,
                             const struct copy_hdf5_field *copyme, FILE *fh);
extern void add_copy_hdf5_field(struct copy_hdf5_field *copyme,
                                const char *name);


#endif	/* HDF5_H */
