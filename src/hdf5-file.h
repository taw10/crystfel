/*
 * hdf5.h
 *
 * Read/write HDF5 data files
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
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
extern void hdfile_close(struct hdfile *f);

#endif	/* HDF5_H */
