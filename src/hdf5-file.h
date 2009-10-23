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

extern int hdf5_write(const char *filename, const uint16_t *data,
                      int width, int height);

extern int hdf5_read(struct image *image, const char *filename);

#endif	/* HDF5_H */
