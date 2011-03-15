/*
 * reflist-utils.h
 *
 * Utilities to complement the core reflist.c
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef REFLIST_UTILS_H
#define REFLIST_UTILS_H


#include "reflist.h"
#include "cell.h"


extern void write_reflections_to_file(FILE *fh, RefList *list, UnitCell *cell);

extern int write_reflist(const char *filename, RefList *list,
                             UnitCell *cell);


#endif	/* REFLIST_UTILS_H */
