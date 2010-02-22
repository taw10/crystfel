/*
 * dirax.h
 *
 * Invoke the DirAx auto-indexing program
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef DIRAX_H
#define DIRAX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "utils.h"

#if HAVE_GLIB

extern void run_dirax(struct image *image);

#else

static void run_dirax(struct image *image)
{
	ERROR("Can't run DirAx without GLib.\n");
}

#endif

#endif	/* DIRAX_H */
