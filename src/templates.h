/*
 * templates.h
 *
 * Handle templates
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef TEMPLATES_H
#define TEMPLATES_H

#include "cell.h"
#include "image.h"

/* An opaque type representing a list of templates */
typedef struct _templatelist TemplateList;

extern TemplateList *generate_templates(UnitCell *cell, struct image params);
extern int try_templates(struct image *image, TemplateList *list);

#endif	/* TEMPLAETS_H */
