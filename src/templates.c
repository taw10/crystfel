/*
 * templates.c
 *
 * Handle templates
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * template_index - Indexing diffraction patterns by template matching
 *
 */


#define _GNU_SOURCE 1
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "templates.h"
#include "cell.h"
#include "image.h"
#include "utils.h"
#include "relrod.h"


struct _templatelist
{
	int		n_templates;
	struct template	*templates;
};


struct template_feature
{
	float		x;
	float		y;
	struct template_feature *next;
};


struct template
{
	float		omega;
	float		tilt;
	struct template_feature *features;
};


static int template_add(TemplateList *list, struct template *template)
{
	if ( list->templates ) {
		list->templates = realloc(list->templates,
		                 (list->n_templates+1)*sizeof(struct template));
	} else {
		assert(list->n_templates == 0);
		list->templates = malloc(sizeof(struct template));
	}

	/* Copy the data */
	list->templates[list->n_templates] = *template;

	list->n_templates++;

	return list->n_templates - 1;
}


static TemplateList *template_list_new()
{
	TemplateList *list;

	list = malloc(sizeof(TemplateList));

	list->n_templates = 0;
	list->templates = NULL;

	return list;
}


TemplateList *generate_templates(UnitCell *cell, struct image params)
{
	TemplateList *list;
	double omega, tilt;

	list = template_list_new();

	omega = deg2rad(40.0);

	//for ( omega=deg2rad(-180); omega<deg2rad(180); omega+=deg2rad(1) ) {

		params.omega = omega;

		for ( tilt=0; tilt<deg2rad(180); tilt+=deg2rad(1) ) {

			struct template t;
			struct template_feature *tfc;
			int nrefl, i;

			t.omega = omega;
			t.tilt = tilt;
			t.features = malloc(sizeof(struct template_feature));
			t.features->next = NULL;
			tfc = t.features;

			params.tilt = tilt;
			get_reflections(&params, cell, 0.01e9);

			nrefl = image_feature_count(params.rflist);
			for ( i=0; i<nrefl; i++ ) {

				struct imagefeature *f;

				f = image_get_feature(params.rflist, i);

				tfc->x = f->x;
				tfc->y = f->y;
				tfc->next = malloc(sizeof(struct template_feature));
				tfc = tfc->next;

			}

			template_add(list, &t);

			image_feature_list_free(params.rflist);

			printf("Generating templates... %+5.2f %+5.2f\r",
			       rad2deg(omega), rad2deg(tilt));

		}
	//}

	return list;
}


int try_template(struct image *image, struct template template)
{
	int fit = 0;
	struct template_feature *f;

	f = template.features;
	while ( f != NULL ) {

		int x, y;

		x = f->x;
		y = f->y;	/* Discards digits after the decimal point */

		fit += image->data[y*image->width+x];

		f = f->next;

	}

	return fit;
}


int try_templates(struct image *image, TemplateList *list)
{
	int i;
	int fit_max = 0;
	int i_max = 0;

	for ( i=0; i<list->n_templates; i++ ) {

		int fit;

		fit = try_template(image, list->templates[i]);
		if ( fit > fit_max ) {
			fit_max = fit;
			i_max = i;
		}

	}
	image->omega = list->templates[i_max].omega;
	image->tilt = list->templates[i_max].tilt;

	return 0;
}
