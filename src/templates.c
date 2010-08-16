/*
 * templates.c
 *
 * Indexing by template matching
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "index.h"
#include "index-priv.h"
#include "symmetry.h"
#include "utils.h"


/* Private data for template indexing */
struct _indexingprivate_template
{
	struct _indexingprivate base;
};


/* Generate templates for the given cell using a representative image */
IndexingPrivate *generate_templates(UnitCell *cell, const char *filename)
{
	struct _indexingprivate_template *priv;
	const char *holo;
	double omega_max, phi_max;

	priv = calloc(1, sizeof(struct _indexingprivate_template));
	priv->base.indm = INDEXING_TEMPLATE;

	/* We can only distinguish orientations within the holohedral cell */
	holo = get_holohedral(cell_get_pointgroup(cell));
	STATUS("%s\n", holo);

	/* These define the orientation in space */
	if ( is_polyhedral(holo) ) {
		ERROR("WARNING: Holohedral point group is polyhedral.\n");
		ERROR("This means I can't properly determine the orientation");
		ERROR(" ranges for template matching.  Expect trouble.\n");
	}
	omega_max = 2.0*M_PI / rotational_order(holo);
	if ( has_bisecting_mirror_or_diad(holo) ) omega_max /= 2.0;
	phi_max = M_PI;
	if ( has_perpendicular_mirror(holo) ) phi_max /= 2.0;

	/* One more axis would define the rotation in the plane of the image */

	STATUS("Orientation ranges: %5.0f -> %5.0f, %5.0f -> %5.0f deg.\n",
	       0.0, rad2deg(omega_max), 0.0, rad2deg(phi_max));

	return (struct _indexingprivate *)priv;
}


void match_templates(struct image *image, IndexingPrivate *ipriv)
{
}
