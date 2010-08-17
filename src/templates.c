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
#include "geometry.h"
#include "hdf5-file.h"


/* Private data for template indexing */
struct _indexingprivate_template
{
	struct _indexingprivate base;
	int n_templates;
	struct template *templates;
};


struct template {
	double omega;
	double phi;
	int n;
	struct reflhit spots;  /* Made into an array by Magic */
};


UnitCell *rotate_cell(UnitCell *in, double omega, double phi)
{
	UnitCell *out;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double xnew, ynew, znew;

	cell_get_reciprocal(in, &asx, &asy, &asz, &bsx, &bsy,
	                        &bsz, &csx, &csy, &csz);

	/* Rotate by "omega" about +z (parallel to c* and c unless triclinic) */
	xnew = asx*cos(omega) + asy*sin(omega);
	ynew = -asx*sin(omega) + asy*cos(omega);
	znew = asz;
	asx = xnew;  asy = ynew;  asz = znew;
	xnew = bsx*cos(omega) + bsy*sin(omega);
	ynew = -bsx*sin(omega) + bsy*cos(omega);
	znew = bsz;
	bsx = xnew;  bsy = ynew;  bsz = znew;
	xnew = csx*cos(omega) + csy*sin(omega);
	ynew = -csx*sin(omega) + csy*cos(omega);
	znew = csz;
	csx = xnew;  csy = ynew;  csz = znew;

	/* Rotate by "phi" about +x (not parallel to anything specific) */
	xnew = asx;
	ynew = asy*cos(phi) + asz*sin(phi);
	znew = -asy*sin(phi) + asz*cos(phi);
	asx = xnew;  asy = ynew;  asz = znew;
	xnew = bsx;
	ynew = bsy*cos(phi) + bsz*sin(phi);
	znew = -bsy*sin(phi) + bsz*cos(phi);
	bsx = xnew;  bsy = ynew;  bsz = znew;
	xnew = csx;
	ynew = csy*cos(phi) + csz*sin(phi);
	znew = -csy*sin(phi) + csz*cos(phi);
	csx = xnew;  csy = ynew;  csz = znew;

	out = cell_new_from_cell(in);
	cell_set_reciprocal(out, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);

	return out;
}


/* Generate templates for the given cell using a representative image */
IndexingPrivate *generate_templates(UnitCell *cell, const char *filename,
                                    struct detector *det)
{
	struct _indexingprivate_template *priv;
	const char *holo;
	double omega_max, phi_max;
	int n_templates;
	const double omega_step = deg2rad(0.5);
	const double phi_step = deg2rad(0.5);
	double omega, phi;
	struct image image;
	struct hdfile *hdfile;

	hdfile = hdfile_open(filename);
	if ( hdfile == NULL ) {
		return NULL;
	} else if ( hdfile_set_image(hdfile, "/data/data0") ) {
		ERROR("Couldn't select path\n");
		return NULL;
	}
	hdf5_read(hdfile, &image, 0);
	hdfile_close(hdfile);
	image.det = det;

	priv = calloc(1, sizeof(struct _indexingprivate_template));
	priv->base.indm = INDEXING_TEMPLATE;

	/* We can only distinguish orientations within the holohedral cell */
	holo = get_holohedral(cell_get_pointgroup(cell));

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

	STATUS("Orientation ranges in %s: %.0f-%.0f, %.0f-%.0f deg.\n",
	       holo, 0.0, rad2deg(omega_max), 0.0, rad2deg(phi_max));

	n_templates = (omega_max * phi_max)/(omega_step * phi_step);
	STATUS("%i templates to be calculated.\n", n_templates);

	for ( omega = 0.0; omega < omega_max; omega += omega_step ) {
	for ( phi = 0.0; phi < phi_max; phi += phi_step ) {

		int n;
		struct reflhit *hits;
		UnitCell *cell_rot;

		cell_rot = rotate_cell(cell, omega, phi);

		hits = find_intersections(&image, cell_rot, 5.0e-3,
		                          3.0/100.0, &n, 0);
		if ( hits == NULL ) {
			ERROR("Template calculation failed.\n");
			return NULL;
		}

		free(cell_rot);

	}
	progress_bar(omega*1000.0, (omega_max-omega_step)*1000.0,
	             "Generating templates");
	}

	priv->n_templates = n_templates;

	return (struct _indexingprivate *)priv;
}


void match_templates(struct image *image, IndexingPrivate *ipriv)
{
}
