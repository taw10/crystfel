/*
 * beam-parameters.h
 *
 * Beam parameters
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef BEAM_PARAMETERS_H
#define BEAM_PARAMETERS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


struct beam_params
{
	double fluence;        /* photons per pulse */
	double beam_radius;    /* metres */
	double photon_energy;  /* eV per photon */

	double dqe;            /* Detector DQE (fraction) */
	double adu_per_photon; /* Detector "gain" */

	double water_radius;   /* metres */
};


extern struct beam_params *get_beam_parameters(const char *filename);


#endif	/* BEAM_PARAMETERS_H */
