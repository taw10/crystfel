/*
 * im-asapo.h
 *
 * ASAP::O data interface
 *
 * Copyright © 2021-2022 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2021-2022 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#ifndef CRYSTFEL_ASAPO_H
#define CRYSTFEL_ASAPO_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

struct im_asapo_params
{
	char *endpoint;
	char *token;
	char *beamtime;
	char *group_id;
	char *source;
	char *stream;
	int wait_for_stream;
	char *output_stream;
};

struct im_asapo;

#if defined(HAVE_ASAPO)

extern struct im_asapo *im_asapo_connect(struct im_asapo_params *params);

extern void im_asapo_shutdown(struct im_asapo *a);

extern void *im_asapo_fetch(struct im_asapo *a, size_t *pdata_size,
                            char **pmeta, char **pfilename, char **pevent,
                            int *pfinished, int *pmessageid);

extern void im_asapo_send(struct im_asapo *a, struct image *image, int hit);

#else /* defined(HAVE_ASAPO) */

static UNUSED struct im_asapo *im_asapo_connect(struct im_asapo_params *params)
{
	ERROR("This installation of CrystFEL was compiled without ASAP::O support.\n");
	return NULL;
}

static UNUSED void im_asapo_shutdown(struct im_asapo *a)
{
}

static UNUSED void *im_asapo_fetch(struct im_asapo *a, size_t *psize,
                                   char **pmeta, char **pfilename, char **pevent,
                                   int *pfinished, int *pmessageid)
{
	*psize = 0;
	*pmeta = NULL;
	*pfilename = NULL;
	*pevent = NULL;
	*pfinished = 1;
	*pmessageid = 0;
	return NULL;
}

static UNUSED void im_asapo_send(struct im_asapo *a, struct image *image, int hit)
{
}

#endif /* defined(HAVE_ASAPO) */

#endif /* CRYSTFEL_ASAPO_H */
