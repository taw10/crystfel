/*
 * zmq.h
 *
 * ZMQ data interface
 *
 * Copyright © 2017-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2018-2021 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
 *   2017      Stijn de Graaf
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


#ifndef CRYSTFEL_ZMQ_H
#define CRYSTFEL_ZMQ_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if defined(HAVE_ZMQ)

extern struct im_zmq *im_zmq_connect(const char *zmq_address,
                                     char **subscriptions,
                                     int n_subscriptions);
extern void im_zmq_shutdown(struct im_zmq *z);
extern void *im_zmq_fetch(struct im_zmq *z, size_t *pdata_size);

#else /* defined(HAVE_ZMQ) */

static UNUSED struct im_zmq *im_zmq_connect(const char *zmq_address,
                                            char *zmq_subscriptions,
                                            int n_subscriptions) { return NULL; }
static UNUSED void im_zmq_shutdown(struct im_zmq *z) { }
static UNUSED void *im_zmq_fetch(struct im_zmq *z, size_t *psize) { *psize = 0; return NULL; }

#endif /* defined(HAVE_ZMQ) */

#endif /* CRYSTFEL_ZMQ_H */
