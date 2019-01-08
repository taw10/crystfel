/*
 * zmq.h
 *
 * ZMQ data interface
 *
 * Copyright © 2017-2019 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2018-2019 Thomas White <taw@physics.org>
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

#include <msgpack.h>

#include "image.h"

extern struct im_zmq *im_zmq_connect(const char *zmq_address);

extern msgpack_object *im_zmq_fetch(struct im_zmq *z);

extern void im_zmq_clean(struct im_zmq *z);

extern void im_zmq_shutdown(struct im_zmq *z);

extern int get_peaks_onda(msgpack_object *obj, struct image *image,
                          int half_pixel_shift);

extern int obj_read(msgpack_object *obj, struct image *image);


#endif /* CRYSTFEL_ZMQ_H */
