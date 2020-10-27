/*
 * zmq.c
 *
 * ZMQ data interface
 *
 * Copyright © 2017-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2018-2020 Thomas White <taw@physics.org>
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>
#include <zmq.h>
#include <msgpack.h>

#include <image.h>
#include <utils.h>

#include "im-zmq.h"

#include "datatemplate_priv.h"


struct im_zmq
{
	void *ctx;
	void *socket;
	zmq_msg_t msg;
	msgpack_unpacked unpacked;
	int unpacked_set;
};


struct im_zmq *im_zmq_connect(const char *zmq_address)
{
	struct im_zmq *z;

	z = malloc(sizeof(struct im_zmq));
	if ( z == NULL ) return NULL;

	z->unpacked_set = 0;

	z->ctx = zmq_ctx_new();
	if ( z->ctx == NULL ) return NULL;

	z->socket = zmq_socket(z->ctx, ZMQ_REQ);
	if ( z->socket == NULL ) return NULL;

	STATUS("Connecting to ZMQ at '%s'\n", zmq_address);
	if ( zmq_connect(z->socket, zmq_address) == -1 ) {
		ERROR("ZMQ connection failed: %s\n", zmq_strerror(errno));
		return NULL;
	}
	STATUS("ZMQ connected.\n");

	return z;
}


msgpack_object *im_zmq_fetch(struct im_zmq *z)
{
	int msg_size;
	int r;

	if ( zmq_send(z->socket, "m", 1, 0) == -1 ) {
		ERROR("ZMQ message send failed: %s\n", zmq_strerror(errno));
		return NULL;
	}

	zmq_msg_init(&z->msg);
	msg_size = zmq_msg_recv(&z->msg, z->socket, 0);
	if ( msg_size == -1 ) {
		ERROR("ZMQ recieve failed: %s\n", zmq_strerror(errno));
		zmq_msg_close(&z->msg);
		return NULL;
	}

	msgpack_unpacked_init(&z->unpacked);
	r = msgpack_unpack_next(&z->unpacked, zmq_msg_data(&z->msg),
	                        msg_size, NULL);
	if ( r != MSGPACK_UNPACK_SUCCESS ) {
		ERROR("Msgpack unpack failed: %i\n", r);
		zmq_msg_close(&z->msg);
		return NULL;
	}
	z->unpacked_set = 1;

	return &z->unpacked.data;
}


/* Clean structures ready for next frame */
void im_zmq_clean(struct im_zmq *z)
{
	if ( z->unpacked_set ) {
		msgpack_unpacked_destroy(&z->unpacked);
		zmq_msg_close(&z->msg);
		z->unpacked_set = 0;
	}
}


void im_zmq_shutdown(struct im_zmq *z)
{
	if ( z == NULL ) return;

	zmq_msg_close(&z->msg);
	zmq_close(z->socket);
	zmq_ctx_destroy(z->ctx);
}
