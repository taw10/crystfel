/*
 * im-asapo.c
 *
 * ASAP::O data interface
 *
 * Copyright © 2021-2023 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2021-2023 Thomas White <taw@physics.org>
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
#include <sys/time.h>
#include <asapo/consumer_c.h>
#include <asapo/producer_c.h>

#include <image.h>
#include <utils.h>
#include <profile.h>

#include "im-asapo.h"

#include "datatemplate_priv.h"


struct im_asapo
{
	char *stream;
	AsapoConsumerHandle consumer;
	AsapoProducerHandle producer;
	AsapoStringHandle group_id;
	int wait_for_stream;
};


static void show_asapo_error(const char *msg, const AsapoErrorHandle err)
{
	char buf[1024];
	char tstr[256];
	time_t t;
	struct tm *tmp;

	t = time(NULL);
	tmp = localtime(&t);
	asapo_error_explain(err, buf, sizeof(buf));
	strftime(tstr, 256, "%d-%m-%y %T %z", tmp);
	ERROR("[%s] %s: %s\n", tstr, msg, buf);
}


static int create_producer(struct im_asapo *a, struct im_asapo_params *params)
{
	char *source;
	AsapoSourceCredentialsHandle cred;
	AsapoErrorHandle err = asapo_new_handle();

	if ( !params->write_output_stream ) {
		a->producer = NULL;
		return 0;
	}

	source = malloc(strlen(params->source)+6);
	if ( source == NULL ) return 1;

	strcpy(source, params->source);
	strcat(source, "_hits");

	cred = asapo_create_source_credentials(kProcessed,
	                                       "auto",        /* instance ID */
	                                       "indexamajig", /* pipeline step */
	                                       params->beamtime,
	                                       "",  /* beamline */
	                                       source,
	                                       params->token);
	STATUS("Writing hits-only output stream as data source '%s'\n", source);
	free(source);

	a->producer = asapo_create_producer(params->endpoint,
	                                    8,      /* Number of sender threads */
	                                    kTcp,
	                                    cred,
	                                    30000,  /* Timeout */
	                                    &err);

	if ( asapo_is_error(err) ) {
		show_asapo_error("Cannot create ASAP::O producer", err);
		asapo_free_handle(&cred);
		asapo_free_handle(&err);
		return 1;
	}

	asapo_producer_set_log_level(a->producer, Debug);

	asapo_free_handle(&err);
	return 0;
}


struct im_asapo *im_asapo_connect(struct im_asapo_params *params)
{
	struct im_asapo *a;
	AsapoSourceCredentialsHandle cred;
	AsapoErrorHandle err = asapo_new_handle();

	if ( params->endpoint == NULL ) {
		ERROR("ASAP::O endpoint not specified.\n");
		return NULL;
	}
	if ( params->beamtime == NULL ) {
		ERROR("ASAP::O beamtime not specified.\n");
		return NULL;
	}
	if ( params->group_id == NULL ) {
		ERROR("ASAP::O consumer group ID not specified.\n");
		return NULL;
	}
	if ( params->source == NULL ) {
		ERROR("ASAP::O data source not specified.\n");
		return NULL;
	}
	if ( params->stream == NULL ) {
		ERROR("ASAP::O stream not specified.\n");
		return NULL;
	}

	a = malloc(sizeof(struct im_asapo));
	if ( a == NULL ) return NULL;

	cred = asapo_create_source_credentials(kProcessed,
	                                       "auto",        /* instance ID */
	                                       "indexamajig", /* pipeline step */
	                                       params->beamtime,
	                                       "",  /* beamline */
	                                       params->source,
	                                       params->token);
	a->consumer = asapo_create_consumer(params->endpoint, "auto", 0, cred, &err);
	if ( asapo_is_error(err) ) {
		show_asapo_error("Cannot create ASAP::O consumer", err);
		asapo_free_handle(&cred);
		free(a);
		return NULL;
	}

	if ( create_producer(a, params) ) return NULL;

	a->stream = strdup(params->stream);
	asapo_consumer_set_timeout(a->consumer, 3000);
	a->group_id = asapo_string_from_c_str(params->group_id);
	a->wait_for_stream = params->wait_for_stream;

	asapo_free_handle(&cred);

	return a;
}


static int stream_empty(struct im_asapo *a)
{
	AsapoErrorHandle err;

	err = asapo_new_handle();
	int64_t size = asapo_consumer_get_current_size(a->consumer, a->stream,
	                                               &err);

	if ( asapo_is_error(err) ) {
		show_asapo_error("Couldn't get stream size", err);
		asapo_free_handle(&err);
		return 0;
	}

	return ( size == 0 );
}


void *im_asapo_fetch(struct im_asapo *a, size_t *pdata_size,
                     char **pmeta, char **pfilename, char **pevent,
                     int *pfinished, int *pmessageid)
{
	void *data_copy;
	AsapoMessageMetaHandle meta;
	AsapoMessageDataHandle data;
	AsapoErrorHandle err;
	uint64_t msg_size;

	*pfinished = 0;

	profile_start("create-handles");
	err = asapo_new_handle();
	meta = asapo_new_handle();
	data = asapo_new_handle();
	profile_end("create-handles");

	profile_start("asapo-get-next");
	asapo_consumer_get_next(a->consumer, a->group_id, &meta, &data,
	                        a->stream, &err);
	profile_end("asapo-get-next");
	if ( asapo_error_get_type(err) == kEndOfStream ) {
		asapo_free_handle(&err);
		asapo_free_handle(&meta);
		asapo_free_handle(&data);
		if ( stream_empty(a) && a->wait_for_stream ) {
			*pfinished = 0;
		} else {
			*pfinished = 1;
		}
		return NULL;
	}

	if ( asapo_is_error(err) ) {
		show_asapo_error("Couldn't get next ASAP::O record", err);
		asapo_free_handle(&err);
		asapo_free_handle(&meta);
		asapo_free_handle(&data);
		return NULL;
	}

	profile_start("get-size");
	msg_size = asapo_message_meta_get_size(meta);
	profile_end("get-size");

	profile_start("malloc-copy");
	data_copy = malloc(msg_size);
	if ( data_copy == NULL ) {
		ERROR("Failed to copy data block.\n");
		asapo_free_handle(&err);
		asapo_free_handle(&meta);
		asapo_free_handle(&data);
		return NULL;
	}
	memcpy(data_copy, asapo_message_data_get_as_chars(data), msg_size);
	profile_end("malloc-copy");

	profile_start("copy-meta");
	*pmeta = strdup(asapo_message_meta_get_metadata(meta));
	*pfilename = strdup(asapo_message_meta_get_name(meta));
	*pevent = strdup("//");
	*pmessageid = asapo_message_meta_get_id(meta);
	profile_end("copy-meta");

	asapo_free_handle(&err);
	asapo_free_handle(&meta);
	asapo_free_handle(&data);

	*pdata_size = msg_size;
	return data_copy;
}


static void send_callback(void *a, AsapoRequestCallbackPayloadHandle payload,
                          AsapoErrorHandle err)
{
	if ( asapo_is_error(err) ) {
		show_asapo_error("ASAP::O send error", err);
	}
}


static void send_real(struct im_asapo *a, struct image *image)
{
	AsapoMessageHeaderHandle header;
	AsapoErrorHandle err;
	char filename[1024];

	snprintf(filename, 1024, "processed/%s_hits/%s-%i.data",
	         a->stream, a->stream, image->serial);

        header = asapo_create_message_header(image->serial,
                                             image->data_block_size,
                                             filename,
                                             image->meta_data,
                                             0,   /* Dataset substream */
                                             0,
                                             0);  /* Auto ID */

	struct timeval tv;
	gettimeofday(&tv,NULL);
	STATUS("sent %s at %lli . %lli\n", filename, tv.tv_sec, tv.tv_usec);

	err = asapo_new_handle();
	asapo_producer_send(a->producer, header, image->data_block,
	                    kTransferData | kStoreInDatabase, a->stream,
	                    send_callback, &err);
	if ( asapo_is_error(err) ) {
		show_asapo_error("Couldn't send ASAP::O message", err);
		asapo_free_handle(&header);
		asapo_free_handle(&err);
		return;
	}

	asapo_free_handle(&header);
	asapo_free_handle(&err);
}


static void send_placeholder(struct im_asapo *a, struct image *image)
{
	AsapoMessageHeaderHandle header;
	AsapoErrorHandle err;
	char filename[1024];

	snprintf(filename, 1024, "processed/%s_hits/%s-%i.placeholder",
	         a->stream, a->stream, image->serial);

	struct timeval tv;
	gettimeofday(&tv,NULL);
	STATUS("sent %s at %lli . %lli\n", filename, tv.tv_sec, tv.tv_usec);

        header = asapo_create_message_header(image->serial,
                                             8,   /* strlen("SKIPPED"+\0) */
                                             filename,
                                             image->meta_data,
                                             0,   /* Dataset substream */
                                             0,
                                             0);  /* Auto ID */

	err = asapo_new_handle();
	asapo_producer_send(a->producer, header, "SKIPPED",
	                    kDefaultIngestMode, a->stream,
	                    send_callback, &err);
	if ( asapo_is_error(err) ) {
		show_asapo_error("Couldn't send ASAP::O message", err);
		asapo_free_handle(&header);
		asapo_free_handle(&err);
		return;
	}

	asapo_free_handle(&header);
	asapo_free_handle(&err);
}


/* Send the image to the output ASAP::O stream, if it's a hit.  Otherwise,
 * send a placeholder */
void im_asapo_send(struct im_asapo *a, struct image *image, int hit)
{
	if ( a == NULL ) return;
	if ( a->producer == NULL ) return;
	profile_start("asapo-send");
	if ( hit ) {
		send_real(a, image);
	} else {
		send_placeholder(a, image);
	}
	profile_end("asapo-send");
}


void im_asapo_shutdown(struct im_asapo *a)
{
	if ( a == NULL ) return;
	asapo_free_handle(&a->consumer);
	asapo_free_handle(&a->group_id);
	free(a);
}
