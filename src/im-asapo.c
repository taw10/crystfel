/*
 * im-asapo.c
 *
 * ASAP::O data interface
 *
 * Copyright © 2021 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2021 Thomas White <taw@physics.org>
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
#include <asapo/consumer_c.h>

#include <image.h>
#include <utils.h>

#include "im-asapo.h"

#include "datatemplate_priv.h"


struct im_asapo
{
	char *stream;
	AsapoConsumerHandle consumer;
	AsapoStringHandle group_id;
};


static void show_asapo_error(const char *msg, const AsapoErrorHandle err)
{
	char buf[1024];
	asapo_error_explain(err, buf, sizeof(buf));
	ERROR("%s: %s\n", msg, buf);
}


char *im_asapo_make_unique_group_id(const char *endpoint,
                                    const char *token)
{
	AsapoConsumerHandle consumer;
	AsapoSourceCredentialsHandle cred;
	AsapoStringHandle group_id;
	AsapoErrorHandle err = asapo_new_handle();

	cred = asapo_create_source_credentials(kProcessed, "", "", "", token);
	consumer = asapo_create_consumer(endpoint, "", 0, cred, &err);
	asapo_free_handle(&cred);
	if ( asapo_is_error(err) ) {
		show_asapo_error("Cannot create temporary ASAP::O consumer", err);
		asapo_free_handle(&consumer);
		return NULL;
	}

	group_id = asapo_consumer_generate_new_group_id(consumer, &err);
	asapo_free_handle(&consumer);
	if ( asapo_is_error(err) ) {
		show_asapo_error("Cannot create ASAP::O group ID", err);
		return NULL;
	}

	return strdup(asapo_string_c_str(group_id));
}


struct im_asapo *im_asapo_connect(const char *endpoint,
                                  const char *token,
                                  const char *beamtime,
                                  const char *group_id,
                                  const char *data_source)
{
	struct im_asapo *a;
	AsapoSourceCredentialsHandle cred;
	AsapoErrorHandle err = asapo_new_handle();

	a = malloc(sizeof(struct im_asapo));
	if ( a == NULL ) return NULL;

	cred = asapo_create_source_credentials(kProcessed, beamtime, "",
	                                       data_source, token);
	a->consumer = asapo_create_consumer(endpoint, "", 0, cred, &err);
	asapo_free_handle(&cred);
	if ( asapo_is_error(err) ) {
		show_asapo_error("Cannot create ASAP::O consumer", err);
		free(a);
		return NULL;
	}

	asapo_consumer_set_timeout(a->consumer, 1000);

	a->group_id = asapo_string_from_c_str(group_id);
	a->stream = NULL;

	return a;
}


static int select_last_stream(struct im_asapo *a)
{
	AsapoStreamInfosHandle si;
	size_t n;
	int i;
	AsapoStreamInfoHandle st;
	AsapoErrorHandle err = asapo_new_handle();

	si = asapo_consumer_get_stream_list(a->consumer, "",
	                                    kAllStreams, &err);

	if ( asapo_is_error(err) ) {
		show_asapo_error("Couldn't get ASAP::O stream list", err);
		asapo_free_handle(&err);
		return 1;
	}

	STATUS("Streams available at start:\n");
	n = asapo_stream_infos_get_size(si);
	for ( i=0; i<n; i++ ) {
		AsapoStreamInfoHandle st = asapo_stream_infos_get_item(si, i);
		STATUS("Stream %i: %s\n", i, asapo_stream_info_get_name(st));
		asapo_free_handle(&st);
	}
	STATUS("End of stream list\n");

	st = asapo_stream_infos_get_item(si, n-1);
	a->stream = strdup(asapo_stream_info_get_name(st));
	asapo_free_handle(&st);
	STATUS("Starting with the last stream: %s\n", a->stream);

	asapo_free_handle(&si);
	asapo_free_handle(&err);
	return 0;
}


static int select_next_stream(struct im_asapo *a)
{
	AsapoStreamInfosHandle si;
	AsapoStreamInfoHandle st;
	AsapoErrorHandle err = asapo_new_handle();
	const char *next_stream;

	si = asapo_consumer_get_stream_list(a->consumer, a->stream,
	                                    kAllStreams, &err);

	if ( asapo_is_error(err) ) {
		show_asapo_error("Couldn't get ASAP::O stream list", err);
		asapo_free_handle(&si);
		asapo_free_handle(&err);
		return 1;
	}

	asapo_free_handle(&err);

	/* Stream list includes the current stream, so we need at least
	 * two entries */
	if ( asapo_stream_infos_get_size(si) < 2 ) {
		STATUS("No newer stream.  Waiting for new data...\n");
		asapo_free_handle(&si);
		return 0;
	}

	/* Stream list includes the current stream, so look at the second one */
	st = asapo_stream_infos_get_item(si, 1);
	next_stream = asapo_stream_info_get_name(st);
	free(a->stream);
	a->stream = strdup(next_stream);
	STATUS("Selecting next stream: %s\n", a->stream);
	asapo_free_handle(&st);
	asapo_free_handle(&si);

	return 0;
}


static void skip_to_stream_end(struct im_asapo *a)
{
	int64_t size;
	AsapoErrorHandle err = asapo_new_handle();

	size = asapo_consumer_get_current_size(a->consumer, a->stream, &err);
	if ( asapo_is_error(err) ) {
		show_asapo_error("Failed to get length of stream", err);
	} else {

		AsapoErrorHandle err = asapo_new_handle();

		asapo_consumer_set_last_read_marker(a->consumer,
		                                    a->group_id, size,
		                                    a->stream, &err);
		if ( asapo_is_error(err) ) {
			show_asapo_error("Failed to skip to end of stream", err);
		} else {
			STATUS("Skipped to end of stream (%lli)\n", size);
		}

		asapo_free_handle(&err);
	}

	asapo_free_handle(&err);
}


void *im_asapo_fetch(struct im_asapo *a, size_t *pdata_size,
                     char **pfilename, char **pevent)
{
	void *data_copy;
	AsapoMessageMetaHandle meta;
	AsapoMessageDataHandle data;
	AsapoErrorHandle err;
	uint64_t msg_size;

	if ( a->stream == NULL ) {
		if ( select_last_stream(a) ) {
			return NULL;
		}
		skip_to_stream_end(a);
	}

	err = asapo_new_handle();
	meta = asapo_new_handle();
	data = asapo_new_handle();

	asapo_consumer_get_next(a->consumer, a->group_id, &meta, &data,
	                        a->stream, &err);
	if ( asapo_error_get_type(err) == kEndOfStream ) {
		select_next_stream(a);
		asapo_free_handle(&err);
		asapo_free_handle(&meta);
		asapo_free_handle(&data);
		return NULL;  /* Sandbox will call try again very soon */
	}

	if ( asapo_is_error(err) ) {
		show_asapo_error("Couldn't get next ASAP::O record", err);
		asapo_free_handle(&err);
		asapo_free_handle(&meta);
		asapo_free_handle(&data);
		return NULL;
	}

	msg_size = asapo_message_meta_get_size(meta);

	data_copy = malloc(msg_size);
	if ( data_copy == NULL ) {
		ERROR("Failed to copy data block.\n");
		asapo_free_handle(&err);
		asapo_free_handle(&meta);
		asapo_free_handle(&data);
		return NULL;
	}
	memcpy(data_copy, asapo_message_data_get_as_chars(data), msg_size);

	*pfilename = strdup(asapo_message_meta_get_id(meta));
	*pevent = strdup("//");

	asapo_free_handle(&err);
	asapo_free_handle(&meta);
	asapo_free_handle(&data);

	*pdata_size = msg_size;
	return data_copy;
}


void im_asapo_shutdown(struct im_asapo *a)
{
	if ( a == NULL ) return;
	asapo_free_handle(&a->consumer);
	asapo_free_handle(&a->group_id);
	free(a);
}
