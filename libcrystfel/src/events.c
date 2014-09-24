/*
 * events.c
 *
 * Event properties
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2014      Valerio Mariani
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

#define _ISOC99_SOURCE
#define _GNU_SOURCE

#include "events.h"

#include <hdf5.h>
#include <string.h>
#include <stdlib.h>


struct event *initialize_event()
{

	struct event *ev;

	ev = malloc(sizeof(struct event));
	ev->path_entries = NULL;
	ev->path_length = 0;

	ev->dim_entries = NULL;
	ev->dim_length = 0;

	return ev;

}


struct event_list *initialize_event_list()
{

	struct event_list *ev_list;

	ev_list = malloc(sizeof(struct event_list));

	ev_list->events = NULL;
	ev_list->num_events = 0;

	return ev_list;

}

struct filename_plus_event *initialize_filename_plus_event()
{

	struct filename_plus_event *fpe;

	fpe = malloc(sizeof(struct filename_plus_event));

	fpe->filename = NULL;
	fpe->ev = NULL;

	return fpe;
}


int event_cmp(struct event *ev1, struct event *ev2)
{

	int pi;
	int di;

	if ( ev1->path_length != ev2->path_length ||
		ev1->dim_length != ev2->dim_length ) {
		return 1;
	}

	for ( pi=0; pi<ev1->path_length; pi++ ) {
		if ( strcmp(ev1->path_entries[pi], ev2->path_entries[pi]) != 0 ) {
			return 1;
		}
	}

	for ( di=0; di<ev1->dim_length; di++ ) {
		if ( ev1->path_entries[di] != ev2->path_entries[di] ) {
			return 1;
		}
	}

	return 0;

}


int add_non_existing_event_to_event_list(struct event_list *ev_list,
                                         struct event *ev)
{

	int evi;
	int found = 0;

	for ( evi=0; evi<ev_list->num_events; evi++ ) {
		if (event_cmp(ev_list->events[evi], ev) == 0 ) {
			found = 1;
			break;
		}
	}

	if ( found == 0) {
		return append_event_to_event_list(ev_list, ev);
	}

	return 0;
}


int append_event_to_event_list(struct event_list *ev_list, struct event *ev)
{

	struct event **new_el;

	new_el = realloc(ev_list->events,
                     (1+ev_list->num_events)*sizeof(struct event*));
	if ( new_el == NULL ) {
        return 1;
	}
	ev_list->events = new_el;
	ev_list->events[ev_list->num_events] = copy_event(ev);
	ev_list->num_events +=1;

	return 0;
}


struct event *copy_event(struct event *ev)
{

	struct event *new_ev;
	int pi, di;

	if ( ev->dim_length == 0 && ev->path_length == 0) {

		new_ev = initialize_event();

	} else {

		new_ev=malloc(sizeof(struct event));

		new_ev->path_entries = malloc(ev->path_length*sizeof(char *));
		new_ev->path_length = ev->path_length;

		new_ev->dim_entries = malloc(ev->dim_length*sizeof(int *));
		new_ev->dim_length = ev->dim_length;

		for ( pi=0; pi<new_ev->path_length; pi++ ) {
			new_ev->path_entries[pi] = strdup(ev->path_entries[pi]);
		}

		for ( di=0; di<new_ev->dim_length; di++ ) {
			new_ev->dim_entries[di] = ev->dim_entries[di];
		}

	}
	return new_ev;
}


struct event_list *copy_event_list(struct event_list *el)
{
	int ei;
	struct event_list *el_copy;
	struct event **events_copy;

	el_copy = malloc(1);
	if ( el_copy == NULL ) {
		return NULL;
	}

	events_copy = malloc(el->num_events);
	if ( events_copy == NULL ) {
		free (el_copy);
		return NULL;
	}
	el_copy->events = events_copy;

	for ( ei=0; ei<el->num_events; ei++ ) {
		el_copy->events[ei]=copy_event(el->events[ei]);
	}

	el_copy->num_events = el->num_events;

	return el_copy;
}



void free_event(struct event *ev)
{
	int pi;

	if ( ev->path_length != 0 ) {
		for ( pi=0; pi<ev->path_length; pi++ ) {
			free(ev->path_entries[pi]);
		}
	}
	free(ev->dim_entries);
	free(ev);
}


void free_event_list(struct event_list *el)
{
	int ei;

	for ( ei=0; ei<el->num_events; ei++ ) {
		free_event(el->events[ei]);
	}
	free(el);
}


void free_filename_plus_event(struct filename_plus_event *fpe)
{

	free(fpe->filename);

	if ( fpe->ev != NULL ) {
		free_event(fpe->ev);
	}
}


char *get_event_string(struct event *ev)
{
	char *ret_string;
	char *new_ret_string;
	int ret_string_len;

	if ( ev == NULL ) return "null event";

	if ( ev->path_length != 0 ) {

		int pi;

		ret_string = strdup(ev->path_entries[0]);
		ret_string_len = strlen(ret_string);

		for ( pi=1; pi<ev->path_length; pi++ ) {

			new_ret_string = realloc(ret_string,
			                 (ret_string_len+1+strlen(ev->path_entries[pi]))
			                 *sizeof(char));
			if ( new_ret_string == NULL ) {
				return NULL;
			}

			ret_string=new_ret_string;
			strncpy(&ret_string[ret_string_len],"/", 1);
			strncpy(&ret_string[ret_string_len+1],ev->path_entries[pi],
			        strlen(ev->path_entries[pi]));

			ret_string_len += 1+strlen(ev->path_entries[pi]);

		}

		new_ret_string = realloc(ret_string,
                             (1+ret_string_len)*sizeof(char));
		if ( new_ret_string == NULL ) {
			return NULL;
		}

		ret_string = new_ret_string;

		strncpy(&ret_string[ret_string_len], "/", 1);
		ret_string_len += 1;

	} else {

		ret_string = strdup("/");
		ret_string_len = strlen(ret_string);

	}

	if ( ev->dim_length !=0 ) {

		char num_buf[64];
		int di;

		for ( di=0; di<ev->dim_length; di++ ) {
			sprintf(num_buf, "%i", ev->dim_entries[di]);

			new_ret_string = realloc(ret_string,
			                         (ret_string_len+1+strlen(num_buf))
			                         *sizeof(char));
			if ( new_ret_string == NULL ) {
				return NULL;
			}

			ret_string=new_ret_string;

			strncpy(&ret_string[ret_string_len],"/", 1);
			strncpy(&ret_string[ret_string_len+1], num_buf,
                    strlen(num_buf));
			ret_string_len += 1+strlen(num_buf);

		}

	} else {

		new_ret_string = realloc(ret_string,
		                         (1+ret_string_len)*sizeof(char));
		if ( new_ret_string == NULL ) {
			return NULL;
		}

		ret_string = new_ret_string;

		strncpy(&ret_string[ret_string_len], "/", 1);
		ret_string_len += 1;

	}

	new_ret_string = realloc(ret_string,
	                         (1+ret_string_len)*sizeof(char));
	if ( new_ret_string == NULL ) {
		return NULL;
	}

	ret_string = new_ret_string;

	strncpy(&ret_string[ret_string_len], "\0", 1);

	return ret_string;
}


struct event *get_event_from_event_string(char *ev_string)
{
	struct event *ev;
	char *ev_sep;
	char buf_path[1024];
	char buf_dim[1024];
	char *sep;
	char *start;

	ev = initialize_event();
	if ( ev == NULL ) {
		return NULL;
	}

	ev_sep = strstr(ev_string, "//");
	if ( ev_sep == NULL ) {
		return NULL;
	}

	strncpy(buf_path, ev_string, ev_sep-ev_string);
	buf_path[ev_sep-ev_string] = '\0';

	strncpy(buf_dim, ev_sep+2, strlen(ev_sep)-2);
	buf_dim[strlen(ev_sep)-2] = '\0';

	if ( strlen(buf_path) !=0 ) {

		do {

			start = buf_path;

			char buf[2014];

			sep = strstr(start, "/");
			if ( sep != NULL ) {

				strncpy(buf, start, sep-start);
				buf[sep-start]='\0';
				push_path_entry_to_event(ev, buf);
				start = sep + 1;

			} else {

				sprintf(buf,"%s",start);
				push_path_entry_to_event(ev, buf);

			}
		} while (sep);

	}


	if ( strlen(buf_dim) !=0 ) {

		start = buf_dim;

		do {

			char buf[2014];
			int buf_int;

			sep = strstr(start, "/");
			if ( sep != NULL ) {
				strncpy(buf, start, sep-start);
				buf[sep-start]='\0';
				buf_int = atoi(buf);
				push_dim_entry_to_event(ev, buf_int);
				start = sep + 1;

			} else {

				sprintf(buf,"%s",start);
				buf_int = atoi(buf);
				push_dim_entry_to_event(ev, buf_int);

			}
		} while (sep);

	}


	return ev;
}


int push_path_entry_to_event(struct event *ev, const char *entry)
{
		char **new_path_entries;

		new_path_entries = realloc(ev->path_entries,
                                 (1+ev->path_length)*sizeof(char *));
		if ( new_path_entries == NULL ) {
			return 1;
		}

		ev->path_entries = new_path_entries;
		ev->path_entries[ev->path_length] = strdup(entry);
		ev->path_length += 1;

		return 0;
}


int push_dim_entry_to_event(struct event *ev, int entry)
{
	int *new_dim_entries;

	new_dim_entries = realloc(ev->dim_entries,
                                 (1+ev->dim_length)*sizeof(int));
	if ( new_dim_entries == NULL ) {
		return 1;
	}

	ev->dim_entries = new_dim_entries;
	ev->dim_entries[ev->dim_length] = entry;
	ev->dim_length += 1;

	return 0;
}


int pop_path_entry_from_event(struct event *ev)
{
	char **new_path_entries;

	if ( ev->path_length == 0 ) {
			return 1;
	}

	free(ev->path_entries[ev->path_length-1]);

	if ( ev->path_length == 1 ) {
		ev->path_entries = NULL;
		ev->path_length = 0;
		return 0;
	}

	new_path_entries = realloc(ev->path_entries,
	                           (ev->path_length-1)*sizeof(char *));

	if ( new_path_entries == NULL) {
		return 1;
	}

	ev->path_entries = new_path_entries;
	ev->path_length = ev->path_length-1;

	return 0;
}


int pop_dim_entry_from_event(struct event *ev)
{
	int *new_dim_entries;

	if ( ev->dim_length == 0 ) {
			return 1;
	}

	if ( ev->dim_length == 1 ) {
		ev->dim_entries = NULL;
		ev->dim_length = 0;
		return 0;
	}

	new_dim_entries = realloc(ev->dim_entries,
                              (ev->dim_length-1)*sizeof(int));

	if ( new_dim_entries == NULL) {
		return 1;
	}

	ev->dim_entries = new_dim_entries;
	ev->dim_length = ev->dim_length-1;

	return 0;
}


char *event_path_placeholder_subst(const char *entry,
                                   const char *data)
{

	char *ph_loc;
	char *full_path;
	int len_head, len_tail;

	full_path = malloc((strlen(data) + strlen(entry)+1)*sizeof(char));
	ph_loc = strstr(data, "%");
	len_head = ph_loc-data;
	len_tail = strlen(ph_loc);

	strncpy(full_path, data, len_head);
	strncpy(full_path+len_head, entry, strlen(entry));
	strncpy(full_path+len_head+strlen(entry), ph_loc+1, len_tail);
	strncpy(&full_path[strlen(data) + strlen(entry)],"\0",1);

	return full_path;
}


char *retrieve_full_path(struct event *ev, const char *data)
{

	int ei ;
	char *return_value;

	return_value = strdup(data);

	for ( ei=0; ei<ev->path_length; ei++ ) {

		char *tmp_subst_data;
		tmp_subst_data = event_path_placeholder_subst(ev->path_entries[ei],
		                                              return_value);

		free(return_value);
		return_value = strdup(tmp_subst_data);
		free(tmp_subst_data);

	}

	return return_value;

}


char *partial_event_substitution(struct event *ev, const char *data)
{
	int ei ;
	char *return_value;
	char *pholder;

	return_value = strdup(data);
	pholder = strstr(return_value,"%");
	ei = 0;

	while( pholder != NULL) {

		char *tmp_subst_data;

		tmp_subst_data = event_path_placeholder_subst(ev->path_entries[ei],
                                                      return_value);
		free(return_value);
		return_value = strdup(tmp_subst_data);
		free(tmp_subst_data);
		pholder = strstr(return_value, "%");
		ei += 1;
	}

	return return_value;
}


struct dim_structure *initialize_dim_structure()
{
	struct dim_structure *hs;
	hs = malloc(sizeof(struct dim_structure));
	if ( hs == NULL ) {
		return NULL;
	}

	hs->dims = NULL;
	hs->num_dims = 0;

	return hs;
}


struct dim_structure *default_dim_structure()
{
	struct dim_structure *hsd;

	hsd = initialize_dim_structure();

	set_dim_structure_entry(hsd, "dim0", "ss");
	set_dim_structure_entry(hsd, "dim1", "fs");

	return hsd;
}


void free_dim_structure(struct dim_structure *hsd)
{
	int di;

	for ( di=0; di<hsd->num_dims; di++ ) {
		free (hsd->dims);
		free (hsd);
	}
}


static int parse_dim_structure_val(const char *val)
{
	if ( strcmp(val,"%") == 0 ) {
		return HYSL_PLACEHOLDER;
	} else if ( strcmp(val,"ss") == 0 ) {
		return HYSL_SS;
	} else if ( strcmp(val,"fs") == 0 ) {
		return HYSL_FS;
	}
	return atoi(val);

}


int set_dim_structure_entry(struct dim_structure *hsd, const char *string_dim,
                            const char *val_string)
{
	int dim_entry;

	dim_entry = atoi(string_dim+3)+1;

	if ( dim_entry > hsd->num_dims ) {

		int di;

		int *new_dims = malloc(dim_entry*sizeof(int));
		if ( new_dims == NULL ) {
			return 0;
		}


		for ( di=0; di<dim_entry; di++ ) {
			new_dims[di] = HYSL_UNDEFINED;
		}

		for ( di=0; di<hsd->num_dims; di++ ) {
			new_dims[di] = hsd->dims[di];
		}

		new_dims[dim_entry-1] = parse_dim_structure_val(val_string);
		if ( hsd->dims == NULL ) {
			free (hsd->dims);
		}
		hsd->dims = new_dims;
		hsd->num_dims = dim_entry;

		return 1;

	}

	hsd->dims[dim_entry] = parse_dim_structure_val(val_string);
	return 1;
}
