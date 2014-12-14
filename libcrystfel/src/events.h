/*
 * events.h
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef EVENTS_H
#define EVENTS_H

struct event
{
	char **path_entries;
	int path_length;
	int *dim_entries;
	int dim_length;
};

struct event_list
{
	struct event **events;
	int num_events;
};

struct filename_plus_event
{
	char *filename;
	struct event *ev;
};

enum
{
	HYSL_UNDEFINED = -99,
	HYSL_PLACEHOLDER = -98,
	HYSL_FS = -1,
	HYSL_SS = -2
};

struct dim_structure
{
	int *dims;
	int num_dims;
};

extern struct event *initialize_event();
extern int push_path_entry_to_event(struct event *ev, const char *entry);
extern int pop_path_entry_from_event(struct event *ev);
extern int push_dim_entry_to_event(struct event *ev, int entry);
extern int pop_dim_entry_from_event(struct event *ev);
extern struct event *copy_event(struct event *ev);
extern void free_event(struct event *ev);
extern char *get_event_string(struct event *ev);
extern struct event *get_event_from_event_string(const char *ev_string);
extern char *event_path_placeholder_subst(const char *ev_name,
                                          const char *data);
extern char *partial_event_substitution(struct event *ev, const char *data);
extern char *retrieve_full_path(struct event *ev, const char *data);


extern struct filename_plus_event *initialize_filename_plus_event();
extern void free_filename_plus_event(struct filename_plus_event *fpe);


extern struct event_list *initialize_event_list();
extern int append_event_to_event_list(struct event_list *ev_list,
                                   struct event *ev);
int add_non_existing_event_to_event_list(struct event_list *ev_list,
                                         struct event *ev);
extern struct event_list *copy_event_list(struct event_list *el);
extern int find_event(struct event *ev, struct event_list *el);
extern void free_event_list(struct event_list *el);


extern struct dim_structure *initialize_dim_structure();
extern struct dim_structure *default_dim_structure();
extern int set_dim_structure_entry(struct dim_structure *hsd,
                                   const char *string_dim,
                                   const char *val_string);
extern void free_dim_structure_entry(struct dim_structure *hsd);
extern void free_dim_structure(struct dim_structure *hsd);

#endif	/* EVENTS_H */
