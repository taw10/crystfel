/*
 * profile.c
 *
 * Simple profiling according to wall clock time
 *
 * Copyright © 2016-2022 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2016-2022 Thomas White <taw@physics.org>
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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>

#include "profile.h"

#ifndef CLOCK_MONOTONIC_RAW
#define CLOCK_MONOTONIC_RAW (CLOCK_MONOTONIC)
#endif

#define MAX_PROFILE_CHILDREN 256

struct _profile_block
{
	char *name;

	time_t start_sec;
	long start_nsec;
	double total_time;

	struct _profile_block *parent;
	struct _profile_block *children[MAX_PROFILE_CHILDREN];
	int n_children;
};


struct _profiledata
{
	struct _profile_block *root;
	struct _profile_block *current;
};


struct _profiledata *pd = NULL;


static struct _profile_block *start_profile_block(const char *name)
{
	struct _profile_block *b;

	b = malloc(sizeof(struct _profile_block));
	if ( b == NULL ) return NULL;

	b->name = strdup(name);
	if ( b->name == NULL ) {
		free(b);
		return NULL;
	}
	b->n_children = 0;

#ifdef HAVE_CLOCK_GETTIME
	struct timespec tp;
	clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
	b->start_sec = tp.tv_sec;
	b->start_nsec = tp.tv_nsec;
#endif

	return b;
}


static void stop_profile_block(struct _profile_block *b)
{
#ifdef HAVE_CLOCK_GETTIME
	struct timespec tp;
	clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
	time_t sec = tp.tv_sec - b->start_sec;
	long nsec = tp.tv_nsec - b->start_nsec;
	b->total_time = sec + nsec*1e-9;
#else
	b->total_time = 0.0;
#endif
}


void profile_init()
{
	if ( pd != NULL ) {
		fprintf(stderr, "Attempted to initialise profiling twice!\n");
		fflush(stderr);
		abort();
	}

	if ( pd == NULL ) {
		pd = malloc(sizeof(struct _profiledata));
		if ( pd == NULL ) return;
	}

	pd->root = start_profile_block("root");
	pd->current = pd->root;
	pd->root->parent = NULL;

#ifndef HAVE_CLOCK_GETTIME
	printf("Profiling disabled because clock_gettime is not available\n");
#endif
}


static char *format_profile_block(struct _profile_block *b)
{
	int i;
	size_t total_len = 0;
	char *subbufs[MAX_PROFILE_CHILDREN];
	char *full_buf;

	total_len = 32 + strlen(b->name);
	for ( i=0; i<b->n_children; i++ ) {
		subbufs[i] = format_profile_block(b->children[i]);
		total_len += 1 + strlen(subbufs[i]);
	}

	full_buf = malloc(total_len);
	snprintf(full_buf, 32, "(%s %.3f", b->name, b->total_time);
	for ( i=0; i<b->n_children; i++ ) {
		strcat(full_buf, " ");
		strcat(full_buf, subbufs[i]);
		free(subbufs[i]);
	}
	strcat(full_buf, ")");

	return full_buf;
}


static void free_profile_block(struct _profile_block *b)
{
	int i;
	for ( i=0; i<b->n_children; i++ ) {
		free_profile_block(b->children[i]);
	}
	free(b->name);
	free(b);
}


void profile_print_and_reset()
{
	char *buf;
	char *buf2;

	if ( pd == NULL ) {
		fprintf(stderr, "Profiling not initialised yet!\n");
		fflush(stderr);
		abort();
	}

	if ( pd->current != pd->root ) {
		fprintf(stderr, "Attempted to finalise profiling while not "
		        "on root block (%s)\n", pd->current->name);
		fflush(stderr);
		abort();
	}

	stop_profile_block(pd->root);

	buf = format_profile_block(pd->root);
	buf2 = malloc(2+strlen(buf));
	strcpy(buf2, buf);
	strcat(buf2, "\n");
	write(STDOUT_FILENO, buf2, strlen(buf2));

	free_profile_block(pd->root);
	pd->root = start_profile_block("root");
	pd->current = pd->root;
	pd->root->parent = NULL;
}


void profile_start(const char *name)
{
	struct _profile_block *b;

	if ( pd == NULL ) return;

	if ( pd->current->n_children >= MAX_PROFILE_CHILDREN ) {
		fprintf(stderr, "Too many profile children "
		                "(opening %s inside %s).\n",
		                pd->current->name, name);
		fflush(stderr);
		abort();
	}

	b = start_profile_block(name);
	b->parent = pd->current;
	pd->current->children[pd->current->n_children++] = b;
	pd->current = b;
}


void profile_end(const char *name)
{
	if ( pd == NULL ) return;

	if ( pd->current == NULL ) {
		fprintf(stderr, "No current profile block!\n");
		fflush(stderr);
		abort();
	}

	if ( strcmp(name, pd->current->name) != 0 ) {
		fprintf(stderr, "Attempt to close wrong profile block (%s) "
		        "current block is %s\n", name, pd->current->name);
		fflush(stderr);
		abort();
	}

	stop_profile_block(pd->current);

	pd->current = pd->current->parent;
}
