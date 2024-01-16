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
#include "utils.h"

#ifndef CLOCK_MONOTONIC_RAW
#define CLOCK_MONOTONIC_RAW (CLOCK_MONOTONIC)
#endif

struct _profile_block
{
	char *name;

	time_t start_sec;
	long start_nsec;
	double total_time;

	struct _profile_block *parent;
	struct _profile_block **children;
	int n_children;
	int max_children;
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

	b = cfmalloc(sizeof(struct _profile_block));
	if ( b == NULL ) return NULL;

	b->name = cfstrdup(name);
	if ( b->name == NULL ) {
		cffree(b);
		return NULL;
	}
	b->n_children = 0;
	b->max_children = 0;
	b->children = NULL;

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
		pd = cfmalloc(sizeof(struct _profiledata));
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
	char **subbufs;
	char *full_buf;

	subbufs = cfmalloc(b->n_children * sizeof(char *));
	if ( subbufs == NULL ) return NULL;

	total_len = 32 + strlen(b->name);
	for ( i=0; i<b->n_children; i++ ) {
		subbufs[i] = format_profile_block(b->children[i]);
		if ( subbufs[i] == NULL ) return NULL;
		total_len += 1 + strlen(subbufs[i]);
	}

	full_buf = cfmalloc(total_len);
	snprintf(full_buf, 32, "(%s %.3f", b->name, b->total_time);
	for ( i=0; i<b->n_children; i++ ) {
		strcat(full_buf, " ");
		strcat(full_buf, subbufs[i]);
		cffree(subbufs[i]);
	}
	strcat(full_buf, ")");

	cffree(subbufs);

	return full_buf;
}


static void free_profile_block(struct _profile_block *b)
{
	int i;
	for ( i=0; i<b->n_children; i++ ) {
		free_profile_block(b->children[i]);
	}
	cffree(b->children);
	cffree(b->name);
	cffree(b);
}


void profile_print_and_reset(int worker_id)
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
	buf2 = cfmalloc(8+strlen(buf));
	size_t len = 8+strlen(buf);
	snprintf(buf2, len, "%i %s\n", worker_id, buf);
	write(STDOUT_FILENO, buf2, strlen(buf2));
	cffree(buf);
	cffree(buf2);

	free_profile_block(pd->root);
	pd->root = start_profile_block("root");
	pd->current = pd->root;
	pd->root->parent = NULL;
}


void profile_start(const char *name)
{
	struct _profile_block *b;

	if ( pd == NULL ) return;

	if ( pd->current->n_children >= pd->current->max_children ) {
		struct _profile_block **nblock;
		int nmax = pd->current->n_children + 64;
		nblock = cfrealloc(pd->current->children, nmax*sizeof(struct _profile_block *));
		if ( nblock == NULL ) {
			fprintf(stderr, "Failed to allocate profiling record. "
			                "Try again without --profile.\n");
			abort();
		}
		pd->current->children = nblock;
		pd->current->max_children = nmax;
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
