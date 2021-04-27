/*
 * thread-pool.c
 *
 * A thread pool implementation
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2019 Thomas White <taw@physics.org>
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

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <assert.h>

#include "utils.h"


/** \file thread-pool.h */

/* --------------------------- Status label stuff --------------------------- */

static int use_status_labels = 0;
static pthread_key_t status_label_key;

struct worker_args
{
	struct task_queue_range *tqr;
	struct task_queue *tq;
	int id;
};


signed int get_status_label()
{
	int *cookie;

	if ( !use_status_labels ) {
		return -1;
	}

	cookie = pthread_getspecific(status_label_key);
	return *cookie;
}


struct task_queue
{
	pthread_mutex_t  lock;

	int              n_started;
	int              n_completed;
	int              max;

	void *(*get_task)(void *);
	void (*finalise)(void *, void *);
	void *queue_args;
	void (*work)(void *, int);
};


static void *task_worker(void *pargsv)
{
	struct worker_args *w = pargsv;
	struct task_queue *q = w->tq;
	int *cookie;

	cookie = malloc(sizeof(int));
	*cookie = w->id;
	pthread_setspecific(status_label_key, cookie);

	free(w);

	do {

		void *task;
		int cookie;

		/* Get a task */
		pthread_mutex_lock(&q->lock);
		if ( (q->max) && (q->n_started >= q->max) ) {
			pthread_mutex_unlock(&q->lock);
			break;
		}
		task = q->get_task(q->queue_args);

		/* No more tasks? */
		if ( task == NULL ) {
			pthread_mutex_unlock(&q->lock);
			break;
		}

		q->n_started++;
		pthread_mutex_unlock(&q->lock);

		cookie = *(int *)pthread_getspecific(status_label_key);
		q->work(task, cookie);

		/* Update totals etc */
		pthread_mutex_lock(&q->lock);
		q->n_completed++;
		if ( q->finalise ) {
			q->finalise(q->queue_args, task);
		}
		pthread_mutex_unlock(&q->lock);

	} while ( 1 );

	free(cookie);

	return NULL;
}


/**
 * \param n_threads The number of threads to run in parallel
 * \param work The function to be called to do the work
 * \param get_task The function which will determine the next unassigned task
 * \param final The function which will be called to clean up after a task
 * \param queue_args A pointer to any data required to determine the next task
 * \param max Stop calling get_task after starting this number of jobs
 * \param cpu_num Ignored
 * \param cpu_groupsize Ignored
 * \param cpu_offset Ignored
 *
 * \p get_task will be called every time a worker is idle.  It returns either
 * NULL, indicating that no further work is available, or a pointer which will
 * be passed to \p work.
 *
 * \p final will be called once per image, and will be given both queue_args
 * and the last task pointer.
 *
 * \p get_task and \p final will be called only under lock, and so do NOT need to
 * be re-entrant or otherwise thread safe.  'work', of course, needs to be
 * thread safe.
 *
 * Work will stop after \p max tasks have been processed whether get_task
 * returned NULL or not.  If \p max is zero, all tasks will be processed.
 *
 * \returns The number of tasks completed.
 **/
int run_threads(int n_threads, TPWorkFunc work,
                TPGetTaskFunc get_task, TPFinalFunc final,
                void *queue_args, int max,
                int cpu_num, int cpu_groupsize, int cpu_offset)
{
	pthread_t *workers;
	int i;
	struct task_queue q;

	pthread_key_create(&status_label_key, NULL);

	workers = malloc(n_threads * sizeof(pthread_t));

	pthread_mutex_init(&q.lock, NULL);
	q.work = work;
	q.get_task = get_task;
	q.finalise = final;
	q.queue_args = queue_args;
	q.n_started = 0;
	q.n_completed = 0;
	q.max = max;

	/* Now it's safe to start using the status labels */
	if ( n_threads > 1 ) use_status_labels = 1;

	/* Start threads */
	for ( i=0; i<n_threads; i++ ) {

		struct worker_args *w;

		w = malloc(sizeof(struct worker_args));

		w->tq = &q;
		w->tqr = NULL;
		w->id = i;

		if ( pthread_create(&workers[i], NULL, task_worker, w) ) {
			/* Not ERROR() here */
			fprintf(stderr, "Couldn't start thread %i\n", i);
			n_threads = i;
			break;
		}

	}

	/* Join threads */
	for ( i=0; i<n_threads; i++ ) {
		pthread_join(workers[i], NULL);
	}

	use_status_labels = 0;

	free(workers);

	return q.n_completed;
}
