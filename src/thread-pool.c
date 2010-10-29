/*
 * thread-pool.c
 *
 * A thread pool implementation
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <assert.h>

#include "utils.h"


/* ---------------------------------- Range --------------------------------- */

enum {
	TASK_READY,
	TASK_RUNNING,
	TASK_FINISHED,
};


struct task_queue_range
{
	pthread_mutex_t  lock;

	int              n_tasks;
	int             *status;
	int              n_done;

	void (*work)(int, void *);
	void *work_args;

	const char      *text;
};


static void *range_worker(void *pargsv)
{
	struct task_queue_range *q = pargsv;

	do {

		int i;
		int found = 0;
		int mytask = -1;

		/* Get a task */
		pthread_mutex_lock(&q->lock);
		for ( i=0; i<q->n_tasks; i++ ) {
			if ( q->status[i] == TASK_READY ) {
				mytask = i;
				found = 1;
				q->status[i] = TASK_RUNNING;
				break;
			}
		}
		pthread_mutex_unlock(&q->lock);

		/* No more tasks? */
		if ( !found ) break;

		q->work(mytask, q->work_args);

		/* Mark this task as done, update totals etc */
		pthread_mutex_lock(&q->lock);
		q->status[mytask] = TASK_FINISHED;
		q->n_done++;
		if ( q->text != NULL ) {
			progress_bar(q->n_done, q->n_tasks, q->text);
		}
		pthread_mutex_unlock(&q->lock);

	} while ( 1 );

	return NULL;
}


void run_thread_range(int n_tasks, int n_threads, const char *text,
                      void (*work)(int, void *), void *work_args)
{
	pthread_t *workers;
	int i;
	struct task_queue_range q;

	/* The nation of CrystFEL prides itself on having 0% unemployment. */
	if ( n_threads > n_tasks ) n_threads = n_tasks;

	workers = malloc(n_threads * sizeof(pthread_t));

	q.status = malloc(n_tasks * sizeof(int));
	pthread_mutex_init(&q.lock, NULL);
	q.n_tasks = n_tasks;
	q.work = work;
	q.work_args = work_args;
	q.n_done = 0;
	q.text = text;

	for ( i=0; i<n_tasks; i++ ) {
		q.status[i] = TASK_READY;
	}

	/* Start threads */
	for ( i=0; i<n_threads; i++ ) {

		if ( pthread_create(&workers[i], NULL, range_worker, &q) ) {
			ERROR("Couldn't start thread %i\n", i);
			n_threads = i;
			break;
		}

	}

	/* Join threads */
	for ( i=0; i<n_threads; i++ ) {
		pthread_join(workers[i], NULL);
	}

	free(q.status);
	free(workers);
}


/* ---------------------------- Custom get_task() --------------------------- */

struct task_queue
{
	pthread_mutex_t  lock;

	int              n_started;
	int              n_completed;
	int              max;
	int              n_cookies;
	int             *cookies;

	void *(*get_task)(void *);
	void (*finalise)(void *, void *);
	void *queue_args;
	void (*work)(void *, int);
};


static void *task_worker(void *pargsv)
{
	struct task_queue *q = pargsv;

	do {

		void *task;
		int i;
		int mycookie = -1;
		int found = 0;

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

		/* Find a cookie */
		for ( i=0; i<q->n_cookies; i++ ) {
			if ( q->cookies[i] == 0 ) {
				mycookie = i;
				found = 1;
				q->cookies[i] = 1;
				break;
			}
		}
		assert(found);

		q->n_started++;
		pthread_mutex_unlock(&q->lock);

		q->work(task, mycookie);

		/* Update totals, release cookie etc */
		pthread_mutex_lock(&q->lock);
		q->n_completed++;
		q->cookies[mycookie] = 0;
		if ( q->finalise ) {
			q->finalise(q->queue_args, task);
		}
		pthread_mutex_unlock(&q->lock);

	} while ( 1 );

	return NULL;
}


int run_threads(int n_threads, void (*work)(void *, int),
                void *(*get_task)(void *), void (*final)(void *, void *),
                void *queue_args, int max)
{
	pthread_t *workers;
	int i;
	struct task_queue q;

	workers = malloc(n_threads * sizeof(pthread_t));

	pthread_mutex_init(&q.lock, NULL);
	q.work = work;
	q.get_task = get_task;
	q.finalise = final;
	q.queue_args = queue_args;
	q.n_started = 0;
	q.n_completed = 0;
	q.max = max;
	q.n_cookies = n_threads;
	q.cookies = malloc(q.n_cookies * sizeof(int));

	for ( i=0; i<n_threads; i++ ) {
		q.cookies[i] = 0;
	}

	/* Start threads */
	for ( i=0; i<n_threads; i++ ) {

		if ( pthread_create(&workers[i], NULL, task_worker, &q) ) {
			ERROR("Couldn't start thread %i\n", i);
			n_threads = i;
			break;
		}

	}

	/* Join threads */
	for ( i=0; i<n_threads; i++ ) {
		pthread_join(workers[i], NULL);
	}

	free(workers);
	free(q.cookies);

	return q.n_completed;
}
