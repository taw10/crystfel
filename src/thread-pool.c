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

#include "utils.h"


struct task_queue
{
	pthread_mutex_t  lock;

	int              n_tasks;
	int             *done;
	int              n_done;

	void (*work)(int, void *);
	void *work_args;

	const char      *text;
};


static void *worker_thread(void *pargsv)
{
	struct task_queue *q = pargsv;

	do {

		int i;
		int found = 0;
		int mytask = -1;

		/* Get a task */
		pthread_mutex_lock(&q->lock);
		for ( i=0; i<q->n_tasks; i++ ) {
			if ( q->done[i] == 0 ) {
				mytask = i;
				found = 1;
			}
		}
		pthread_mutex_unlock(&q->lock);

		/* No more tasks? */
		if ( !found ) break;

		q->work(mytask, q->work_args);

		/* Mark this task as done, update totals etc */
		pthread_mutex_lock(&q->lock);
		q->done[mytask] = 1;
		q->n_done++;
		progress_bar(q->n_done, q->n_tasks, q->text);
		pthread_mutex_unlock(&q->lock);

	} while ( 1 );

	return NULL;
}


void munch_threads(int n_tasks, int n_threads, const char *text,
                   void (*work)(int, void *), void *work_args)
{
	pthread_t *workers;
	int i;
	struct task_queue q;

	/* The nation of CrystFEL prides itself on having 0% unemployment. */
	if ( n_threads > n_tasks ) n_threads = n_tasks;

	workers = malloc(n_threads * sizeof(pthread_t));

	q.done = malloc(n_tasks * sizeof(int));
	pthread_mutex_init(&q.lock, NULL);
	q.n_tasks = n_tasks;
	q.work = work;
	q.work_args = work_args;
	q.n_done = 0;
	q.text = text;

	for ( i=0; i<n_tasks; i++ ) {
		q.done[i] = 0;
	}

	/* Start threads */
	for ( i=0; i<n_threads; i++ ) {

		if ( pthread_create(&workers[i], NULL, worker_thread, &q) ) {
			ERROR("Couldn't start thread %i\n", i);
			n_threads = i;
			break;
		}

	}

	/* Join threads */
	for ( i=0; i<n_threads; i++ ) {
		pthread_join(workers[i], NULL);
	}

	free(q.done);
	free(workers);
}
