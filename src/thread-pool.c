/*
 * thread-pool.c
 *
 * A thread pool implementation
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
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


static int use_status_labels = 0;
static pthread_key_t status_label_key;
pthread_mutex_t stderr_lock = PTHREAD_MUTEX_INITIALIZER;

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
	struct worker_args *w = pargsv;
	struct task_queue_range *q = w->tqr;
	int *cookie;

	cookie = malloc(sizeof(int));
	*cookie = w->id;
	pthread_setspecific(status_label_key, cookie);

	free(w);

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

	free(cookie);

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

	pthread_key_create(&status_label_key, NULL);

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

	/* Now it's safe to start using the status labels */
	if ( n_threads > 1 ) use_status_labels = 1;

	/* Start threads */
	for ( i=0; i<n_threads; i++ ) {

		struct worker_args *w;

		w = malloc(sizeof(struct worker_args));

		w->tqr = &q;
		w->tq = NULL;
		w->id = i;

		if ( pthread_create(&workers[i], NULL, range_worker, w) ) {
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


int run_threads(int n_threads, void (*work)(void *, int),
                void *(*get_task)(void *), void (*final)(void *, void *),
                void *queue_args, int max)
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
