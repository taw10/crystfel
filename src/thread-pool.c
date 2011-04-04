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

#ifdef HAVE_CPU_AFFINITY
#include <sched.h>
#endif


#include "utils.h"


/**
 * SECTION:thread-pool
 * @short_description: The thread pool
 * @title: The thread pool
 * @section_id:
 * @see_also:
 * @include: "thread-pool.h"
 * @Image:
 *
 * The thread pool helps when running many tasks in parallel.  It takes care of
 * starting and stopping threads, and presents a relatively simple interface to
 * the individual programs.
 */

/* ------------------------------ CPU affinity ------------------------------ */

#ifdef HAVE_CPU_AFFINITY

static void set_affinity(int n, int cpu_num, int cpu_groupsize, int cpu_offset)
{
	cpu_set_t c;
	int group;
	int n_cpu_groups;
	int i;

	if ( cpu_num == 0 ) return;

	CPU_ZERO(&c);

	/* Work out which group this thread belongs to */
	group = (n / cpu_groupsize) + cpu_offset;

	/* Work out which CPUs should be used for this group */
	n_cpu_groups = cpu_num / cpu_groupsize;
	group = group % n_cpu_groups;

	/* Set flags */
	for ( i=0; i<cpu_groupsize; i++ ) {

		int cpu = cpu_groupsize*group + i;

		CPU_SET(cpu, &c);

	}

	if ( sched_setaffinity(0, sizeof(cpu_set_t), &c) ) {

		/* Cannot use ERROR() just yet */
		fprintf(stderr, "%i: Failed to set CPU affinity.\n", n);

	}
}

#else /* HAVE_CPU_AFFINITY */

static void set_affinity(int n, int cpu_num, int cpu_groupsize, int cpu_offset)
{
	/* Do absolutely nothing */
}

#endif /* HAVE_CPU_AFFINITY */


/* --------------------------- Status label stuff --------------------------- */

static int use_status_labels = 0;
static pthread_key_t status_label_key;
pthread_mutex_t stderr_lock = PTHREAD_MUTEX_INITIALIZER;

struct worker_args
{
	struct task_queue_range *tqr;
	struct task_queue *tq;
	int id;
	int cpu_num;
	int cpu_groupsize;
	int cpu_offset;
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

	set_affinity(w->id, w->cpu_num, w->cpu_groupsize, w->cpu_offset);

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
 * run_threads:
 * @n_threads: The number of threads to run in parallel
 * @work: The function to be called to do the work
 * @get_task: The function which will determine the next unassigned task
 * @final: The function which will be called to clean up after a task
 * @queue_args: A pointer to any data required to determine the next task
 * @max: Stop calling get_task() after starting this number of jobs
 * @cpu_num: The number of CPUs in the system
 * @cpu_groupsize: The group size into which the CPUs are grouped
 * @cpu_offset: The CPU group number at which to start pinning threads
 *
 * Returns: The number of tasks completed.
 **/
int run_threads(int n_threads, void (*work)(void *, int),
                void *(*get_task)(void *), void (*final)(void *, void *),
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
		w->cpu_num = cpu_num;
		w->cpu_groupsize = cpu_groupsize;
		w->cpu_offset = cpu_offset;

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
