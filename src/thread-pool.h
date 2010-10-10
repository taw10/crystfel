/*
 * thread-pool.h
 *
 * A thread pool implementation
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


/* work() will be called with a number and work_args.  The number will be
 * unique and in the range 0..n_tasks.  A progress bar will be shown using
 * "text" and the progress through the tasks. */
extern void run_thread_range(int n_tasks, int n_threads, const char *text,
                             void (*work)(int, void *), void *work_args);


/* get_task() will be called every time a worker is idle.  It returns either
 * NULL, indicating that no further work is available, or a pointer which will
 * be passed to work().  Work will stop after 'max' tasks have been processed.
 * get_task() does not need to be re-entrant.
 * Returns: the number of tasks processed. */
extern int run_threads(int n_threads, void (*work)(void *),
                       void *(*get_task)(void *), void *queue_args, int max);


#endif	/* THREAD_POOL_H */
