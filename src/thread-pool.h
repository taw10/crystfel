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


#include <pthread.h>

extern pthread_mutex_t stderr_lock;
extern signed int get_status_label(void);


/**
 * TPGetTaskFunc:
 * @qargs: The queue_args pointer which was given to run_threads().
 * Returns: A pointer which will be passed to the worker function.
 *
 * This function is called, non-reentrantly, to get a new work item to give to
 * your work function.  The stuff you need to generate the new work item should
 * have been stored in @qargs which was passed to run_threads().
 *
 **/
typedef void *(*TPGetTaskFunc)(void *qargs);


/**
 * TPWorkFunc:
 * @work: The queue_args pointer which was given to run_threads().
 * @cookie: A small integral number which is guaranteed to be unique among all
 * currently running threads.
 *
 * This function is called, reentrantly, for each work item.
 *
 **/
typedef void (*TPWorkFunc)(void *work, int cookie);


/**
 * TPFinalFunc:
 * @qargs: The queue_args pointer which was given to run_threads().
 * @work: The pointer which was returned by your get_task function.
 *
 * This function is called, non-reentrantly, after each work item has been
 * completed.  A typical use might be to update some counters inside @qargs
 * according to fields withing @work which were filled by your 'work' function.
 *
 **/
typedef void (*TPFinalFunc)(void *qargs, void *work);


extern int run_threads(int n_threads, TPWorkFunc work,
                       TPGetTaskFunc get_task, TPFinalFunc final,
                       void *queue_args, int max,
                       int cpu_num, int cpu_groupsize, int cpu_offset);


#endif	/* THREAD_POOL_H */
