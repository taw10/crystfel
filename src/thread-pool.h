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

extern int run_threads(int n_threads, void (*work)(void *, int),
                       void *(*get_task)(void *), void (*final)(void *, void *),
                       void *queue_args, int max,
                       int cpu_num, int cpu_groupsize, int cpu_offset);


#endif	/* THREAD_POOL_H */
