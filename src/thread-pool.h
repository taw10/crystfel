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


extern void munch_threads(int n_tasks, int n_threads, const char *text,
                          void (*work)(int, void *), void *work_args);


#endif	/* THREAD_POOL_H */
