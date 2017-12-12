/*
 * im-sandbox.h
 *
 * Sandbox for indexing
 *
 * Copyright © 2012-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2017 Thomas White <taw@physics.org>
 *   2011      Richard Kirian
 *   2012      Lorenzo Galli
 *   2012      Chunhong Yoon
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

#ifndef IM_SANDBOX_H
#define IM_SANDBOX_H

#include <semaphore.h>

struct sb_shm;

#include "index.h"
#include "stream.h"
#include "cell.h"
#include "process_image.h"

/* Length of event queue */
#define QUEUE_SIZE (256)

/* Maximum length of an event ID including serial number */
#define MAX_EV_LEN (1024)

/* Maximum number of workers */
#define MAX_NUM_WORKERS (1024)

struct sb_shm
{
	pthread_mutex_t term_lock;

	pthread_mutex_t queue_lock;
	int n_events;
	char queue[QUEUE_SIZE][MAX_EV_LEN];
	int no_more;
	char last_ev[MAX_NUM_WORKERS][MAX_EV_LEN];
	int pings[MAX_NUM_WORKERS];

	pthread_mutex_t totals_lock;
	int n_processed;
	int n_hadcrystals;
	int n_crystals;
};

extern void create_sandbox(struct index_args *iargs, int n_proc, char *prefix,
                           int config_basename, FILE *fh,  Stream *stream,
                           const char *tempdir, int serial_offset);

#endif /* IM_SANDBOX_H */
