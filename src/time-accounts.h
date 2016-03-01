/*
 * time-accounts.h
 *
 * Simple profiling according to wall clock time
 *
 * Copyright © 2016 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2016 Thomas White <taw@physics.org>
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

#ifndef TIME_ACCOUNTS_H
#define TIME_ACCOUNTS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

enum timeaccount
{
	TACC_NOTHING,
	TACC_STREAMREAD,
	TACC_SELECT,
	TACC_SIGNALS,
	TACC_QUEUETOPUP,
	TACC_STATUS,
	TACC_ENDCHECK,
	TACC_WAKEUP,
	TACC_WAITPID,
	TACC_HDF5OPEN,
	TACC_HDF5READ,
	TACC_FILTER,
	TACC_RESRANGE,
	TACC_PEAKSEARCH,
	TACC_INDEXING,
	TACC_PREDPARAMS,
	TACC_INTEGRATION,
	TACC_TOTALS,
	TACC_WRITESTREAM,
	TACC_CLEANUP,
	TACC_EVENTWAIT,
	TACC_FINALCLEANUP,
};

typedef struct _timeaccounts TimeAccounts;

extern TimeAccounts *time_accounts_init(void);
extern void time_accounts_free(TimeAccounts *accs);

extern void time_accounts_set(TimeAccounts *accs, enum timeaccount new_acc);

extern void time_accounts_print(TimeAccounts *accs);

#endif	/* TIME_ACCOUNTS_H */
