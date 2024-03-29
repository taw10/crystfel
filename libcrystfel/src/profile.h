/*
 * profile.h
 *
 * Simple profiling according to wall clock time
 *
 * Copyright © 2016-2022 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2016-2022 Thomas White <taw@physics.org>
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

#ifndef PROFILE_H
#define PROFILE_H

/**
 * \file profile.h
 * Simple wall-clock profiling
 */

extern void profile_init();
extern void profile_print_and_reset(int worker_id);
extern void profile_start(const char *name);
extern void profile_end(const char *name);

#endif	/* PROFILE_H */
