/*
 * fromfile.h
 *
 * Perform indexing from solution file
 *
 *
 * Authors:
 *   2020 Pascal Hogan-Lamarre <pascal.hogan.lamarre@mail.utoronto.ca>
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

#ifndef FROMFILE_H
#define FROMFILE_H

struct fromfile_keys;
struct fromfile_entries;
struct fromfile_private;

#include "image.h"
#include "cell.h"
#include "uthash.h"

extern void print_struct(struct fromfile_entries *sol_hash);

extern void full_print_struct(struct fromfile_entries *sol_hash);

extern void *fromfile_prepare(char *solution_filename, UnitCell *cell);

extern int fromfile_index(struct image *image, void *mpriv, int crystal_number);


#endif	/* FROMFILE_H */