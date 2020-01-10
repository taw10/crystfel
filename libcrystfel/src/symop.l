/*
 * symop.l
 *
 * Lexical scanner for symmetry operations
 *
 * Copyright © 2019-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019 Thomas White <taw@physics.org>
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

%{
  #define YYDEBUG 1
  #include "rational.h"
  #include "symop-parse.h"
%}

%option prefix="symop"
%option noyywrap nounput noinput

%%

[,]        { return COMMA; }
[0-9]+     { symoplval.n = atoi(yytext); return NUMBER; }
[/]        { return DIVIDE; }
[+]        { return PLUS; }
[-]        { return MINUS; }
[ahx]      { return H; }
[bky]      { return K; }
[clz]      { return L; }
[(]        { return OPENB; }
[)]        { return CLOSEB; }
[;]        { return SEMICOLON; }

%%
