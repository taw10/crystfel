/*
 * symop.y
 *
 * Parser for symmetry operations
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
  #include <stdio.h>

  #include "rational.h"
  #include "symmetry.h"

  #include "symop-parse.h"
  #include "symop-lex.h"

  void symoperror(void *scanner, RationalMatrix *m, SymOpList *list, const char *s);
%}

%define api.prefix {symop}
%define api.pure full
%parse-param {void *scanner} {RationalMatrix *m} {SymOpList *list}
%lex-param {void *scanner}

%code requires {
  #include "symmetry.h"
}

%union {
  RationalMatrix *m;  /* Full rational matrix */
  Rational rv[3];     /* Rational vector, e.g. '1/2h+3k' */
  Rational r;         /* Rational number */
  int n;              /* Just a number */
}

%token SEMICOLON
%token COMMA
%token NUMBER
%token OPENB CLOSEB
%token H K L

%left PLUS MINUS
%left DIVIDE
%precedence MUL
%precedence NEG

%type <m> symop
%type <rv> axexpr
%type <rv> part
%type <n> NUMBER
%type <r> fraction

%{
static int try_add_symop(void *scanner, SymOpList *list, RationalMatrix *m, int complain)
{
	if ( list == NULL ) {
		/* Only complain if this isn't the only operation provided */
		if ( complain ) {
			yyerror(scanner, m, list, "Must be a single symmetry operation");
		}
		return 1;
	} else {
		IntegerMatrix *im;
		im = intmat_from_rtnl_mtx(m);
		if ( im == NULL ) {
			yyerror(scanner, m, list, "Symmetry operations must all be integer");
			return 1;
		} else {
			add_symop(list, im);
		}
	}
	return 0;
}
%}

%%

symoplist:
  symop                       { try_add_symop(scanner, list, m, 0); }
| symoplist SEMICOLON symop   { if ( try_add_symop(scanner, list, m, 1) ) YYERROR; }
;

symop:
  axexpr COMMA axexpr COMMA axexpr { rtnl_mtx_set(m, 0, 0, $1[0]);
                                     rtnl_mtx_set(m, 1, 0, $1[1]);
                                     rtnl_mtx_set(m, 2, 0, $1[2]);
                                     rtnl_mtx_set(m, 0, 1, $3[0]);
                                     rtnl_mtx_set(m, 1, 1, $3[1]);
                                     rtnl_mtx_set(m, 2, 1, $3[2]);
                                     rtnl_mtx_set(m, 0, 2, $5[0]);
                                     rtnl_mtx_set(m, 1, 2, $5[1]);
                                     rtnl_mtx_set(m, 2, 2, $5[2]);
                                   }
;

axexpr:
  part                      { int i;  for ( i=0; i<3; i++ ) $$[i] = $1[i]; }
| axexpr PLUS axexpr        { int i;  for ( i=0; i<3; i++ ) $$[i] = rtnl_add($1[i], $3[i]); }
| axexpr MINUS axexpr       { int i;  for ( i=0; i<3; i++ ) $$[i] = rtnl_sub($1[i], $3[i]); }
| MINUS axexpr %prec NEG    { int i;  for ( i=0; i<3; i++ ) $$[i] = rtnl_sub(rtnl_zero(), $2[i]); }
| OPENB axexpr CLOSEB       { int i;  for ( i=0; i<3; i++ ) $$[i] = $2[i]; }
| axexpr DIVIDE NUMBER      { int i;  for ( i=0; i<3; i++ ) $$[i] = rtnl_div($1[i], rtnl($3, 1)); }
| NUMBER axexpr %prec MUL   { int i;  for ( i=0; i<3; i++ ) $$[i] = rtnl_mul($2[i], rtnl($1, 1)); }
| fraction axexpr %prec MUL { int i;  for ( i=0; i<3; i++ ) $$[i] = rtnl_mul($2[i], $1); }
;

part:
  H  { $$[0] = rtnl(1, 1); $$[1] = rtnl_zero(); $$[2] = rtnl_zero(); }
| K  { $$[1] = rtnl(1, 1); $$[0] = rtnl_zero(); $$[2] = rtnl_zero(); }
| L  { $$[2] = rtnl(1, 1); $$[0] = rtnl_zero(); $$[1] = rtnl_zero(); }
;

fraction:
  NUMBER DIVIDE NUMBER     { $$ = rtnl($1, $3); }
;

%%

void symoperror(void *scanner, RationalMatrix *m, SymOpList *list, const char *s) {
	printf("Error: %s\n", s);
}
