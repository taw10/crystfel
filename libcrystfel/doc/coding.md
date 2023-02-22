\page coding CrystFEL coding standards

### Licensing

CrystFEL is distributed under the terms of the GNU General Public License
version 3 or higher.  Contributions are very welcome provided they also use this
license.  If your code is not already licensed compatibly, or if the license if
not clear, we will ask you to re-license it.

Whenever you edit a source file, don't forget to update the copyright dates at
the top.  Add your name and email address if they're not there already.  Be sure
to add your name to the 'AUTHORS' file in the top level folder, as well.

### Scope

The remainder of these rules apply to C code in libcrystfel and the core
CrystFEL programs.  There are currently no specific guidelines for Python,
Perl, shell script, build scripts or other parts of the codebase.

### Indentation

*Indentation* is done with *tabs* and *alignment* is done with spaces.
For example:

    int function(int a, int b)
    {
    <-tab-->int p;  /* <--- Tab character used to indent code inside function */
            char *str;

    <-tab-->do_something(a, "A long string which takes up a lot of space",
    <-tab-->.............str, &p);   /* <--- spaces used to align with bracket */
    }

**Rationale:** Using tab characters makes it easy to align code correctly,
because you can't slip out of alignment by one character.  It also makes the
code look neat whatever width you configure your editor to display tabs as.

### Wrap width

Code should fit into 80 columns (counting tab characters as 8 columns) wherever
possible, with exceptions only in cases where not doing so would cause line
breaks at ugly positions, such as straight after an *opening* bracket.  The
absolute maximum allowable line length is 120 columns, with no exceptions
whatsoever.

For example, this is preferred because it fits into 80 columns (just!):

    very_long_variable_name = very_long_function_name(even_longer_long_variablename,
                                                      var2, var3);

However, if the variable names are even longer, then the following is preferred,
even though it goes over 80 columns:

    very_very_long_variable_name = very_long_function_name(even_longer_long_variablename,
                                                           var2, var3);

This is preferred over the following type of thing, which is just plain ugly:

    very_very_long_variable_name = very_long_function_name(
                                                   even_longer_long_variablename,
                                                           var2, var3);

However, see the point below regarding variable names.

**Rationale:** Aside from ensuring it's always possible to display two parts of
the code side-by-side with a reasonable font size, this is not really about how
the code looks.  Rather, it is to discourage excessive levels of nesting and
encourage smaller, more easily understood functions.  I don't think I've yet
seen any examples of code where the intelligibility could not be improved while
simultaneously reducing the number of levels of indentation.

### Variable, function and type names

Shorter variable and function names, with explanatory comments, are preferred
over long variable names:

    double wavelength_in_angstrom_units;  /* <--- not preferred */

Preferred:

    /* The wavelength in Angstrom units */
    double wl;

***Rationale:*** Shorter variable names make it easier to see the overall
logical or mathematical structure.

"Snake case" is preferred over "camel case":

    int model_option;   /* <--- Preferred */
    int modelOption;    /* <--- Discouraged */

Capitalisation is used for complex type names:

    UnitCell *cell;

### Nested loops

When performing a two or three dimensional iteration which could be considered
as one larger iteration, for example over image coordinates or Miller indices,
it is acceptable to indent as follows:

    for ( h=-10; h<+10; h++ ) {
    for ( k=-10; k<+10; k++ ) {
    for ( l=-10; l<+10; l++ ) {

            /* Do stuff */

    }
    }
    }

In this case, there must be no lines at all, not even blank ones, between each
of the "for" statements and also between each of the final closing braces.

***Rationale:*** Large multi-dimensional loops are common in scientific code,
with more than three levels not at all uncommon.  Any reasonable limit on code
width would be overshot by indenting each level separately.

### Comments

C-style comments are preferred over C++-style, even for one-line comments:

    /* Preferred type of comment */

    // Discouraged type of comment

***Rationale:***  CrystFEL is written in C, not C++.

Multi-line comments should have stars down one side:

    /* This is a multiple-line comment.
     * Here is the second line of the comment */

It's also acceptable to put the closing slash on its own line:

    /* Here is another multiple-line comment.
     * Here is the second line of the comment
     */

### Bracket positions

The opening brace of a function should be on its own line:

    void somefunction(int someparam)
    {
            /* Here is the body of the function */
    }

The opening brace of an if statement should be on the same line as the 'if',
unless the condition is very long, especially if it spans multiple lines:

            if ( a < b ) {
                    do_something(a);
            } else {
                    do_other_something(a);
            }

            if ( a_very_long_condition == structure->wibble.value
              && other_very_long_condition )
            {
                    do_something_completely_different(someparam);
            }

***Rationale:*** This keeps the appearance of 'if' statements compact instead of
spread out, while allowing some visual separation between long conditions and
the conditional code.

### Space around parantheses

Parentheses should have spaces after them in 'if' statements, but not in
function calls. Function arguments should have spaces after the comma.  There
should be no space between the function name and the opening bracket.  That
means:

    if ( something < 3 ) {
            do_something(a, b, c);
    }

or:

    if ( h>3 && k<3 && l==-1 ) {
            do_something(h, k, l);
    }

instead of:

    if (something<3) {
            do_something (a,b,c);
    }

***Rationale:*** I find this guideline helps to encourage good grouping of
logical statements.  Some flexibility is allowed here, though.

### Whitespace

No trailing whitespace at all, even in empty lines.

***Rationale:*** Invisible characters do nothing other than generate VCS
conflicts.

### Cleverness

Transparent, easily understood solutions are preferred over faster ones, except
where you can demonstrate that the code is performance-critical and the benefit
is significant.  Even in that case, copious comments should be provided.

Use of undefined behaviour, even if "it always works", is absolutely forbidden.

### Git/VCS usage

This style of commit message is preferred:

> Strip out libfrosticle references and add new function model

This style of commit message is discouraged:

> Stripped out libfrosticle references, and added new function model

**Rationale:** this encourages you to think in terms of small, self-contained
changes to the code, rather than successive versions with many different changes
combined together.  It also matches the conventions used by most other projects.
