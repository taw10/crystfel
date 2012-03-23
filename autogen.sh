#!/bin/sh

libtoolize --force --copy \
&& aclocal -I m4 --force \
&& gtkdocize --copy \
&& autoheader --force \
&& automake --add-missing --copy --force \
&& autoconf --force
