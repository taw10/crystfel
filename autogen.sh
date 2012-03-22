#!/bin/sh

   aclocal -I m4 \
&& libtoolize --force \
&& gtkdocize --copy \
&& autoheader --force \
&& automake --add-missing --copy --force \
&& autoconf --force
