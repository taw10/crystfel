#!/bin/sh

   aclocal -I m4 \
&& libtoolize --force --copy \
&& gtkdocize --copy \
&& autoheader --force \
&& automake --add-missing --copy --force \
&& autoconf --force
