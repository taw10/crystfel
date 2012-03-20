#!/bin/sh

   aclocal -I m4 \
&& libtoolize \
&& gtkdocize --copy \
&& autoheader \
&& automake -ac \
&& autoconf
