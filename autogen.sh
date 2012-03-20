#!/bin/sh

gtkdocize --copy
libtoolize
aclocal -I m4
autoconf
autoheader
automake -ac
