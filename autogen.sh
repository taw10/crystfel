#!/bin/sh

aclocal -I m4
autoconf
autoheader
automake -ac
gtkdocize --copy || exit 1
