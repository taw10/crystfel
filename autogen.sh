#!/bin/sh

aclocal -I m4
gtkdocize --copy
autoconf
autoheader
automake -ac
