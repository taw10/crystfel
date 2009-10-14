dnl Process this file with autoconf to produce a configure script.

AC_INIT(template_index, 0.1.0, thomas.white@desy.de)
VERSION=AC_PACKAGE_VERSION

AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE(template_index, "$VERSION")

AC_PROG_CC
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S

AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h unistd.h])
AC_C_CONST
AC_FUNC_MALLOC
AC_TYPE_SIGNAL
AC_CHECK_FUNCS([strdup])

LIBS="$LIBS -lhdf5"

AC_OUTPUT(Makefile src/Makefile)
