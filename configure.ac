dnl Process this file with autoconf to produce a configure script.

AC_INIT(pattern_sim, 0.1.0, thomas.white@desy.de)
VERSION=AC_PACKAGE_VERSION

AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE(pattern_sim, "$VERSION")

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

AC_ARG_WITH(hdf5,
[AS_HELP_STRING([--with-hdf5], [specify location of HDF5 library])],
[HDF5_CFLAGS="-I$withval/include"
 HDF5_LIBS="-L$withval/lib -lhdf5"],
[HDF5_LIBS="-lhdf5"])

CFLAGS="$CFLAGS $HDF5_CFLAGS"
LIBS="$LIBS $HDF5_LIBS -lm -lz -lgsl -lgslcblas"

AC_OUTPUT(Makefile src/Makefile)
