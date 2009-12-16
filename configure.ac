dnl Process this file with autoconf to produce a configure script.

AC_INIT(crystfel, 0.1.0, taw@physics.org)
VERSION=AC_PACKAGE_VERSION

AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE(crystfel, "$VERSION")

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

havegtk=false
AM_PATH_GTK_2_0(2.0.0, havegtk=true,
AC_MSG_WARN([GTK not found.  hdfsee will not be built.]))

AM_CONDITIONAL([HAVE_GTK], test x$havegtk = xtrue)

CFLAGS="$CFLAGS $HDF5_CFLAGS $GTK_CFLAGS"
LIBS="$LIBS $HDF5_LIBS -lm -lz -lgsl -lgslcblas $GTK_LIBS -lgthread-2.0"

AC_OUTPUT(Makefile src/Makefile data/Makefile)
