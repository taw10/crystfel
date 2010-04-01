dnl Process this file with autoconf to produce a configure script.

AC_INIT(crystfel, 0.1.0, taw@physics.org)
AM_SILENT_RULES([yes])
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


AC_ARG_WITH(gsl,
[AS_HELP_STRING([--with-gsl], [specify location of GSL])],
[GSL_CFLAGS="-I$withval/include"
 GSL_LIBS="-L$withval/lib -lgsl -lgslcblas"],
[GSL_LIBS="-lgsl -lgslcblas"])


AC_ARG_ENABLE(opencl, AS_HELP_STRING([--enable-opencl], [Enable the use of OpenCL]))
AC_MSG_CHECKING([whether to use OpenCL])
AS_IF([test "x$enable_opencl" == "xyes"],
[
   AC_MSG_RESULT([yes])
   AC_CHECK_HEADERS([/System/Library/Frameworks/OpenCL.framework/Headers/cl.h],
   [
      OPENCL_CFLAGS="-I/System/Library/Frameworks/OpenCL.framework/Headers"
      OPENCL_LIBS="-framework OpenCL"
   ],
   [
      OPENCL_CFLAGS="-I/usr/include/CL"
      OPENCL_LIBS="-lOpenCL"
   ])
   AC_DEFINE([HAVE_OPENCL], [1], [Define to 1 if OpenCL is available])
   have_opencl=true
],
[
   AC_MSG_RESULT([no])
])
AM_CONDITIONAL([HAVE_OPENCL], test x$have_opencl = xtrue)


AC_ARG_ENABLE(gtk, AS_HELP_STRING([--disable-gtk], [Disable GTK+/GLib]))
havegtk=false
AC_MSG_CHECKING([whether to use GTK])
AS_IF([test "x$enable_gtk" != "xno"],
[
   AC_MSG_RESULT([yes])
   AM_PATH_GTK_2_0(2.0.0,
   [
      havegtk=true
      GTK_LIBS="$GTK_LIBS"
   ],
   [
      AC_MSG_WARN([GTK not found.  hdfsee will not be built.])
   ])
],
[
   AC_MSG_RESULT([no])
])
AM_CONDITIONAL([HAVE_GTK], test x$havegtk = xtrue)


AC_CHECK_FUNCS(forkpty,
   AC_DEFINE([HAVE_FORKPTY_BSD], [1],
             [Define to 1 if a BSD-style forkpty is available]),
   AC_CHECK_LIB(util,forkpty,
   [
      AC_DEFINE([HAVE_FORKPTY_LINUX], [1],
                [Define to 1 if a Linux-style forkpty is available])
      LIBS="$LIBS -lutil"
   ])
)


CFLAGS="$CFLAGS $HDF5_CFLAGS $GTK_CFLAGS $GSL_CFLAGS $OPENCL_CFLAGS -pthread"
LIBS="$LIBS $HDF5_LIBS -lm -lz $GSL_LIBS $GTK_LIBS $OPENCL_LIBS -pthread"

AC_OUTPUT(Makefile src/Makefile data/Makefile)
