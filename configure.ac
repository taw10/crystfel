dnl Process this file with autoconf to produce a configure script.

AC_INIT(crystfel, 0.2.0, taw@physics.org)
AC_CONFIG_AUX_DIR([build-aux])
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([subdir-objects])

AC_PROG_CC
gl_EARLY
AM_PROG_CC_C_O
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_RANLIB
AM_SILENT_RULES([yes])

AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h stdlib.h unistd.h])
AC_C_CONST
AC_FUNC_MALLOC
AC_TYPE_SIGNAL

gl_INIT


AC_ARG_WITH(hdf5,
[AS_HELP_STRING([--with-hdf5], [specify location of HDF5 library])],
[HDF5_CFLAGS="-I$withval/include"
 HDF5_LIBS="-L$withval/lib -lhdf5"],
[HDF5_LIBS="-lhdf5"])


PKG_CHECK_MODULES([GSL], [gsl], [],
[
   AC_ERROR([GSL is required])
])
PKG_CHECK_MODULES([GSL_version], [gsl >= 1.13],
[
   gl_IGNORE_UNUSED_LIBRARIES
])


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
      AC_DEFINE([HAVE_GTK], [1], [Define to 1 if GTK and GDK are available])
   ],
   [
      AC_MSG_WARN([GTK not found.  hdfsee will not be built.])
   ])
],
[
   AC_MSG_RESULT([no])
])
AM_CONDITIONAL([HAVE_GTK], test "x$havegtk" = "xtrue")


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


AC_MSG_CHECKING([whether to use PNG])
AC_ARG_ENABLE([png],
    AS_HELP_STRING([--disable-png], [Disable the use of libPNG]))
AS_IF([test "x$enable_png" != "xno"],
[
   AC_MSG_RESULT([yes])
   PKG_CHECK_MODULES([libPNG], [libpng >= 1.2.0],
   [
     AC_DEFINE([HAVE_LIBPNG], [1], [Define to 1 if libPNG is available])
   ])
],[
   AC_MSG_RESULT([no])
])


AC_MSG_CHECKING([whether to use TIFF])
AC_ARG_ENABLE([tiff],
    AS_HELP_STRING([--disable-tiff], [Disable the use of libTIFF]))
AS_IF([test "x$enable_tiff" != "xno"],
[
   AC_ARG_WITH(libtiff,
   [AS_HELP_STRING([--with-libtiff], [specify location of libTIFF library])],
   [LIBTIFF_CFLAGS="-I$withval/include"
    LIBTIFF_LIBS="-L$withval/lib -ltiff"],
   [LIBTIFF_LIBS="-ltiff"])
   AC_MSG_RESULT([yes])
   AC_DEFINE([HAVE_TIFF], [1], [Define to 1 if libTIFF is available])
],
[
   AC_MSG_RESULT([no])
])


AC_MSG_CHECKING([whether to use Cairo])
AC_ARG_ENABLE([cairo],
              AS_HELP_STRING([--disable-cairo], [Disable the use of Cairo]))
AS_IF([test "x$enable_cairo" != "xno"],
[
   AC_MSG_RESULT([yes])
   PKG_CHECK_MODULES([Cairo], [cairo >= 1.2.0],
   [
     havecairo=true
     AC_DEFINE([HAVE_CAIRO], [1], [Define to 1 if Cairo is available])
   ])
],
[
   AC_MSG_RESULT([no])
])
AM_CONDITIONAL([HAVE_CAIRO], test x$havecairo = xtrue)


PKG_CHECK_MODULES([GDK_pixbuf], [gdk-pixbuf], [],
[
   PKG_CHECK_MODULES([GDK_pixbuf_2], [gdk-pixbuf-2.0], [],
   [
      AC_ERROR([gdk-pixbuf is required.])
   ])
])


CFLAGS="$CFLAGS $HDF5_CFLAGS $GSL_CFLAGS $OPENCL_CFLAGS -pthread $GTK_CFLAGS"
CFLAGS="$CFLAGS $LIBTIFF_CFLAGS -D_GNU_SOURCE $libPNG_CFLAGS $Cairo_CFLAGS"
CFLAGS="$CFLAGS $GDK_pixbuf_CFLAGS $GDK_pixbuf_2_CFLAGS"
LIBS="$LIBS $HDF5_LIBS -lm -lz $GSL_LIBS $GTK_LIBS $OPENCL_LIBS -pthread"
LIBS="$LIBS $LIBTIFF_LIBS $libPNG_LIBS $Cairo_LIBS $GDK_pixbuf_LIBS"
LIBS="$LIBS $GDK_pixbuf_2_LIBS $LDFLAGS"


AC_CONFIG_FILES(Makefile lib/Makefile)
AC_OUTPUT
