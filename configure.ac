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


AC_ARG_WITH(gsl,
[AS_HELP_STRING([--with-gsl],
                [specify location of GSL (instead of using pkg-config)])],
[GSL_CFLAGS="-I$withval/include"
 GSL_LIBS="-L$withval/lib -lgsl -lgslcblas -lm"],
[PKG_CHECK_MODULES([GSL], [gsl])])


AC_ARG_ENABLE(gsl-fudge, AS_HELP_STRING([--enable-gsl-fudge],
   [Use this option if you get lots of linker errors concerning GSL and CBLAS]))
AC_MSG_CHECKING([whether to use the GSL CBLAS fudge])
AS_IF([test "x$enable_gsl_fudge" == "xyes"],
[
   AC_MSG_RESULT([yes])
   AC_DEFINE([GSL_FUDGE], [1], [Define to 1 if the GSL fudge is needed])
],
[
   AC_MSG_RESULT([no])
])


AC_ARG_WITH(opencl, AS_HELP_STRING([--with-opencl], [Use OpenCL]))
AC_MSG_CHECKING([whether to use OpenCL])
AS_IF([test "x$with_opencl" == "xyes"],
[
   AC_MSG_RESULT([yes])

   AC_MSG_CHECKING([cl.h])
   AS_IF([test -f /System/Library/Frameworks/OpenCL.framework/Headers/cl.h],
   [
      AC_MSG_RESULT([/System/Library/Frameworks/OpenCL.framework/Headers/cl.h])
      OPENCL_CFLAGS="-I/System/Library/Frameworks/OpenCL.framework/Headers"
      OPENCL_LIBS="-framework OpenCL"
   ], [
      AS_IF([test -f /usr/local/cuda/cuda/include/CL/cl.h],
      [
         AC_MSG_RESULT([/usr/local/cuda/cuda/include/CL/cl.h])
         OPENCL_CFLAGS="-I/usr/local/cuda/cuda/include"
         OPENCL_LIBS="-lOpenCL"
         cl_cl_h=true
      ], [
         AS_IF([test -f /usr/local/cuda/include/CL/cl.h],
         [
            AC_MSG_RESULT([/usr/local/cuda/include/CL/cl.h])
            OPENCL_CFLAGS="-I/usr/local/cuda/include"
            OPENCL_LIBS="-lOpenCL"
            cl_cl_h=true
         ], [
            AC_MSG_RESULT([not found, assuming /usr/include/CL/cl.h])
            OPENCL_CFLAGS=""
            OPENCL_LIBS="-lOpenCL"
            cl_cl_h=true
          ])
       ])
   ])
   AC_DEFINE([HAVE_OPENCL], [1], [Define to 1 if OpenCL is available])
   have_opencl=true
], [
   AC_MSG_RESULT([no])
])
AS_IF([test x$cl_cl_h = xtrue],
[
   AC_DEFINE([HAVE_CL_CL_H], [1], [Define to 1 if CL/cl.h should be used])
])


AC_ARG_ENABLE(gtk, AS_HELP_STRING([--disable-gtk], [Disable GTK+/GLib]))
have_gtk=false
AC_MSG_CHECKING([whether to use GTK])
AS_IF([test "x$enable_gtk" != "xno"],
[
   AC_MSG_RESULT([yes])
   AM_PATH_GTK_2_0(2.0.0,
   [
      have_gtk=true
      AC_DEFINE([HAVE_GTK], [1], [Define to 1 if GTK and GDK are available])
   ],
   [
      AC_MSG_WARN([GTK not found.  hdfsee will not be built.])
   ])
],
[
   AC_MSG_RESULT([no])
])


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
      have_libpng=true
   ], [
      have_libpng=false
   ])
],[
   AC_MSG_RESULT([no])
])


AC_MSG_CHECKING([whether to use TIFF])
AC_ARG_ENABLE([tiff],
    AS_HELP_STRING([--disable-tiff], [Disable the use of libTIFF]))
AS_IF([test "x$enable_tiff" != "xno"],
[
   AC_MSG_RESULT([yes])
   AC_ARG_WITH(libtiff,
   [AS_HELP_STRING([--with-libtiff], [specify location of libTIFF library])],
   [
      LIBTIFF_CFLAGS="-I$withval/include"
      LIBTIFF_LIBS="-L$withval/lib -ltiff"
      have_libtiff=true
   ], [
      AC_CHECK_LIB([tiff], [TIFFOpen], [
         have_libtiff=true
         LIBTIFF_LIBS="-ltiff"
      ])
   ])
], [
   AC_MSG_RESULT([no])
])
AS_IF([test x$have_libtiff = xtrue],
[
   AC_DEFINE([HAVE_TIFF], [1], [Define to 1 if libTIFF is available])
])


AC_MSG_CHECKING([whether to use Cairo])
AC_ARG_ENABLE([cairo],
              AS_HELP_STRING([--disable-cairo], [Disable the use of Cairo]))
AS_IF([test "x$enable_cairo" != "xno"],
[
   AC_MSG_RESULT([yes])
   PKG_CHECK_MODULES([Cairo], [cairo >= 1.2.0],
   [
      have_cairo=true
      AC_DEFINE([HAVE_CAIRO], [1], [Define to 1 if Cairo is available])
   ], [
      have_cairo=false
   ])
], [
   AC_MSG_RESULT([no])
])


AC_MSG_CHECKING([whether to use GDK-pixbuf])
AC_ARG_ENABLE([gdk-pixbuf], AS_HELP_STRING([--disable-gdk-pixbuf],
                                           [Disable the use of gtk-pixbuf]))
AS_IF([test "x$enable_gdk_pixbuf" != "xno"],
[
   AC_MSG_RESULT([yes])
   PKG_CHECK_MODULES([GDK_pixbuf], [gdk-pixbuf], [],
   [
      PKG_CHECK_MODULES([GDK_pixbuf_2], [gdk-pixbuf-2.0],
      [
         AC_DEFINE([HAVE_GDK_PIXBUF], [1],
                   [Define to 1 if gdk-pixbuf is available])
         have_gdk_pixbuf=true
      ], [
         have_gdk_pixbuf=false
      ])
   ])
], [
   AC_MSG_RESULT([no])
])


AC_CHECK_LIB([pthread], [sched_setaffinity], [
   AC_DEFINE([HAVE_CPU_AFFINITY], [1],
             [Define to 1 if sched_setaffinity is available.])
])


dnl Conditionals...
AM_CONDITIONAL([HAVE_OPENCL], test x$have_opencl = xtrue)

AM_CONDITIONAL([BUILD_HDFSEE], test x$have_gdk_pixbuf = xtrue \
                                -a "x$have_gtk" = "xtrue" \
                                -a "x$have_cairo" = "xtrue")

AM_CONDITIONAL([BUILD_CUBEIT], test x$have_cairo = xtrue \
                                 -a x$have_libpng = xtrue)


gl_IGNORE_UNUSED_LIBRARIES


CFLAGS="$CFLAGS $HDF5_CFLAGS $GSL_CFLAGS $OPENCL_CFLAGS -pthread $GTK_CFLAGS"
CFLAGS="$CFLAGS $LIBTIFF_CFLAGS -D_GNU_SOURCE $libPNG_CFLAGS $Cairo_CFLAGS"
CFLAGS="$CFLAGS $GDK_pixbuf_CFLAGS $GDK_pixbuf_2_CFLAGS"
LIBS="$LIBS $HDF5_LIBS -lm -lz $GSL_LIBS $GTK_LIBS $OPENCL_LIBS -pthread -lrt"
LIBS="$LIBS $LIBTIFF_LIBS $libPNG_LIBS $Cairo_LIBS $GDK_pixbuf_LIBS"
LIBS="$LIBS $GDK_pixbuf_2_LIBS $LDFLAGS"


GTK_DOC_CHECK([1.11],[--flavour no-tmpl])
AC_CONFIG_MACRO_DIR(m4)


AC_CONFIG_FILES(Makefile lib/Makefile doc/reference/Makefile)
AC_OUTPUT
