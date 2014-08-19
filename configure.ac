dnl Process this file with autoconf to produce a configure script.

AC_INIT([crystfel],[0.5.3a],[taw@physics.org])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([subdir-objects])

AC_PROG_CC
gl_EARLY
AM_PROG_CC_C_O
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL
AM_SILENT_RULES([yes])

AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h stdlib.h unistd.h])
AC_C_CONST
AC_FUNC_MALLOC

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


AC_ARG_WITH(opencl, AS_HELP_STRING([--with-opencl], [Use OpenCL]))
AC_ARG_WITH(opencl_dir, AS_HELP_STRING([--with-opencl-dir],
            [Specify location of OpenCL headers and libraries]))
AC_MSG_CHECKING([whether to use OpenCL])
AS_IF([test "x$with_opencl_dir" != "x"], [with_opencl=yes])
AS_IF([test "x$with_opencl" == "xyes"],
[
   AC_MSG_RESULT([yes])
   AS_IF([test "x$with_opencl_dir" != "x"],
   [
      dnl Path manually given
      AC_MSG_CHECKING([cl.h])
      AS_IF([test -f ${with_opencl_dir}/include/CL/cl.h],
      [
         AC_MSG_RESULT([${with_opencl_dir}/include/CL/cl.h])
         OPENCL_CFLAGS="-I${with_opencl_dir}/include"
         OPENCL_LIBS="-lOpenCL -L${with_opencl_dir}/lib64 -L${with_opencl_dir}/lib"
         cl_cl_h=true
         have_opencl=true
      ], [
         AS_IF([test -f ${with_opencl_dir}/include/cl.h],
         [
            AC_MSG_RESULT([${with_opencl_dir}/include/cl.h])
            OPENCL_CFLAGS="-I${with_opencl_dir}/include/cl.h/include"
            OPENCL_LIBS="-lOpenCL -L${with_opencl_dir}/lib64 -L${with_opencl_dir}/lib"
            cl_cl_h=false
            have_opencl=true
         ], [
            AC_MSG_RESULT([not found.])
            AC_MSG_ERROR([Could not find cl.h])
         ])
      ])

   ], [

      dnl Attempt to automatically find the headers
      AC_MSG_CHECKING([cl.h])
      AS_IF([test -f /System/Library/Frameworks/OpenCL.framework/Headers/cl.h],
      [
         AC_MSG_RESULT([/System/Library/Frameworks/OpenCL.framework/Headers/cl.h])
         OPENCL_CFLAGS="-I/System/Library/Frameworks/OpenCL.framework/Headers"
         OPENCL_LIBS="-framework OpenCL"
      ], [
         AS_IF([test -f /opt/cuda/cuda/include/CL/cl.h],
         [
            AC_MSG_RESULT([/opt/cuda/cuda/include/CL/cl.h])
            OPENCL_CFLAGS="-I/opt/cuda/cuda/include"
            OPENCL_LIBS="-lOpenCL"
            cl_cl_h=true
         ], [
            AS_IF([test -f /opt/cuda/include/CL/cl.h],
            [
               AC_MSG_RESULT([/opt/cuda/include/CL/cl.h])
               OPENCL_CFLAGS="-I/opt/cuda/include"
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
      have_opencl=true
   ])
], [
   AC_MSG_RESULT([no])
])
AS_IF([test x$cl_cl_h = xtrue],
[
   AC_DEFINE([HAVE_CL_CL_H], [1], [Define to 1 if CL/cl.h should be used])
])
AS_IF([test x$have_opencl = xtrue],
[
   AC_DEFINE([HAVE_OPENCL], [1], [Define to 1 if OpenCL is available])
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


AC_MSG_CHECKING([whether to use Pango])
AC_ARG_ENABLE([pango],
              AS_HELP_STRING([--disable-pango], [Disable the use of Pango]))
AS_IF([test "x$enable_pango" != "xno"],
[
   AC_MSG_RESULT([yes])
   PKG_CHECK_MODULES([Pango], [pango >= 1.0.0],
   [
      have_pango=true
      AC_DEFINE([HAVE_PANGO], [1], [Define to 1 if Pango is available])
   ], [
      have_pango=false
   ])
   PKG_CHECK_MODULES([PangoCairo], [pangocairo],
   [
      have_pangocairo=true
      AC_DEFINE([HAVE_PANGOCAIRO], [1],
                [Define to 1 if PangoCairo is available])
   ], [
      have_pangocairo=false
   ])
], [
   AC_MSG_RESULT([no])
])


have_gdk_pixbuf=false
AC_MSG_CHECKING([whether to use GDK-pixbuf])
AC_ARG_ENABLE([gdk-pixbuf], AS_HELP_STRING([--disable-gdk-pixbuf],
                                           [Disable the use of gdk-pixbuf]))
AS_IF([test "x$enable_gdk_pixbuf" != "xno"],
[
   AC_MSG_RESULT([yes])

   PKG_CHECK_MODULES([GDK_pixbuf_2], [gdk-pixbuf-2.0],
   [
      have_gdk_pixbuf=true
   ], [
      PKG_CHECK_MODULES([GDK_pixbuf], [gdk-pixbuf], [],
      [
         have_gdk_pixbuf=true
      ])
   ])
], [
   AC_MSG_RESULT([no])
])

AS_IF([test "$have_gdk_pixbuf" != "false"],
[
  AC_DEFINE([HAVE_GDK_PIXBUF], [1], [Define to 1 if gdk-pixbuf is available])
])


AC_CHECK_LIB([pthread], [sched_setaffinity], [
   AC_DEFINE([HAVE_CPU_AFFINITY], [1],
             [Define to 1 if sched_setaffinity is available.])
])


LIBS_SAVE=$LIBS
LIBS=$HDF5_LIBS
AC_CHECK_LIB([hdf5], [H5Fcreate], [],
[
   AC_MSG_ERROR([HDF5 not found!])
])
LIBS=$LIBS_SAVE


AC_CHECK_LIB([rt], [clock_gettime],
[
   AC_DEFINE([HAVE_CLOCK_GETTIME], [1],
             [Define to 1 if clock_gettime is available.])
   TIMER_LIBS="-lrt"
], [
   AC_MSG_WARN([Falling back on gettimeofday().])
   AC_MSG_WARN([Time intervals may not be correct in certain circumstances.])
])


PKG_CHECK_MODULES([FFTW], [fftw3],
[
   have_fftw=true
   AC_DEFINE([HAVE_FFTW], [1], [Define to 1 if FFTW is available])
], [
   AC_MSG_WARN([ReAx indexing will not be available.])
   have_fftw=false
])

AX_WITH_CURSES

dnl Conditionals...
AM_CONDITIONAL([HAVE_OPENCL], test x$have_opencl = xtrue)

AM_CONDITIONAL([BUILD_HDFSEE], test x$have_gdk_pixbuf = xtrue \
                                -a "x$have_gtk" = "xtrue" \
                                -a "x$have_cairo" = "xtrue")

AM_CONDITIONAL([BUILD_EXPLORER], test "x$have_gtk" = "xtrue" \
                                  -a "x$have_cairo" = "xtrue")

AM_CONDITIONAL([HAVE_CAIRO], test x$have_cairo = xtrue)

AM_CONDITIONAL([HAVE_PANGOCAIRO], test x$have_pangocairo = xtrue)

AM_CONDITIONAL([HAVE_PANGO], test x$have_pango = xtrue)

AM_CONDITIONAL([HAVE_FFTW], test x$have_fftw = xtrue)

gl_IGNORE_UNUSED_LIBRARIES


MAIN_CFLAGS="$CFLAGS $HDF5_CFLAGS $GSL_CFLAGS $OPENCL_CFLAGS -pthread"
MAIN_CFLAGS="$MAIN_CFLAGS $GTK_CFLAGS $LIBTIFF_CFLAGS -D_GNU_SOURCE "
MAIN_CFLAGS="$MAIN_CFLAGS $libPNG_CFLAGS $Cairo_CFLAGS $GDK_pixbuf_CFLAGS"
MAIN_CFLAGS="$MAIN_CFLAGS $GDK_pixbuf_2_CFLAGS $FFTW_CFLAGS $Pango_CFLAGS"
MAIN_CFLAGS="$MAIN_CFLAGS $PangoCairo_CFLAGS"
AC_SUBST([MAIN_CFLAGS])

LIBCRYSTFEL_CFLAGS="$CFLAGS $HDF5_CFLAGS $GSL_CFLAGS $FFTW_CFLAGS"
AC_SUBST([LIBCRYSTFEL_CFLAGS])

MAIN_LIBS="$LIBS $HDF5_LIBS -lm -lz $GSL_LIBS $GTK_LIBS $OPENCL_LIBS -pthread"
MAIN_LIBS="$MAIN_LIBS $LIBTIFF_LIBS $libPNG_LIBS $Cairo_LIBS $GDK_pixbuf_LIBS"
MAIN_LIBS="$MAIN_LIBS $GDK_pixbuf_2_LIBS $TIMER_LIBS $FFTW_LIBS"
MAIN_LIBS="$MAIN_LIBS $Pango_LIBS $PangoCairo_LIBS $LDFLAGS"
AC_SUBST([MAIN_LIBS])

LIBCRYSTFEL_LIBS="$LIBS $HDF5_LIBS $GSL_LIBS $FFTW_LIBS $CURSES_LIB $LDFLAGS"
AC_SUBST([LIBCRYSTFEL_LIBS])

GTK_DOC_CHECK([1.9],[--flavour no-tmpl])
AC_CONFIG_MACRO_DIR(m4)

HTML_DIR=${docdir}/reference/libcrystfel


AC_CONFIG_FILES(Makefile lib/Makefile doc/reference/libcrystfel/Makefile libcrystfel/Makefile)
AC_CONFIG_FILES(libcrystfel/crystfel.pc)
AC_OUTPUT
