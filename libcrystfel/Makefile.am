lib_LTLIBRARIES = libcrystfel.la
libcrystfel_la_LIBADD = ../lib/libgnu.la @LIBCRYSTFEL_LIBS@ $(PTY_LIB)
libcrystfel_la_LDFLAGS = -version-info 5:0:0

libcrystfel_la_SOURCES = src/reflist.c src/utils.c src/cell.c src/detector.c \
                         src/thread-pool.c src/image.c src/hdf5-file.c \
                         src/geometry.c src/statistics.c \
                         src/symmetry.c src/stream.c src/peaks.c \
                         src/reflist-utils.c src/filters.c \
                         src/render.c src/index.c src/dirax.c src/mosflm.c \
                         src/cell-utils.c src/integer_matrix.c src/crystal.c \
                         src/grainspotter.c src/xds.c src/integration.c \
                         src/histogram.c src/events.c

if HAVE_FFTW
libcrystfel_la_SOURCES += src/reax.c
endif

libcrystfel_la_includedir=$(includedir)/crystfel/

libcrystfel_la_include_HEADERS = ${top_srcdir}/version.h \
                                 src/hdf5-file.h src/reflist.h src/symmetry.h \
                                 src/cell.h src/reflist-utils.h \
                                 src/thread-pool.h  src/statistics.h \
                                 src/utils.h src/detector.h \
                                 src/geometry.h src/peaks.h src/stream.h \
                                 src/render.h src/index.h src/image.h \
                                 src/filters.h src/dirax.h src/mosflm.h \
                                 src/reax.h src/cell-utils.h \
                                 src/integer_matrix.h src/crystal.h \
                                 src/grainspotter.h src/xds.h \
                                 src/integration.h src/histogram.h \
                                 src/events.h

AM_CPPFLAGS = -DDATADIR=\""$(datadir)"\" -I$(top_builddir)/lib -Wall
AM_CPPFLAGS += -I$(top_srcdir)/lib @LIBCRYSTFEL_CFLAGS@

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = crystfel.pc
