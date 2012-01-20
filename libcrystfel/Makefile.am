lib_LTLIBRARIES = libcrystfel.la
libcrystfel_la_LIBADD = @LIBCRYSTFEL_LIBS@

libcrystfel_la_SOURCES = src/reflist.c src/utils.c src/cell.c src/detector.c \
                         src/thread-pool.c src/image.c src/hdf5-file.c \
                         src/beam-parameters.c src/geometry.c src/statistics.c \
                         src/symmetry.c src/stream.c src/peaks.c \
                         src/reflist-utils.c src/filters.c \
                         src/render.c src/index.c src/dirax.c src/mosflm.c \
                         src/reax.c

libcrystfel_la_includedir=$(includedir)/crystfel/

libcrystfel_la_include_HEADERS = src/beam-parameters.h src/hdf5-file.h \
                                 src/reflist.h src/symmetry.h src/cell.h \
                                 src/reflist-utils.h src/thread-pool.h \
                                 src/statistics.h src/utils.h src/detector.h \
                                 src/geometry.h src/peaks.h src/stream.h \
                                 src/render.h src/index.h src/image.h \
                                 src/filters.h

INCLUDES = "-I$(top_srcdir)/data"
AM_CPPFLAGS = -DDATADIR=\""$(datadir)"\" -I$(top_builddir)/lib
AM_CPPFLAGS += -I$(top_srcdir)/lib @LIBCRYSTFEL_CFLAGS@
