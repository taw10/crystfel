lib_LTLIBRARIES        = libcrystfel.la
libcrystfel_la_SOURCES = src/reflist.c src/utils.c src/cell.c src/detector.c \
                         src/thread-pool.c src/image.c src/hdf5-file.c \
                         src/beam-parameters.c src/geometry.c src/statistics.c \
                         src/symmetry.c src/stream.c src/peaks.c \
                         src/reflist-utils.c src/filters.c src/diffraction.c \
                         src/diffraction-gpu.c src/cl-utils.c

INCLUDES = "-I$(top_srcdir)/data"
AM_CPPFLAGS = -DDATADIR=\""$(datadir)"\" -I$(top_builddir)/lib
AM_CPPFLAGS += -I$(top_srcdir)/lib
