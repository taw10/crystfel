bin_PROGRAMS = pattern_sim process_hkl get_hkl indexamajig compare_hkl \
               powder_plot

if HAVE_GTK
bin_PROGRAMS += hdfsee
endif

AM_CFLAGS = -Wall
AM_CPPFLAGS = -DDATADIR=\""$(datadir)"\"

pattern_sim_SOURCES = pattern_sim.c diffraction.c utils.c image.c cell.c \
                      hdf5-file.c detector.c sfac.c intensities.c \
                      reflections.c
if HAVE_OPENCL
pattern_sim_SOURCES += diffraction-gpu.c cl-utils.c
endif

pattern_sim_LDADD = @LIBS@

process_hkl_SOURCES = process_hkl.c sfac.c statistics.c cell.c utils.c \
                      reflections.c
process_hkl_LDADD = @LIBS@

indexamajig_SOURCES = indexamajig.c hdf5-file.c utils.c cell.c image.c \
                      intensities.c peaks.c index.c filters.c \
                      diffraction.c detector.c sfac.c
indexamajig_LDADD = @LIBS@
if HAVE_OPENCL
indexamajig_SOURCES += diffraction-gpu.c cl-utils.c
endif

if HAVE_GLIB
indexamajig_SOURCES += dirax.c
endif

if HAVE_GTK
hdfsee_SOURCES = hdfsee.c displaywindow.c render.c hdf5-file.c utils.c image.c \
                 filters.c
hdfsee_LDADD = @LIBS@
endif

get_hkl_SOURCES = get_hkl.c sfac.c cell.c utils.c reflections.c
get_hkl_LDADD = @LIBS@

compare_hkl_SOURCES = compare_hkl.c sfac.c cell.c utils.c reflections.c
compare_hkl_LDADD = @LIBS@

powder_plot_SOURCES = powder_plot.c cell.c utils.c image.c hdf5-file.c \
                      detector.c
powder_plot_LDADD = @LIBS@

if HAVE_GLIB
powder_plot_SOURCES += dirax.c
endif

INCLUDES = "-I$(top_srcdir)/data"
