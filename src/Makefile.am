bin_PROGRAMS = pattern_sim process_hkl get_hkl indexamajig compare_hkl \
               powder_plot render_hkl calibrate_detector facetron cubeit \
               reintegrate estimate_background check_hkl

if HAVE_GTK
bin_PROGRAMS += hdfsee
endif

AM_CFLAGS = -Wall
AM_CPPFLAGS = -DDATADIR=\""$(datadir)"\" -I$(top_builddir)/lib -I$(top_srcdir)/lib
LDADD = lib/libgnu.a
ACLOCAL_AMFLAGS = -I m4

pattern_sim_SOURCES = pattern_sim.c diffraction.c utils.c image.c cell.c \
                      hdf5-file.c detector.c sfac.c peaks.c reflections.c \
                      beam-parameters.c symmetry.c
if HAVE_OPENCL
pattern_sim_SOURCES += diffraction-gpu.c cl-utils.c
endif

pattern_sim_LDADD = @LIBS@

process_hkl_SOURCES = process_hkl.c sfac.c statistics.c cell.c utils.c \
                      reflections.c symmetry.c stream.c beam-parameters.c
process_hkl_LDADD = @LIBS@

indexamajig_SOURCES = indexamajig.c hdf5-file.c utils.c cell.c image.c \
                      peaks.c index.c filters.c diffraction.c detector.c \
                      sfac.c dirax.c mosflm.c reflections.c templates.c \
                      symmetry.c geometry.c thread-pool.c beam-parameters.c
indexamajig_LDADD = @LIBS@
if HAVE_OPENCL
indexamajig_SOURCES += diffraction-gpu.c cl-utils.c
endif


if HAVE_GTK
hdfsee_SOURCES = hdfsee.c displaywindow.c render.c hdf5-file.c utils.c image.c \
                 filters.c
hdfsee_LDADD = @LIBS@
endif

get_hkl_SOURCES = get_hkl.c sfac.c cell.c utils.c reflections.c symmetry.c \
                  beam-parameters.c
get_hkl_LDADD = @LIBS@

compare_hkl_SOURCES = compare_hkl.c sfac.c cell.c utils.c reflections.c \
                      statistics.c symmetry.c
compare_hkl_LDADD = @LIBS@

check_hkl_SOURCES = check_hkl.c sfac.c cell.c utils.c reflections.c \
                    statistics.c symmetry.c
check_hkl_LDADD = @LIBS@

powder_plot_SOURCES = powder_plot.c cell.c utils.c image.c hdf5-file.c \
                      detector.c
powder_plot_LDADD = @LIBS@

render_hkl_SOURCES = render_hkl.c cell.c reflections.c utils.c povray.c \
                     symmetry.c render.c hdf5-file.c image.c filters.c
render_hkl_LDADD = @LIBS@

calibrate_detector_SOURCES = calibrate_detector.c utils.c hdf5-file.c image.c \
                             filters.c peaks.c detector.c cell.c thread-pool.c
calibrate_detector_LDADD = @LIBS@

facetron_SOURCES = facetron.c cell.c hdf5-file.c utils.c detector.c peaks.c \
                   image.c geometry.c reflections.c stream.c thread-pool.c \
                   beam-parameters.c symmetry.c post-refinement.c hrs-scaling.c
facetron_LDADD = @LIBS@

if HAVE_CAIRO
cubeit_SOURCES = cubeit.c cell.c hdf5-file.c utils.c detector.c render.c \
                 filters.c image.c symmetry.c stream.c thread-pool.c
cubeit_LDADD = @LIBS@
endif

reintegrate_SOURCES = reintegrate.c cell.c hdf5-file.c utils.c detector.c \
                      peaks.c image.c stream.c \
                      index.c dirax.c mosflm.c templates.c geometry.c \
                      symmetry.c thread-pool.c
reintegrate_LDADD = @LIBS@

estimate_background_SOURCES = estimate_background.c stream.c utils.c cell.c
estimate_background_LDADD = @LIBS@

INCLUDES = "-I$(top_srcdir)/data"

EXTRA_DIST = cell.h hdf5-file.h image.h  utils.h diffraction.h detector.h \
             sfac.h reflections.h list_tmp.h statistics.h displaywindow.h \
             render.h hdfsee.h dirax.h mosflm.h peaks.h index.h  filters.h \
             diffraction-gpu.h cl-utils.h symmetry.h \
             povray.h index-priv.h geometry.h templates.h render_hkl.h \
             stream.h thread-pool.h beam-parameters.h post-refinement.h \
             hrs-scaling.h
