EXTRA_DIST = configure m4/gnulib-cache.m4
SUBDIRS = lib
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = src/pattern_sim src/process_hkl src/get_hkl src/indexamajig \
               src/compare_hkl src/powder_plot src/render_hkl \
               src/calibrate_detector src/partialator src/reintegrate \
               src/estimate_background src/check_hkl

if HAVE_GTK
bin_PROGRAMS += src/hdfsee
endif

if HAVE_CAIRO
bin_PROGRAMS += src/cubeit
endif

AM_CFLAGS = -Wall
AM_CPPFLAGS = -DDATADIR=\""$(datadir)"\" -I$(top_builddir)/lib -I$(top_srcdir)/lib
LDADD = $(top_builddir)/lib/libgnu.a @IGNORE_UNUSED_LIBRARIES_CFLAGS@

src_pattern_sim_SOURCES = src/pattern_sim.c src/diffraction.c src/utils.c \
                          src/image.c src/cell.c src/hdf5-file.c \
                          src/detector.c src/sfac.c src/peaks.c \
                          src/reflections.c src/beam-parameters.c \
                          src/symmetry.c src/thread-pool.c src/reflist.c
if HAVE_OPENCL
src_pattern_sim_SOURCES += src/diffraction-gpu.c src/cl-utils.c
endif

src_process_hkl_SOURCES = src/process_hkl.c src/sfac.c src/statistics.c \
                          src/cell.c src/utils.c src/reflections.c \
                          src/symmetry.c src/stream.c src/beam-parameters.c \
                          src/thread-pool.c

src_indexamajig_SOURCES = src/indexamajig.c src/hdf5-file.c src/utils.c \
                          src/cell.c src/image.c src/peaks.c src/index.c \
                          src/filters.c src/diffraction.c src/detector.c \
                          src/sfac.c src/dirax.c src/mosflm.c \
                          src/reflections.c src/templates.c src/symmetry.c \
                          src/geometry.c src/thread-pool.c \
                          src/beam-parameters.c src/reflist.c
if HAVE_OPENCL
src_indexamajig_SOURCES += src/diffraction-gpu.c src/cl-utils.c
endif


if HAVE_GTK
src_hdfsee_SOURCES = src/hdfsee.c src/displaywindow.c src/render.c \
                     src/hdf5-file.c src/utils.c src/image.c src/filters.c \
                     src/thread-pool.c
endif

src_get_hkl_SOURCES = src/get_hkl.c src/sfac.c src/cell.c src/utils.c \
                      src/reflections.c src/symmetry.c src/beam-parameters.c \
                      src/thread-pool.c

src_compare_hkl_SOURCES = src/compare_hkl.c src/sfac.c src/cell.c src/utils.c \
                          src/reflections.c src/statistics.c src/symmetry.c \
                          src/thread-pool.c

src_check_hkl_SOURCES = src/check_hkl.c src/sfac.c src/cell.c src/utils.c \
                        src/reflections.c src/statistics.c src/symmetry.c \
                        src/thread-pool.c

src_powder_plot_SOURCES = src/powder_plot.c src/cell.c src/utils.c src/image.c \
                          src/hdf5-file.c src/detector.c src/thread-pool.c

src_render_hkl_SOURCES = src/render_hkl.c src/cell.c src/reflections.c \
                         src/utils.c  src/povray.c src/symmetry.c src/render.c \
                         src/hdf5-file.c src/image.c src/filters.c \
                         src/thread-pool.c

src_calibrate_detector_SOURCES = src/calibrate_detector.c src/utils.c \
                                 src/hdf5-file.c src/image.c src/filters.c \
                                 src/peaks.c src/detector.c src/cell.c \
                                 src/thread-pool.c src/reflist.c

src_partialator_SOURCES = src/partialator.c src/cell.c src/hdf5-file.c \
                          src/utils.c src/detector.c src/peaks.c src/image.c \
                          src/geometry.c src/reflections.c src/stream.c \
                          src/thread-pool.c src/beam-parameters.c \
                          src/symmetry.c src/post-refinement.c \
                          src/hrs-scaling.c src/reflist.c

if HAVE_CAIRO
src_cubeit_SOURCES = src/cubeit.c src/cell.c src/hdf5-file.c src/utils.c \
                     src/detector.c src/render.c src/filters.c src/image.c \
                     src/symmetry.c src/stream.c src/thread-pool.c
endif

src_reintegrate_SOURCES = src/reintegrate.c src/cell.c src/hdf5-file.c \
                          src/utils.c src/detector.c src/peaks.c src/image.c \
                          src/stream.c src/index.c src/dirax.c src/mosflm.c \
                          src/templates.c src/geometry.c src/symmetry.c \
                          src/thread-pool.c src/reflist.c

src_estimate_background_SOURCES = src/estimate_background.c src/stream.c \
                                  src/utils.c src/cell.c src/thread-pool.c

INCLUDES = "-I$(top_srcdir)/data"

EXTRA_DIST += src/cell.h src/hdf5-file.h src/image.h src/utils.h \
              src/diffraction.h src/detector.h src/sfac.h src/reflections.h \
              src/list_tmp.h src/statistics.h src/displaywindow.h \
              src/render.h src/hdfsee.h src/dirax.h src/mosflm.h src/peaks.h \
              src/index.h src/filters.h src/diffraction-gpu.h src/cl-utils.h \
              src/symmetry.h src/povray.h src/index-priv.h src/geometry.h \
              src/templates.h src/render_hkl.h src/stream.h src/thread-pool.h \
              src/beam-parameters.h src/post-refinement.h src/hrs-scaling.h

hdfseedir = $(datadir)/hdfsee
hdfsee_DATA = data/displaywindow.ui

crystfeldir = $(datadir)/crystfel
crystfel_DATA = data/diffraction.cl data/defs.h data/sfac/Ca.nff \
                data/sfac/C.nff data/sfac/Fe.nff data/sfac/H.nff \
                data/sfac/Mg.nff data/sfac/N.nff data/sfac/O.nff \
                data/sfac/P.nff data/sfac/S.nff data/sfac/f0_WaasKirf.dat

EXTRA_DIST += data/sfac/Ca.nff data/sfac/C.nff data/sfac/Fe.nff \
              data/sfac/H.nff data/sfac/Mg.nff data/sfac/N.nff data/sfac/O.nff \
              data/sfac/P.nff data/sfac/S.nff data/sfac/f0_WaasKirf.dat \
              data/defs.h data/diffraction.cl data/displaywindow.ui

crystfeldocdir = $(docdir)
crystfeldoc_DATA = doc/geometry.txt doc/indexamajig.txt doc/pattern_sim.txt \
                   doc/process_hkl.txt doc/symmetry.txt \
                   doc/twin-calculator.pdf \
                   doc/0-INDEX doc/examples/lcls-dec.geom \
                   doc/examples/lcls-dec.beam \
                   doc/examples/lcls-june-r0013-r0128.geom \
                   doc/examples/lcls-june.beam \
                   doc/examples/lcls-xpp-estimate.geom \
                   doc/examples/lcls-xpp-estimate.beam \
                   doc/examples/simple.geom

EXTRA_DIST += doc/geometry.txt doc/indexamajig.txt doc/pattern_sim.txt \
              doc/process_hkl.txt doc/symmetry.txt doc/twin-calculator.pdf \
              doc/0-INDEX doc/examples/lcls-dec.geom \
              doc/examples/lcls-june-r0013-r0128.geom \
              doc/examples/lcls-xpp-estimate.geom doc/examples/simple.geom

scriptdir = $(docdir)/scripts
script_DATA = scripts/check-near-bragg scripts/double-hit \
              scripts/random-image scripts/stream-split scripts/wibbletron \
              scripts/zone-axis scripts/create-mtz scripts/frequency \
              scripts/i0-analysis scripts/mtz2hkl scripts/sequence-image \
              scripts/test-facetron scripts/zone-axes \
              scripts/alternate-stream scripts/README scripts/indexing-rate \
              scripts/observed-peak-rate scripts/cell-please \
              scripts/check-hkl-completeness.gp scripts/check-hkl-quality.gp \
              scripts/compare-hkl.gp

EXTRA_DIST += scripts/check-near-bragg scripts/double-hit scripts/random-image \
              scripts/stream-split scripts/wibbletron scripts/zone-axis \
              scripts/create-mtz scripts/frequency scripts/i0-analysis \
              scripts/mtz2hkl scripts/sequence-image scripts/test-facetron \
              scripts/zone-axes scripts/alternate-stream scripts/README \
              scripts/indexing-rate scripts/observed-peak-rate \
              scripts/cell-please scripts/check-hkl-completeness.gp \
              scripts/check-hkl-quality.gp scripts/compare-hkl.gp
