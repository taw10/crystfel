EXTRA_DIST = configure m4/gnulib-cache.m4
SUBDIRS = lib doc/reference
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = src/pattern_sim src/process_hkl src/get_hkl src/indexamajig \
               src/compare_hkl src/powder_plot src/render_hkl \
               src/calibrate_detector src/partialator \
               src/check_hkl src/sum_stack

noinst_PROGRAMS = tests/list_check

TESTS = tests/list_check

if BUILD_HDFSEE
bin_PROGRAMS += src/hdfsee
endif

if BUILD_CUBEIT
bin_PROGRAMS += src/cubeit
endif

AM_CFLAGS = -Wall
AM_CPPFLAGS = -DDATADIR=\""$(datadir)"\" -I$(top_builddir)/lib -I$(top_srcdir)/lib
LDADD = $(top_builddir)/lib/libgnu.a @IGNORE_UNUSED_LIBRARIES_CFLAGS@

src_pattern_sim_SOURCES = src/pattern_sim.c src/diffraction.c src/utils.c \
                          src/image.c src/cell.c src/hdf5-file.c \
                          src/detector.c src/peaks.c \
                          src/reflist-utils.c src/beam-parameters.c \
                          src/symmetry.c src/thread-pool.c src/reflist.c
if HAVE_OPENCL
src_pattern_sim_SOURCES += src/diffraction-gpu.c src/cl-utils.c

TESTS += tests/gpu_sim_check
noinst_PROGRAMS += tests/gpu_sim_check

tests_gpu_sim_check_SOURCES = tests/gpu_sim_check.c src/utils.c \
                              src/diffraction.c src/diffraction-gpu.c \
                              src/thread-pool.c src/cell.c src/cl-utils.c \
                              src/symmetry.c src/detector.c src/hdf5-file.c \
                              src/image.c
endif

src_process_hkl_SOURCES = src/process_hkl.c src/statistics.c \
                          src/cell.c src/utils.c \
                          src/symmetry.c src/stream.c src/beam-parameters.c \
                          src/thread-pool.c src/image.c src/detector.c \
                          src/hdf5-file.c src/reflist.c src/reflist-utils.c

src_indexamajig_SOURCES = src/indexamajig.c src/hdf5-file.c src/utils.c \
                          src/cell.c src/image.c src/peaks.c src/index.c \
                          src/filters.c src/diffraction.c src/detector.c \
                          src/dirax.c src/mosflm.c \
                          src/reflist-utils.c src/symmetry.c \
                          src/geometry.c src/thread-pool.c \
                          src/beam-parameters.c src/reflist.c src/stream.c

if HAVE_OPENCL
src_indexamajig_SOURCES += src/diffraction-gpu.c src/cl-utils.c
endif


if BUILD_HDFSEE
src_hdfsee_SOURCES = src/hdfsee.c src/dw-hdfsee.c src/render.c \
                     src/hdf5-file.c src/utils.c src/image.c src/filters.c \
                     src/thread-pool.c src/detector.c
endif

src_get_hkl_SOURCES = src/get_hkl.c src/cell.c src/utils.c \
                      src/reflist-utils.c src/symmetry.c src/beam-parameters.c \
                      src/thread-pool.c src/reflist.c

src_compare_hkl_SOURCES = src/compare_hkl.c src/cell.c src/utils.c \
                          src/reflist-utils.c src/statistics.c src/symmetry.c \
                          src/thread-pool.c src/reflist.c

src_check_hkl_SOURCES = src/check_hkl.c src/cell.c src/utils.c \
                        src/reflist-utils.c src/statistics.c src/symmetry.c \
                        src/thread-pool.c src/reflist.c

src_powder_plot_SOURCES = src/powder_plot.c src/cell.c src/utils.c src/image.c \
                          src/hdf5-file.c src/detector.c src/thread-pool.c \
                          src/beam-parameters.c

src_render_hkl_SOURCES = src/render_hkl.c src/cell.c src/reflist-utils.c \
                         src/utils.c src/povray.c src/symmetry.c src/render.c \
                         src/hdf5-file.c src/image.c src/filters.c \
                         src/thread-pool.c src/reflist.c

src_sum_stack_SOURCES = src/sum_stack.c src/utils.c src/hdf5-file.c \
                        src/image.c src/filters.c src/peaks.c src/detector.c \
                        src/cell.c src/thread-pool.c src/reflist.c

src_calibrate_detector_SOURCES = src/calibrate_detector.c src/utils.c \
                                 src/hdf5-file.c src/image.c src/thread-pool.c \
                                 src/detector.c

src_partialator_SOURCES = src/partialator.c src/cell.c src/hdf5-file.c \
                          src/utils.c src/detector.c src/peaks.c src/image.c \
                          src/geometry.c src/reflist-utils.c src/stream.c \
                          src/thread-pool.c src/beam-parameters.c \
                          src/symmetry.c src/post-refinement.c \
                          src/hrs-scaling.c src/reflist.c

if BUILD_CUBEIT
src_cubeit_SOURCES = src/cubeit.c src/cell.c src/hdf5-file.c src/utils.c \
                     src/detector.c src/render.c src/filters.c src/image.c \
                     src/symmetry.c src/stream.c src/thread-pool.c \
                     src/reflist.c src/reflist-utils.c
endif

tests_list_check_SOURCES = tests/list_check.c src/reflist.c src/thread-pool.c \
                           src/utils.c

INCLUDES = "-I$(top_srcdir)/data"

EXTRA_DIST += src/cell.h src/hdf5-file.h src/image.h src/utils.h \
              src/diffraction.h src/detector.h src/reflist-utils.h \
              src/list_tmp.h src/statistics.h src/dw-hdfsee.h \
              src/render.h src/hdfsee.h src/dirax.h src/mosflm.h src/peaks.h \
              src/index.h src/filters.h src/diffraction-gpu.h src/cl-utils.h \
              src/symmetry.h src/povray.h src/index-priv.h src/geometry.h \
              src/render_hkl.h src/stream.h src/thread-pool.h \
              src/beam-parameters.h src/post-refinement.h src/hrs-scaling.h \
              src/reflist.h

crystfeldir = $(datadir)/crystfel
crystfel_DATA = data/diffraction.cl data/defs.h data/hdfsee.ui

EXTRA_DIST += $(crystfel_DATA)

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
                   doc/examples/simple.geom \
                   doc/examples/cspad-feb2011.geom \
                   doc/examples/lcls-cxi-9keV.beam

EXTRA_DIST += doc/geometry.txt doc/indexamajig.txt doc/pattern_sim.txt \
              doc/process_hkl.txt doc/symmetry.txt doc/twin-calculator.pdf \
              doc/0-INDEX doc/examples/lcls-dec.geom \
              doc/examples/lcls-june-r0013-r0128.geom \
              doc/examples/lcls-xpp-estimate.geom doc/examples/simple.geom \
              doc/examples/lcls-dec.beam doc/examples/lcls-june.beam \
              doc/examples/lcls-xpp-estimate.beam \
              doc/examples/cspad-feb2011.geom \
              doc/examples/lcls-cxi-9keV.beam

scriptdir = $(docdir)/scripts
script_DATA = scripts/check-near-bragg scripts/double-hit \
              scripts/random-image scripts/stream-split scripts/wibbletron \
              scripts/zone-axis scripts/create-mtz scripts/frequency \
              scripts/i0-analysis scripts/mtz2hkl scripts/sequence-image \
              scripts/test-facetron scripts/zone-axes \
              scripts/alternate-stream scripts/README scripts/indexing-rate \
              scripts/observed-peak-rate scripts/cell-please \
              scripts/check-hkl-completeness.gp scripts/check-hkl-quality.gp \
              scripts/compare-hkl.gp scripts/check-peak-detection

EXTRA_DIST += scripts/check-near-bragg scripts/double-hit scripts/random-image \
              scripts/stream-split scripts/wibbletron scripts/zone-axis \
              scripts/create-mtz scripts/frequency scripts/i0-analysis \
              scripts/mtz2hkl scripts/sequence-image scripts/test-facetron \
              scripts/zone-axes scripts/alternate-stream scripts/README \
              scripts/indexing-rate scripts/observed-peak-rate \
              scripts/cell-please scripts/check-hkl-completeness.gp \
              scripts/check-hkl-quality.gp scripts/compare-hkl.gp \
              scripts/check-peak-detection


DISTCHECK_CONFIGURE_FLAGS=--enable-gtk-doc
