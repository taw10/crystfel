EXTRA_DIST = configure m4/gnulib-cache.m4
SUBDIRS = lib doc/reference
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = src/pattern_sim src/process_hkl src/get_hkl src/indexamajig \
               src/compare_hkl src/powder_plot src/render_hkl \
               src/calibrate_detector src/partialator \
               src/check_hkl src/sum_stack src/partial_sim \
               contrib/alter_stream

noinst_PROGRAMS = tests/list_check tests/integration_check \
                  tests/pr_gradient_check tests/symmetry_check

TESTS = tests/list_check tests/first_merge_check tests/second_merge_check \
        tests/third_merge_check tests/fourth_merge_check \
        tests/integration_check tests/pr_gradient_check \
        tests/symmetry_check

EXTRA_DIST += tests/first_merge_check tests/second_merge_check \
              tests/third_merge_check tests/fourth_merge_check \
              tests/first_merge_check.hkl tests/third_merge_check.hkl \
              tests/fourth_merge_check.hkl

if BUILD_HDFSEE
bin_PROGRAMS += src/hdfsee
endif

if BUILD_CUBEIT
bin_PROGRAMS += contrib/cubeit
endif

AM_CFLAGS = -Wall
AM_CPPFLAGS = -DDATADIR=\""$(datadir)"\" -I$(top_builddir)/lib
AM_CPPFLAGS += -I$(top_srcdir)/lib
LDADD = $(top_builddir)/lib/libgnu.a @IGNORE_UNUSED_LIBRARIES_CFLAGS@

src_partial_sim_SOURCES = src/partial_sim.c src/cell.c src/detector.c \
                          src/beam-parameters.c src/thread-pool.c src/utils.c \
                          src/reflist-utils.c src/reflist.c src/symmetry.c \
                          src/hdf5-file.c src/image.c src/geometry.c \
                          src/peaks.c src/stream.c

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
                          src/beam-parameters.c src/reflist.c src/stream.c \
                          src/reax.c

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
                          src/reflist-utils.c src/beam-parameters.c src/stream.c \
                          src/reflist.c src/symmetry.c

src_render_hkl_SOURCES = src/render_hkl.c src/cell.c src/reflist-utils.c \
                         src/utils.c src/povray.c src/symmetry.c src/render.c \
                         src/hdf5-file.c src/image.c src/filters.c \
                         src/thread-pool.c src/reflist.c src/detector.c

src_sum_stack_SOURCES = src/sum_stack.c src/utils.c src/hdf5-file.c \
                        src/image.c src/filters.c src/peaks.c src/detector.c \
                        src/cell.c src/thread-pool.c src/reflist.c

src_calibrate_detector_SOURCES = src/calibrate_detector.c src/utils.c \
                                 src/hdf5-file.c src/image.c src/thread-pool.c \
                                 src/detector.c src/stream.c src/cell.c \
                                 src/reflist-utils.c src/reflist.c \
                                 src/symmetry.c src/peaks.c

src_partialator_SOURCES = src/partialator.c src/cell.c src/hdf5-file.c \
                          src/utils.c src/detector.c src/peaks.c src/image.c \
                          src/geometry.c src/reflist-utils.c src/stream.c \
                          src/thread-pool.c src/beam-parameters.c \
                          src/symmetry.c src/post-refinement.c \
                          src/hrs-scaling.c src/reflist.c

if HAVE_CAIRO
src_partialator_SOURCES += src/scaling-report.c
endif

if BUILD_CUBEIT
contrib_cubeit_SOURCES = contrib/cubeit.c src/cell.c src/hdf5-file.c \
                     src/utils.c src/detector.c src/render.c src/filters.c \
                     src/image.c src/symmetry.c src/stream.c src/thread-pool.c \
                     src/reflist.c src/reflist-utils.c
endif

contrib_alter_stream_SOURCES = contrib/alter_stream.c src/stream.c src/image.c \
                               src/thread-pool.c src/detector.c src/utils.c \
                               src/cell.c src/reflist-utils.c src/reflist.c \
                               src/symmetry.c src/hdf5-file.c \
                               src/beam-parameters.c src/geometry.c src/peaks.c

tests_list_check_SOURCES = tests/list_check.c src/reflist.c src/thread-pool.c \
                           src/utils.c

tests_integration_check_SOURCES = tests/integration_check.c src/peaks.c \
                                  src/reflist.c src/cell.c src/detector.c \
                                  src/thread-pool.c src/utils.c src/image.c \
                                  src/hdf5-file.c

tests_symmetry_check_SOURCES = tests/symmetry_check.c src/symmetry.c \
                               src/utils.c src/thread-pool.c


tests_pr_gradient_check_SOURCES = tests/pr_gradient_check.c src/detector.c \
                                  src/cell.c src/geometry.c src/reflist.c \
                                  src/thread-pool.c src/utils.c src/peaks.c \
                                  src/symmetry.c src/image.c src/hdf5-file.c \
                                  src/post-refinement.c

INCLUDES = "-I$(top_srcdir)/data"

EXTRA_DIST += src/cell.h src/hdf5-file.h src/image.h src/utils.h \
              src/diffraction.h src/detector.h src/reflist-utils.h \
              src/list_tmp.h src/statistics.h src/dw-hdfsee.h \
              src/render.h src/hdfsee.h src/dirax.h src/mosflm.h src/peaks.h \
              src/index.h src/filters.h src/diffraction-gpu.h src/cl-utils.h \
              src/symmetry.h src/povray.h src/index-priv.h src/geometry.h \
              src/render_hkl.h src/stream.h src/thread-pool.h \
              src/beam-parameters.h src/post-refinement.h src/hrs-scaling.h \
              src/reflist.h src/scaling-report.h src/reax.h

crystfeldir = $(datadir)/crystfel
crystfel_DATA = data/diffraction.cl data/defs.h data/hdfsee.ui

EXTRA_DIST += $(crystfel_DATA)

man_MANS = doc/man/indexamajig.1 doc/man/process_hkl.1 doc/man/pattern_sim.1 \
           doc/man/crystfel_geometry.1 doc/man/powder_plot.1

crystfeldocdir = $(docdir)
crystfeldoc_DATA = doc/twin-calculator.pdf doc/examples/lcls-dec.geom \
                   doc/examples/lcls-dec.beam \
                   doc/examples/lcls-june-r0013-r0128.geom \
                   doc/examples/lcls-june.beam \
                   doc/examples/lcls-xpp-estimate.geom \
                   doc/examples/lcls-xpp-estimate.beam \
                   doc/examples/simple.geom \
                   doc/examples/cspad-feb2011.geom \
                   doc/examples/lcls-cxi-9keV.beam

EXTRA_DIST += doc/twin-calculator.pdf doc/examples/lcls-dec.geom \
              doc/examples/lcls-june-r0013-r0128.geom \
              doc/examples/lcls-xpp-estimate.geom doc/examples/simple.geom \
              doc/examples/lcls-dec.beam doc/examples/lcls-june.beam \
              doc/examples/lcls-xpp-estimate.beam \
              doc/examples/cspad-feb2011.geom \
              doc/examples/lcls-cxi-9keV.beam

scriptdir = $(docdir)/scripts
script_DATA = scripts/check-near-bragg  scripts/random-image \
              scripts/wibbletron scripts/create-mtz scripts/mtz2hkl \
              scripts/sequence-image scripts/zone-axes \
              scripts/alternate-stream scripts/README scripts/cell-please \
              scripts/check-hkl-completeness.gp scripts/check-hkl-quality.gp \
              scripts/compare-hkl.gp scripts/check-peak-detection \
              scripts/indexed-filenames scripts/gen-sfs

EXTRA_DIST += scripts/check-near-bragg  scripts/random-image \
              scripts/wibbletron scripts/create-mtz scripts/mtz2hkl \
              scripts/sequence-image scripts/zone-axes \
              scripts/alternate-stream scripts/README scripts/cell-please \
              scripts/check-hkl-completeness.gp scripts/check-hkl-quality.gp \
              scripts/compare-hkl.gp scripts/check-peak-detection \
              scripts/indexed-filenames scripts/gen-sfs


DISTCHECK_CONFIGURE_FLAGS=--enable-gtk-doc
