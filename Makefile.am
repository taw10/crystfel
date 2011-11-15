EXTRA_DIST = configure m4/gnulib-cache.m4
SUBDIRS = lib doc/reference libcrystfel
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = src/pattern_sim src/process_hkl src/get_hkl src/indexamajig \
               src/compare_hkl src/powder_plot src/render_hkl \
               src/calibrate_detector src/partialator \
               src/check_hkl src/partial_sim

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

AM_CFLAGS = -Wall
AM_CPPFLAGS = -DDATADIR=\""$(datadir)"\" -I$(top_builddir)/lib
AM_CPPFLAGS += -I$(top_srcdir)/lib
LDADD = $(top_builddir)/lib/libgnu.a @IGNORE_UNUSED_LIBRARIES_CFLAGS@
LDADD += $(top_builddir)/libcrystfel/libcrystfel.la

src_partial_sim_SOURCES = src/partial_sim.c

src_pattern_sim_SOURCES = src/pattern_sim.c
if HAVE_OPENCL
TESTS += tests/gpu_sim_check
noinst_PROGRAMS += tests/gpu_sim_check

tests_gpu_sim_check_SOURCES = tests/gpu_sim_check.c
endif

src_process_hkl_SOURCES = src/process_hkl.c

src_indexamajig_SOURCES = src/indexamajig.c

if BUILD_HDFSEE
src_hdfsee_SOURCES = src/hdfsee.c src/dw-hdfsee.c
endif

src_get_hkl_SOURCES = src/get_hkl.c

src_compare_hkl_SOURCES = src/compare_hkl.c

src_check_hkl_SOURCES = src/check_hkl.c

src_powder_plot_SOURCES = src/powder_plot.c

src_render_hkl_SOURCES = src/render_hkl.c src/povray.c

src_calibrate_detector_SOURCES = src/calibrate_detector.c

src_partialator_SOURCES = src/partialator.c src/post-refinement.c \
                          src/hrs-scaling.c

if HAVE_CAIRO
src_partialator_SOURCES += src/scaling-report.c
endif

tests_list_check_SOURCES = tests/list_check.c

tests_integration_check_SOURCES = tests/integration_check.c

tests_symmetry_check_SOURCES = tests/symmetry_check.c


tests_pr_gradient_check_SOURCES = tests/pr_gradient_check.c \
                                  src/post-refinement.c

INCLUDES = -I$(top_srcdir)/libcrystfel/src -I$(top_srcdir)/data

EXTRA_DIST += src/dw-hdfsee.h src/hdfsee.h src/povray.h \
              src/render_hkl.h src/post-refinement.h src/hrs-scaling.h \
              src/scaling-report.h

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
