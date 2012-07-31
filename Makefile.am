EXTRA_DIST = configure m4/gnulib-cache.m4
SUBDIRS = lib doc/reference libcrystfel
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = src/pattern_sim src/process_hkl src/get_hkl src/indexamajig \
               src/compare_hkl src/powder_plot \
               src/partialator src/check_hkl src/partial_sim

noinst_PROGRAMS = tests/list_check tests/integration_check \
                  tests/pr_gradient_check tests/symmetry_check

TESTS = tests/list_check tests/first_merge_check tests/second_merge_check \
        tests/third_merge_check tests/fourth_merge_check \
        tests/integration_check tests/pr_gradient_check \
        tests/symmetry_check

EXTRA_DIST += tests/first_merge_check tests/second_merge_check \
              tests/third_merge_check tests/fourth_merge_check

if BUILD_HDFSEE
bin_PROGRAMS += src/hdfsee
endif

if HAVE_CAIRO
bin_PROGRAMS += src/render_hkl
endif

AM_CFLAGS = -Wall @MAIN_CFLAGS@
AM_CPPFLAGS = -DDATADIR=\""$(datadir)"\" -I$(top_builddir)/lib
AM_CPPFLAGS += -I$(top_srcdir)/lib
LDADD = $(top_builddir)/lib/libgnu.la @IGNORE_UNUSED_LIBRARIES_CFLAGS@
LDADD += $(top_builddir)/libcrystfel/libcrystfel.la @MAIN_LIBS@

src_partial_sim_SOURCES = src/partial_sim.c

src_pattern_sim_SOURCES = src/pattern_sim.c src/diffraction.c

if HAVE_OPENCL
src_pattern_sim_SOURCES += src/cl-utils.c src/diffraction-gpu.c
TESTS += tests/gpu_sim_check
noinst_PROGRAMS += tests/gpu_sim_check
tests_gpu_sim_check_SOURCES = tests/gpu_sim_check.c src/diffraction.c \
                              src/diffraction-gpu.c src/cl-utils.c
endif

src_process_hkl_SOURCES = src/process_hkl.c

src_indexamajig_SOURCES = src/indexamajig.c src/im-sandbox.c

if BUILD_HDFSEE
src_hdfsee_SOURCES = src/hdfsee.c src/dw-hdfsee.c src/hdfsee-render.c
endif

src_get_hkl_SOURCES = src/get_hkl.c

src_compare_hkl_SOURCES = src/compare_hkl.c

src_check_hkl_SOURCES = src/check_hkl.c

src_powder_plot_SOURCES = src/powder_plot.c

if HAVE_CAIRO
src_render_hkl_SOURCES = src/render_hkl.c
endif

src_partialator_SOURCES = src/partialator.c src/post-refinement.c \
                          src/hrs-scaling.c

if HAVE_CAIRO
if HAVE_PANGO
if HAVE_PANGOCAIRO
src_partialator_SOURCES += src/scaling-report.c
endif
endif
endif

tests_list_check_SOURCES = tests/list_check.c

tests_integration_check_SOURCES = tests/integration_check.c

tests_symmetry_check_SOURCES = tests/symmetry_check.c


tests_pr_gradient_check_SOURCES = tests/pr_gradient_check.c \
                                  src/post-refinement.c

INCLUDES = -I$(top_srcdir)/libcrystfel/src -I$(top_srcdir)/data

EXTRA_DIST += src/dw-hdfsee.h src/hdfsee.h src/render_hkl.h \
              src/post-refinement.h src/hrs-scaling.h src/scaling-report.h \
              src/cl-utils.h src/hdfsee-render.h src/diffraction.h \
              src/diffraction-gpu.h src/pattern_sim.h src/list_tmp.h \
              src/im-sandbox.h

crystfeldir = $(datadir)/crystfel
crystfel_DATA = data/diffraction.cl data/hdfsee.ui

EXTRA_DIST += $(crystfel_DATA)

man_MANS = doc/man/crystfel.7 doc/man/check_hkl.1 doc/man/compare_hkl.1 \
           doc/man/crystfel_geometry.5 doc/man/get_hkl.1 doc/man/hdfsee.1 \
           doc/man/indexamajig.1 doc/man/partialator.1 doc/man/partial_sim.1 \
           doc/man/pattern_sim.1 doc/man/powder_plot.1 doc/man/process_hkl.1 \
           doc/man/render_hkl.1


EXTRA_DIST += $(man_MANS)

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
              doc/examples/lcls-cxi-9keV.beam \
              doc/reference/xml/coding-standards.xml

scriptdir = $(docdir)/scripts
script_DATA = scripts/check-near-bragg  scripts/random-image \
              scripts/wibbletron scripts/create-mtz scripts/mtz2hkl \
              scripts/sequence-image scripts/zone-axes \
              scripts/alternate-stream scripts/README scripts/cell-please \
              scripts/check-hkl-completeness.gp scripts/check-hkl-quality.gp \
              scripts/compare-hkl.gp scripts/check-peak-detection \
              scripts/indexed-filenames scripts/gen-sfs scripts/cif2hkl

EXTRA_DIST += scripts/check-near-bragg  scripts/random-image \
              scripts/wibbletron scripts/create-mtz scripts/mtz2hkl \
              scripts/sequence-image scripts/zone-axes \
              scripts/alternate-stream scripts/README scripts/cell-please \
              scripts/check-hkl-completeness.gp scripts/check-hkl-quality.gp \
              scripts/compare-hkl.gp scripts/check-peak-detection \
              scripts/indexed-filenames scripts/gen-sfs scripts/cif2hkl


DISTCHECK_CONFIGURE_FLAGS=--enable-gtk-doc
