docdir = ${datadir}/doc/${PACKAGE}
EXTRA_DIST = configure m4/gnulib-cache.m4 version.sh version.h.in
SUBDIRS = lib doc/reference/libcrystfel libcrystfel
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = src/pattern_sim src/process_hkl src/get_hkl src/indexamajig \
               src/compare_hkl src/partialator src/check_hkl src/partial_sim \
               src/ambigator src/geoptimiser src/whirligig src/list_events

noinst_PROGRAMS = tests/list_check tests/integration_check \
                  tests/pr_p_gradient_check tests/symmetry_check \
                  tests/centering_check tests/transformation_check \
                  tests/cell_check tests/ring_check \
                  tests/prof2d_check tests/ambi_check \
		  tests/prediction_gradient_check

MERGE_CHECKS = tests/first_merge_check tests/second_merge_check \
               tests/third_merge_check tests/fourth_merge_check

PARTIAL_CHECKS = tests/partialator_merge_check_1 \
                 tests/partialator_merge_check_2 \
                 tests/partialator_merge_check_3 \
                 tests/pr_p_gradient_check

TESTS = tests/list_check $(MERGE_CHECKS) $(PARTIAL_CHECKS) \
        tests/integration_check \
        tests/symmetry_check tests/centering_check tests/transformation_check \
        tests/cell_check tests/ring_check tests/prof2d_check tests/ambi_check \
        tests/prediction_gradient_check

EXTRA_DIST += $(MERGE_CHECKS) $(PARTIAL_CHECKS)
EXTRA_DIST += relnotes-0.6.0

if BUILD_HDFSEE
bin_PROGRAMS += src/hdfsee
endif

if BUILD_EXPLORER
bin_PROGRAMS += src/cell_explorer
endif

if HAVE_CAIRO
bin_PROGRAMS += src/render_hkl
endif

AM_CFLAGS = -Wall @MAIN_CFLAGS@
AM_CPPFLAGS = -DDATADIR=\""$(datadir)"\" -I$(top_builddir)/lib
AM_CPPFLAGS += -I$(top_srcdir)/lib
LDADD = $(top_builddir)/lib/libgnu.la @IGNORE_UNUSED_LIBRARIES_CFLAGS@
LDADD += $(top_builddir)/libcrystfel/libcrystfel.la @MAIN_LIBS@

tests_integration_check_LDADD = $(top_builddir)/lib/libgnu.la
tests_integration_check_LDADD += $(top_builddir)/libcrystfel/libcrystfel.la
tests_integration_check_LDADD += @LIBCRYSTFEL_LIBS@

tests_prof2d_check_LDADD = $(top_builddir)/lib/libgnu.la
tests_prof2d_check_LDADD += $(top_builddir)/libcrystfel/libcrystfel.la
tests_prof2d_check_LDADD += @LIBCRYSTFEL_LIBS@

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

src_list_events_SOURCES = src/list_events.c

src_indexamajig_SOURCES = src/indexamajig.c src/im-sandbox.c \
                          src/process_image.c src/predict-refine.c

if BUILD_HDFSEE
src_hdfsee_SOURCES = src/hdfsee.c src/dw-hdfsee.c src/hdfsee-render.c
endif

if BUILD_EXPLORER
src_cell_explorer_SOURCES = src/cell_explorer.c src/multihistogram.c
endif

src_get_hkl_SOURCES = src/get_hkl.c

src_compare_hkl_SOURCES = src/compare_hkl.c

src_check_hkl_SOURCES = src/check_hkl.c

if HAVE_CAIRO
src_render_hkl_SOURCES = src/render_hkl.c
endif

src_partialator_SOURCES = src/partialator.c src/post-refinement.c \
                          src/hrs-scaling.c src/rejection.c

src_ambigator_SOURCES = src/ambigator.c

src_geoptimiser_SOURCES = src/geoptimiser.c src/hdfsee-render.c

src_whirligig_SOURCES = src/whirligig.c

if HAVE_CAIRO
if HAVE_PANGO
if HAVE_PANGOCAIRO
src_partialator_SOURCES += src/scaling-report.c
endif
endif
endif

tests_list_check_SOURCES = tests/list_check.c

tests_integration_check_SOURCES = tests/integration_check.c

tests_prof2d_check_SOURCES = tests/prof2d_check.c

tests_symmetry_check_SOURCES = tests/symmetry_check.c

tests_ambi_check_SOURCES = tests/ambi_check.c

tests_pr_p_gradient_check_SOURCES = tests/pr_p_gradient_check.c \
                                  src/post-refinement.c

tests_prediction_gradient_check_SOURCES = tests/prediction_gradient_check.c

tests_centering_check_SOURCES = tests/centering_check.c

tests_transformation_check_SOURCES = tests/transformation_check.c

tests_ring_check_SOURCES = tests/ring_check.c

tests_cell_check_SOURCES = tests/cell_check.c

INCLUDES = -I$(top_srcdir)/libcrystfel/src -I$(top_srcdir)/data

EXTRA_DIST += src/dw-hdfsee.h src/hdfsee.h src/render_hkl.h \
              src/post-refinement.h src/hrs-scaling.h src/scaling-report.h \
              src/cl-utils.h src/hdfsee-render.h src/diffraction.h \
              src/diffraction-gpu.h src/pattern_sim.h src/list_tmp.h \
              src/im-sandbox.h src/process_image.h src/multihistogram.h \
              src/rejection.h src/predict-refine.h

crystfeldir = $(datadir)/crystfel
crystfel_DATA = data/diffraction.cl data/hdfsee.ui

EXTRA_DIST += $(crystfel_DATA)

man_MANS = doc/man/crystfel.7 doc/man/check_hkl.1 doc/man/compare_hkl.1 \
           doc/man/crystfel_geometry.5 doc/man/get_hkl.1 doc/man/hdfsee.1 \
           doc/man/indexamajig.1 doc/man/partialator.1 doc/man/partial_sim.1 \
           doc/man/pattern_sim.1 doc/man/process_hkl.1 doc/man/render_hkl.1 \
           doc/man/cell_explorer.1 doc/man/ambigator.1 doc/man/geoptimiser.1 \
           doc/man/whirligig.1


EXTRA_DIST += $(man_MANS)

crystfeldocdir = $(docdir)
crystfeldoc_DATA = doc/twin-calculator.pdf doc/examples/lcls-dec.geom \
                   doc/examples/lcls-june-r0013-r0128.geom \
                   doc/examples/lcls-xpp-estimate.geom \
                   doc/examples/simple.geom \
                   doc/examples/cspad-single.geom \
                   doc/examples/cspad-cxiformat.geom \
                   doc/examples/pilatus.geom \
                   doc/examples/cell-example.cell \
                   doc/hitrate.html doc/hitrate.png

EXTRA_DIST += $(crystfeldoc_DATA) doc/twin-calculator.odt \
              doc/reference/libcrystfel/xml/coding-standards.xml

scriptdir = $(docdir)/scripts
script_DATA = scripts/alternate-stream scripts/cell-please \
              scripts/check-hkl-completeness.gp scripts/check-hkl-quality.gp \
              scripts/check-near-bragg scripts/check-peak-detection \
              scripts/cif2hkl scripts/compare-hkl.gp scripts/create-mtz \
              scripts/create-xscale scripts/eV-to-A scripts/find-filename \
              scripts/gen-sfs scripts/gen-sfs-ano scripts/hkl2hkl \
              scripts/indexed-filenames scripts/mtz2hkl scripts/peak-intensity \
              scripts/random-image scripts/README scripts/sequence-image \
              scripts/split-indexed scripts/stream_grep scripts/zone-axes \
              scripts/Rsplit_surface scripts/Rsplit_surface.py \
              scripts/clean-stream.py scripts/fg-graph scripts/truncate-stream \
              scripts/gen-sfs-expand scripts/add-beam-params \
              scripts/find-pairs scripts/plot-cc-and-scale.R \
              scripts/ave-resolution scripts/crystal-frame-number \
              scripts/plot-radius-resolution scripts/plot-predict-refine \
              scripts/detector-shift

EXTRA_DIST += $(script_DATA)


DISTCHECK_CONFIGURE_FLAGS=--enable-gtk-doc

mostlyclean-local: mostlyclean-version
.PHONY: mostlyclean-version
mostlyclean-version:
	rm -f $(top_builddir)/version.h

all: update_version_header
.PHONY: version_header
update_version_header:
	$(top_srcdir)/version.sh $(top_srcdir)
