# libcrystfel

zlibdep = dependency('zlib', required : true)
if cc.has_function('gzbuffer',
                   prefix: '#include <zlib.h>',
                   dependencies: zlibdep)
  conf_data.set10('HAVE_GZBUFFER', 1)
endif

fftwdep = dependency('fftw3', required : false)
if fftwdep.found()
  conf_data.set10('HAVE_FFTW', 1)
endif

ncursesdep = dependency('ncurses', required: false)
if ncursesdep.found()
  conf_data.set10('HAVE_CURSES', 1)
endif

xgandalfdep = dependency('xgandalf', required: false)
if xgandalfdep.found()
  conf_data.set10('HAVE_XGANDALF', 1)
endif

pinkindexerdep = dependency('pinkIndexer', required: false)
if pinkindexerdep.found()
  conf_data.set10('HAVE_PINKINDEXER', 1)
endif

fdipdep = dependency('fdip', required: false)
if fdipdep.found()
  conf_data.set10('HAVE_FDIP', 1)
endif

libcrystfel_versionc = vcs_tag(input: 'src/libcrystfel-version.c.in',
                               output: 'libcrystfel-version.c')


libcrystfel_includes = include_directories('src')

# Find forkpty()
utildep = cc.find_library('util', required : true)
if cc.has_function('forkpty', dependencies : utildep, prefix : '#include <pty.h>')
  conf_data.set10('HAVE_FORKPTY_PTY_H', 1)
elif cc.has_function('forkpty', dependencies : utildep, prefix : '#include <util.h>')
  conf_data.set10('HAVE_FORKPTY_UTIL_H', 1)
else
  error('Couldn\'t find forkpty()')
endif


# Symmetry operation parser Flex/Bison stuff
flex = find_program('flex')
bison = find_program('bison')

flex_gen = generator(flex,
                     output : ['@BASENAME@-lex.c', '@BASENAME@-lex.h'],
                     arguments : ['--outfile=@OUTPUT0@',
                                  '--header-file=@OUTPUT1@',
                                  '@INPUT@'])

bison_gen = generator(bison,
                      output : ['@BASENAME@-parse.c', '@BASENAME@-parse.h'],
                      arguments : ['--output=@OUTPUT0@',
                                   '--defines=@OUTPUT1@',
                                   '--report=all',
                                   '@INPUT@'])

symop_parse_ch = bison_gen.process('src/symop.y')
symop_lex_ch = flex_gen.process('src/symop.l')


libcrystfel_sources = ['src/image.c',
                       'src/cell.c',
                       'src/index.c',
                       'src/spectrum.c',
                       'src/cell-utils.c',
                       'src/integer_matrix.c',
                       'src/stream.c',
                       'src/crystal.c',
                       'src/integration.c',
                       'src/symmetry.c',
                       'src/peakfinder8.c',
                       'src/thread-pool.c',
                       'src/peaks.c',
                       'src/utils.c',
                       'src/predict-refine.c',
                       'src/filters.c',
                       'src/rational.c',
                       'src/geometry.c',
                       'src/reflist.c',
                       'src/reflist-utils.c',
                       'src/datatemplate.c',
                       'src/colscale.c',
                       'src/detgeom.c',
                       'src/fom.c',
                       'src/image-cbf.c',
                       'src/image-hdf5.c',
                       'src/indexers/dirax.c',
                       'src/indexers/felix.c',
                       'src/indexers/mosflm.c',
                       'src/indexers/taketwo.c',
                       'src/indexers/xds.c',
                       'src/indexers/asdf.c',
                       'src/indexers/xgandalf.c',
                       'src/indexers/pinkindexer.c',
                       symop_lex_ch,
                       symop_parse_ch,
                      ]

if msgpackdep.found()
  libcrystfel_sources += 'src/image-msgpack.c'
endif

configure_file(input : 'config.h.in',
               output : 'config.h',
               configuration : conf_data)

libcrystfel_conf_inc = include_directories('.')

libcrystfel = library('crystfel', [libcrystfel_sources, libcrystfel_versionc],
                      include_directories : [libcrystfel_includes, libcrystfel_conf_inc],
                      dependencies : [mdep, utildep, fftwdep, gsldep, zlibdep,
                                      hdf5dep, pthreaddep, ncursesdep,
                                      xgandalfdep, pinkindexerdep, fdipdep],
                      install : true)

libcrystfeldep = declare_dependency(include_directories : libcrystfel_includes,
                                    link_with : libcrystfel,
                                    dependencies: gsldep)


install_headers(['src/reflist.h',
                 'src/symmetry.h',
                 'src/cell.h',
                 'src/reflist-utils.h',
                 'src/thread-pool.h',
                 'src/utils.h',
                 'src/geometry.h',
                 'src/peakfinder8.h',
                 'src/peaks.h',
                 'src/stream.h',
                 'src/index.h',
                 'src/image.h',
                 'src/filters.h',
                 'src/cell-utils.h',
                 'src/integer_matrix.h',
                 'src/crystal.h',
                 'src/predict-refine.h',
                 'src/integration.h',
                 'src/peakfinder8.h',
                 'src/rational.h',
                 'src/spectrum.h',
                 'src/datatemplate.h',
                 'src/colscale.h',
                 'src/detgeom.h',
                 'src/image-msgpack.h',
                 'src/fom.h'],
                subdir: 'crystfel')

# API documentation (Doxygen)
doxygen = find_program('doxygen')

index_md = files('doc/index.md')
coding_md = files('doc/coding.md')

doc_data = configuration_data()
doc_data.set('DOX_TOP', join_paths(meson.current_source_dir(), 'src'))
doc_data.set('INDEX_MD', join_paths(meson.current_source_dir(), 'doc/index.md'))
doc_data.set('CODING_MD', join_paths(meson.current_source_dir(), 'doc/coding.md'))
doc_data.set('OUTPUT_DIR', join_paths(meson.current_build_dir(), 'docs'))
doc_data.set('VERSION', meson.project_version())

doxyfile = configure_file(input: 'doc/Doxyfile.in',
                          output: 'Doxyfile',
                          configuration: doc_data,
                          install: false)

api_docs = run_target('api-docs', command : [doxygen, doxyfile])


# pkg-config file
pkg = import('pkgconfig')
pkg.generate(libcrystfel,
             description: 'Data processing for serial crystallography (shared library)')
