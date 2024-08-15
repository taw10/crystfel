# libcrystfel

zlibdep = dependency('zlib', required: false)
if zlibdep.found()
  if cc.has_function('gzbuffer',
                     prefix: '#include <zlib.h>',
                     dependencies: zlibdep)
    conf_data.set10('HAVE_GZBUFFER', true)
  endif
else
  zlibproj = subproject('zlib')
  zlibdep = zlibproj.get_variable('zlib_dep')
  # Cannot check for gzbuffer in this case, because Zlib won't be built
  # until compile time.  However, we already know that the wrap version is
  # new enough
  conf_data.set10('HAVE_GZBUFFER', true)
endif

if zlibdep.found()
  conf_data.set10('HAVE_ZLIB', true)
endif

fftwdep = dependency('fftw3', required: false)
if fftwdep.found()
  conf_data.set10('HAVE_FFTW', true)
endif

seedeedep = dependency('seedee', required: false)
if seedeedep.found()
  cjsondep = dependency('libcjson',
                        required: true,
                        fallback: ['cjson', 'libcjson_dep'])
  if cjsondep.found() and seedeedep.found()
    conf_data.set10('HAVE_SEEDEE', true)
  endif
else
  # Default 'not found' dependency
  cjsondep = dependency('', required : false)
endif

xgandalfdep = dependency('xgandalf',
                         required: false,
                         fallback: ['xgandalf', 'xgandalf_dep'])
if xgandalfdep.found()
  conf_data.set10('HAVE_XGANDALF', true)
endif

pinkindexerdep = dependency('pinkIndexer',
                            required: false,
                            fallback: ['pinkindexer', 'pinkindexer_dep'])
if pinkindexerdep.found()
  conf_data.set10('HAVE_PINKINDEXER', true)
endif

fdipdep = dependency('fdip',
                     required: false,
                     fallback: ['fdip', 'fdip_dep'])
if fdipdep.found()
  conf_data.set10('HAVE_FDIP', true)
endif

ccp4dep = dependency('libccp4c',
                     required: false,
                     fallback: ['libccp4c', 'libccp4c_dep'])
if ccp4dep.found()
  conf_data.set10('HAVE_LIBCCP4', true)
endif

msgpackdep = dependency('msgpack-c', required: false)
if not msgpackdep.found()
  msgpackdep = dependency('msgpack', required: false)
endif
if msgpackdep.found()
  conf_data.set10('HAVE_MSGPACK', true)
endif


libcrystfel_versionc = vcs_tag(input: 'src/libcrystfel-version.c.in',
                               output: 'libcrystfel-version.c')


libcrystfel_includes = include_directories('src')

# Find forkpty()
utildep = cc.find_library('util', required: true)
if cc.has_function('forkpty', dependencies: utildep, prefix: '#include <pty.h>')
  conf_data.set10('HAVE_FORKPTY_PTY_H', true)
elif cc.has_function('forkpty', dependencies: utildep, prefix: '#include <util.h>')
  conf_data.set10('HAVE_FORKPTY_UTIL_H', true)
else
  error('Couldn\'t find forkpty()')
endif


# Symmetry operation parser Flex/Bison stuff
flex = find_program('flex')
bison = find_program('bison')

flex_gen = generator(flex,
                     output: ['@BASENAME@-lex.c', '@BASENAME@-lex.h'],
                     arguments: ['--outfile=@OUTPUT0@',
                                 '--header-file=@OUTPUT1@',
                                 '@INPUT@'])

bison_gen = generator(bison,
                      output: ['@BASENAME@-parse.c', '@BASENAME@-parse.h'],
                      arguments: ['--output=@OUTPUT0@',
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
                       'src/profile.c',
                       'src/crystfel-mille.c',
                       'src/image-cbf.c',
                       'src/image-hdf5.c',
                       'src/image-msgpack.c',
                       'src/image-seedee.c',
                       'src/indexers/dirax.c',
                       'src/indexers/felix.c',
                       'src/indexers/mosflm.c',
                       'src/indexers/taketwo.c',
                       'src/indexers/xds.c',
                       'src/indexers/asdf.c',
                       'src/indexers/xgandalf.c',
                       'src/indexers/pinkindexer.c',
                       'src/indexers/fromfile.c',
                       symop_lex_ch,
                       symop_parse_ch]

configure_file(input: 'libcrystfel-config.h.meson.in',
               output: 'libcrystfel-config.h',
               configuration: conf_data)

libcrystfel_conf_inc = include_directories('.')

libcrystfel = library('crystfel', [libcrystfel_sources, libcrystfel_versionc],
                      soversion: libcrystfel_api_version,
                      include_directories: [libcrystfel_includes,
                                            libcrystfel_conf_inc],
                      dependencies: [mdep, utildep, fftwdep, gsldep, zlibdep,
                                     hdf5dep, pthreaddep,
                                     xgandalfdep, pinkindexerdep, fdipdep,
                                     ccp4dep, msgpackdep, seedeedep, cjsondep],
                      install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib:$ORIGIN',
                      install: true)

libcrystfeldep = declare_dependency(include_directories: libcrystfel_includes,
                                    link_with: libcrystfel,
                                    dependencies: gsldep)


install_headers(['src/reflist.h',
                 'src/symmetry.h',
                 'src/cell.h',
                 'src/reflist-utils.h',
                 'src/thread-pool.h',
                 'src/utils.h',
                 'src/geometry.h',
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
                 'src/rational.h',
                 'src/spectrum.h',
                 'src/datatemplate.h',
                 'src/colscale.h',
                 'src/detgeom.h',
                 'src/fom.h',
                 'src/crystfel-mille.h'],
                subdir: 'crystfel')

# API documentation (Doxygen)
doxygen = find_program('doxygen', required: false)

if doxygen.found()
  vers = vcs_tag(input: 'doc/Doxyfile-version.in',
                 output: 'Doxyfile-version')

  doc_data = configuration_data()
  doc_data.set('DOX_TOP', join_paths(meson.current_source_dir(), 'src'))
  doc_data.set('INDEX_MD', join_paths(meson.current_source_dir(), 'doc/index.md'))
  doc_data.set('CODING_MD', join_paths(meson.current_source_dir(), 'doc/coding.md'))
  doc_data.set('OUTPUT_DIR', join_paths(meson.current_build_dir(), 'docs'))
  doc_data.set('DT_PRIV', join_paths(meson.current_source_dir(), 'src/datatemplate_priv.h'))
  doc_data.set('UTHASH', join_paths(meson.current_source_dir(), 'src/uthash.h'))
  doc_data.set('VERSION', meson.project_version())
  doc_data.set('API_REVISION', libcrystfel_api_version)
  doc_data.set('LOGO_FILE', join_paths(meson.current_source_dir(), 'doc/crystfel-logo.png'))

  doxyfile = configure_file(input: 'doc/Doxyfile.in',
                            output: 'Doxyfile',
                            configuration: doc_data,
                            install: false)

  api_docs = run_target('api-docs',
                        command: [doxygen, doxyfile],
                        depends: vers)
endif


# pkg-config file
pkg = import('pkgconfig')
pkg.generate(libcrystfel,
             description: 'Data processing for serial crystallography (shared library)')
