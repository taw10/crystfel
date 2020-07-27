# libcrystfel

zlibdep = dependency('zlib', required : true)

fftwdep = dependency('fftw3', required : false)
if fftwdep.found()
  conf_data.set10('HAVE_FFTW', 1)
endif

# FIXME: Curses, Xgandalf, pinkindexer, FDIP

libcrystfel_includes = include_directories('src')

# Find forkpty()
utildep = cc.find_library('util', required : true)
if cc.has_function('forkpty', dependencies : utildep, prefix : '#include <pty.h>')
  add_project_arguments('-DHAVE_FORKPTY_PTY_H', language : 'c')
elif cc.has_function('forkpty', dependencies : utildep, prefix : '#include <util.h>')
  add_project_arguments('-DHAVE_FORKPTY_UTIL_H', language : 'c')
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
                       'src/mosflm.c',
                       'src/taketwo.c',
                       'src/dirax.c',
                       'src/peakfinder8.c',
                       'src/thread-pool.c',
                       'src/peaks.c',
                       'src/utils.c',
                       'src/felix.c',
                       'src/predict-refine.c',
                       'src/xds.c',
                       'src/filters.c',
                       'src/rational.c',
                       'src/xgandalf.c',
                       'src/pinkindexer.c',
                       'src/geometry.c',
                       'src/reflist.c',
                       'src/reflist-utils.c',
		       'src/datatemplate.c',
		       'src/colscale.c',
		       'src/detgeom.c',
		       'src/image-cbf.c',
		       'src/image-hdf5.c',
		       'src/image-msgpack.c',
                       symop_lex_ch,
                       symop_parse_ch,
                      ]

if fftwdep.found()
  libcrystfel_sources += 'src/asdf.c'
endif

configure_file(input : 'config.h.in',
               output : 'config.h',
               configuration : conf_data)

libcrystfel_conf_inc = include_directories('.')

libcrystfel = library('crystfel', libcrystfel_sources,
                   include_directories : [libcrystfel_includes, libcrystfel_conf_inc],
                   dependencies : [mdep, utildep, fftwdep, gsldep, zlibdep,
                                   hdf5dep, pthreaddep],
                   install : true)

libcrystfeldep = declare_dependency(include_directories : libcrystfel_includes,
                                    link_with : libcrystfel)


# FIXME: pkgconfig file
# FIXME: API documents (Doxygen)
# FIXME: install headers
