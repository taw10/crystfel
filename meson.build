# Meson file for CrystFEL
project('crystfel', 'c',
        version: '0.10.2',
        license: 'GPL3+',
        meson_version: '>=0.55.0',
        default_options: ['buildtype=debugoptimized'])

libcrystfel_api_version = 16
crystfel_rpath = '$ORIGIN/../lib64/:$ORIGIN/../lib:$ORIGIN/../../asapo-libraries/lib:$ORIGIN/../../asapo-libraries/lib64'

add_project_arguments('-DHAVE_CONFIG_H', language: 'c')

versionc = vcs_tag(input: 'src/version.c.in',
                   output: 'version.c')


conf_data = configuration_data()

# Mandatory dependencies
cc = meson.get_compiler('c')
mdep = cc.find_library('m', required: true)
gsldep = dependency('gsl', required: true)
pthreaddep = dependency('threads', required: true)

# Find HDF5 using inbuilt Meson methods.  Requires Meson >= 0.50.0
hdf5dep = dependency('hdf5', language: 'c', required: get_option('hdf5'))
if hdf5dep.found()
  conf_data.set10('HAVE_HDF5', true)
  # Request HDF5 1.10-style API
  add_project_arguments('-DH5_USE_110_API', language: 'c')
  add_project_arguments('-DH5O_GET_INFO_BY_IDX1', language: 'c')
endif

fftwdep = dependency('fftw3', required: false)
if fftwdep.found()
  conf_data.set10('HAVE_FFTW', true)
endif

gnome = import('gnome')

gtkdep = dependency('gtk+-3.0', required: false)
if gtkdep.found()
  add_project_arguments('-Wno-deprecated-declarations', language: 'c')
  conf_data.set10('HAVE_GTK', true)
endif

cairodep = dependency('cairo', required: false)
if cairodep.found()
  conf_data.set10('HAVE_CAIRO', true)
endif

zmqdep = dependency('libzmq', required: false)
if zmqdep.found()
  conf_data.set10('HAVE_ZMQ', true)
endif

# There are two separate dependencies for ASAP::O (consumer and producer)
# We need both of them.
asapoproddep = dependency('libasapo-producer', required: false)
if asapoproddep.found()
  asapodep = dependency('libasapo-consumer', required: true)
  if asapodep.found()
    conf_data.set10('HAVE_ASAPO', true)
  endif
else
  asapodep = dependency('', required: false)
endif

if cc.has_function('clock_gettime', prefix: '#include <time.h>')
  conf_data.set10('HAVE_CLOCK_GETTIME', true)
endif

if cc.has_function('sched_setaffinity',
                   prefix: '#include <sched.h>',
                   args: '-D_GNU_SOURCE')
  conf_data.set10('HAVE_SCHED_SETAFFINITY', true)
endif

# ************************ libcrystfel (subdir) ************************

subdir('libcrystfel')

# ************************ The programs ************************

configure_file(input: 'config.h.in',
               output: 'config.h',
               configuration: conf_data)
conf_inc = include_directories('.')

# cell_tool
executable('cell_tool',
           ['src/cell_tool.c', versionc],
           dependencies: [mdep, libcrystfeldep],
           install: true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# process_hkl
process_hkl = executable('process_hkl',
                         ['src/process_hkl.c', versionc],
                         dependencies: [mdep, libcrystfeldep],
                         install: true,
                         install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# list_events
executable('list_events',
           ['src/list_events.c', versionc],
           dependencies: [mdep, libcrystfeldep],
           install: true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# get_hkl
executable('get_hkl',
           ['src/get_hkl.c', versionc],
           dependencies: [mdep, libcrystfeldep, gsldep],
           install: true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# compare_hkl
executable('compare_hkl',
           ['src/compare_hkl.c', versionc],
           dependencies: [mdep, libcrystfeldep, gsldep],
           install: true,
           install_rpath: crystfel_rpath)

# check_hkl
executable('check_hkl',
           ['src/check_hkl.c', versionc],
           dependencies: [mdep, libcrystfeldep, gsldep],
           install: true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# partialator
partialator = executable('partialator',
                         ['src/partialator.c',
                          'src/post-refinement.c',
                          'src/merge.c',
                          'src/rejection.c',
                          'src/scaling.c',
                          versionc],
                         dependencies: [mdep, libcrystfeldep, gsldep, pthreaddep],
                         install: true,
                         install_rpath: crystfel_rpath)

# ambigator
executable('ambigator',
           ['src/ambigator.c', versionc],
           dependencies: [mdep, libcrystfeldep, gsldep, hdf5dep],
           install: true,
           install_rpath: crystfel_rpath)

# whirligig
executable('whirligig',
           ['src/whirligig.c', versionc],
           dependencies: [mdep, libcrystfeldep],
           install: true,
           install_rpath: crystfel_rpath)

# indexamajig
indexamajig_sources = ['src/indexamajig.c', 'src/im-sandbox.c',
                       'src/process_image.c',
                       versionc]
if zmqdep.found()
  indexamajig_sources += ['src/im-zmq.c']
endif

if asapodep.found()
  indexamajig_sources += ['src/im-asapo.c']
endif

indexamajig = executable('indexamajig', indexamajig_sources,
                         dependencies: [mdep, libcrystfeldep, gsldep,
                                        pthreaddep, zmqdep, asapodep, asapoproddep],
                         install: true,
                         install_rpath: crystfel_rpath)

# make_pixelmap
if hdf5dep.found()
  executable('make_pixelmap',
             ['src/make_pixelmap.c', versionc],
             dependencies: [mdep, libcrystfeldep, hdf5dep],
             install: true,
             install_rpath: crystfel_rpath)
endif


# CrystFEL GUI
if gtkdep.found()

  gui_sources = ['src/crystfel_gui.c',
                 'src/crystfelimageview.c',
                 'src/crystfelindexingopts.c',
                 'src/crystfelmergeopts.c',
                 'src/crystfelsymmetryselector.c',
                 'src/crystfelcolourscale.c',
                 'src/crystfelfomgraph.c',
                 'src/gtk-util-routines.c',
                 'src/gui_import.c',
                 'src/gui_peaksearch.c',
                 'src/gui_index.c',
                 'src/gui_merge.c',
                 'src/gui_fom.c',
                 'src/gui_export.c',
                 'src/gui_ambi.c',
                 'src/gui_backend_local.c',
                 'src/gui_backend_slurm.c',
                 'src/gui_project.c',
                 'src/gtkmultifilechooserbutton.c',
                 versionc]

  gresources = gnome.compile_resources('crystfel-gui-resources',
                                       'data/crystfel.gresource.xml',
                                       source_dir: 'data', c_name: 'crystfel')

  executable('crystfel', gui_sources, gresources,
             dependencies: [mdep, libcrystfeldep, gtkdep, gsldep],
             install: true,
             install_rpath: crystfel_rpath)

  install_data(['data/crystfel.desktop'],
               install_dir: get_option('datadir') / 'applications')

endif

# render_hkl
executable('render_hkl',
           ['src/render_hkl.c', versionc],
           dependencies: [mdep, libcrystfeldep, cairodep, gsldep],
           install: true,
           install_rpath: crystfel_rpath)

# cell_explorer
if gtkdep.found()
  executable('cell_explorer',
             ['src/cell_explorer.c', 'src/multihistogram.c', versionc],
             dependencies: [mdep, libcrystfeldep, gtkdep, gsldep],
             install: true,
             install_rpath: crystfel_rpath)

  install_data(['data/cell_explorer.desktop'],
               install_dir: get_option('datadir') / 'applications')
endif


# ************************ Misc resources ************************

install_data(['data/crystfel.svg'],
             install_dir: get_option('datadir')+'/icons/hicolor/scalable/apps')


# ************************ Scripts ************************

install_data(['scripts/detector-shift',
              'scripts/peakogram-stream'],
             install_dir: 'bin')


# ************************ Testing ************************

subdir('tests')


# ************************ Manual pages ************************
install_man(['doc/man/ambigator.1',
             'doc/man/cell_explorer.1',
             'doc/man/cell_tool.1',
             'doc/man/check_hkl.1',
             'doc/man/compare_hkl.1',
             'doc/man/crystfel.7',
             'doc/man/crystfel_geometry.5',
             'doc/man/get_hkl.1',
             'doc/man/indexamajig.1',
             'doc/man/list_events.1',
             'doc/man/list_events.1',
             'doc/man/make_pixelmap.1',
             'doc/man/partialator.1',
             'doc/man/process_hkl.1',
             'doc/man/render_hkl.1',
             'doc/man/whirligig.1'])
