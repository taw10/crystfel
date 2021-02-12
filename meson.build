# Meson file for CrystFEL
project('crystfel', 'c',
        version : '0.9.1',
        license : 'GPL3+',
        default_options : ['buildtype=debugoptimized'])

add_project_arguments('-DHAVE_CONFIG_H', language : 'c')

versionc = vcs_tag(input: 'src/version.c.in',
                   output: 'version.c')


conf_data = configuration_data()

# Mandatory dependencies
cc = meson.get_compiler('c')
mdep = cc.find_library('m', required : true)
gsldep = dependency('gsl', required : true)
zlibdep = dependency('zlib', required : true)
pthreaddep = dependency('threads', required : true)

# Try via pkg-config first:
slurmdep = dependency('slurm', required : false)
if not slurmdep.found()
  slurmdep = cc.find_library('slurm', required : false)
endif

if slurmdep.found()
  conf_data.set10('HAVE_SLURM', 1)
endif

# HDF5 >= 1.8.21 and >= 1.10.2 has a pkg-config file:
hdf5dep = dependency('hdf5', language : 'c', required : false)

# Otherwise, try finding it ourselves:
if not hdf5dep.found()
  hdf5dep = cc.find_library('hdf5', required : false)
endif

if not hdf5dep.found()
  error('Couldn\'t find HDF5')
endif

# Request HDF5 1.10-style API
add_project_arguments('-DH5_USE_110_API', language : 'c')
add_project_arguments('-DH5O_GET_INFO_BY_IDX1', language : 'c')

# "Optional" dependencies
fftwdep = dependency('fftw3', required : false)
if fftwdep.found()
  conf_data.set10('HAVE_FFTW', 1)
endif

gtkdep = dependency('gtk+-3.0', required : false)
if gtkdep.found()
  add_project_arguments('-Wno-deprecated-declarations', language : 'c')
  conf_data.set10('HAVE_GTK', 1)
endif

gdkdep = dependency('gdk-3.0', required : false)
if gdkdep.found()
  conf_data.set10('HAVE_GDK', 1)
endif

gdkpixbufdep = dependency('gdk-pixbuf-2.0', required : false)
if gdkpixbufdep.found()
  conf_data.set10('HAVE_GDKPIXBUF', 1)
endif

cairodep = dependency('cairo', required : false)
if cairodep.found()
  conf_data.set10('HAVE_CAIRO', 1)
endif

opencldep = dependency('OpenCL', required : false)
if opencldep.found()
  conf_data.set10('HAVE_OPENCL', 1)
endif
if cc.has_header('CL/cl.h')
  conf_data.set10('HAVE_CL_CL_H', 1)
endif

zmqdep = dependency('libzmq', required : false)
if zmqdep.found()
  conf_data.set10('HAVE_ZMQ', 1)
endif

msgpackdep = dependency('msgpack', required : false)
if msgpackdep.found()
  conf_data.set10('HAVE_MSGPACK', 1)
endif

if cc.has_function('clock_gettime', prefix : '#include <time.h>')
  conf_data.set10('HAVE_CLOCK_GETTIME', 1)
endif

ccp4dep = dependency('libccp4c', required: false)
if ccp4dep.found()
  conf_data.set10('HAVE_LIBCCP4', 1)
endif

# ************************ libcrystfel (subdir) ************************

subdir('libcrystfel')

# ************************ The programs ************************

configure_file(input : 'config.h.in',
               output : 'config.h',
               configuration : conf_data)
conf_inc = include_directories('.')

# cell_tool
executable('cell_tool',
           ['src/cell_tool.c', versionc],
           dependencies : [mdep, libcrystfeldep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# partial_sim
executable('partial_sim',
           ['src/partial_sim.c', versionc],
           dependencies : [mdep, libcrystfeldep, gsldep, pthreaddep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# process_hkl
process_hkl = executable('process_hkl',
                         ['src/process_hkl.c', versionc],
                         dependencies : [mdep, libcrystfeldep],
                         install : true,
                         install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# list_events
executable('list_events',
           ['src/list_events.c', versionc],
           dependencies : [mdep, libcrystfeldep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# get_hkl
executable('get_hkl',
           ['src/get_hkl.c', versionc],
           dependencies : [mdep, libcrystfeldep, gsldep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# compare_hkl
executable('compare_hkl',
           ['src/compare_hkl.c', versionc],
           dependencies : [mdep, libcrystfeldep, gsldep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# check_hkl
executable('check_hkl',
           ['src/check_hkl.c', versionc],
           dependencies : [mdep, libcrystfeldep, gsldep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# partialator
partialator = executable('partialator',
                         ['src/partialator.c',
                          'src/post-refinement.c',
                          'src/merge.c',
                          'src/rejection.c',
                          'src/scaling.c',
                          versionc],
                         dependencies : [mdep, libcrystfeldep, gsldep, pthreaddep],
                         install : true,
                         install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# ambigator
executable('ambigator',
           ['src/ambigator.c', versionc],
           dependencies : [mdep, libcrystfeldep, gsldep, hdf5dep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# whirligig
executable('whirligig',
           ['src/whirligig.c', versionc],
           dependencies : [mdep, libcrystfeldep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# indexamajig
indexamajig_sources = ['src/indexamajig.c', 'src/im-sandbox.c',
                       'src/process_image.c', 'src/time-accounts.c',
                       versionc]
if zmqdep.found() and msgpackdep.found()
  indexamajig_sources += ['src/im-zmq.c']
endif

executable('indexamajig', indexamajig_sources,
           dependencies : [mdep, libcrystfeldep, gsldep, pthreaddep, zmqdep,
                           msgpackdep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# make_pixelmap
executable('make_pixelmap',
           ['src/make_pixelmap.c', versionc],
           dependencies : [mdep, libcrystfeldep, hdf5dep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# geoptimiser
# FIXME: restore
#executable('geoptimiser',
#           ['src/geoptimiser.c', 'src/hdfsee-render.c', versionc],
#           dependencies : [mdep, libcrystfeldep, gsldep, gdkpixbufdep, gdkdep],
#          install : true,
#          install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# CrystFEL GUI
if gtkdep.found()

  gui_sources = ['src/crystfel_gui.c',
                 'src/crystfelimageview.c',
                 'src/crystfelindexingopts.c',
                 'src/crystfelmergeopts.c',
                 'src/crystfelsymmetryselector.c',
                 'src/gtk-util-routines.c',
                 'src/gui_peaksearch.c',
                 'src/gui_index.c',
                 'src/gui_merge.c',
                 'src/gui_fom.c',
                 'src/gui_export.c',
                 'src/gui_backend_local.c',
                 'src/gui_project.c',
                 versionc]

  if slurmdep.found()
    gui_sources += 'src/gui_backend_slurm.c'
  endif


  executable('crystfel', gui_sources,
             dependencies : [mdep, libcrystfeldep, gtkdep, gsldep, slurmdep,
	                     ccp4dep],
             install : true,
             install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

endif

# render_hkl
executable('render_hkl',
           ['src/render_hkl.c', versionc],
           dependencies : [mdep, libcrystfeldep, cairodep, gsldep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')

# cell_explorer
if gtkdep.found()
  executable('cell_explorer',
             ['src/cell_explorer.c', 'src/multihistogram.c', versionc],
             dependencies : [mdep, libcrystfeldep, gtkdep, gsldep],
             install : true,
             install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')
endif

# pattern_sim
pattern_sim_sources = ['src/pattern_sim.c', 'src/diffraction.c', versionc]
if opencldep.found()
  pattern_sim_sources += ['src/diffraction-gpu.c', 'src/cl-utils.c']
endif
executable('pattern_sim', pattern_sim_sources,
           dependencies : [mdep, libcrystfeldep, gsldep, hdf5dep, opencldep],
           install : true,
           install_rpath: '$ORIGIN/../lib64/:$ORIGIN/../lib')
# For testing (see tests/meson.build)
simulation_bits = files(['src/diffraction.c', 'src/diffraction-gpu.c', 'src/cl-utils.c'])


# ************************ Misc resources ************************

install_data(['data/crystfel.svg'],
             install_dir : get_option('datadir')+'/icons/hicolor/scalable/apps')


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
             'doc/man/geoptimiser.1',
             'doc/man/get_hkl.1',
             'doc/man/indexamajig.1',
             'doc/man/list_events.1',
             'doc/man/list_events.1',
             'doc/man/make_pixelmap.1',
             'doc/man/partialator.1',
             'doc/man/partial_sim.1',
             'doc/man/pattern_sim.1',
             'doc/man/process_hkl.1',
             'doc/man/render_hkl.1',
             'doc/man/whirligig.1'])
