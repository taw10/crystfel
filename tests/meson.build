# CrystFEL unit tests

# Test wrapper scripts for process_hkl
process_hkl_tests = ['process_hkl_check_1',
                     'process_hkl_check_2',
                     'process_hkl_check_3',
                     'process_hkl_check_4']

foreach name : process_hkl_tests
  exe = find_program(name)
  test(name, exe, args : [process_hkl.full_path()])
endforeach


# Test wrapper scripts for partialator
partialator_tests = ['partialator_merge_check_1',
                     'partialator_merge_check_2',
                     'partialator_merge_check_3']

foreach name : partialator_tests
  exe = find_program(name)
  test(name, exe, args : [partialator.full_path()])
endforeach


# Test of waiting for files
if hdf5dep.found()
  test('file-wait',
       find_program('file-wait'),
       args : [indexamajig.full_path(),
               files('wavelength_geom.h5'),
               files('wavelength_geom1.geom')])
endif


test('indexamajig-missing-file',
     find_program('indexamajig-missing-file'),
     args : [indexamajig.full_path(),
             files('wavelength_geom1.geom')])


# Easy unit tests of libcrystfel functions
simple_tests = ['ambi_check',
                'cell_check',
                'centering_check',
                'list_check',
                'ring_check',
                'symmetry_check',
                'transformation_check',
                'rational_check',
                'spectrum_check',
                'cellcompare_check',
                'evparse1',
                'evparse2',
                'evparse3',
                'evparse4',
                'evparse5',
                'evparse6',
                'evparse7',
                'symop_parse']

foreach name : simple_tests
  exe = executable(name, ''.join([name, '.c']),
                   dependencies : [libcrystfeldep, mdep, gsldep],
		   include_directories: conf_inc)

  test(name, exe, timeout : 60)
endforeach


# Less easy unit tests for libcrystfel functions
exe = executable('stream_roundtrip',
                 ['stream_roundtrip.c'],
                 dependencies : [libcrystfeldep])
test('stream_roundtrip', exe,
     args: [files('stream_roundtrip.geom')])

exe = executable('stream_read',
                 ['stream_read.c'],
                 dependencies : [libcrystfeldep])
test_stream = files('test.stream')
test('stream_read',
     exe,
     args: [test_stream])

if hdf5dep.found()

  exe = executable('header_read_single',
                   ['header_read.c'],
                   dependencies : [libcrystfeldep, hdf5dep])
  test('header_read_single', exe,
       args: [files('header_read_single.geom'), files('header_read.h5'), '//', '/bigint_single', 'header_read_single.stream'])

  exe = executable('header_read_arr',
                   ['header_read.c'],
                   dependencies : [libcrystfeldep, hdf5dep])
  test('header_read_arr', exe,
       args: [files('header_read_arr.geom'), files('header_read.h5'), '//1', '/bigint_arr', 'header_read_arr.stream'])

endif

exe = executable('integration_check',
                 ['integration_check.c',
                  'histogram.c'],
                 dependencies : [libcrystfeldep, mdep, gsldep])
test('integration_check', exe)

exe = executable('prof2d_check',
                 ['prof2d_check.c',
                  'histogram.c'],
                 dependencies : [libcrystfeldep, mdep, gsldep])
test('prof2d_check', exe)


# Refinement gradient checks, part 1: panel translations
panel_gradient_tests = [['gradient_panel_x', 'cnx', 'GPARAM_DET_TX'],
                        ['gradient_panel_y', 'cny', 'GPARAM_DET_TY'],
                        ['gradient_panel_z', 'cnz', 'GPARAM_DET_TZ']]

foreach name : panel_gradient_tests
  exe = executable(name[0],
                   ['gradient_check.c',
                    'gradient_check_utils.c'],
                   c_args : ['-DTHING_TO_MOVE='+name[1],
                             '-DTEST_GPARAM='+name[2],
                             '-DTRANSLATE_PANEL=1'],
                   dependencies : [libcrystfeldep, mdep, gsldep])
  test(name[0], exe)
endforeach


# Refinement gradient checks, part 2: panel rotations
panel_gradient_tests = [['gradient_panel_rx', '-DROTATE_PANEL_X', 'GPARAM_DET_RX'],
                        ['gradient_panel_ry', '-DROTATE_PANEL_Y', 'GPARAM_DET_RY'],
                        ['gradient_panel_rz', '-DROTATE_PANEL_Z', 'GPARAM_DET_RZ']]

foreach name : panel_gradient_tests
  exe = executable(name[0],
                   ['gradient_check.c',
                    'gradient_check_utils.c'],
                   c_args : [name[1],
                             '-DTEST_GPARAM='+name[2]],
                   dependencies : [libcrystfeldep, mdep, gsldep])
  test(name[0], exe)
endforeach


# Refinement gradient checks, part 3: diffraction physics
panel_gradient_tests = [['gradient_cell_asx', 'asx', 'GPARAM_ASX'],
                        ['gradient_cell_asy', 'asy', 'GPARAM_ASY'],
                        ['gradient_cell_asz', 'asz', 'GPARAM_ASZ'],
                        ['gradient_cell_bsx', 'bsx', 'GPARAM_BSX'],
                        ['gradient_cell_bsy', 'bsy', 'GPARAM_BSY'],
                        ['gradient_cell_bsz', 'bsz', 'GPARAM_BSZ'],
                        ['gradient_cell_csx', 'csx', 'GPARAM_CSX'],
                        ['gradient_cell_csy', 'csy', 'GPARAM_CSY'],
                        ['gradient_cell_csz', 'csz', 'GPARAM_CSZ']]

foreach name : panel_gradient_tests
  exe = executable(name[0],
                   ['gradient_check.c',
                    'gradient_check_utils.c'],
                   c_args : ['-DTHING_TO_MOVE='+name[1],
                             '-DTEST_GPARAM='+name[2],
                             '-DCHANGE_CELL=1'],
                   dependencies : [libcrystfeldep, mdep, gsldep])
  test(name[0], exe)
endforeach


# Event enumeration tests
if hdf5dep.found()
  ev_enum_tests = ['ev_enum1',
                   'ev_enum2',
                   'ev_enum3']

  foreach name : ev_enum_tests
    exe = executable(name, ''.join([name, '.c']),
                     dependencies : [libcrystfeldep, hdf5dep])
    h5 = files(''.join([name, '.h5']))
    geom = files(''.join([name, '.geom']))
    test(name, exe, args : [h5, geom])
  endforeach
endif


# Wavelength tests
if hdf5dep.found()
  wavelength_tests = [['wavelength_geom1', '1e-10'],
                      ['wavelength_geom2', '1.3776e-10'],
                      ['wavelength_geom3', '1.3776e-10'],
                      ['wavelength_geom4', '1.9687e-12'],
                      ['wavelength_geom5', '1.9687e-12'],
                      ['wavelength_geom6', '1.3776e-10'],
                      ['wavelength_geom7', '1.3776e-10'],
                      ['wavelength_geom8', '1.9687e-12'],
                      ['wavelength_geom9', '1.3776e-10'],
                      ['wavelength_geom10', '1.3776e-10'],
                      ['wavelength_geom11', '1.125e-10'],
                      ['wavelength_geom12', '1.125e-10']]

  exe = executable('wavelength_geom', 'wavelength_geom.c',
                   dependencies : [libcrystfeldep, hdf5dep])

  h5 = files('wavelength_geom.h5')

  foreach p : wavelength_tests
    geom = files(''.join([p[0], '.geom']))
    expected_wl = p[1]
    test(p[0], exe, args : [h5, geom, expected_wl])
  endforeach
endif


test('geom_roundtrip',
     find_program('geom_roundtrip'),
     args: [adjust_detector.full_path()])
