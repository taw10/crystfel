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
test('file-wait',
     find_program('file-wait'),
     args : [indexamajig.full_path(),
             files('wavelength_geom.h5'),
             files('wavelength_geom1.geom')])


# Easy unit tests of libcrystfel functions
simple_tests = ['ambi_check',
                'cell_check',
                'centering_check',
                'list_check',
                'prediction_gradient_check',
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
                   dependencies : [libcrystfeldep, mdep, gsldep, hdf5dep],
		   include_directories: conf_inc)

  test(name, exe, timeout : 60)
endforeach


# Less easy unit tests for libcrystfel functions
exe = executable('stream_read',
                 ['stream_read.c'],
                 dependencies : [libcrystfeldep])
test_stream = files('test.stream')
test('stream_read',
     exe,
     args: [test_stream])

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

if opencldep.found()
  exe = executable('gpu_sim_check',
                   ['gpu_sim_check.c',
                    simulation_bits],
                   dependencies : [libcrystfeldep, mdep, gsldep, opencldep],
                   include_directories: conf_inc)
  geom = files('gpu_sim_check.geom')
  test('gpu_sim_check', exe, args: [geom])
endif

# Event enumeration tests
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


# Wavelength tests
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
