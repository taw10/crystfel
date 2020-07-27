# CrystFEL unit tests

exe = executable('ambi_check', 'ambi_check.c',
                        dependencies : libcrystfeldep)
test('ambi_check', exe)

exe = executable('cell_check', 'cell_check.c',
                        dependencies : [libcrystfeldep, mdep, gsldep])
test('cell_check', exe)


exe = executable('centering_check', 'centering_check.c',
                 dependencies : [libcrystfeldep, mdep, gsldep])
test('centering_check', exe)

exe = executable('integration_check',
                 ['integration_check.c',
                  'histogram.c'],
                 dependencies : [libcrystfeldep, mdep, gsldep])
test('integration_check', exe)
