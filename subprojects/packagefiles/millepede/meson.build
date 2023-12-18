# Meson file for Millepede-II
project('millepede', ['c', 'fortran'],
        version: 'V04-13-06',
        license: 'LGPL2+',
        default_options: ['buildtype=debugoptimized'])

pede_c_args = []

zlibdep = dependency('zlib', required: false)
if zlibdep.found()
  pede_c_args += '-DUSE_ZLIB'
endif

pede = executable('pede',
           ['mpdef.f90',
            'mpdalc.f90',
            'mpmod.f90',
            'mpmon.f90',
            'mpbits.f90',
            'mpqldec.f90',
            'mptest1.f90',
            'mptest2.f90',
            'mille.f90',
            'mpnum.f90',
            'mptext.f90',
            'mphistab.f90',
            'minresDataModule.f90',
            'minresModule.f90',
            'minresqlpDataModule.f90',
            'minresqlpBlasModule.f90',
            'minresqlpModule.f90',
            'randoms.f90',
            'vertpr.f90',
            'linesrch.f90',
            'Dbandmatrix.f90',
            'pede.f90',
            'readc.c'],
           fortran_args: ['-DREAD_C_FILES', '-cpp'],
           c_args: pede_c_args,
           build_by_default: true,
           dependencies: [zlibdep],
           install: true)
