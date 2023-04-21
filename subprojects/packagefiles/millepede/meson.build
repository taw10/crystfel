# Meson file for Millepede-II
project('millepede', ['c', 'cpp','fortran'],
        version: '2.0.0',
        license: 'LGPL2+',
        default_options: ['buildtype=debugoptimized'])


mille = library('mille', ['Mille.cc'],
                install: true)

millepede_dep = declare_dependency(include_directories: '.',
                                   link_with: mille)

executable('pede',
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
           install: true)
