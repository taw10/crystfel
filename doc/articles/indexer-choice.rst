======================================
Which indexing method(s) should I use?
======================================

The short answer
================

For most applications, the best choice of indexing method is ``xgandalf``.

In most cases, it's not necessary to use several algorithms in sequence
(e.g. ``--indexing=xgandalf,mosflm,asdf``).

Read on for recommendations in specific cases.


Unknown crystal lattice parameters
==================================

Most of the indexing algorithms run fine without prior information about the
lattice parameters.  The indexing algorithms will generate a set of estimated
lattice parameters for each crystal, and you can use ``cell_explorer`` to plot
histograms to determine the most "popular" parameters.

Problems arise because there are an infinite number of ways to represent any
given lattice.  Different indexing algorithms have different "opinions"
about which representation to use when not constrained by prior information.

In this case, ``mosflm`` is a good choice of indexing method because it can
detect centered lattices instead of producing a primitive unit cell every time.
This makes it easier and quicker to see the correct symmetry (including indexing
ambiguities) and get to the crystallographically conventional representation of
the structure.

Other indexing methods, particularly ``xgandalf``, are fine as long as you use
proper crystallographic knowledge to determine the symmetry for merging and
structure solution.

It's particularly advisable to avoid using multiple indexing methods together
when the lattice parameters are unknown, because the "preferences" of different
algorithms will be mixed up together.


Fastest indexing (e.g. real-time analysis pipelines)
====================================================

Read the `document about data processing speed <speed.rst>`_.  If you have a
GPU available, use ``ffbidx``.  Otherwise, the best choice is probably ``asdf``
or ``xgandalf`` with ``--xgandalf-fast``, but there are many other
considerations for increasing processing speed.


Wide-bandwidth X-ray diffraction (Laue/pink beam)
=================================================

Use ``pinkindexer``, especially if the (FWHM) bandwidth is larger than about 3%
of the wavelength.  For smaller bandwidths, try with ``xgandalf``.

For larger bandwidths, it's also advisable to use options ``--no-refine`` and
``--no-check-peaks``, because these parts of the program assume that the
bandwidth is small.


Electron diffraction data
=========================

Use ``pinkindexer``, but first read the `document about electrons <electrons.rst>`_.


Small unit cell
===============

For very small unit cells (<20 Angstrom axis length), use ``smallcell``.  For
larger (but still "small") cells, researchers have also had success using
``xgandalf``, and it's also worth trying with ``taketwo``.


Multiple overlapping diffraction patterns
=========================================

Try ``felix``, which is specifically designed for heavily overlapping
diffraction patterns.  Note, however, that it's unfortunately difficult to get
hold of a copy of the Felix executable.


CrystFEL was installed minimally
================================

The ``taketwo`` algorithm has no dependencies, neither compile-time nor
run-time, so is always available.  Note that ``taketwo`` always needs prior
lattice parameters.

Even if you have no system-wide access to the computer system, you can still
make a user-local installation of ``dirax``, ``mosflm``, ``xds`` or ``felix``.
These methods have no compile-time requirements, and only need the
corresponding programs to be available at run-time.
