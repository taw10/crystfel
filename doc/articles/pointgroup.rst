===============================================
How to choose the right point group for merging
===============================================

A common question from our users is how to choose the correct symmetry for
merging, i.e. the correct ``-y`` option.  It's actually not that difficult, but
it does touch on several areas of crystallography theory. This document aims
to be a gentle introduction to the process, introducing the concepts step by
step. For a somewhat terse explanation, see section 6 of the following paper:

T. A. White, A. Barty, F. Stellato et al
"Crystallographic data processing for free-electron laser sources"
Acta Cryst. D69 (2013), p1231–1240.
`doi:10.1107/S0907444913013620 <http://dx.doi.org/10.1107/S0907444913013620>`_

Another useful article is the following:

M. Nespolo, M. I. Aroyo and B. Souvignier
"Crystallographic shelves: space-group hierarchy explained"
J. Applied Cryst. 51 (2018) p1481-1491
`doi:10.1107/S1600576718012724 <https://doi.org/10.1107/S1600576718012724>`_

Step 1: Temporarily forget about space groups
=============================================

To merge your reflection data, CrystFEL needs to know which reflections are
symmetrically equivalent.  This information is given by the *point group*.
The 230 space groups can be classified into 17 categories, each corresponding
to a single point group.  If you know the *space* group for your crystals in
advance, that's a big advantage.  However, for now you only need to know the
*point* group.

If you're working on an unknown structure, don't get ahead of yourself!
Many crystallographic data processing programs start suggesting possible space
groups very early in the process, such as when the patterns are indexed.
The space group is **only a hypothesis until the structure is solved**, so you
always need to take these early suggestions with a pinch of salt.  CrystFEL's
design philosophy is not to deal with space group determination at all.
CrystFEL will never ask you to tell it the space group of your structure ahead
of time, nor will it suggest one automatically for your structure [#f1]_.

The following table shows the point group corresponding to each of the space
groups.  To keep things simple, the table only contains the `Sohncke space
groups <https://dictionary.iucr.org/Sohncke_groups>`_, which are the ones
relevant to biological structures.  The point groups are given in exactly the
form you will type them into CrystFEL:

===========    ============
Point group    Space groups
===========    ============
``1``          P1
``2``          P2, P2\ :sub:`1`, C2 (pay special attention to step 3 below)
``222``        P222, P222\ :sub:`1`, P2\ :sub:`1`\ 2\ :sub:`1`\ 2, P2\ :sub:`1`\ 2\ :sub:`1`\ 2\ :sub:`1`, C222\ :sub:`1`, C222, F222, I222, I2\ :sub:`1`\ 2\ :sub:`1`\ 2\ :sub:`1`
``4``          P4 P4\ :sub:`1`, P4\ :sub:`2`, P4\ :sub:`3`, I4, I4\ :sub:`1`
``422``        P422, P42\ :sub:`1`\ 2, P4\ :sub:`1`\ 22, P4\ :sub:`1`\ 2\ :sub:`1`\ 2, P4\ :sub:`2`\ 22, P4\ :sub:`2`\ 2\ :sub:`1`\ 2, P4\ :sub:`3`\ 22, P4\ :sub:`3`\ 2\ :sub:`1`\ 2, I4222, I4\ :sub:`1`\ 22
``32_R``       R32 (rhombohedral axes, pay special attention to step 6)
``3_R``        R3 (rhombohedral axes, pay special attention to step 6)
``3_H``        H3 (hexagonal axes, pay special attention to step 6), P3, P3\ :sub:`1`, P3\ :sub:`2`
``321_H``      H32 (hexagonal axes, pay special attention to step 6), P321, P3\ :sub:`1`\ 21, P3\ :sub:`2`\ 21
``312_H``      P312, P3\ :sub:`1`\ 12, P3\ :sub:`2`\ 12
``6``          P6, P6\ :sub:`1`, P6\ :sub:`2`, P6\ :sub:`3`, P6\ :sub:`4`, P6\ :sub:`5`
``622``        P622, P6\ :sub:`1`\ 22, P6\ :sub:`2`\ 22, P6\ :sub:`3`\ 22, P6\ :sub:`4`\ 22, P6\ :sub:`5`\ 22
``23``         P23, F23, I23, P2\ :sub:`1`\ 3, I2\ :sub:`1`\ 3
``432``        P432, P4\ :sub:`2`\ 32, F432, F4\ :sub:`1`\ 32, I432, P4\ :sub:`3`\ 32, P4\ :sub:`1`\ 32, I4\ :sub:`1`\ 32
===========    ============

Notice that, in most cases, the correct point group can easily be recognised
from the space group, without memorizing the entire table.

If you are in the fortunate situation of knowing the space group for your
sample before processing the data, look up the point group in the table above
and keep it in mind as you read the next sections.  If you can't find your
space group in the table (for example, *A112*), your source of information is
using a non-standard setting.  Everything should become clear in step 3.

If you don't know the space group, no problem: we will work everything out in
the steps below.


Step 2: Determine the apparent symmetry
=======================================

The orientation of each crystal in your dataset was determined by the indexing
procedure inside ``indexamajig``.  There's a choice of indexing algorithms
which work in many different ways, but they all share one thing in common: they
only look at the positions of the Bragg peaks, not the intensities.

As you should know from basic diffraction theory, the positions of the Bragg
peaks are determined by the translational symmetry of the structure (the
*lattice*), whereas the intensities are determined by the contents of the
unit cell.

This leads to a problem for some symmetry classes.  If the overall crystal
structure, taking into account the unit cell contents, has lower symmetry than
the lattice, there will be an *indexing ambiguity*.  In these cases, the Bragg
peak positions don't provide enough information to correctly determine the
orientation of the crystal.  The results will be an equal mixture of correctly
indexed patterns, and ones where the Miller indices for the reflections are
wrong.  But, we're getting ahead of ourselves...

If we ignore the possibility of indexing ambiguities, then the symmetry of the
intensity dataset we will get by merging the diffraction patterns will be
whatever symmetry the indexing algorithm can discern, which is the symmetry
of the crystal lattice.

The following table shows the possible cases and the symmetry of the lattice
in each case.  Use the furthest down row that is compatible with your data, for
example if the axis lengths are all equal (*a=b=c*) and the angles are all 90°,
you should use ``432``, even though ``32_R``, ``422``, ``222``, ``2`` and ``1``
would all fit.

For this step, what matters are the *approximate* symmetries of the lattice.
You should consider an angle to be equal to 90° if it's within about 1° of that
value, and axis lengths to be equal if they're within about 1% of the same
length.  If ``indexamajig`` gets confused between the axes (shown by double
peaks in the ``cell_explorer`` histograms), then they should be considered
equal.  Conversely, if ``indexamajig`` was able to tell the axes apart (clear,
single peaks for each axis length, with significantly different lengths and
angles), then you can consider them distinct.

The centering of the cell (P, A, B, C, I, F, R or H) is irrelevant at this
step, unless you have "H-centering", which is a special case that we will come
to later.

=================================== =======================
Unit cell parameters                Symmetry of lattice
=================================== =======================
No restrictions                     ``1``
alpha=beta=90°                      ``2``
H-centering and a=b and gamma=120°  ``321_H``
a=b and gamma=120°                  ``622``
a=b=c                               ``32_R``
All angles 90°                      ``222``
All angles 90° and a=b              ``422``
All angles 90° and a=b=c            ``432``
=================================== =======================

Perhaps your cell parameters resemble one of the cases, but with the axes
"re-named".  For example, you might have beta=gamma=90°, alpha≠90°, and all
axes different.  This matches the second row in the table, but with the axes
*a,b,c* re-labelled as *b,c,a*.  We can deal with that, as described in the
next step.


Step 3: Make sure the "unique axis" is correct
==============================================

Let's say your point group is *2*.  In this case, there is a single twofold
axis of rotational symmetry.  The symmetry axis can be along the *a*, *b* or
*c* direction of the lattice - these letters are just the names we use to refer
to the axes.  In theory, you can define the unit cell any way you like, and
CrystFEL will be able to cope (with one exception, mentioned below).  However,
some possibilities are more "conventional" than others, and it can help to
avoid problems if you follow the established conventions.  For example, not all
software can handle all of the possibilities smoothly.  It's also easier to
compare structures when they're described in the same way.

You can tell the direction of the twofold rotation axis, because it has to be
along the axis perpendicular to the angle that isn't 90°.  For example, the
following cell parameters show that the twofold rotation axis is along *b*.
We refer to *b* as the *unique axis*:

a=34 Å, b=123 Å, c=44 Å, alpha=gamma=90°, beta=131°

The following cell has *unique axis a*:

a=92 Å, b=74 Å, c=34 Å, alpha=128°, beta=gamma=90°

However, *a* as the unique axis is a very unconventional situation.  It would
make things easier for yourself to change your target unit cell to make *b* or
*c* the unique axis, and re-run ``indexamajig`` [#f2]_.

**If you've been told that the space group is simply "P2", check carefully to
make sure which convention is being used, because unique axis b or c are
considered equally acceptable.**

If your non-90° angle is very close to 90°, then you should instead be using
point group *222*.  As mentioned above, what matters are the *approximate*
symmetries that can be discerned by the indexing algorithm.

Other types of unit cell have a 'unique' axis, as well.  For example, a
tetragonal cell has all angles 90°, two axes with the same length and one
different.  The different length axis could be labelled as *a*, *b* or *c*.
However, in this case, anything other than unique axis *c* is highly
unconventional.  Nevertheless, check carefully here as well.

When you tell ``process_hkl`` or ``partialator`` the symmetry, you'll need to
tell it the unique axis.  By default, CrystFEL programs assume that the unique
axis is *c*.   If you have anything else, append ``_uaa``, ``_uab`` or ``_uac``
to the point group symbol (from the tables above) to indicate which is the
'unique' axis.  For the first example from above, we would use ``2_uab``:

a=34 Å, b=123 Å, c=44 Å, alpha=gamma=90°, beta=131°

For the tetragonal unit cell parameters shown below, we would use ``422``,
which is a synonym for ``422_uac`` since the unique axis is assumed to be *c*:

a=123 Å, b=123 Å, c=44 Å, alpha=beta=gamma=90°


Step 4: Worry about indexing ambiguities
========================================

At this point, you're in a position to merge your data.  If your prior
information about the point group from step 1 agreed with what you determined
in step 2, then everything is OK and you're finished already!  Simply give the
point group symbol to ``partialator`` or ``process_hkl`` with the ``-y``
argument (or via the CrystFEL GUI).  For example: ``-y 4/mmm``.

However, maybe something is still not right.  Perhaps the structure solution
software is complaining about "twinned data", strange statistics or "poor"
L-test results.  Or, perhaps your prior information about the structure doesn't
match the point group you determined in the previous steps.  In this case, you
may be facing an indexing ambiguity, where the true symmetry is lower than what
can be distinguished by the indexing algorithm.

An *indexing ambiguity* is when the positions of the Bragg peaks do *not* give
sufficient information to uniquely identify the orientation of the crystal.
Instead, there are a small number (usually 2) of possible orientations which
give the *same Bragg peak positions*.  The correct orientation can be
determined by looking at the peak intensities, so it requires a separate
processing step after indexing and integration.

Indexing ambiguities can be resolved in CrystFEL using ``ambigator``.  This
program takes a stream (from ``indexamajig``), works out the correct indexing
assignments, and writes a new stream with the incorrectly assigned reflections
re-labelled with their correct indices.  Here, "correct" means "consistent with
the other patterns in the dataset" - you should keep in mind that the indexing
ambiguity allows separately-processed datasets to have inconsistent labels.

The mechanics of running ``ambigator`` will be described in a separate
document.  However, you will need to know the "real" and "apparent" point
groups.  The apparent point group is the one we already determined.  The real
point group is so far unknown (unless you have prior information!), but there
are a small number of possibilities.  Here they are:

============================  ======================================================
Apparent point group          Real point group
============================  ======================================================
``422``                       ``4``
``32_R`` (rhombohedral axes)  ``3_R`` (rhombohedral axes)
``432``                       ``23``
``622``                       ``3_H`` (hexagonal axes) - double ambiguity, see below
``622``                       ``6``
``622``                       ``312_H`` (hexagonal axes)
``622``                       ``321_H`` (hexagonal axes)
============================  ======================================================

Notice that structures with hexagonal lattices (apparent point group *622*) are
particularly problematic, with quite a large number of real point groups giving
the apparent *622* symmetry.  One of those cases, point group ``3_H`` exhibits
a *double ambiguity* where there are four indexing possibilities for each
pattern, not just the usual two.


Step 5: Add an inversion center to merge Friedel pairs
======================================================

Remember that the point group tells CrystFEL which reflections to consider
as symmetrically equivalent.  The point group you have, at this point, will
*not* include an inversion center, i.e. it will *not* consider reflections
h,k,l and -h,-k,-l as equivalent.  This means that the merging process will
preserve any anomalous signal present in your data.

If you don't expect (or want) an anomalous signal, you can get better results
by merging Friedel pairs of reflections.  This doubles the number of
measurements per symmetrically unique reflection, which can make a large
improvement!  To do this, simply add the missing inversion center to the point
group.  This will change the point group symbol in a way that's not immediately
logical.  The following table shows the results of adding an inversion symmetry
to each of the point groups, so you just have to look up your case.

===========    =================================
Point group    Point group with inversion center
===========    =================================
``1``          ``-1``
``2``          ``2/m``
``222``        ``mmm``
``4``          ``4/m``
``422``        ``4/mmm``
``3_R``        ``-3_R``
``32_R``       ``-3m_R``
``3_H``        ``-3_H``
``321_H``      ``-3m1_H``
``312_H``      ``-31m_H``
``6``          ``6/m``
``622``        ``6/mmm``
``23``         ``m-3``
``432``        ``m-3m``
===========    =================================

The point group symbols in the table above look quite strange.  If you need to
look up one of these symbols in a crystallographic textbook, you just need to
know that the minus signs are supposed to indicate a "bar" over the following
digit.  However, there's usually no need to worry about that.

If you've added a unique axis suffix, add the same suffix to your new point
group.  For example, ``622_uab`` goes to ``6/mmm_uab`` (although, either of
these cases would be considered very unconventional).


Step 6: Extra information about "H cells"
=========================================

A rhombohedral unit cell (all axes the same length, all angles the same but not
equal to 90°) can be represented in two ways.  The first way is using the axes
exactly as just described.  In this case, we talk about "rhombohedral axes" and
use space group symbols *R3* and *R32*.  The second way is to embed the
rhombohedral cell inside a hexagonal unit cell (a=b≠c, alpha=beta=90°,
gamma=120°) while having multiple lattice points (i.e. extra copies of the
crystal structure) within the unit cell.  In this case, we talk about
"hexagonal axes" and use space group symbols *H3* and *H32*.

You will find both representations in space group tables - for example
`here, in the International Tables Volume A <https://it.iucr.org/Ac/ch2o3v0001/sgtable2o3o155/>`_.
Rhombohedral axes are easier to think about, but hexagonal axes are commonly
used for protein structures.  If you've downloaded a rhombohedral structure
from the PDB, it's probably (but not always!) using hexagonal axes.

Different software packages use different conventions for labelling these
cells.  For example, you might also encounter *R3:h* and *R3:r* for hexagonal
and rhombohedral axes respectively.  Unfortunately, sometimes you might even
encounter programs which use *R3* to refer to *hexagonal* axes, and *H3* for
*rhombohedral* axes!  However, you can always tell the difference by looking
at the unit cell parameters.  For some more discussion, including a useful
diagram, see `this classic article
<http://www.phenix-online.org/phenixwebsite_static/mainsite/files/newsletter/CCN_2011_01.pdf#page=12>`_.

The most important thing to keep in mind is that representing the unit cell in
a different way will never change any of the physical properties.  If the
symmetry is *R3* or *H3*, there's an indexing ambiguity, and if it's *R32* or
*H32* then there's no ambiguity. The *R3* and *H3* cases are the same thing, as
are the *R32* and *H32* cases. In both cases, the number of symmetry
equivalents for each reflection is the same.  If there's a strange accidental
indexing ambiguity for one version (see step 7), the same accidental indexing
ambiguity applies to the other version as well.

However, you need to tell CrystFEL which representation you're using.  For all
trigonal point groups - that is, anything with a rhombohedral lattice, or a
hexagonal lattice but no sixfold symmetry - you will need to append either
``_H`` or ``_R`` to the space group symbol.  For example, for point group
*3* on rhombohedral axes, use ``3_R``.  For hexagonal axes, use ``3_H``.

You *cannot* use the unique axis and axis definition suffixes together, for
example ``321_H_uab``.  Always use unique axis *c* for trigonal cells on
hexagonal axes.

There's a further complication.  There are actually two ways that the
rhombohedral cell can be "embedded" into the hexagonal cell.  The two ways are
called *obverse* and *reverse*.  The *International Tables* uses the *obverse*
representation [#f3]_, and so does all the software that I know about.
This complication affects the point group symbol that you must use for space
group *R32*/*H32* (it makes no difference for *R3*/*H3*).  Here are all the
cases for *R32*/*H32*:

============   =========  ================================  ==================
Axes           Setting    Point group as given to CrystFEL  Comment
============   =========  ================================  ==================
Rhombohedral   n/a        ``32_R``
Hexagonal      Obverse    ``321_H``
Hexagonal      Reverse    ``312_H``                         Don't use this one
============   =========  ================================  ==================

Just "for fun", here's the same table for *R3*/*H3*:

============   =========  ================================  ==================
Axes           Setting    Point group as given to CrystFEL  Comment
============   =========  ================================  ==================
Rhombohedral   n/a        ``3_R``
Hexagonal      Obverse    ``3_H``
Hexagonal      Reverse    ``3_H``                           Same as for obverse
============   =========  ================================  ==================

As you can see, your life will be much easier if you just use rhombohedral axes
all the time.  However, due to the prevalence of hexagonal axes in deposited
structures, this is likely to mean that you have to convert from one
representation to the other.  Converting atomic locations (i.e. a structural
model) is outside the scope of CrystFEL, but CrystFEL *can* convert just the
unit cell parameters.  For example, given an "H-centered" unit cell file::

  CrystFEL unit cell file version 1.0

  lattice_type = hexagonal
  centering = H
  unique_axis = c

  a = 66.2 A
  b = 66.2 A
  c = 150.2 A

  al = 90.0 deg
  be = 90.0 deg
  ga = 120.0 deg

CrystFEL's ``cell_tool`` can calculate the rhombohedral representation::

  $ cell_tool -p example.cell --uncenter
  Input unit cell: cell-example.cell
  ------------------> The input unit cell:
  hexagonal H, unique axis c, right handed.
  a      b      c            alpha   beta  gamma
   66.20  66.20 150.20 A     90.00  90.00 120.00 deg
  ------------------> The primitive unit cell:
  rhombohedral R, right handed.                                <<-----------
  a      b      c            alpha   beta  gamma               <<-----------  Look here!
   62.99  62.99  62.99 A     63.40  63.40  63.40 deg           <<-----------
  ------------------> The centering transformation:
  [    1    0    1 ]
  [   -1    1    1 ]
  [    0   -1    1 ]
  ------------------> The un-centering transformation:
  [  2/3 -1/3 -1/3 ]
  [  1/3  1/3 -2/3 ]
  [  1/3  1/3  1/3 ]



Step 7: "It still isn't working!"
=================================

The ambiguities described in step 4 are the most common cases, but there are
more possibilities.  Sometimes, the lattice parameters "accidentally" give rise
to indexing ambiguities.  As noted above, it's the *apparent* symmetries of the
lattice that matter here.  For example, unless the indexing is *very* accurate
(within 1/20 of a degree), the following unit cell will need to be merged with
point group *222* (or *mmm* to merge Friedel pairs), even though it is
technically monoclinic:

a=63 Å, b=82 Å, c=95 Å, alpha=gamma=90°, beta=90.04°

In this case, there will be an indexing ambiguity, because the true symmetry
is *2* (unique axis *b*), but the apparent symmetry is *222*.

Things can get even more complicated than this, and some very "interesting"
ambiguities have turned up over the years.  CrystFEL's ``cell_tool`` utility
can analyse your unit cell and spot possible ambiguities.  See `the manual
<https://desy.de/~twhite/crystfel/manual-cell_tool.html>`_ for details.

Crystal structures seem to have a way of finding new ways to cause trouble.
So, if things are still not working, or if you're just confused, we're happy to
help.  Just send an email!  See the `contact <https://desy.de/~twhite/crystfel/contact.html>`_
page on the CrystFEL website for details.

**Good luck, and may all your indexing be unambiguous!**


.. rubric:: Footnotes

.. [#f1] There are a couple of small exceptions here, when the data is exported
   to XScale or MTZ format.  These formats *require* a space group to be
   nominated, because of the aforementioned reliance on early space group
   nomination.  Here, CrystFEL chooses the lowest-symmetry space group that
   reflects the point symmetry according to which the merging was performed.
   The "downstream" structure solution software should be clever enough to
   assign the correct space group, regardless of what's in the data file.

.. [#f2] It's also possible to change the indexing assignments in the stream
   without re-running indexing, but this could be considered "advanced" usage.
   As mentioned above, it's also possible to continue using the non-standard
   setting, at least as far as CrystFEL is concerned.  However, in that case
   you can expect to have difficulty with other software or when depositing the
   structure.

.. [#f3] If you're interested, this is made explicit in section 2.1.3.6.6 of
   International Tables Volume A (2016 edition), which you can read
   `here <https://it.iucr.org/Ac/ch2o1v0001/sec2o1o3/>`_ (subscription
   required).
