=====================================
How to increase data processing speed
=====================================

You want ``indexamajig`` to run faster?  You're probably already using ``-j``,
causing it to divide its work between parallel processes.  Maybe you're even
already using a compute cluster with a batch system, via the GUI or the
``turbo-index-slurm`` and ``turbo-index-lsf`` scripts.  But you want even more
speed?  Here are some tips for getting things to run as fast as possible:


Compile CrystFEL and dependencies with optimisations
====================================================

Note that CMake's default is to compile *without* optimisations.  You need to
add the option ``-DCMAKE_BUILD_TYPE=Release`` (or ``RelWithDebInfo``) to your
CMake invokation to tell it to enable optimisations.  In CrystFEL, it's
particularly important to do this for the HDF5 compression plugins (this makes
a factor of 3 difference in decompression speed!), XGandalf and PinkIndexer.
If you compile CrystFEL itself using CMake (not recommended - use Meson
instead), then you should add the option here as well.


Tune or avoid compression
=========================

Data compression always trades speed for disk space.  For the highest speed,
disable it altogether.  Obviously, there needs to be a trade-off with available
disk space.

When compressing data in HDF5, pay careful attention to the
`chunk size <https://support.hdfgroup.org/HDF5/doc/Advanced/Chunking/>`_.
A badly selected chunk size can cause a very large slowdown.


Bin the pixel data
==================

If you're using a high-resolution detector such as an Eiger 16M, consider
whether you really need the full resolution or not.  Most experiments don't
need anything close to 16 megapixel resolution.  If not, bin the detector
frames down to 4M or even 1M.  This makes a huge difference because the peak
search algorithm must look at all pixels, so binning your data from 16M to 4M
can make it four times faster.  Note that the peak search is one of the only
processing stages which needs to be done on every single frame, hit or non-hit!


Avoid x/y bad regions
=====================

For a similar reason, avoid defining bad regions in x/y coordinates.  If you
can, define them in fs/ss coordinates instead, or use in-band bad pixel flags
(i.e. set the bad pixel values to NaN).  If you specify bad regions in x/y
coordinates, CrystFEL has to figure out which detector pixels fall into the
specified area in the lab coordinate system, for which it (currently) uses a
slow brute-force algorithm.


Avoid bad pixel masks
=====================

In many cases, e.g. Pilatus and Eiger detectors, the bad pixel information is
included in the image data itself, so there's no need to put the information in
a separate file.  Bad pixels have a special flag value, usually 65535.  With
recent versions, you can tell CrystFEL to take note of these values using
``flag_morethan = 65535`` in the geometry file.


Skip non-hits
=============

Use ``--min-peaks``, so that only plausible hits get processed.  At the same
time, add ``--no-non-hits-in-stream`` so that time isn't wasted recording
information about non-hits.


Choose the fastest peak search algorithms
=========================================

If the background is low and/or smooth, you can use the faster ``zaef`` peak
search algorithm instead of ``peakfinder8`` without compromising on the
results.


Choose the fastest indexing algorithms
======================================

Don't use PinkIndexer, unless you really need it (wide bandwidth or electron
diffraction data).  PinkIndexer is a very general and accurate indexing
algorithm, but these advantages must be "paid for" in speed.  Prefer DirAx,
TakeTwo, Mosflm and XGandalf.


Try less hard to index each frame
=================================

The default behaviour is to try very hard to index each frame: all indexing
methods will be tried up to six times, deleting the weakest peaks after each
unsuccessful attempt, and trying again with the leftover peaks after a
successful attempt.  If you enable a large number of indexing methods, this can
add up to over 30 attempts to index each frame!  The options ``--no-retry``
and ``--no-multi`` will disable this behaviour.  In addition, you should
reduce the number of indexing methods in operation: ``xgandalf`` alone is a
good choice.

Of course, doing the above will probably decrease the fraction of indexed
frames somewhat, but the trade-off might be positive for your data.


Integrate to lower resolution
=============================

Restrict the resolution of data for integration by setting
``indexamajig --push-res``.  This affects data quality, so you will need to
try different values to find the best one.  Start with ``--push-res=1.5``,
which will cause spots to be integrated up to 1.5 nm :sup:`-1` higher than the
conservatively-estimated resolution of each diffraction pattern.  This gives a
reasonable balance between integrating weak "invisible" high-resolution data,
and not including too much "junk" data.  If the metrics for the final merged
data suggest that there might be more information at higher resolution, use a
larger value.

Normally, we recommend limiting the resolution only at the merging stage
(``partialator --push-res``), because this gives you the most flexibility - you
can set any ``--push-res`` value without re-integrating the entire dataset.  If
you limit the resolution at the integration stage, the number of reflections to
be integrated will be much smaller, which can lead to a significant speed
improvement.  However, the ``--push-res`` value that you use for merging must
be smaller than the value used for integration.


Use a static detector geometry
==============================

CrystFEL geometry files allow some aspects of the geometry to come from the
data files, such as the panel z-positions ("clen"/camera length) and overall
detector shifts.  If you can instead give fixed numerical values for
everything, then some parts of CrystFEL can prepare calculations in advance.
In some cases, this can make a significant speed improvement.
