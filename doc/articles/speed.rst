=====================================
How to increase data processing speed
=====================================

You want ``indexamajig`` to run faster?  You're probably already using ``-j``,
causing it to divide its work between parallel processes.  Maybe you're even
already using a compute cluster with a batch system, via the GUI or the
``turbo-index-slurm`` and ``turbo-index-lsf`` scripts.  But you want even more
speed?  Here are some tips for getting things to run as fast as possible:

* If you're using a high-resolution detector such as an Eiger 16M, consider
  whether you really need the full resolution or not.  Most experiments don't
  need anything close to 16 megapixel resolution.  If not, bin the detector
  frames down to 4M or even 1M.  This makes a huge difference because the peak
  search algorithm must look at all pixels, so binning your data from 16M to 4M
  can make it four times faster.  Note that the peak search is one of the only
  processing stages which needs to be done on every single frame.

* For a similar reason, avoid defining bad regions in x/y coordinates.  If you
  can, define them in fs/ss coordinates instead, or use in-band bad pixel
  flags (i.e. set the bad pixel values to NaN).  If you specify bad regions in
  x/y coordinates, CrystFEL has to figure out which detector pixels fall into
  the specified area in the lab coordinate system, and currently this uses a
  slow method.

* Use ``--min-peaks``, so that only plausible hits get processed.  At the same
  time, add ``--no-non-hits-in-stream`` so that time and disk space isn't used
  for recording information about non-hits.

* Avoid HDF5 files which use compression, as well as the "virtual data set"
  feature.  Obviously, this may involve a trade-off with disk space and data
  organisation.

* Don't use PinkIndexer, unless you really need it (wide bandwidth or electron
  diffraction data).  PinkIndexer is a very general and accurate indexing
  algorithm, but its generality makes it very slow.  Prefer DirAx, TakeTwo,
  Mosflm and XGandalf.

* Use ``-no-retry`` and ``--no-multi``.  This will probably decrease fraction
  of indexed frames somewhat, but the trade-off might be positive for your data.

* Restrict the resolution of data for integration by setting ``indexamajig
  --push-res``.  This affects data quality, so you will need to try different
  values to find the best one.  Normally, we recommend limiting the resolution
  only at the merging stage (``partialator --push-res``), because this gives
  you the most flexibility - you can set any ``--push-res`` value without
  re-integrating the entire dataset.  If you limit the resolution at the
  integration stage, the number of reflections to be integrated will be much
  smaller, which can lead to a significant speed improvement.  However, the
  ``--push-res`` value that you use for merging must be smaller than the value
  used for integration.  If you don't know where to start with this option, try
  ``--push-res=1.5``.

