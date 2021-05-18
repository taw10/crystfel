==================================================
Processing electron diffraction data with CrystFEL
==================================================

Electron diffraction data is fully supported in recent versions of CrystFEL.
For some background, see these papers:

* R. Bücker, P. Hogan-Lamarre, P. Mehrabi, E. C. Schulz et al. "Serial protein
  crystallography in an electron microscope". Nature Communications 11 (2020)
  996.
  `doi:10.1038/s41467-020-14793-0 <https://doi.org/10.1038/s41467-020-14793-0>`_
* Robert Bücker, Pascal Hogan-Lamarre, R. J. Dwayne Miller. "Serial Electron
  Diffraction Data Processing with diffractem and CrystFEL"
  https://arxiv.org/abs/2011.02977

Here are some tips to get you started with your own data.


Setting the wavelength
======================

You can set the accelerating voltage in the geometry file::

  electron_voltage 200 kV

This should be the accelerating voltage, not the relativistically-corrected
energy of the accelerated electrons (the difference is small).  Alternatively,
give the wavelength directly::

  wavelength 2.51e-12 m

You can also use ``V`` instead of ``kV`` and ``A`` (Angstroms) instead of
``m``.

Obviously, do not use ``photon_energy`` - that's only for electromagnetic
radiation.


Lens distortions and variable direct beam position
==================================================

Depending on the method of data acquisition, the central beam may vary from
frame to frame quite a lot.  These shifts will need to be determined and stored
in the image headers as part of a pre-processing step such as that performed by
`diffractem <https://github.com/robertbuecker/diffractem>`_.  Specify the
location of the shifts in the geometry file, for example::

  detector_shift_x = /%/shots/shift_x_mm mm
  detector_shift_y = /%/shots/shift_y_mm mm

The location (``/%/shots/shift_x_mm``) will, of course, depend on your
pre-processing.  These are the correct locations for data from ``diffractem``.

Note that the coordinates are shifts of the *detector* (not the *beam*) in the
*laboratory coordinate system* (see ``man crystfel_geometry``).  Shifts in
"detector coordinates" (i.e. ``detector_shift_fs``) would only make sense for
single-panel detectors and are therefore not implemented.

Correcting elliptical distortions of the projector lens can be done quite well
by setting the panel direction vectors (see ``man crystfel_geometry`` again) to
be non-orthogonal.  For example::

  p0/fs = +0.9934x -0.0092y
  p0/ss = -0.0092x +1.0067y

This is also covered by ``diffractem``.


Indexing patterns
=================

The only indexing method suitable for electrons is PinkIndexer, but more
methods might be added in the future.  As well as selecting the indexing method
using ``--indexing=pinkindexer``, you will need to configure PinkIndexer's
assumption about the size of a reflection in reciprocal space using
``--pinkIndexer-reflection-radius=0.003``::

  indexamajig -g sample.geom -i input-files.lst -o sample.stream \
              --indexing=pinkindexer --pinkIndexer-reflection-radius=0.003 \
              --peaks=zaef --threshold=300 --min-squared-gradient=50000

That's it!  You will have to tune the peak search parameters as usual, of
course, but this should already be enough to get started with electron
diffraction data.

The default parameters give quick and reasonably good indexing.  You can get
more accurate results by increasing the granularity of refinement within
Pinkindexer using ``--pinkIndexer-refinement-type=N``, with ``N`` from 0 (least
precise and fastest) to 5 (most precise and slowest).  The default value is 1.

For speed, you might want to add ``--max-indexer-threads=N``, where ``N`` is
the number of threads that PinkIndexer should use within itself.  You should
also use ``-j M`` to specify the number of frames that ``indexamajig`` should
process in parallel.  The sum of ``N`` and ``M`` should be roughly equal to the
number of available CPU cores.  It's also helpful to disable the
default``indexamajig`` behaviour of making multiple attempts to index each
pattern (``--no-retry``).

If the camera length is specified in the geometry file by referencing values
from image file headers, then you will also need to give an estimate of the
camera length using ``--camera-length-estimate``.  The same applies if the
wavelength is specified indirectly, in which case you will need
``--wavelength-estimate``.  However, neither of these are usually necessary for
electron diffraction data.


Saving and re-loading indexing results
======================================

Indexing patterns with PinkIndexer can take a long time, which makes it painful
to try (for example) different integration parameters.  To mitigate this, you
can store the results from one run of ``indexamajig`` and re-run the process
using the previous indexing results instead of repeating the PinkIndexer step.

First, create a *solution file* from the stream using the ``stream2sol.py``
program from the CrystFEL scripts folder::

  stream2sol.py -i sample.stream -o sample.sol

Then, run ``indexamajig`` using ``--indexing=file`` and giving the filename::

  indexamajig -g sample.geom -i input-files.lst -o output.stream \
              --indexing=file --fromfile-input-file=sample.sol \
              --no-refine --no-check-peaks --no-retry \
              --int-radius=10,12,14

Note the options ``--no-refine --no-check-peaks --no-retry`` are added to skip
the checks and refinement stages normally performed after indexing - these are
not necessary when replaying solutions like this.
