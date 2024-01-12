% indexamajig(1)

NAME
====

indexamajig - bulk indexing and data reduction program


SYNOPSIS
========

indexamajig
-i _filename_ -o _output.stream_ -g _detector.geom_ --peaks=_method_ --indexing=_method_ ...


DESCRIPTION
===========

**indexamajig** takes a list of diffraction snapshots from crystals in random
orientations and attempts to find peaks, index and integrate each one.  The
input is a list of diffraction image files and some auxiliary files and
parameters.  The output is a long text file ('stream') containing the results
from each image in turn.

For minimal basic use, you need to provide the list of diffraction patterns,
the method which will be used to index, a file describing the geometry of the
detector, and a file which contains the unit cell which will be used for the
indexing.  Here is what the minimal use might look like on the command line:

    indexamajig
       -i mypatterns.lst
       -g mygeometry.geom
       --indexing=xgandalf
       --peaks=hdf5
       -p mycell.pdb
       -o test.stream

More typical use includes all the above, but might also include extra
parameters to modify the behaviour. For example, you'll probably want to run
more than one indexing job at a time (-j <n>).

See **man crystfel_geometry** for information about how to create a file
describing the detector geometry and beam characteristics.


DIFFRACTION PATTERN LIST
========================

Indexamajig requires an input file with a list of diffraction patterns to
process. In its simplest form, this is just a text files containing a list of
image filenames. The image files might be in some folder a long way from the
current directory, so you might want to specify a full pathname to be added in
front of each filename. The geometry file includes a description of the data
layout within the files. Indexamajig uses this description to determine the
number of diffraction patterns stored in each file, and tries to process them
all.  You can also specify explicity which framess you would like to process
by putting a string describing the frames after the file name(s) in this list.


PEAK DETECTION
==============

You can control the peak detection on the command line.  First, you can choose
the peak detection method using **--peaks=method**.  **--peaks=hdf5** or
**--peaks=cxi** will take the peak locations from the input file.  See the
documentation for peak_list and peak_list_type in **crystfel_geometry**(5) for
details.

If you use **--peaks=zaef**, indexamajig will use a simple gradient search
after Zaefferer (2000).  You can control the overall threshold and minimum
squared gradient for finding a peak using **--threshold** and
**--min-squared-gradient**.  The threshold has arbitrary units matching the
pixel values in the data, and the minimum gradient has the equivalent squared
units.  Peaks will be rejected if the 'foot point' is further away from the
'summit' of the peak by more than the inner integration radius (see below).
They will also be rejected if the peak is closer than twice the inner
integration radius from another peak.

If you instead use **--peaks=peakfinder8**, indexamajig will use the
"peakfinder8" peak finding algorithm describerd in Barty et al. (2014). Pixels
above a radius-dependent intensity threshold are considered as candidate peaks
(although the user sets an absolute minimum threshold for candidate peaks).
Peaks are then only accepted if their signal to noise level over the local
background is sufficiently high. Peaks can include multiple pixels and the user
can reject a peak if it includes too many or too few. The distance of a peak
from the center of the detector can also be used as a filtering criterion. Note
that the peakfinder8 will not report more than 2048 peaks for each panel: any
additional peak is ignored.

If you instead use **--peaks=peakfinder9**, indexamajig will use the
"peakfinder9" peak finding algorithm described in the master thesis "Real-time
image analysis and data compression in high throughput X-ray diffraction
experiments" by Gevorkov. Other than peakFinder8, peakFinder9 uses local
background estimation based on border pixels in a specified radius
(**--local-bg-radius**). For being fast and precise, a hierarchy of
conditions is used. First condition is only useful for speed consideration, it
demands that a pixel that is the biggest pixel in a peak must be larger than
every border pixel by a constant value (**--min-peak-over-neighbour**).
Second condition ensures, that the pixel passing the previous condition is the
highest pixel in the peak. It assumes, that peaks rise monotonically towards
the biggest pixel. Third condition ensures, that the biggest pixel in the peak
is significantly over the noise level (**--min-snr-biggest-pix**) by
computing the local statistics from the border pixels in a specified radius.
Fourth condition sums up all pixels belonging to the peak
(**--min-snr-peak-pix**) and demands that the whole peak must be
significantly over the noise level (**--min-snr**). Only if all conditions
are passed, the peak is accepted.


INDEXING METHODS
================

You can choose between a variety of indexing methods.  You can choose more than
one method, in which case each method will be tried in turn until one of them
reports that the pattern has been successfully indexed.  Choose from:

**dirax**
: Invoke DirAx.  See Duisenberg, J. Applied Crystallography 25 (1992) p92,
: https://doi.org/10.1107/S0021889891010634.

**mosflm**
: Invoke Mosflm.  See Powell, Acta Crystallographica D55 (1999) p1690,
: https://doi.org/10.1107/S0907444999009506.

**asdf**
: This is a implementation of the dirax algorithm, with some very small changes
: such as using a 1D Fourier transform for finding the lattice repeats.  This
: algorithm is implemented natively within CrystFEL meaning that no external
: software is required.

**felix**
: Invoke Felix, which will use your cell parameters to find multiple crystals in
: each pattern.
: The Felix indexer has been developed by Soeren Schmidt <ssch@fysik.dtu.dk>. To
: use this option, 'Felix' must be in your shell's search path. This can be a
: link to the latest version of Felix. If you see the Felix version information
: when you run Felix on the command line, things are set up correctly.

**xds**
: Invoke XDS, and use its REFIDX procedure to attempt to index the pattern.

**taketwo**
: Use the TakeTwo algorithm.  See Ginn et al., Acta Crystallographica D72
(2016), p956, https://doi.org/10.1107/S2059798316010706.

**xgandalf**
: Use the eXtended GrAdieNt Descent Algorithm for Lattice Finding.
: See Gevorkov et al., Acta Crystallographica A75 (2019) p694,
: https://doi.org/10.1107/S2053273319010593.

**pinkIndexer**
: Use the pinkIndexer algorithm.  See Gevorkov et al., Acta Crystallographica
: A76 (2020) p121, https://doi.org/10.1107/S2053273319015559.

**file**
: See **Re-playing old indexing and injecting external results** below.

Most of the indexing methods require some extra software to be installed,
either at the time of compiling CrystFEL or afterwards.  CrystFEL is
distributed with a script (scripts/install-indexers) which can help you to
quickly install all the required programs.

If you don't specify any indexing methods, indexamajig will try to
automatically determine which indexing methods are available.  You can also
specify indexing method none, in which case no indexing will be done.  This is
useful if you just want to check that the peak detection is working properly.


### Prior unit cell information

You can add one or more of the following to the above indexing methods, to
control what information should be provided to them.  Note that indexamajig
performs a series of checks on the indexing results, including checking that
the result is consistent with the target unit cell parameters.  To get
completely "raw" indexing, you need to disable these checks (see below) and not
provide prior information.

**-latt**
: Provide the Bravais lattice type (e.g. the knowledge that the lattice is
: tetragonal primitive), as prior information to the indexing engine.

**-nolatt**
: The opposite of -latt: do not provide Bravais lattice type information to the
: indexing engine.

**-cell**
: Provide your unit cell parameters as prior information to the indexing
: engine.

**-nocell**
: The opposite of -cell: do not provide unit cell parameters as prior information
: to the indexing engine.

Example: **--indexing=mosflm-cell-latt** means to use Mosflm for indexing, and
provide it with unit cell parameters and Bravais lattice type information.

Usually, you do not need to explicitly specify anything more than the indexing
method itself (e.g. mosflm or asdf).  The default behaviour for all indexing
methods is to make the maximum possible use of prior information such as the
lattice type and cell parameters.  If you do not provide this information, for
example if you do not give any unit cell file or if the unit cell file does not
contain cell parameters (only lattice type information), the indexing methods
you give will be modified accordingly.  If you only specify the indexing
methods themselves, in most cases indexamajig will do what you want and
intuitively expect!  However, the options are available if you need finer
control.


### Post-indexing stages

The indexing results from the indexing engine will be put through a number of
refinement and checking stages.  See the options **--no-check-cell**,
**--no-multi**, **--no-retry** and **--no-refine** below for more details.


### Re-playing old indexing and injecting external results

With **--indexing=file**, indexamajig will read indexing results from a text
file.  Use **--fromfile-input-file** to specify the filename.  Each line of the
file represents one diffraction pattern in one image, and should be formatted
as follows:

    filename frameID asx asy asz bsx bsy bsz csx csy csz xshift yshift latt_cen

To describe multiple overlapping diffraction patterns in one frame, simply
use the same _filename_ and _frameID_ in multiple lines.

The lattice type and centering information are contained in _latt\_cen_.
This should consist of two characters, e.g. **tI** or **aP**.  The first letter
represents the lattice type (**a**, **m**, **o**, **t**, **c**, **h**, **r**
for, respectively, triclinic, monoclinic, orthorhombic, tetragonal, cubic,
hexagonal or rhombohedral).  The centering symbol (**P**, **A**, **B**, **C**,
**I**, **F**, **H** or **R**) follows.

The vector components _asx_,_asy_,... are the reciprocal lattice basis vectors
in reciprocal nanometres.  The _xshift_ and _yshift_ values are the offsets of
the detector position in millimetres.  Usually, these offsets should be zero.

This method can be used to inject indexing results from an external indexing
program.  It can also be used to re-play previous indexing results, which is
useful for changing integration parameters while using slower indexing methods
(such as **pinkIndexer**).  To generate an indexing result from a previous
**indexamajig** run, use the script **stream2sol** in the CrystFEL **scripts**
folder.


REFLECTION INTEGRATION
======================

If the pattern could be successfully indexed, peaks will be predicted in the
pattern and their intensities measured.  You have a choice of integration
methods, and you specify the method using **--integration**.  Choose from:

**rings**
: Use three concentric rings to determine the peak, buffer and background
: estimation regions.  The radius of the smallest circle sets the peak region.
: The radius of the middle and outer circles describe an annulus from which the
: background will be estimated.  You can set the radii of the rings using
: **--int-radius** (see below).  The default behaviour with rings is not to
: center the peak boxes first.  Use rings-cen if you want to use centering.

**prof2d**
: Integrate the peaks using 2D profile fitting with a planar background, close to
: the method described by Rossmann (1979) J. Appl. Cryst. 12 p225.  The default
: behaviour with prof2d is to center the peak first - use prof2d-nocen to skip
: this step.

You can add one or more of the following to the above integration methods:

**-cen**
: Center the peak boxes iteratively on the actual peak locations.  The opposite
: is -nocen, which is the default.

**-sat**
: Normally, reflections which contain one or more pixels above max_adu (defined
: in the detector geometry file) will not be integrated and written to the
: stream.  Using this option skips this check, and allows saturated reflections
: to be passed to the later merging stages.  The opposite is -nosat, which is the
: default for all integration methods.  However, note that the saturation check
: will only be done if max_adu is set in the geometry file.  Usually, it's better
: to exclude saturated reflections at the merging stage.  See the documentation
: for max_adu in crystfel_geometry(5).

**-grad**
: Fit the background around the reflection using gradients in two dimensions.
: This was the default until version 0.6.1.  Without the option (or with its
: opposite, -nograd, which is the default), the background will be considered to
: have the same value across the entire integration box, which gives better
: results in most cases.


BASIC OPTIONS
-------------

**-i filename**, **--input=filename**
: Read the list of images to process from filename.  **--input=-** means to
: read from stdin.  There is no default.

**-o filename**, **--output=filename**
: Write the output data stream to filename.

**-g filename**, **--geometry=filename**
: Read the detector geometry description from _filename_.  See **man
: crystfel_geometry** for more information.

: **--zmq-input=address**
: Receive data over ZeroMQ from address.  The options **--input** and
: **--zmq-input** are mutually exclusive - you must specify exactly one of
: them.  Example: **--zmq-input=tcp://127.0.0.1:5002**.
: If you use this option, you should also use either **--zmq-subscribe** to add
: a ZeroMQ subscription, or **--zmq-request** to define how to request data.

**--zmq-subscribe=tag**
: Subscribe to ZeroMQ message type tag.  You can use this option multiple times
: to add multiple subscriptions.  This option and **--zmq-request** are mutually
: exclusive.

**--zmq-request=msg**
: Request new data over ZeroMQ by sending string msg.  This will cause
: indexamajig's ZeroMQ socket to use REQ mode instead of SUB.  This option and
: **--zmq-subscribe** are mutually exclusive.

**--asapo-endpoint=endpoint**
: Receive data via the specified ASAP::O endpoint.  This option and **--zmq-input**
: are mutually exclusive.

**--asapo-token=token**
: Authentication token for ASAP::O data.

**--asapo-beamtime=beamtime**
: Beamtime ID for ASAP::O data.

**--asapo-source=source**
: Data source for ASAP::O data.

**--asapo-group=group**
: Consumer group name for ASAP::O data.  Concurrent **indexamajig** processes
: working on the same data should use the same value for this, to have the
: data shared between them.

**--asapo-stream=stream**
: Stream name for ASAP::O data.

**--asapo-output-stream=stream**
: Send an output stream via ASAP::O.  For non-hits, a small placeholder will be
: sent.

**--asapo-wait-for-stream**
: If the ASAP::O stream does not exist, wait for it to be appear.  Without this
: option, indexamajig will exit immediately if the stream is not found.

**--no-data-timeout**
: Shut down the entire indexamajig process if the specified number of seconds
: elapse without any data being seen.  This currently applies to ASAP::O data
: only, but might be extended to other streaming systems in future.  The
: default is 60 seconds.

**--asapo-consumer-timeout**
: Set the timeout used for "get next" calls from ASAP::O, in ms.  The default
: is 500 ms.

**--data-format=format**
: Specify the data format for data received over ZeroMQ or ASAP::O.  Possible
: values in this version are msgpack, hdf5 and seedee.

**--basename**
: Remove the directory parts of the filenames taken from the input file.  If
: **--prefix** or -x is also given, the directory parts of the filename will be
: removed before adding the prefix.

**-x prefix**, **--prefix=prefix**
: Prefix the filenames from the input file with prefix.  If **--basename** is also
: given, the filenames will be prefixed after removing the directory parts of the
: filenames.

**-j n**
: Run n analyses in parallel.  Default: 1.  See also **--max-indexer-threads.**
: Tip: use -j `nproc` (note the backticks) to use as many processes as there
: are available CPUs.

**--cpu-pin**
: Pin worker processes to CPUs.  Usually this is not needed or desirable, but in
: some cases it dramatically improves performance.

**--no-check-prefix**
: Don't attempt to correct the prefix (see **--prefix**) if it doesn't look correct.

**--highres=n**
: Mark all pixels on the detector higher than n Angstroms as bad.  This might be
: useful when you have noisy patterns and don't expect any signal above a certain
: resolution.

**--profile**
: Display timing data for performance monitoring.

**--temp-dir=path**
: Put the temporary folder under path.

**--wait-for-file=n**
: Wait at most n seconds for each image file in the input list to be created
: before trying to process it.  This is useful for some automated processing
: pipelines.  It obviously only really works for single-frame files.  If a file
: exists but is not readable when this option is set non-zero, a second attempt
: will be made after ten seconds.  This is to allow for incompletely written
: files.  A value of -1 means to wait forever.  The default value is
: **--wait-for-file=0**.

**--no-image-data**
: Do not load the actual image data (or bad pixel masks), only the metadata.
: This allows you to check if patterns can be indexed, without high data
: bandwidth requirements.  Obviously, any feature requiring the image data,
: especially peak search procedures and integration, cannot be used in this case.
: Therefore, you'll need to get the peaks from somewhere else (see
: **--peaks=msgpack** or **--peaks=hdf5**).


PEAK SEARCH OPTIONS
-------------------

**--peaks=method**
: Find peaks in the images using method.  See the second titled PEAK DETECTION
: (above) for more information.

**--peak-radius=inner,middle,outer**
: Set the inner, middle and outer radii for three-ring integration during the
: peak search.  See the section about PEAK INTEGRATION, above, for details of how
: to determine
: these.  The default is to use the same values as for **--int-radius**.

**--min-peaks=n**
: Do not try to index frames with fewer than n peaks.  These frames will still be
: described in the output stream.  To exclude them, use **--no-non-hits-in-stream**.
: The default is **--min-peaks=0**, which means that all frames will be considered
: hits, even if they have no peaks at all.

**--median-filter=n**
: Apply a median filter with box "radius" n to the image.  The median of the
: values from a (n+1)x(n+1) square centered on the pixel will be subtracted from
: each pixel.  This might help with peak detection if the background is high
: and/or noisy.  The unfiltered image will be used for the final integration of
: the peaks.  If you also use **--noise-filter**, the median filter will be applied
: first.

**--filter-noise**
: Apply a noise filter to the image with checks 3x3 squares of pixels and sets
: all of them to zero if any of the nine pixels have a negative value.  This
: filter may help with peak detection under certain circumstances.  The
: unfiltered image will be used for the final integration of the peaks, because
: the filter is destroys a lot of information from the pattern.  If you also use
: **--median-filter**, the median filter will be applied first.

**--threshold=thres**
: Set the overall threshold for peak detection using **--peaks=zaef** or
: **--peaks=peakfinder8** to thres, which has the same units as the detector data.
: The default is **--threshold=800**.

**--min-squared-gradient=grad**
: Set the square of the gradient threshold for peak detection using **--peaks=zaef**
: to grad, which has units of "squared detector units per pixel".  The default is
: **--min-squared-gradient=100000**.  **--min-sq-gradient** and **--min-gradient** are
: synonyms for this option, however the latter should not be used to avoid
: confusion.

**--min-snr=snr**
: Set the minimum I/sigma(I) for peak detection when using **--peaks=zaef**,
: **--peaks=peakfinder8** or **--peaks=peakfinder9.**  The default is **--min-snr=5**.

**--min-snr-biggest-pix=<n>**
: (peakFinder9 only) min snr of the biggest pixel in the peak, given as a factor
: of the standard deviation. Default is 7.0.

**--min-snr-peak-pix=<n>**
: (peakFinder9 only) min snr of a peak pixel, given as a factor of the standard
: deviation. Should be smaller or equal to sig_fac_biggest_pix. Default is 6.0.

**--min-sig=<n>**
: (peakFinder9 only) minimum standard deviation of the background. Prevents
: finding of peaks in erroneous or highly shadowed unmasked regions. Default is
: 11.0.

**--min-peak-over-neighbour=<n>**
: (peakFinder9 only) just for speed. Biggest pixel must be n higher than the
: pixels in window_radius distance to be a candidate for the biggest pixel in a
: peak. Should be chosen as a small positive number, a few times smaller than the
: weakest expected peak. The default is -INFINITY, which turns off the speedup
: and searches with maximum precision.

**--min-pix-count=cnt**
: Accepts peaks only if they include more than cnt pixels, when using
: --peaks=peakfinder8.  The default is **--min-pix-count=2**.

**--max-pix-count=cnt**
: Accepts peaks only if they include less than cnt pixels, when using
: **--peaks=peakfinder8**.  The default is **--max-pix-count=200**.

**--local-bg-radius=r**
: Radius (in pixels) used for the estimation of the local background when using
: **--peaks=peakfinder8** or **--peaks=peakfinder9**.  The default is
: **--local-bg-radius=3**.

**--min-res=px**
: Only accept peaks if they lay at more than px pixels from the center of the
: detector when using **--peaks=peakfinder8**.  The default is **--min-res=0**.

**--max-res=px**
: Only accept peaks if they lay at less than px pixels from the center of the
: detector when using **--peaks=peakfinder8**.  The default is **--max-res=1200**.

**--no-use-saturated**
: Normally, peaks which contain one or more pixels above max_adu (defined in the
: detector geometry file) will be used for indexing (but not used in the final
: integration - see the section on peak integration above).  Using this option
: causes saturated peaks to be ignored completely.  The opposite is
: **--use-saturated**, which is the default.

**--no-revalidate**
: When using **--peaks=hdf5**, **--peaks=cxi** or **--peaks=msgpack**, the peaks will be put
: through some of the same checks as if you were using **--peaks=zaef**.  These
: checks reject peaks which are too close to panel edges, are saturated (unless
: you use **--use-saturated**), have other nearby peaks (closer than twice the inner
: integration radius, see **--int-radius**), or have any part in a bad region.  Using
: this option skips this validation step, and uses the peaks directly.

**--no-half-pixel-shift**
: CrystFEL considers all peak locations to be distances from the corner of the
: detector panel, in pixel units, consistent with its description of detector
: geometry (see 'man crystfel_geometry').  The software which generates the image
: data files, including Cheetah, may instead consider the peak locations to be
: pixel indices in the data array.  Therefore, the peak coordinates from
: **--peaks=cxi** or **--peaks=hdf5** will by default have 0.5 added to them.  This
: option disables this half-pixel offset.

**--check-hdf5-snr**
: With this option with **--peaks=hdf5** (or cxi or msgpack), the peaks will
: additionally be checked to see that they satisfy the minimum SNR specified with
: **--min-snr**.

**--peakfinder8-fast**
: (peakfinder8 only) Increase speed by restricting the number of sampling
: points used for the background statistics calculation.


INDEXING OPTIONS
----------------

**--indexing=method**
: Index the patterns using method.  See the section titled INDEXING METHODS
: (above) for more information.  The default is to automatically detect which
: indexing methods to use.

**-p unitcell.cell**, **-p unitcell.pdb**, **--pdb=unitcell.pdb**
: Specify the name of the file containing unit cell information, in PDB or
: CrystFEL format.

**--tolerance=tol**
: Set the tolerances for unit cell comparison.  tol takes the form a,b,c,ang.  a,
: b and c are the tolerances, in percent, for the respective reciprocal space
: axes, and ang is the tolerance in degrees for the reciprocal space angles.  If
: the unit cell is centered, the tolerances are applied to the corresponding
: primitive unit cell.  The default is **--tolerance=5,5,5,1.5**.

**--no-check-cell**
: Do not check the cell parameters against the reference unit cell (given with
: -p).  If you've used older versions of CrystFEL, this replaces putting "-raw"
: in the indexing method.

**--no-check-peaks**
: Do not check that most of the peaks can be accounted for by the indexing
: solution.  The thresholds for a successful result are influenced by option
: --multi.

**--multi**
: Enable the "subtract and retry" method, where after a successful indexing
: attempt the spots accounted for by the indexing solution are removed before
: trying to index again in the hope of finding a second lattice.  This doesn't
: have anything to do with the multi-lattice indexing algorithms such as Felix.
: This option also adjusts the thresholds for identifying successful indexing
: results (see **--no-check-peaks**).

**--no-retry**
: Disable retry indexing.  After an unsuccessful indexing attempt, indexamajig
: would normally remove the 10% weakest peaks and try again.  This option
: disables that, which makes things much faster but decreases the indexing
: success rate.

**--no-refine**
: Skip the prediction refinement step.  Usually this will decrease the quality of
: the results and allow false solutions to get through, but occasionally it might
: be necessary.

**--mille**
: Write detector calibration data in Millepede-II format.

**--mille-dir=dirname**
: Write the Millepede-II data into _dirname_.

**--max-mille-level=n**
: Write the Millepede-II data up to a maximum hierarchy depth _n_.  If _n_ is
: 0, only the overall detector position can be refined.  Larger numbers allow
: finer-grained refinement.  Set this to a large number (9 is high enough!)
: to disable the limit, although this can make the Millepede-II files quite
: large.

**--wavelength-estimate=m** **--camera-length-estimate=m**
: Some indexing algorithms need to know the camera length or the wavelength of
: the incident radiation in advance, e.g. to prepare an internal look-up table.
: However, if these values are taken from image headers, then they are not
: available at start-up.  In this case, you will be prompted to add one of these
: options to give approximate values (in metres).  A warning will be generated if
: the actual value differs from this value by more than 10%.

**--max-indexer-threads=n**
: Some indexing algorithms (e.g. pinkIndexer) can use multiple threads for faster
: calculations.  This is in addition to the frame-based parallelism already
: available in indexamajig (see -j).  This option sets the maximum number of
: threads that each indexing engine is allowed to use.  Default: 1.

**--taketwo-member-threshold=n**
: Minimum number of vectors in the network before the pattern is considered
: indexed.  Default 20.

**--taketwo-len-tolerance=n**
: The length tolerance in reciprocal Angstroms.  Default 0.001.

**--taketwo-angle-tolerance=n**
: The angle tolerance in degrees.  Default 0.6.

**--taketwo-trace-tolerance=n**
: The rotation matrix trace tolerance in degrees.  Default 3.

**--felix-domega=n**
: Low-level parameter for the Felix indexing algorithm.

**--felix-fraction-max-visits=n**
: Low-level parameter for the Felix indexing algorithm.

**--felix-max-internal-angle=n**
: Low-level parameter for the Felix indexing algorithm.

**--felix-max-uniqueness=n**
: Low-level parameter for the Felix indexing algorithm.

**--felix-min-completeness=n**
: Low-level parameter for the Felix indexing algorithm.

**--felix-min-visits=n**
: Low-level parameter for the Felix indexing algorithm.

**--felix-num-voxels=n**
: Low-level parameter for the Felix indexing algorithm.

**--felix-sigma=n**
: Low-level parameter for the Felix indexing algorithm.

**--felix-tthrange-max=n**
: Low-level parameter for the Felix indexing algorithm.

**--felix-tthrange-min=n**
: Low-level parameter for the Felix indexing algorithm.

**--xgandalf-sampling-pitch=n**
: Selects how dense the reciprocal space is sampled. [0-4]: extremelyLoose to
: extremelyDense. [5-7]: standardWithSeondaryMillerIndices to
: extremelyDenseWithSeondaryMillerIndices. Default is 6
: (denseWithSeondaryMillerIndices).

**--xgandalf-grad-desc-iterations**
: Selects how many gradient descent iterations are performed. [0-5]: veryFew to
: extremelyMany. Default is 4 (manyMany).

**--xgandalf-tolerance**
: relative tolerance of the lattice vectors. Default is 0.02.

**--xgandalf-no-deviation-from-provided-cell**
: If a prior unit cell was provided, and this flag is set, the found unit cell
: will have exactly the same size as the provided one.

**--xgandalf-min-lattice-vector-length** **--xgandalf-min-lattice-vector-length**
: Minimum and maximum possible lattice vector lengths (unit is A). Used for
: fitting without prior lattice as starting point for gradient descent, so the
: final minimum lattice vector length can be smaller/highier as min/max. Note:
: This is valid for the uncentered cell, i.e. the P-cell! Default is 30A and 250A
: respectively.

**--xgandalf-max-peaks**
: Maximum number of peaks used for indexing. For refinement all peaks are used.
: Peaks are selected by increasing radius. Limits the maximum execution time for
: Patterns with a huge amount of peaks - either real ones or false positives.
: Default is 250.

**--xgandalf-fast-execution**
: Shortcut to set **--xgandalf-sampling-pitch=2 --xgandalf-grad-desc-iterations=3**.

**--pinkIndexer-considered-peaks-count**
: Selects how many peaks are considered for indexing. [0-4] (veryFew to
: manyMany). Default is 4 (manyMany).

**--pinkIndexer-angle-resolution**
: Selects how dense the orientation angles of the sample lattice are sampled.
: [0-4] (extremelyLoose to extremelyDense). Default is 2 (normal).

**--pinkIndexer-refinement-type**
: Selects the refinement type. 0 = none, 1 = fixedLatticeParameters, 2 =
: variableLatticeParameters, 3 = firstFixedThenVariableLatticeParameters, 4 =
: firstFixedThenVariableLatticeParametersMultiSeed, 5 =
: firstFixedThenVariableLatticeParametersCenterAdjustmentMultiSeed.

**--pinkIndexer-tolerance**
: Selects the tolerance of the pinkIndexer (relative tolerance of the lattice
: vectors). Default is 0.06. For bad geometrys or cell parameters use a high
: tolerance. For a well known geometry and cell use a small tolerance. Only
: important for refinement and indexed/not indexed identificaton. Too small
: tolerance will lead to refining to only a fraction of the peaks and possibly
: discarding of correctly indexed images. Too high tolerance will lead to bad
: Fitting in presence of multiples or noise and can mark wrongly-indexed patterns
: as indexed.

**--pinkIndexer-reflection-radius**
: Sets radius of the reflections in reciprocal space in 1/A. Default is 2%% of a*
: (which works quiet well for X-rays). Should be chosen much bigger for electrons
: (~0.002).

**--pinkIndexer-max-resolution-for-indexing**
: Sets the maximum resolution in 1/A used for indexing. Peaks at high resolution
: Don't add much information, but they add a lot of computation time. Default is
: infinity. Does not influence the refinement.

**--pinkIndexer-max-refinement-disbalance**
: Indexing solutions are dismissed if the refinement refined very well to one
: side of the detector and very badly to the other side. Allowed values range
: from 0 (no disbalance) to 2 (extreme disbalance), default 0.4. Disbalance after
: Refinement usually appears for bad geometries or bad prior unit cell
: parameters.

**--asdf-fast**
: This enables a faster mode of operation for asdf indexing, which is around 3
: times faster but only about 7% less successful.


INTEGRATION OPTIONS
-------------------

**--integration=method**
: Integrate the reflections using method.  See the section titled PEAK
: INTEGRATION (above) for more information.  The default is
: **--integration=rings-nocen**.

**--fix-profile-radius=n**, **--fix-divergence=n**
: Fix the beam and crystal paramters to the given values.  The profile radius is
: given in m^-1 and the divergence in radians (full angle).  The default is to
: set the divergence to zero, and then to automatically determine the profile
: radius.
: You do not have to use all three of these options together.  For example, if
: the automatic profile radius determination is not working well for your data
: set, you could fix that alone and continue using the default values for the
: other parameters (which might be automatically determined in future versions of
: CrystFEL, but are not currently).

**--int-radius=inner,middle,outer**
: Set the inner, middle and outer radii for three-ring integration.  See the
: section about PEAK INTEGRATION, above, for details of how to determine
: these.  The defaults are probably not appropriate for your situation.
: The default is **--int-radius=4,5,7**.

**--int-diag=condition**
: Show detailed information about reflection integration when condition is met.
: The condition can be all, none, a set of Miller indices separated by commas,
: random, implausible or negative.  random means to show information about a
: random 1% of the peaks.  negative means to show peaks with intensities which
: are negative by more than 3 sigma.  implausible means to show peaks with
: intensities which are negative by more than 5 sigma.  strong means to show
: peaks with intensities which are positive by more than 3 sigma  The default is
: **--int-diag=none**.

**--push-res=n**
: Integrate n nm^-1 higher than the apparent resolution limit of each individual
: crystal.  n can be negative to integrate lower than the apparent resolution
: limit.  The default is --push-res=infinity, which means that no cutoff is
: applied.  Note that you can also apply this cutoff at the merging stage using
: **process_hkl/partialator --push-res**, which is usually better: reflections
: which are thrown away at the integration stage cannot be brought back later.
: However, applying a resolution cutoff during integration will make the stream
: file significantly smaller and faster to merge.

**--overpredict**
: Over-predict reflections.  This is needed to provide a buffer zone when using
: post-refinement, but makes it difficult to judge the accuracy of the
: predictions because there are so many reflections.  It will also reduce the
: quality of the merged data if you merge without partiality estimation.

**--cell-parameters-only**
: Do not predict reflections at all.  Use this option if you're not at all
: interested in the integrated reflection intensities or even the positions of
: the reflections.  You will still get unit cell parameters, and the process will
: be much faster, especially for large unit cells.

OUTPUT OPTIONS
--------------

**--no-non-hits-in-stream**
: Completely exclude 'non-hit' frames in the stream.  When this option is given,
: frames with fewer than the number of peaks given to **--min-peaks** will not
: have chunks written to the stream at all.  Note that the default value for
: **--min-peaks** is zero, which means all frames will be written to the stream,
: even if they have no peaks at all.

**--copy-header=header**
: Copy the information from header in the image file into the output stream.  For
: HDF5 files, header is interpreted as a path within the file.  This option is
: sometimes useful to allow data to be separated after indexing according to some
: condition such the presence of an optical pump pulse.  You can give this option
: as many times as you need to copy multiple bits of information.
: The old option **--copy-hdf5-field** is a synonym for this option.

**--no-peaks-in-stream**
: Do not record peak search results in the stream.  You won't be able to check
: that the peak detection was any good, but the stream will be around 30%
: smaller.

**--no-refls-in-stream**
: Do not record integrated reflections in the stream.  The resulting output won't
: be usable for merging, but will be a lot smaller.  This option might be useful
: if you're only interested in things like unit cell parameters and orientations.

**--serial-offset=n**
: Start the serial numbers in the stream at n instead of 1.  Use this if you are
: splitting an indexing job up into several smaller ones, so that the streams can
: be concatenated into a single one with consistent numbering.  This is important
: if you use **whirligig**.

**--harvest-file=fn**
: Write a list of parameters to fn, in JSON format.  This is intended to be used
: for harvesting data into a database system.  This option has no effect if
: --serial-offset is set to a number larger than 1, to avoid the file being
: overwritten multiple times in a batch system.


HISTORICAL OPTIONS
------------------

**--no-sat-corr**
: This option is here for historical purposes only, to disable a correction which
: is done if certain extra information is included in the image data file.


IDENTIFYING SINGLE PATTERNS IN THE INPUT FILE
=============================================

By default indexamajig processes all diffraction patterns ("events") in each of
the data files listed in the input list. It is however, possible, to only
process single events in a multi-event file, by adding in the list an event
description string after the data filename. The event description always
includes a first section with alphanumeric strings separated by forward slashes
("/") and a second section with integer numbers also separated by forward
slashes. The two sections are in turn separated by a double forward slash
('//'). Any of the two sections can be empty, but the double forward slash
separator must always be present.  Indexamajig matches the strings and the
numbers in the event description with the event placeholders ('%') present
respectively in the 'data' and 'dim' properties defined in the geometry file,
and tries to retrieve the full HDF path to the event data and the the its
location in a multi-dimensional data space. Consider the following examples:

### Example 1

The 'data' and 'dim' properties have been defined like this in the geometry file:

    data = /data/%/rawdata
    dim0 = ss
    dim1 = fs

The event list contains the following line:

    filename.h5  event1//

This identifies an event in the 2-dimensional data block located at
/data/event1/rawdata in the HDF5 file called filename.h5.

### Example 2

The 'data' and 'dim' properties have been defined like this in the geometry
file:

    data = /data/rawdata
    dim0 = %
    dim1 = ss
    dim2 = fs

The event list contains the following line:

    filename.h5  //3

This identifies an event in the 3-dimensional data block located at
/data/rawdata in the HDF5 file called filename.h5, specifically the
2-dimensional data slice defined by the value 3 of the first axis of the data
space.

Indexamajig tries to match the alphanumerical strings to the placeholders in
the 'dim' property defined in the geometry file. The first string is matched to
the first placeholder, the second to the second placeholder, and so on. A
similar strategy is followed to match integer numbers to the placeholders in
the 'dim' property defined in the geometry file. For a full explanation of how
the internal layout of the data file can be  described in the geometry file,
please see **man crystfel_geometry**.

You can use **list_events** to prepare a list of each event in one or more
input files.  Note that you only need to do this if you need to perform some
sorting or filtering on this list.  If you want to process every event in a
file, simply specify the filename in the input file.


AUTHOR
=======

This page was written by Thomas White, Yaroslav Gevorkov and Valerio Mariani.


REPORTING BUGS
==============

Report bugs to <taw@physics.org>, or visit <http://www.desy.de/~twhite/crystfel>.


COPYRIGHT AND DISCLAIMER
========================

Copyright © 2012-2023 Deutsches Elektronen-Synchrotron DESY, a research centre
of the Helmholtz Association.

indexamajig, and this manual, are part of CrystFEL.

CrystFEL is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CrystFEL is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CrystFEL.  If not, see <http://www.gnu.org/licenses/>.


SEE ALSO
========

**crystfel**(7), **crystfel_geometry**(5), **cell_explorer**(1),
**process_hkl**(1), **partialator**(1), **list_events**(1), **whirligig**(1)
