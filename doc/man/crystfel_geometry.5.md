% crystfel_geometry(5)

INTRODUCTION
============

A CrystFEL "geometry file", usually named _something.geom_, is used by programs
in the CrystFEL suite to get information about:

* The physical position of the detector, including all sub-detectors if
  applicable.

* The layout of data in the file, especially if a "container format" such as
  HDF5 is used.

* The values of parameters such as the incident radiation wavelength, or a
  description of where to get these values in the data files.


COORDINATE SYSTEM
=================

It's important to distinguish between the **data coordinate system** and the
**laboratory coordinate system**.

The **data coordinate system** consists of the **fast scan** and **slow scan**
directions.  **Fast scan** refers to the direction whose coordinate changes
most quickly as the bytes in the data file are moved through.  **Slow scan**
refers to the other dimension of the 2D data array.  Arrays in the data file
can have more than two dimensions - see secion **DATA DIMENSIONS** below.

CrystFEL's **laboratory coordinate system** defined as follows:

* +z is the beam direction, and points along the beam (i.e. away from the source)

* +y points towards the zenith (ceiling).

* +x completes the right-handed coordinate system.

The CrystFEL GUI shows +x horizontally (left to right), and +y vertically
(bottom to top).  This means that the GUI shows images from the "into the beam"
perspective.


PANELS AND GEOMETRY FILE SYNTAX
===============================

CrystFEL's representation of a detector is broken down into one or more
**panels**, each of which has its own position, size and various other
parameters.

Lines in a CrystFEL geometry file have the following general form:

    parameter = value

Many parameters, however, are specified for a specific panel, as follows:

    panel_name/parameter = value

Most parameters can theoretically be specified on a per-panel level, but in
practice are the same for all panels.  For example, the pixel size is usually
the same for all panels.  In this case, specify the parameter once, without
a panel name, at the top of the geometry file.  The value will then be used
for all panels which are *first mentioned* later in the file.

The panel names can be anything of your choosing, except that the names must
not start with **bad**, **group** or **rigid_group**.  These are reserved for
other purposes (see below).

You can also add comments, for example:

    ; This is a comment
    mypanel/min_fs = 34    ; This is also a comment


SPECIFYING THE WAVELENGTH AND BEAM PARAMETERS
=============================================

To specify the incident beam wavelength, you need to use one of the following
forms (only one):

**wavelength** = _nnn_ [**m**|**A**]
: Specifies a wavelength directly, in meters or Angstroms according to the
: suffix.

**photon_energy** = _nnn_ [**eV**|**keV**]
: Specifies the energy per photon of electromagnetic radiation (e.g. X-rays).
: If no units suffix is given, **eV** (electron volts) will be assumed.

**electron_voltage** = _nnn_ [**V**|**kV**]
: Specifies the accelerating voltage in an electron microscope.  This should be
: the accelerating voltage, not the relativistically-corrected energy of the
: accelerated electrons (the difference is small).

For all of these, a data location in the input files can be given instead of
a literal number.  For example, with data in HDF5 format:

    wavelength = /data/LCLS/wavelength m

If there are multiple frames per input file, the program will do what you
expect depending on the type of `/data/LCLS/wavelength`.  For example, a scalar
value in the metadata will be applied to all frames, or an array of values can
be used to provide a separate wavelength for each frame.

You can also specify the radiation bandwidth, as follows:

**bandwidth** = _bw_
: The bandwidth of the radiation, expressed as a fraction of the wavelength.
: The bandwidth will be interpreted as the standard deviation of a Gaussian
: spectrum, and used for calculating reflection positions.


PHYSICAL PANEL LOCATIONS
========================

For each panel, the physical location is controlled by the **fs**, **ss**,
**corner_x**, **corner_y** and **clen** parameters, for example:

    q3a15/fs = -0.0010058292x +0.9999995232y
    q3a15/ss = -0.9999995232x -0.0010058292y
    q3a15/corner_x = 575.475
    q3a15/corner_y = -221.866
    q3a15/coffset = 0.01


**fs**, **ss**
: The vectors in the lab coordinate system of the fast and slow scan directions
: of the panel data, measured in pixels.  Inclusion of a component in the z
: direction means that the panel is not perpendicular to the X-ray beam.

**corner_x**, **corner_y**
: The position in x and y (lab coordinate system) of the corner of this panel.
: The corner is defined as the first point in the panel to appear in the image
: data. The units are pixel widths of the current panel.  This should be the
: location of the very corner of the panel (not the center of the first pixel).

**coffset**
: The offset of the panel along the z-direction from the position given by
: **clen**.

The overall detector position in the z-direction is given by **clen**, which
can only be specified once in the geometry file (not for each panel):

**clen** = _nnn_ [mm|m]
: The overall z-position ("camera length") for the detector, or a data file
: location, and units (m or mm).  If no units are given, **m** will be assumed
: if the value is given as a literal number, or **mm** if it's a data file
: header location.  This discrepancy is for historical reasons, and you should
: always specify the units.  Like when specifying the wavelength (see above),
: CrystFEL should do what you expect with multi-frame data files.

Getting the **clen** value from the image headers gives the illusion of
avoiding the need to create a new geometry file every time you change the
detector position.  However, in practice this doesn't work very well because
the detector movement direction is usually not exactly parallel to the beam
axis.  That means that the beam center position varies with the camera length,
and you would have to prepare a new geometry file for each position anyway.
For the best results, simply make sure that the experimental geometry stays as
static as possible.

### Per-frame beam center position

You can specify an overall detector shift.  The most common use of this is for
serial scanning electron diffraction experiments, where the beam center moves
from frame to frame.  You should avoid using this feature, especially with a
metadata location, because it limits CrystFEL's ability to pre-calculate
certain data structures.  This feature is likely to be removed in future
versions.

**detector_shift_**[x,y] = _nnn_ [m|mm]
: These specify that the entire detector should be shifted by this amount in the
: x and y directions.  The units should be specified as m or mm.  If units are
: not specified, the value will be taken as metres.  nnn can be a file metadata
: location (e.g. an HDF5 path).


PANEL DATA LOCATIONS
====================

For each panel, you have to specify where to find the data in the file.  Don't
forget that many of these parameters will be the same for each panel, so you
can set them once at the top (see the section **Panels and geometry file
syntax**).

**data** = _/location/of/data_
: The location in the data file of the array that contains the panel's data.
: The interpretation of the value depends on the file type.  If you're using
: HDF5 files, it will be a path such as `/data/run_1/imagedata`.  The default
: value is `/data/data`.

**min_fs**, **max_fs**, **min_ss**, **max_ss** = _nnn_
: The range of pixels in the data block that correspond to this panel.  Often,
: multiple panels are grouped together into one "slab".  The pixel ranges are
: in the *data coordinate system*, and are specified *inclusively*.

### Multiple frames per file

The best performance is achieved when each file on disk contains a large number
of images, rather than just one.  In this case, you have to additionally
specify the data layout to CrystFEL.

Consider a file format where each frame has its own data array under a separate
name, for example an HDF5 file with the following layout:

    /data/run_1/image_1/data
    /data/run_1/image_2/data
    /data/run_1/image_3/data
    /data/run_1/image_4/data
    /data/run_1/image_5/data
    ...

In this case, you can use **%** as a placeholder in the data location.
For example:

    data = /data/run_1/%/data

For HDF5 files, the **%** must be a whole name at a certain hierarchy level,
i.e. the following is not allowed:

    data = /data/run_1/image_%/data  ;; This won't work

Next, consider a file format where the frames and detector panels are grouped
into a four-dimensional array.  The first two dimensions are the image data
axes, the third dimension is the panel number, and the fourth dimension is the
frame number.  In this case, you can specify the data location as follows:

    data = /data/run_1/image4Darray
    dim0 = %
    dim2 = ss
    dim3 = fs
    min_fs = 0
    max_fs = 255
    min_ss = 0
    max_ss = 255

    panel1/dim1 = 0
    panel1/corner_x = ....
    ...

    panel2/dim1 = 1
    panel2/corner_x = ....
    ...

    panel3/dim1 = 2
    panel3/corner_x = ....
    ...

**dim** values can be a literal number, a placeholder (**%**), or **fs** or
**ss**.  Note that, in this example, the common parameter values have been
placed at the top, avoiding some repetition.

CrystFEL assumes that the data block defined by the 'data' property has a
dimensionality corresponding to the axis with the highest value of n defined by
the 'dim' property.  That is, if the geometry file specifies dim0, dim1 and
dim2, then the data block is expected to be three-dimensional.  The size of the
data block along each of those axes comes from the image metadata (e.g. the
array sizes in the HDF5 file).

The lowest number of n corresponds to the most slowly-changing array index as
the data block is traversed.  The default values are dim0=ss and dim1=fs.  The
value of n corresponding to fs must not be lower than the value assigned to ss,
i.e. "fast scan is always fast scan".


PEAK LISTS
==========

It's possible to include lists of peak positions into the data file.  If the
data is pre-processed using a "hit finding" procedure, usually a peak search
will already have been performed.  It makes sense to re-use these peak search
results, instead of performing a new peak search inside CrystFEL.

In this case, you need to specify the location of the peak list in the data
file, and the format of the peak list.

**peak_list** = _loc_
: Peak list location in the data files.

**peak_list_type** = _type_
: Specify the layout of the peak list.  Allowed values are **cxi**, **list3**
: and **auto**.

The possible list types are:

**list3**
: The peak list is a two dimensional array whose size in the first dimension
: equals the number of peaks and whose size in the second dimension is exactly
: three.  The first two columns contain the fast scan and slow scan coordinates,
: the third contains the intensities.  This is the correct option for
: "single-frame" HDF5 files as written by older versions of Cheetah.

**cxi**
: The peak list is an HDF5 group containing four separate HDF5 datasets: nPeaks,
: peakXPosRaw, peakYPosRaw and peakTotalIntensity.  See the specification for the
: CXI file format at http://www.cxidb.org/ for more details.  This is the correct
: option for "multi-event" HDF5 files as output by recent versions of Cheetah.

**auto**
: CrystFEL will decide between the above options based on the file extension.

### Important note about coordinate conventions

Note that CrystFEL considers all peak locations to be distances from the corner
of the detector panel, in pixel units, consistent with its description of
detector geometry (see the section about **corner_x** above).  The software
which generates the HDF5 or CXI files, including Cheetah, may instead consider
the peak locations to be pixel indices in the data array.  In the former case,
a peak position (0,0) corresponds to the very corner of the detector panel.  In
the latter case, position (0,0) corresponds to the center of the first pixel,
and the very corner would be (-0.5,-0.5).

To compensate for this discrepancy, CrystFEL will, by default, add 0.5 to all
peak coordinates. See the **indexamajig** option **--no-half-pixel-shift** if
this isn't what you want.


PIXEL SIZE
==========

You will need to specify the size of the pixels, of course.  Use one of the
following:

**pixel_pitch** = _pixelSize_
: The width of the pixels, in meters.

**res** = _pixelsPerMeter_
: The resolution, in pixels per metre, i.e. one divided by the pixel size in
: metres.

These values effectively give the scale factor between the length of the
**fs,ss** vectors and physical space.  If the **fs** and **ss** vectors have
different magnitudes, the pixels will not be square.  This is allowed, but
comes with a possibility of strange problems, because many algorithms assume
square pixels.


DETECTOR GAIN
=============

CrystFEL needs to know the gain of the detector, in order to determine how
many photons correspond to a particular signal level and hence calculate error
estimates on the intensity values.  These gain values are **not** used to
correct the pixel values for different gains among the panels.

Use one of the following:

**adu_per_photon**
: The number of detector intensity units which will arise from one quantum of
: intensity (one X-ray photon, or one electron in an electron microscope).

**adu_per_eV**
: The number of detector intensity units which will arise from a 1 eV photon
: of electromagnetic radiation.  This will be scaled by the photon energy
: (see **photon_energy**) to calculate the intensity per photon at the
: wavelength used by the experiment.  This option should only be used for
: electromagnetic radiation.


DETECTOR SATURATION
===================

You can specify the saturation value in the geometry file, which will allow
**indexamajig** to avoid integrating saturated reflections.  However, usually
it's best to include all reflections at this stage, and exclude the saturated
reflections at the merging stage (see **process_hkl** and **partialator**
options **--max-adu**).

**max_adu**
: The saturation value for the panel.  A warning will be displayed if you use
: this option, because it's better to exclude saturated reflections at the
: merging stage.

Some combinations of detectors and processing methods result in the saturation
level varying pixel-to-pixel.  For this case, you can provide a per-pixel map
of saturation values.  Note that **both** the map values and the **max_adu**
values will both be honoured.

**saturation_map**
: This specifies the location of the per-pixel saturation map in the data file.

**saturation_map_file**
: Specifies that the saturation map should come from the file named here,
: instead of the file being processed.  This can be an absolute filename or
: relative to the working directory.


BAD REGIONS
===========

"Bad region" refers to any set of pixels that should be completely ignored by
CrystFEL.  There are multiple ways to mark pixels as bad.

### Marking a whole panel

To flag all pixels in one panel as bad, simply set the **no_index** parameter:

**no_index**
: If set to **true** or any numerical value other than 0, indicates that the
: panel should be ignored.  The slightly misleading name is for historical
: reasons.

### Marking pixels at the panel edges

With many detectors, the pixels at the edge of the detector panels behave
differently and should be masked out.

**mask_edge_pixels** = _n_
: Mark a border of _n_ pixels around the edge of the panel as bad.

### Marking pixels according to value

Many data files contain information about bad pixels encoded in the pixel
values, for example a value of 65535 often indicates a bad pixel.

**flag_lessthan** = _n_
: Mark pixels as bad if their value is less than _n_.

**flag_morethan** = _n_
: Mark pixels as bad if their value is more than _n_.

**flag_equal** = _n_
: Mark pixels as bad if their value exactly _n_.

Note carefully that the inequalities are strict, not inclusive: "less than",
not "less than or equal to".

Note also that **flag_equal** will be difficult to use for data in
floating-point format.  With floating-point data, you should use
**flag_lessthan** and **flag_morethan**.

### Marking pixels in rectangles

You can specify a range of pixels to ignore in the *data coordinate system* or
the *laboratory coordinate system*.

To mask pixels in the *data coordinate system*, use the following syntax:

    badregionB/min_fs = 128
    badregionB/max_fs = 160
    badregionB/min_ss = 256
    badregionB/max_ss = 512
    badregionB/panel = q0a1

A bad region is distinguished from a panel because it starts with **bad**.
Apart from that, the region can use any name of your choice.

The pixel ranges are specified *inclusively*.  The *panel* name has to be
specified, because the pixel range alone might not be unique (see section
**Multiple frames per file**).  Bad regions specified in this way therefore
cannot stretch across multiple panels.

To mask pixels in the *laboratory coordinate system*, use the following syntax:

    badregionA/min_x = -20.0
    badregionA/max_x = +20.0
    badregionA/min_y = -100.0
    badregionA/max_y = +100.0

In this case, the panel name is not required, and the bad region can span
multiple panels.  However, bad regions specified in laboratory coordinates take
longer to process (when loading images) than regions specified in fs/ss (image
data) coordinates.  You should therefore use fs/ss coordinates unless the
convenience of x/y coordinates outweighs the speed reduction.

### Providing a separate bad pixel mask

You can provide an array, separate to the image data array, containing
information about the bad pixels.  Up to 8 such masks can be provided for each
detector panel.  Specify the mask location using the following directives,
where you should substitute **N** for a number between 0 and 7 inclusive:

**maskN_data** = _location_
: The location (inside the image data file) of the mask array.  Placeholders
: ('%') in the location will be substituted with the same values as used for the
: placeholders in the image data, although there may be fewer of them for the
: masks than for the image data.

**maskN_file** = _filename_
: Filename to use for the mask data, if not the same as the image data.
: The filename : may be specified as an absolute filename, or relative to the
: working directory.

**maskN_goodbits** = _bitmask_
: Bit mask for good pixels (see below).

**maskN_badbits** = _bitmask_
: Bit mask for bad pixels (see below).

A pixel will be considered *bad* unless all of the bits which are set in
**goodbits** are set.  A pixel will *also* be considered bad if *any* of the
bits which are set in **badbits** are set.  In pseudocode, where **&** is a
bitwise "and", the algorithm is:

    if (mask_value & mask_goodbits) != mask_goodbits:
        mark_pixel_as_bad

    if (mask_value & mask_badbits) != 0:
        mark_pixel_as_bad

Example:

    mask2_data = /data/bad_pixel_map
    mask2_file = /home/myself/mybadpixels.h5
    mask2_goodbits = 0x00
    mask2_badbits = 0xff

There are some older mask directives which are still understood by this version
of CrystFEL.  They are synonyms of the new directives as follows:

    mask       ----->   mask0_data
    mask_file  ----->   mask0_file
    mask_good  ----->   mask0_goodbits
    mask_bad   ----->   mask0_badbits


DETECTOR HIERARCHY
==================

Detector panels can be combined into **groups**.  Certain operations,
especially detector geometry refinement, use these groups to conveniently move
panels.  Groups are specified as follows:

    group_abc = panel1,panel2
    group_def = panel3,panel4

This creates a group called **abc**, containing panels **panel1** and
**panel2**, and a group **def** containing **panel3** and **panel4**.

Groups can themselves be combined into higher-level groups, for example:

    group_all = abc,def

This defines a group called *all** which contains both of the groups created
above.

The highest-level group should always be called **all**.

The group definitions must come **after** the panel definitions.  A good way
to go is to put all the group definitions at the very end of the geometry file.

If the detector consists of only one panel, CrystFEL will automatically create
the **all** group containing it.

The **group** system replaces the **rigid_group** system used in older versions
of CrystFEL.  If the geometry file contains any **rigid_group** lines, they
will be ignored in this version.


EXAMPLES
========

For examples, look in the **examples** folder, which can be found online at
https://gitlab.desy.de/thomas.white/crystfel/-/tree/master/doc/examples


AUTHOR
======

This page was written by Thomas White and Valerio Mariani.


REPORTING BUGS
==============

Report bugs to <taw@physics.org>, or visit <http://www.desy.de/~twhite/crystfel>.


COPYRIGHT AND DISCLAIMER
========================

Copyright © 2023 Deutsches Elektronen-Synchrotron DESY, a research centre of
the Helmholtz Association.

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

**crystfel**(7), **indexamajig**(1), **adjust_detector**(1),
**align_detector**(1)
