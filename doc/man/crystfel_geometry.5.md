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

data
: The location in the HDF5 file of the data block that contains the panel's data.
: The default value is /data/data.  If the HDF5 file contains multiple events,
: and each event is stored in a different data block, the variable part of the
: path can be represented using the % character placeholder.
:
: Example:
:
:     data = /data/%/rawdata
:
: The CrystFEL programs will look for the first event at
: /data/event1_name/rawdata, for the second at /data/event2_name/rawdata, etc.,
: where event_name and event2_name are simply whatever the program could find in
: the HDF5 file which matched the pattern you gave.

dimn
: Information about the layout of the data block identified by the 'data'
: property. n is an integer number identifying an axis in a multidimensional HDF5
: data block. The property value defines the kind of information encoded by the
: axis. Possible values are:
: % - event placeholder,the axis encodes events
: ss - the axis encoding the slow scan index
: fs - the axis encodes the fast scan index
: number -  the index in this dimension should be fixed at number.
:
: CrystFEL assumes that the data block defined by the 'data' property has a
: dimensionality corresponding to the axis with the highest value of n defined by
: the 'dim' property.  That is, if the geometry file specifies dim0, dim1 and
: dim2, then the data block is expected to be three-dimensional.  The size of the
: data block along each of those axes comes from the image metadata (e.g. the
: array sizes in the HDF5 file).
:
: The lowest number of n corresponds to the most slowly-changing array index as
: the data block is traversed.  The default values are dim0=ss and dim1=fs.  The
: value of n corresponding to fs must not be lower than the value assigned to ss,
: i.e. "fast scan is always fast scan".
:
: Example:
:
:     dim0 = %
:     dim1 = 4
:     dim2 = ss
:     dim3 = fs
:
: The above snippet specifies that the data block is 4-dimensional. The first
: axis represents the event number, the index in the second axis is always 4, and
: the remaining two axes are the image coordinates.

min_fs, min_ss, max_fs, max_ss
: The range of pixels in the data block specified by the 'data' property that
: corresponds to the panel, in fast scan/slow scan coordinates, specified
: inclusively.


PEAK LISTS
==========

peak_list = loc
: This gives the location of the peak list in the data files, for peak detection
: methods hdf5 and cxi (see man indexamajig).

peak_list_type = layout
: Specify the layout of the peak list.  Allowed values are cxi, list3 and auto.
:
: list3 expects the peak list to be a two dimensional array whose size in the
: first dimension equals the number of peaks and whose size in the second
: dimension is exactly three.  The first two columns contain the fast scan and
: slow scan coordinates, the third contains the intensities.  This is the correct
: option for "single-frame" HDF5 files as written by older versions of Cheetah.
:
: cxi expects the peak list to be a group containing four separate HDF5 datasets:
: nPeaks, peakXPosRaw, peakYPosRaw and peakTotalIntensity.  See the specification
: for the CXI file format at http://www.cxidb.org/ for more details.  This is the
: correct option for "multi-event" HDF5 files as output by recent versions of
: Cheetah.
:
: auto tells CrystFEL to decide between the above options based on the file extension.
:
: Note that CrystFEL considers all peak locations to be distances from the corner
: of the detector panel, in pixel units, consistent with its description of
: detector geometry (see 'man crystfel_geometry').  The software which generates
: the HDF5 or CXI files, including Cheetah, may instead consider the peak
: locations to be pixel indices in the data array.  To compensate for this
: discrepancy, CrystFEL will, by default, add 0.5 to all peak coordinates. Use
: --no-half-pixel-shift if this isn't what you want.


DETECTOR RESPONSE PROPERTIES
============================

adu_per_eV, adu_per_photon
: The number of detector intensity units (ADU) which will arise from either one
: electron-Volt of photon energy, or one photon.  This is used to estimate
: Poisson errors.  Note that setting different values for this parameter for
: different panels does not result in the intensities being scaled accordingly
: when integrating data.  You should only specify one out of adu_per_eV and
: adu_per_photon.

res
: The resolution (in pixels per metre) for this panel.  This is one over the
: pixel size in metres.

max_adu
: The saturation value for the panel.  You can use this to exclude saturated
: peaks from the peak search or to avoid integrating saturated reflections.
: However, usually it's best to include saturated peaks, and exclude saturated
: reflections with the --max-adu option of process_hkl and partialator.
: Therefore you should avoid setting this parameter - a warning will be displayed
: if you do.

saturation_map
: This specifies the location of the per-pixel saturation map in the HDF5 file.
: This works just like mask in that it can come from the current file or a
: separate one (see saturation_map_file).  Reflections will be rejected if they
: contain any pixel above the per-pixel values, in addition to the other checks
: (see max_adu).

saturation_map_file
: Specifies that the saturation map should come from the HDF5 file named here,
: instead of the HDF5 file being processed.  It can be an absolute filename or
: relative to the working directory.


BAD REGIONS
===========

Bad regions will be completely ignored by CrystFEL.  You can specify the pixels
to exclude in pixel units, either in the lab coordinate system (see above) or
in fast scan/slow scan coordinates (mixtures are not allowed).   In the latter
case, the range of pixels is specified inclusively.  Bad regions are
distinguished from normal panels by the fact that they begin with the three
letters "bad".

If you specify a bad region in fs/ss (image data) coordinates, you must also
specify which panel name you are referring to.

Note that bad regions specified in x/y (lab frame) coordinates take longer to
process (when loading images) than regions specified in fs/ss (image data)
coordinates.  You should use fs/ss coordinates unless the convenience of x/y
coordinates outweighs the speed reduction.

no_index
: Set this to 1 or "true" to ignore this panel completely.

flag_lessthan, flag_morethan, flag_equal
: Mark pixels as "bad" if their values are respectively less than, more than or
: equal to the given value.  Note carefully that the inequalities are strict, not
: inclusive: "less than", not "less than or equal to".

mask_edge_pixels
: Mark the specified number of pixels, at the edge of the panel, as "bad".

maskN_data, maskN_file, maskN_goodbits, maskN_badbits
: These specify the parameters for bad pixel mask number N.  You can have up to 8
: bad pixel masks, numbered from 0 to 7 inclusive.  Placeholders ('%') in the
: location (maskN_data) will be substituted with the same values as used for the
: placeholders in the image data, although there may be fewer of them for the
: masks than for the image data.
: 
: You can optionally give a filename for each mask with maskN_file.  The filename
: may be specified as an absolute filename, or relative to the working directory.
: If you don't specify a filename, the mask will be read from the same file as
: the image data.
: 
: A pixel will be considered bad unless all of the bits which are set in goodbits
: are set.  A pixel will also be considered bad if any of the bits which are set
: in badbits are set.  Note that pixels can additionally be marked as bad via
: other mechanisms as well (e.g. no_index or bad).

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

Examples:

    badregionA/min_x = -20.0
    badregionA/max_x = +20.0
    badregionA/min_y = -100.0
    badregionA/max_y = +100.0

    badregionB/min_fs = 128
    badregionB/max_fs = 160
    badregionB/min_ss = 256
    badregionB/max_ss = 512
    badregionB/panel = q0a1


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
