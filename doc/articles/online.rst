==========================================
Real-time data processing with indexamajig
==========================================

Instead of reading from files, indexamajig can receive its data from an
"online" streaming system.  In this version, two streaming systems are
implemented in CrystFEL: ZeroMQ and ASAP::O.

When using streamed data, most things work as they normally would.  However,
there are some limitations and "gotchas".  Read this document closely for
details.

Streamed data can be in HDF5, MsgPack or Seedee format.  CBF (as well as
gzipped CBF) may be added in the future.  To specify the format of the data,
use ``--data-format=hdf5``, ``--data-format=seedee`` or
``--data-format=msgpack``.

You can use all of the usual peak search methods for streaming data, but
note that HDF5/CXI peak lists are not currently supported for streamed data.
This may also be added in future.

The option ``--no-image-data`` will be honoured, if given.  This makes it
possible to quickly check streaming data for "indexability", which is more
likely to be useful with streamed data than with files.  You will be able
to do almost all of the usual downstream analysis operations on the resulting
stream, except that attempting to merge it using partialator or process_hkl
will result in zeroes everywhere.


Streaming data over ZeroMQ
==========================

To tell indexamajig to receive data over a ZeroMQ socket, use ``--zmq-input``
instead of ``--input`` or ``-i``.  An error will be generated if you use
``--zmq-input`` and ``--input``  or ``-i`` simultaneously.

Indexamajig can use either a SUB (subscriber) or a REQ (request) socket.  The
SUB socket type can be used for receiving data from OnDA/OM via the same
mechanism that the OnDA/OM GUI uses.  In this case, you will also need to
specify which message prefixes to subscribe to using ``--zmq-subscribe``::

  indexamajig --zmq-input=tcp://127.0.0.1:5002 \
              --zmq-subscribe=ondaframedata \
              -o output.stream -g Eiger.geom ...

You can use ``--zmq-subscribe`` multiple times to subscribe to multiple message
prefixes.

Note that this mode of operation does not combine well with multi-threading
in indexamajig - all worker processes will receive the same data!  For anything
more than "taking a peek" at the data, use the REQ socket instead by using
``--zmq-request`` instead of ``--zmq-subscribe``.  The argument to this option
is the string which should be sent in the request message::

  indexamajig --zmq-input=tcp://127.0.0.1:5002 \
              --zmq-request=next \
              -o output.stream -g Eiger.geom ...

Because they represent completely different modes of operation, the two options
``--zmq-request`` and ``--zmq-subscribe`` are mutually exclusive.


Streaming data using ASAP::O
============================

To tell indexamajig to receive data via ASAP::O, use ``--asapo-endpoint``.
You must additionally specify ``--asapo-token`` (giving the ASAP::O
authentication token), ``--asapo-beamtime`` and ``--asapo-source`` (specifying
the ASAP::O beamtime ID and data source, respectively), and ``--asapo-group``
with the ASAP::O consumer group ID.  If you run multiple copies of
``indexamajig`` on the same data stream, you should make sure that they all use
the same consumer group ID.

Tip: Since the ASAP::O token is a long text string, put it in a separate file
and use ``cat`` in backticks, as follows::

   indexamajig \
       --asapo-endpoint=my-endpoint.facility.de:8400 \
       --asapo-token=`cat /path/to/asapo-token.txt` \
       --asapo-beamtime=mybeamtime1234 \
       --asapo-source=eiger \
       --asapo-group=online

ASAP::O will remember the last retrieved frame in the stream for the given
consumer group ID.  Therefore, if you run ``indexamajig`` again, it will start
from where it left off.  To start again from the beginning of a stream, you
will need to reset the 'last read' location separately.  Instructions and tools
for this are pending...


MsgPack data format
===================

For data in MessagePack format, the following assumptions are made:

* The data consists of either a single MsgPack 'map' object, or an array of
  maps.
  If there are multiple map objects in the array, only the first one will be
  used.  The others will be ignored.
* The image data is given as a two-dimensional array (i.e. no 3D arrays with
  'panel number' as one of the dimensions).
* The image data itself is given as a MsgPack 'map' object representing a
  serialised NumPy array.  That is, it should contain ``type``, ``data`` and
  ``shape`` keys.
* The data ``type`` field should contain either ``<i4`` (if the data is in
  little-endian 32-bit signed integer format) or ``<f4`` for 32-bit (IEEE754
  single precision) floating-point.
* The data ``shape`` field should be a 1D array of two values.  The first
  element is the slow-scan size, the second is that fast-scan size.
* The data array within the NumPy map should be in a binary object called
  ``data``.

Note that *all* of these assumptions are 'open for negotiation' and will be
relaxed in future CrystFEL versions, as new online data formats arise.

You can specify which map objects to look at in the geometry file.  The
following example will get the incident photon energy (in eV) and detector
distance (in mm) by looking up the ``beam_energy`` and ``detector_distance``
keys, respectively, in the MsgPack map object.  It will then look up
``detector_data`` to find the image data itself.  See the next section for an
explanation of ``peak_list``::

  photon_energy = beam_energy eV
  adu_per_photon = 1
  clen = detector_distance mm
  res = 5814.0
  peak_list = peak_list
  
  thepanel/data = detector_data
  thepanel/min_fs = 0
  thepanel/max_fs = 2067
  thepanel/min_ss = 0
  thepanel/max_ss = 2161
  thepanel/corner_x = -1034
  thepanel/corner_y = -1081
  thepanel/fs = x
  thepanel/ss = y


You can use ``--peaks=msgpack`` to get the peak locations from the MsgPack
data.  In this case, the ``peak_list`` directive in the geometry file specifies
the key for the peak information in the MsgPack map object. The peak
information itself is expected to be a map object with three keys: ``fs``,
``ss`` and ``intensity``.  Each of these keys should correspond to an array
containing (respectively) the fast scan and slow scan coordinates of each peak,
and their intensities.  Obviously, the three arrays must have equal sizes.

Note that there is no way, in this structure, to communicate which detector
panel contains a peak, in the case where different detector panels cover the
same pixel ranges (in this case, the pixel data would from multiple data
blocks).  In practice, this means that the ``data`` directives for all panels
need to be the same when using ``--peaks=msgpack``.

Note also that the options ``--no-revalidate`` and ``--check-hdf5-snr`` apply
to the peak lists from ``--peaks=msgpack``.
