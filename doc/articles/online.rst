==========================================
Real-time data processing with indexamajig
==========================================

Instead of reading from files, indexamajig can get its data over a ZeroMQ
socket.  Use ``--zmq-input`` instead of ``--input`` or ``-i``.  An error will
be generated if you use ``--zmq-input`` and ``--input``  or ``-i``
simultaneously.

Indexamajig assumes that the socket is a *pub/sub* socket.  You will also need
to specify which message prefixes to subscribe to.  Use ``--zmq-subscribe`` for
this::

  indexamajig --zmq-input=tcp://127.0.0.1:5002 \
              --zmq-subscribe=ondaframedata \
              -o output.stream -g Eiger.geom ...

You can use ``--zmq-subscribe`` multiple times to subscribe to multiple message
prefixes.

The option ``--no-image-data`` will be honoured, if given.  This makes it
possible to quickly check streaming data for "indexability".


Data format
===========

In this version, CrystFEL makes the following assumptions about the data
received via ZeroMQ:

* The data is serialised using MessagePack
* The data consists of either a single MsgPack 'map' object, or an array of
  maps.
  If there are multiple map objects in the array, only the first one will be
  used.  The others will be ignored.
* The image data is given as a two-dimensional array (i.e. no 3D arrays with
  'panel number' as one of the dimensions).
* The image data itself is given as a MsgPack 'map' object with "type"=>"<i4".

Note that *all* of the above assumptions are 'open for negotiation' and will be
relaxed in future CrystFEL versions, as new online data formats arise.  Anyone
interested in streaming CBF files over ZeroMQ?

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


Peak lists
==========

You can use all of the usual peak search methods for streaming data.

In addition, you can use ``--peaks=msgpack`` to get the peak locations from
the MsgPack data.  In this case, the ``peak_list`` directive in the geometry
file specifies the key for the peak information in the MsgPack map object.
The peak information itself is expected to be a map object with three keys:
``fs``, ``ss`` and ``intensity``.  Each of these keys should correspond to an
array containing (respectively) the fast scan and slow scan coordinates of each
peak, and their intensities.  Obviously, the three arrays must have equal sizes.

Note that there is no way, in this structure, to communicate which detector
panel contains a peak, in the case where different detector panels cover the
same pixel ranges (in this case, the pixel data would from multiple data
blocks).  In practice, this means that the ``data`` directives for all panels
need to be the same when using ``--peaks=msgpack``.
