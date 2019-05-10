\mainpage libcrystfel index page

Version
=======
This documentation is for libcrystfel from CrystFEL ${CRYSTFEL_VERSION}, API
revision ${CRYSTFEL_API_VERSION}.

Abstract
========
This is the internal documentation for CrystFEL.  Unless you are looking at
the code, writing new programs or fixing bugs, you should not need to read
this.  You might use the information here when reading the code or to better
understand how the software works, or refer to it when creating a new
program within the suite.

Coding standards
================
Please see the \ref coding "section on coding standards" for CrystFEL's coding
style rules (including libcrystfel and the core CrystFEL programs).

API documentation
=================

* \ref image.h "The image structure and image data file handling"
* Handling reflection data:
   * \ref reflist.h "Reflection list structure"
   * \ref reflist-utils.h "Reflection list utility functions"
* Unit cells:
   * \ref cell.h "Unit cell structure"
   * \ref cell-utils.h "Unit cell utility functions"
* \ref crystal.h "Crystal structure"
* \ref events.h "Frame descriptors for multi-event files"
* \ref geometry.h "Geometry of diffraction (prediction/partiality calculations)"
* Peak search
   * \ref peaks.h "Main peak search functions"
   * \ref peakfinder8.h "The peakfinder8 algorithm"
* \ref filters.h "Image (noise) filters"
* \ref symmetry.h "Point group symmetry"
* Mathematical constructions:
   * \ref integer_matrix.h "Integer matrices"
   * \ref rational.h "Rational numbers (including rational matrices)"
* \ref index.h "Top-level indexing system"
   * \ref xgandalf.h "XGANDALF indexer interface"
   * \ref xds.h "XDS indexer inderface"
   * \ref mosflm.h "MOSFLM indexer interface"
   * \ref dirax.h "DirAx indexer interface"
   * \ref taketwo.h "TakeTwo indexing algorithm"
   * \ref asdf.h "ASDF indexing algorithm"
   * \ref felix.h "Felix indexer interface"
* \ref predict-refine.h "Prediction refinement"
* \ref integration.h "Integration of reflections"
* \ref detector.h "Detector geometry descriptions"
* \ref spectrum.h "Radiation spectrum object"
* \ref hdf5-file.h "HDF5 file interface"
* \ref stream.h "Stream format for indexing/integration results"
* \ref render.h "Miscellaneous rendering functions (colour scale)"
* \ref thread-pool.h "Thread pool"
* \ref utils.h "Miscellaneous utility functions"
