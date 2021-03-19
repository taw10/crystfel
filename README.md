CrystFEL - Data processing for serial crystallography
=====================================================

Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
                      a research centre of the Helmholtz Association.

See AUTHORS as well as individual source code files for full details of contributors.

Introduction
------------

CrystFEL is a suite of programs for processing (and simulating) Bragg
diffraction data from "serial crystallography" experiments, often (but not
always) performed using an X-ray Free-Electron Laser.  Compared to rotation data,
some of the particular characteristics of such data which call for a
specialised software suite are:

* The sliced, rather than integrated, measurement of intensity data.  Many, if
  not all reflections are partially integrated.

* Many patterns (thousands) are required - high throughput is needed.

* The crystal orientations in each pattern are random and uncorrelated.

* Merging into lower symmetry point groups may require the resolution of
  indexing ambiguities.


Getting started
---------------

The best way to get started, after installation, is to run command ```crystfel```
to start the graphical user interface.


Installation
------------

CrystFEL installation is supported on GNU/Linux and Mac OS X.  The terse
installation instructions below should be enough if you're experienced with
installing software from source.  More detailed installation information is
available [on the website](https://www.desy.de/~twhite/crystfel/install.html).

Here are the mandatory dependencies - you cannot install CrystFEL without these:

* Either [CMake](https://cmake.org/) 3.12 or later or [Meson](https://mesonbuild.com/) (Meson is preferred)
* [HDF5](https://www.hdfgroup.org/downloads/hdf5/) 1.8.0 or later (1.10.0 or later is required for many recent data formats)
* [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
* [Bison](https://www.gnu.org/software/bison/) 2.6 or later
* [Flex](https://www.gnu.org/software/flex/)
* [Zlib](https://www.zlib.net/) (1.2.3.5 or later preferred for better decompression speed)

The following dependencies are "optional", in the sense that you can install
CrystFEL without them.  However, a CrystFEL installation without these will lack
important features such as the graphical user interface:

* GTK3 or later
* Cairo
* Pango
* gdk-pixbuf
* NCurses
* libPNG
* [libccp4](ftp://ftp.ccp4.ac.uk/opensource/)

Note that all of the dependencies mentioned above (including libccp4) should be
available from your Linux distribution's package manager, or from
[Homebrew](https://brew.sh/) on Mac OS.  You should not need to download and
install any of them separately from source, and we emphatically recommend
against trying to do so!

Note that using the libraries from the full CCP4 suite is not recommended.  CCP4
includes so many other libraries that it becomes very difficult to link using
the correct versions of everything.

Processing data relies on indexing algorithms.  The more of the following are
installed, the better your experience will be:

* [FFTW3](http://fftw.org/)
* [XGandalf](https://stash.desy.de/users/gevorkov/repos/xgandalf)
* [PinkIndexer](https://stash.desy.de/users/gevorkov/repos/pinkindexer)
* [Mosflm](https://www.mrc-lmb.cam.ac.uk/mosflm/mosflm/)
* [DirAx](http://www.crystal.chem.uu.nl/distr/dirax/)
* [XDS](http://xds.mpimf-heidelberg.mpg.de/)

Installation follows the normal CMake procedure:

```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ sudo make install
```

Or, with Meson:

```
$ meson . build
$ ninja -C build
$ sudo ninja -C build install
```


Licence
-------

CrystFEL is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CrystFEL is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CrystFEL.  If not, see <http://www.gnu.org/licenses/>.


Funding acknowledgements
------------------------

Development of CrystFEL is primarily funded by the Helmholtz Association.

Partial funding for CrystFEL has previously been provided by:

- "X-Probe", a project of the European Union's 2020 Research and Innovation
  Program Under the Marie Skłodowska-Curie grant agreement 637295 (2015-2018).

- The BMBF German-Russian Cooperation "SyncFELMed", grant 05K14CHA (2014-2017).

- BioStruct-X, a project funded by the Seventh Framework Programme (FP7) of the
  European Commission (2011-2016).

