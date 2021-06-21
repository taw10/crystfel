Installation instructions
=========================

CrystFEL installation is supported on GNU/Linux (including on Windows via
[Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/))
and Mac OS X.  A native installation on Windows is reported as possible, but is
not supported.


Supported installations at facilities
-------------------------------------

Before starting, check the [list of facility installations](https://www.desy.de/~twhite/crystfel/facilities.html)
to find out if there's already a CrystFEL installation available at your site.

If you want to use CrystFEL on a facility computer system, we recommend leaving
the installation to facility IT staff.  This means that a single installation
can be shared and maintained.  In addition, installing software from source, on
a system where you don't have administrative access, can be quite difficult.
Please feel free to tell the IT staff to get in touch with us for assistance,
and also to make sure that the installation is documented on the list of
supported facility installations.


Dependencies
------------

Here are the mandatory dependencies - you cannot install CrystFEL without these:

* Either [CMake](https://cmake.org/) 3.12 or later or [Meson](https://mesonbuild.com/) 0.50.0 or later (Meson is preferred)
* [HDF5](https://www.hdfgroup.org/downloads/hdf5/) 1.8.0 or later (1.10.0 or later is required for many recent data formats which use the 'virtual data set' feature)
* [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
* [Bison](https://www.gnu.org/software/bison/) 2.6 or later
* [Flex](https://www.gnu.org/software/flex/)
* [Zlib](https://www.zlib.net/) (1.2.3.5 or later preferred for better decompression speed)

The following dependencies are "optional", in the sense that you can install
CrystFEL without them.  However, a CrystFEL installation without these will lack
important features such as the graphical user interface.  The following list is
roughly in order of importance:

* GTK3 (required for GUI)
* Cairo (required for GUI)
* Pango (required for GUI)
* gdk-pixbuf (required for GUI)
* [libccp4](ftp://ftp.ccp4.ac.uk/opensource/) (required for MTZ import/export)
* [FFTW3](http://fftw.org/) (required for `asdf` indexing)
* [FDIP](https://stash.desy.de/users/gevorkov/repos/fastdiffractionimageprocessing/) (for `peakFinder9` peak search algorithm)
* NCurses (for integration diagnostics: `indexamajig --int-diag`)

Apart from FDIP, all of the dependencies mentioned above (including libccp4 -
yes, really!) should be available from your Linux distribution's package
manager, or from [Homebrew](https://brew.sh/) on Mac OS.  You should not need
to download and install any of them separately from source, and we emphatically
recommend against trying to do so.

Processing data relies on indexing 'engines'.  By default, you will have access
to the [TakeTwo](https://journals.iucr.org/d/issues/2016/08/00/rr5128/) indexing
algorithm.  Provided that FFTW is available, you will also have access to the
`asdf` algorithm.  The more of the following are additionally installed, the
better your experience will be:

* [XGandalf](https://stash.desy.de/users/gevorkov/repos/xgandalf)
* [PinkIndexer](https://stash.desy.de/users/gevorkov/repos/pinkindexer)
* [Mosflm](https://www.mrc-lmb.cam.ac.uk/mosflm/mosflm/)
* [DirAx](http://www.crystal.chem.uu.nl/distr/dirax/)
* [XDS](http://xds.mpimf-heidelberg.mpg.de/)

The indexing engines have their own installation instructions.  XGandalf and
PinkIndexer need to be installed before compiling CrystFEL, but the others can
be installed at any time afterwards.

Note that you only need the Mosflm binary, not the full `iMosflm` user interface.
[Download it here](https://www.mrc-lmb.cam.ac.uk/mosflm/mosflm/ver730/pre-built/mosflm-linux-64-noX11.zip).

You probably have CCP4 installed already, with its own copy of `mosflm`.
However, you should keep your CCP4 installation separate from CrystFEL.  The
reason is that CCP4 comes with its own private copies of all its dependencies,
and it becomes very difficult to make the build process for CrystFEL (or
anything else) find the right versions of everything.  Do not 'source' the CCP4
setup file before trying to install CrystFEL, and make sure that the setup file
is not automatically referenced in your shell setup files (`.bashrc` and
others).

Here are the commands to install all the basic dependencies (including the
optional dependencies, but not the indexing engines) on CentOS and Fedora 22 or
later (for CentOS, replace `dnf` with `yum`):
```
$ sudo dnf group install 'Development Tools'
$ sudo dnf install hdf5-devel gsl-devel gtk3-devel cairo-devel pango-devel gdk-pixbuf2-devel meson cmake
$ sudo dnf install fftw-devel ncurses-devel zeromq-devel msgpack-devel libccp4-devel
```

Here are the commands for Ubuntu and Debian:
```
$ sudo apt install build-essential
$ sudo apt install libhdf5-dev libgsl-dev gtk3-devel cairo-devel pango-devel gdk-pixbuf2-devel meson cmake
$ sudo apt install fftw-devel ncurses-devel libmsgpack-dev libzmq3-dev libccp4-dev
```

For Mac OS X, first install [Homebrew](https://brew.sh/), which will also cause
[Xcode](https://developer.apple.com/xcode/) to be installed.  Then:
```
$ brew install gsl hdf5 flex bison argp-standalone pkg-config doxygen gtk+3 cairo pango gdk-pixbuf fftw meson cmake
$ export PATH="$(brew --prefix)/opt/bison/bin:$(brew --prefix)/opt/flex/bin:$PATH"
$ export LDFLAGS="-L$(brew --prefix)/opt/bison/lib -L$(brew --prefix)/opt/flex/lib -L$(brew --prefix)/opt/argp-standalone/lib -largp $LDFLAGS"
$ export CFLAGS="-I$(brew --prefix)/opt/flex/include -I$(brew --prefix)/opt/argp-standalone/include/ $CFLAGS"
```
The `export` commands ensure that the libraries installed via Homebrew can be
found by CrystFEL's build system.


Building and installing using Meson (recommended)
-------------------------------------------------

The actual compilation and installation step is very easy.  Here are the
commands:
```
$ meson build
$ ninja -C build
$ sudo ninja -C build install
```


Building and installing using CMake (deprecated)
------------------------------------------------

Compiling and installing using CMake is also simple:
```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ sudo make install
```


Starting up
-----------

Run `indexamajig` for a basic check that the installation has succeeded.  It
should complain that `You need to provide the input filename (use -i)`.
Alternatively, run `crystfel` to start the graphical user interface, provided
that GTK3 was available during installation.

Refer to the tutorial to see where to go from here!


Notes about strange problems
----------------------------
The following problems are usually only encountered when installing dependencies
manually, in non-standard locations and with multiple conflicting parallel
installations of certain dependencies.

* If you get an error about an incorrect number of parameters for `getpgrp`,
  check that HDF5 was installed using CMake, not autotools (`./configure`).

  The compiler flags generated by `h5cc -show` include all kinds of stuff that
  shouldn't be there, such as `-D_BSD_SOURCE` which changes some function
  prototypes (e.g. `getpgrp`).  Meson uses this method as a fallback if
  `pkg-config` fails to locate HDF5, breaking the CrystFEL build because
  CrystFEL calls `getpgrp` and relies on a certain function signature.
  Therefore, you need to make sure that HDF5 is discovered using `pkg-config`,
  not `h5cc`.  However, HDF5 only installs its `pkg-config` file when built
  using CMake, not autotools.  This problem does not occur when building
  CrystFEL with CMake, because CMake extracts the include paths from the `h5cc`
  output.

* The following error means that CrystFEL was accidentally compiled with headers
  from a new HDF5 version (1.10 or higher) but linked with an older one (1.8 or
  lower):
    ```
    libcrystfel/libcrystfel.so.0.9.1: undefined reference to `H5P_CLS_LINK_CREATE_ID_g'
    libcrystfel/libcrystfel.so.0.9.1: undefined reference to `H5Oget_info_by_idx1'
    ```
  This can happen if there are headers in the include path corresponding to an
  HDF5 installation different to the one found during the configuration step of
  CrystFEL's build process.

  The solution is to make sure that HDF5 is found consistently: make sure that
  the `pkg-config` file (`libhdf5.pc`) exists and is found correctly (set
  `PKG_CONFIG_PATH`).  If not using `pkg-config`, make sure that `PATH` is set
  so that the location of the correct version of `h5cc` comes first.

* If, after fixing the above problem, you still get an error about
  `H5Oget_info_by_idx` (note: no `1` at the end), there could be some residual
  `ccache` files causing problems.  The (somewhat drastic) fix is simply to
  `rm -rf ~/.ccache`.

* If the build process complains about FFTW and `-fPIC`, make sure that your
  installation of FFTW was configured with `./configure --enable-shared`:
    ```
    /usr/bin/ld: /usr/lib/libfftw3.a(mapflags.o): relocation R_X86_64_32 against `.rodata' can not be used when making a shared object; recompile with -fPIC
    /usr/lib/libfftw3.a: could not read symbols: Bad value
    ```
  This is because CrystFEL uses the FFTW shared library via `libcrystfel`, which
  is itself a shared library.
