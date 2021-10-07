Installation instructions
=========================

CrystFEL installation is easiest on GNU/Linux.  Installation on Mac OS X is
supported, but more difficult because you have to get all the dependencies
from a third-party repository.  Installation in Windows via
[Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/))
is also possible.


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

* Either [CMake](https://cmake.org/) 3.12 or later or [Meson](https://mesonbuild.com/) 0.55.0 or later (Meson is preferred)
* [HDF5](https://www.hdfgroup.org/downloads/hdf5/) 1.8.0 or later (1.10.0 or later is required for many recent data formats which use the 'virtual data set' feature)
* [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
* [Bison](https://www.gnu.org/software/bison/) 2.6 or later
* [Flex](https://www.gnu.org/software/flex/)
* [Zlib](https://www.zlib.net/) (1.2.3.5 or later preferred for better decompression speed)

The following dependencies are "optional", in the sense that you can install
CrystFEL without them.  However, a CrystFEL installation without these will lack
important features such as the graphical user interface.  The following list is
roughly in order of importance:

* [GTK](https://gtk.org/) version between 3.12 and 3.24 (inclusive) - required for GUI
* [Cairo](https://www.cairographics.org/) 1.2 or later (required for GUI)
* [Pango](https://pango.gnome.org/) 1.0 or later, including [PangoCairo](https://docs.gtk.org/PangoCairo/) (required for GUI)
* [gdk-pixbuf](https://docs.gtk.org/gdk-pixbuf/) 2.0 or later (required for GUI)
* [libccp4](ftp://ftp.ccp4.ac.uk/opensource/) \[\*\] (required for MTZ import/export)
* [XGandalf](https://stash.desy.de/users/gevorkov/repos/xgandalf) \[\*\] (for `xgandalf` indexing)
* [PinkIndexer](https://stash.desy.de/users/gevorkov/repos/pinkindexer) \[\*\] (for indexing electron or wide bandwidth diffraction patterns)
* [FFTW](http://fftw.org/) 3.0 or later (required for `asdf` indexing)
* [FDIP](https://stash.desy.de/users/gevorkov/repos/fastdiffractionimageprocessing/) \[\*\] (for `peakFinder9` peak search algorithm)
* NCurses (for integration diagnostics: `indexamajig --int-diag`)

Most of the dependencies mentioned above should be available from your Linux
distribution's package manager, or from [Homebrew](https://brew.sh/) on Mac OS.
We emphatically recommend against trying to install GTK, Cairo, Pango or
gdk-pixbuf from source.

If you compile using Meson, dependencies marked with \[\*\] above will be
downloaded and compiled automatically if they are not available on the system.
If you don't want this, add option `--wrap-mode=nofallback` when invoking
Meson.  See the Meson manual for other possibilities, such as using
locally-provided files instead of downloading them.

If you install libraries in a non-system location, set `PKG_CONFIG_PATH` so
that they can be found.  For example:
```
$ export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/home/user/xgandalf/lib64/pkgconfig
```
Tip: build using Meson wherever possible, because it will take care of several
of the dependencies automatically and you will never have to worry about this.

We also do not recommend using dependencies from Conda/Anaconda.  Do not
activate any Conda environment before compiling CrystFEL, not even the "base"
environment.  Don't even "source" the Conda setup file before installing
CrystFEL - keep it completely separate.  A Conda recipe for CrystFEL might be
coming soon, though, if development resources allow for it.

If OpenCL headers and corresponding GPU drivers are available on your system,
GPU-accelerated diffraction calculation will be enabled for `pattern_sim`.

Processing data relies on indexing 'engines'.  By default, you will have access
to the [TakeTwo](https://journals.iucr.org/d/issues/2016/08/00/rr5128/)
indexing algorithm.  Some of the optional dependencies listed above, if found,
will add more indexing algorithms.  In addition, the following programs will be
used for indexing if they are available:

* [Mosflm](https://www.mrc-lmb.cam.ac.uk/mosflm/mosflm/)
* [DirAx](http://www.crystal.chem.uu.nl/distr/dirax/)
* [XDS](http://xds.mpimf-heidelberg.mpg.de/)

You can install these at any time before or after installing CrystFEL.  Note
that you only need the Mosflm binary, not the full `iMosflm` user interface.
[Download it here](https://www.mrc-lmb.cam.ac.uk/mosflm/mosflm/ver740/pre-built/mosflm-linux-64-noX11.zip).

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
$ sudo dnf install hdf5-devel gsl-devel gtk3-devel cairo-devel pango-devel gdk-pixbuf2-devel meson cmake gcc-c++
$ sudo dnf install fftw-devel ncurses-devel zeromq-devel msgpack-devel libccp4-devel
```

Here are the commands for Ubuntu and Debian:
```
$ sudo apt install build-essential
$ sudo apt install libhdf5-dev libgsl-dev gtk3-devel cairo-devel pango-devel gdk-pixbuf2-devel meson cmake gcc-c++
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

Run `indexamajig` for a basic check that the installation has succeeded.  A
healthy newborn CrystFEL should complain that `You need to provide the input
filename (use -i)` when run with no other command-line options.

Alternatively, run `crystfel` to start the graphical user interface, provided
that the dependencies for the GUI were met (see above).

Refer to the tutorial to see where to go from here!


Notes about strange problems
----------------------------
The following problems are usually only encountered when installing dependencies
manually, in non-standard locations and with multiple conflicting parallel
installations of certain dependencies.

* **Problem**: Linker error about HDF5 and `fPIC`:
    ```
    /usr/bin/ld: <install path>/lib/libhdf5.a(H5.o): relocation R_X86_64_32 against `.rodata' can not be used when making a shared object; recompile with -fPIC
    ```

  **Explanation**: The dependency libraries, including HDF5, GSL and FFTW,
  must be built such that they can be used from within CrystFEL's shared
  library *libcrystfel*.  In particular, they must be compiled into relocatable
  code (position-independent code, hence 'PIC').  This is usually (but not
  always!) the case when the libraries are built as shared objects ("`.so`"),
  but not for static-linking libraries ("`.a`").  HDF5's `h5cc` tool prefers to
  use static linking by default, and the preference gets picked up and used by
  Meson for CrystFEL.  However, because *libcrystfel* is itself a shared
  library, this will only work if the HDF5 static-linking library was compiled
  with `-fPIC`.

  **Solution**: When building HDF5, disable the static linking libraries
  altogether, to force the use of shared libraries: `./configure
  --enable-shared --disable-static`.  Alternatively, add the `-fPIC` option so
  that the static libraries can be used within *libcrystfel*: `./configure
  H5_CFLAGS=-fPIC`.  Compiling HDF5 using CMake also avoids the problem.

* **Problem**: Linker error about FFTW and `fPIC`:
    ```
    /usr/bin/ld: /usr/lib/libfftw3.a(mapflags.o): relocation R_X86_64_32 against `.rodata' can not be used when making a shared object; recompile with -fPIC
    /usr/lib/libfftw3.a: could not read symbols: Bad value
    ```

  **Explanation**: As above.

  **Solution**: When building FFTW: `./configure --enabled-shared`.

* **Problem**: Compiler error about `getpgrp`:
    ```
    ../libcrystfel/src/utils.c: In function ‘progress_bar’:
    ../libcrystfel/src/utils.c:383:2: error: too few arguments to function ‘getpgrp’
      if ( tcgetpgrp(STDERR_FILENO) != getpgrp() ) return;
      ^
    ````

  **Explanation**: This one is is quite complicated, and ultimately due to a
  misfeature of the HDF5 build process.  The compiler flags for HDF5, set via
  the `CFLAGS` environment variable while compiling HDF5, are tracked through
  into the output of the `h5cc` compiler driver program.  From here, they get
  picked up by Meson and added to the compiler flags for CrystFEL.

  It's not obvious from the HDF5 documentation, but `H5_CFLAGS` should be used
  instead of `CFLAGS` for controlling how HDF5 itself is compiled.  `H5_CFLAGS`
  are private to HDF5, whereas `CFLAGS` get added to `h5cc`.  Several Linux
  distributions, including Arch and CentOS, apparently set some problematic
  flags in `CFLAGS`.  One such flag is `-D_BSD_SOURCE`, which changes the
  prototypes of certain system calls including `getpgrp`, breaking CrystFEL.

  **Solution**: Three options:

  Option 1: Do not use your distribution's HDF5 package.  Instead, compile
  HDF5 yourself.  If you compile HDF5 with autotools (`./configure` et al.),
  take care to use `H5_CFLAGS` instead of `CFLAGS` to add any compiler flags.
  Ensure that the `h5cc` corresponding to the correct version of HDF5 is found,
  by making sure that its location comes first in `PATH` when compiling
  CrystFEL.

  Option 2: Compile HDF5 using CMake instead of autotools.  In this case,
  HDF5 will install a `pkg-config` file which Meson will use in preference to
  `h5cc`.

  Option 3: Compile CrystFEL using CMake instead of Meson.  CMake extracts
  the include paths from the `h5cc` output and ignores the others.  However,
  note that CrystFEL's CMake build process is deprecated and will eventually be
  removed.


* **Problem**: Linker error about missing HDF5 symbols:
    ```
    libcrystfel/libcrystfel.so.0.9.1: undefined reference to `H5P_CLS_LINK_CREATE_ID_g'
    libcrystfel/libcrystfel.so.0.9.1: undefined reference to `H5Oget_info_by_idx1'
    ```

  **Explanation**: CrystFEL was accidentally compiled with headers from a new
  HDF5 version (1.10 or higher) but linked with an older one (1.8 or lower).
  This can happen if there are headers in the include path corresponding to an
  HDF5 installation different to the one found during the configuration step of
  CrystFEL's build process.

  **Solution**: Ensure that HDF5 is found consistently:  Make sure that
  the `pkg-config` file (`libhdf5.pc` or `libhdf5-<version>.pc`) exists and is
  found correctly (set `PKG_CONFIG_PATH`).  Make sure that `PATH` is set so
  that the location of the correct version of `h5cc` comes first.

* **Problem**: After fixing the above problem, there is still a complaint about
  HDF5: `H5Oget_info_by_idx` (note: no `1` at the end).

  **Explanation**: This change to the HDF5 headers seems to confuse `ccache`.

  **Solution**: The (somewhat drastic) fix is simply to `rm -rf ~/.ccache`.

* **Problem**: After all of the above, HDF5 is *still* not found correctly.

  **Explanation**: The `pkg-config` files for HDF5, if they exist, might
  contain the version number.  Example: `libhdf5-1.12.0.pc`.  This makes them
  difficult to find using `pkg-config`, because you have to know the version
  number in advance!  For example, `pkg-config --cflags libhdf5-1.12.0`.

  **Solution**: Upgrade to Meson 0.53.0 or later, which is aware of this quirk
  and looks for any HDF5 `pkg-config` file.

