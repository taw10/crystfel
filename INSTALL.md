Installation instructions
=========================

CrystFEL installation is supported on GNU/Linux and Mac OS X.  Installation in
Windows via [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/)
is also possible.


Supported installations at facilities
-------------------------------------

Before starting, check the [list of facility installations](https://www.desy.de/~twhite/crystfel/facilities.html)
to find out if there's already a CrystFEL installation available at your site.

If you want to use CrystFEL on a facility computer system, we recommend leaving
the installation to facility IT staff.  This means that a single installation
can be shared and maintained centrally.  In addition, installing software from
source, on a system where you don't have administrative access, can be quite
difficult. Please feel free to tell the IT staff to get in touch with us for
assistance, and also to make sure that the installation is documented on the
list of supported facility installations.


Installation via container registry
-----------------------------------

The easiest way to get started is to download our container image and run it
using a virtualization tool of your choice (e.g. Docker, Podman,
[Singularity/Apptainer](http://apptainer.org/)).  For example, using Apptainer:
```
$ apptainer pull docker://gitlab.desy.de:5555/thomas.white/crystfel/crystfel:latest
$ apptainer run -B /path/to/data crystfel_latest.sif
```

After the second command, you are working inside the container and should be
able to run all the CrystFEL commands. Start by running ```crystfel``` to run
the GUI.

By default, only a few directories will be accessible inside the container,
including your home directory. The argument ```--bind /path/to/data``` tells
Apptainer to additionally make the given path available. You will probably need
to use this to access your data.

"Singularity" changed name to "Apptainer" in 2021. If you're using a slightly
older version, simply replace ```apptainer``` with ```singularity``` in the
commands above. Don't worry, it's the same software! Further documentation is
available on [their website](http://apptainer.org/).


Installation via package manager
--------------------------------

### Nix

CrystFEL is available through <a href="https://nixos.org/">NixOS</a> since
22.05, for x86_64, Darwin and aarch64.  Two packages are
available.  Package `crystfel` contains all tools including the GUI,
whereas `crystfel-headless` excludes the GUI, making it easier to
install and more suitable for backend data processing servers.

To install via NixOS, simply add the package globally to your
`environment.systemPackages`, or enter a temporary shell via `nix shell
nixpkgs#crystfel` to have all CrystFEL tools in your `PATH`.

### Homebrew

To install the development version of CrystFEL using
[Homebrew](https://brew.sh/), download the latest version of
[crystfel.rb](https://gitlab.desy.de/thomas.white/crystfel/-/blob/master/crystfel.rb)
from the Git repository, and use `brew install`:
```
$ wget https://gitlab.desy.de/thomas.white/crystfel/-/raw/master/crystfel.rb
$ brew install -v -s ./crystfel.rb --HEAD
```
Omit `--HEAD` to install the latest stable version.

Note that the Homebrew formulae included in the "tarball" release files always
refer to the *previous* release.  This is because the formula has to contain
the SHA256 hash of the downloaded tarball, which is not known until after the
tarball has been created!

We hope to contribute the Homebrew formula to the main repository soon, so that
you can install CrystFEL with a simple `brew install crystfel`.  See
[issue #94](https://gitlab.desy.de/thomas.white/crystfel/-/issues/94) for some
discussion.


Installation from source
------------------------

CrystFEL is compiled and installed using  [Meson](https://mesonbuild.com/).
First, download the latest Git version or unpack the package from the
[downloads page](https://desy.de/~twhite/crystfel/download.html):
```
$ git clone https://gitlab.desy.de/thomas.white/crystfel.git
$ cd crystfel
    or
$ tar -xzf crystfel-0.11.0.tar.gz    # or whichever other version
$ cd crystfel-0.11.0
```
Then, simply:
```
$ meson setup build
$ meson compile -C build
$ meson install -C build
```
The `meson setup` command will report if dependencies are missing - see the
next section for details.  If necessary, the `meson install` command will ask
for your password to gain administrative privileges.

Run `indexamajig` for a basic check that the installation has succeeded.  A
healthy newborn CrystFEL should complain that `You need to provide the input
filename (use -i)` when run with no other command-line options.

Alternatively, run `crystfel` to start the graphical user interface, provided
that the dependencies for the GUI were met (see above).

Refer to the [tutorial](doc/articles/tutorial.rst) to see where to go from
here!


Dependencies
------------

There are very few mandatory dependencies for CrystFEL.  They are:

* [Meson](https://mesonbuild.com/) 0.55.0 or later
* [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
* [Bison](https://www.gnu.org/software/bison/) 2.6 or later
* [Flex](https://www.gnu.org/software/flex/)

The following dependencies are "optional", in the sense that you can install
CrystFEL without them.  However, a CrystFEL installation without these will lack
important features such as the graphical user interface.  The following list is
roughly in order of importance:

* [HDF5](https://www.hdfgroup.org/downloads/hdf5/) 1.8.0 or later (required for HDF5 file read/write.  Version 1.10.0 or later is required for many recent data formats which use the 'virtual data set' feature)
* [GTK](https://gtk.org/) 3.12 or later, but not 4.x (required for GUI)
* [Cairo](https://www.cairographics.org/) 1.2 or later (required for GUI)
* [Pango](https://pango.gnome.org/) 1.0 or later, including [PangoCairo](https://docs.gtk.org/PangoCairo/) (required for GUI)
* [gdk-pixbuf](https://docs.gtk.org/gdk-pixbuf/) 2.0 or later (required for GUI)
* [libccp4](ftp://ftp.ccp4.ac.uk/opensource/) \[\*\] (required for MTZ import/export)
* [XGandalf](https://stash.desy.de/users/gevorkov/repos/xgandalf) \[\*\] (for `xgandalf` indexing)
* [Zlib](https://www.zlib.net/) \[\*\] (required for reading gzipped CBF files.  Version 1.2.3.5 or later preferred for better decompression speed)
* [PinkIndexer](https://stash.desy.de/users/gevorkov/repos/pinkindexer) \[\*\] (for indexing electron or wide bandwidth diffraction patterns)
* [FFTW](http://fftw.org/) 3.0 or later (required for `asdf` indexing)
* [FDIP](https://stash.desy.de/users/gevorkov/repos/fastdiffractionimageprocessing/) \[\*\] (for `peakFinder9` peak search algorithm)
* [libZMQ](https://github.com/zeromq/libzmq/) (for online data streaming)
* [libasapo-consumer](https://stash.desy.de/projects/ASAPO/repos/asapo/browse) (for online and offline data streaming via DESY's ASAP::O framework)
* [msgpack-c](https://github.com/msgpack/msgpack-c) (for streaming data in MsgPack format)
* [Seedee](https://gitlab.desy.de/fs-sc/seedee) (for streaming data serialised with Seedee)
* [cJSON](https://github.com/DaveGamble/cJSON/) \[\*\] (extra dependency if Seedee is found)
* [Pandoc](https://pandoc.org/) (to convert documentation to `man` format)

Most of the dependencies mentioned above should be available from your Linux
distribution's package manager, or from [Homebrew](https://brew.sh/) on Mac OS.
We emphatically recommend against trying to install GTK, Cairo, Pango or
gdk-pixbuf from source.

Dependencies marked with \[\*\] above will be downloaded and compiled
automatically by Meson if they are not available on the system.
If you don't want this, add option `--wrap-mode=nofallback` when invoking
Meson.  See the Meson manual for other possibilities, such as using
locally-provided files instead of downloading them.

If you install libraries in a non-system location, set `PKG_CONFIG_PATH` so
that they can be found.  For example:
```
$ export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/home/user/xgandalf/lib64/pkgconfig
```
If the libraries are automatically installed by Meson (see above about \[\*\]),
you will not have to worry about this.

We also do not recommend using dependencies from Conda/Anaconda.  Do not
activate any Conda environment before compiling CrystFEL, not even the "base"
environment.  Don't even "source" the Conda setup file before installing
CrystFEL - keep it completely separate.  A Conda recipe for CrystFEL might be
coming soon, though, if development resources allow for it.


Installing the indexing engines
-------------------------------

Processing data relies on indexing 'engines'.  At the absolute minimum, you
will have access to the [TakeTwo](https://journals.iucr.org/d/issues/2016/08/00/rr5128/)
indexing algorithm, since this is built into CrystFEL and does not have any
other dependencies.  [Xgandalf](https://journals.iucr.org/a/issues/2019/05/00/ae5071/),
[PinkIndexer](https://journals.iucr.org/a/issues/2020/02/00/ae5078/) and
`asdf` will very likely be available as well - they have some dependencies
(see above), but the required packages are quite common.  Still more indexing
methods can be enabled by installing one or more of the following external
programs:

* [Mosflm](https://www.mrc-lmb.cam.ac.uk/mosflm/mosflm/)
* [DirAx](http://www.crystal.chem.uu.nl/distr/dirax/)
* [XDS](http://xds.mpimf-heidelberg.mpg.de/)
* [Felix](https://doi.org/10.1107/S1600576717007506)

You can install these at any time before or after installing CrystFEL.  Note
that you only need the Mosflm binary, not the full `iMosflm` user interface.
[Download it here](https://www.mrc-lmb.cam.ac.uk/mosflm/mosflm/ver740/pre-built/mosflm-linux-64-noX11.zip).

You probably have CCP4 installed already, which comes with its own copy of `mosflm`.
However, you should keep your CCP4 installation separate from CrystFEL.  The
reason is that CCP4 comes with its own private copies of all its dependencies,
and it becomes very difficult to make the build process for CrystFEL (or
anything else) find the right versions of everything.  Do not 'source' the CCP4
setup file before trying to install CrystFEL, and make sure that the setup file
is not automatically referenced in your shell setup files (`.bashrc` and
others).

The script `install-indexers`, found in the `scripts` directory of the CrystFEL
source code, can help you to install Mosflm, DirAx and XDS.  Run
`scripts/install-indexers --help` for information.


Finding syminfo.lib
-------------------

The CCP4 libraries need to refer to symmetry information in a file called
`syminfo.lib`.  Unfortunately, they can only find this file if either `CLIBD`
or `SYMINFO` are set in the environment.  The CrystFEL GUI (`crystfel`),
`get_hkl` and `mosflm` therefore all need this variable to be set.  You could
set this variable in your shell setup file (e.g. as part of a `module load`),
or create wrapper scripts for these three programs.  The `install-indexers`
script, which assists with installing indexing engines, creates such a wrapper
for Mosflm already.


Compiling without HDF5
----------------------

When building with Meson, it's possible to compile without HDF5.  This is
useful if you only need to use other file formats, e.g. CBF.  You may have
noticed that most compilation problems (see *Notes about strange problems*
below) are related to HDF5.  To build without HDF5, set up your build directory
as follows, replacing the `meson build` step:
```
meson build -Dhdf5=disabled
```


Fedora 22 or later
------------------

All dependencies can be taken from the Fedora repositories:

```
$ sudo dnf group install 'Development Tools'
$ sudo dnf install hdf5-devel gsl-devel gtk3-devel cairo-devel pango-devel gdk-pixbuf2-devel meson gcc-c++ fftw-devel zeromq-devel msgpack-devel flex bison gcc-gfortran pandoc
$ cd /home/user/downloads/crystfel
$ meson build
$ ninja -C build
$ sudo ninja -C build install
```

Up to Fedora 32 (inclusive), you can also install `libccp4-devel` via `dnf`.
This package was removed starting from Fedora 33, but CrystFEL's build system
will take care of it.


Ubuntu 18.04 LTS ("Bionic") and 20.04 LTS ("Focal")
---------------------------------------------------

The Meson version in these Ubuntu releases is slightly too old for CrystFEL.
Install version 0.55.0 or later from the
[Meson website](https://mesonbuild.com/Getting-meson.html).
You don't need to "install" it, but do remember the location where you unpacked
it.  You will need to additionally install `python3`, if it's not already
present.

All other dependencies are available in the Ubuntu repositories.  Make sure
that the "universe" repository is enabled -
[instructions here](https://help.ubuntu.com/community/Repositories/Ubuntu).
```
$ apt install -y build-essential libhdf5-dev libgsl-dev \
                 libgtk-3-dev libcairo2-dev libpango1.0-dev \
                 libgdk-pixbuf2.0-dev libfftw3-dev \
                 git flex bison libzmq3-dev libmsgpack-dev \
                 libeigen3-dev libccp4-dev \
                 ninja-build python3
$ cd /home/user/downloads/crystfel
$ /home/user/downloads/meson/meson.py build
$ ninja -C build
$ sudo ninja -C build install
$ sudo ldconfig
```

CrystFEL will install the remaining dependencies (Xgandalf, PinkIndexer, FDIP)
automatically as part of its build process.


Ubuntu 16.04 LTS ("Xenial")
--------------------------

The Meson version in Ubuntu 16.04 is too old for CrystFEL, and unfortunately
the Python version is also "end of life" and too old for recent versions of
Meson.  However, the last release of Meson before it became incompatible with
older Python versions, 0.56.2, is new enough to build CrystFEL.

Ninja is not available in this Ubuntu release, but fortunately it's an easy
download.

The Eigen version is too old for XGandalf and PinkIndexer, so you should not
install Eigen from Ubuntu (`libeigen3-dev`), instead let Meson take care of
it during the build process.

Finally, the dynamic library loader has to be made aware that libraries have
been installed under `/usr/local`.

```
$ apt install -y build-essential libhdf5-dev libgsl-dev \
                 libgtk-3-dev libcairo2-dev libpango1.0-dev \
                 libgdk-pixbuf2.0-dev libfftw3-dev \
                 git flex bison libzmq3-dev libmsgpack-dev \
                 libccp4-dev wget python3 unzip

$ wget https://github.com/mesonbuild/meson/releases/download/0.56.2/meson-0.56.2.tar.gz
$ tar -xzf meson-0.56.2.tar.gz
$ export MESON=`pwd`/meson-0.56.2/meson.py

$ wget https://github.com/ninja-build/ninja/releases/download/v1.10.2/ninja-linux.zip
$ unzip ninja-linux.zip
$ export NINJA=`pwd`/ninja

$ cd /home/user/downloads/crystfel
$ $MESON build
$ $NINJA -C build
$ sudo $NINJA -C build install

$ sudo echo "/usr/local/lib/x86_64-linux-gnu/" > /etc/ld.so.conf.d/local.conf
$ sudo ldconfig
```

Debian 11 (Bullseye)
--------------------

All dependencies can be taken from the Debain repositories:

```
$ sudo apt install -y build-essential libhdf5-dev libgsl-dev \
                      libgtk-3-dev libcairo2-dev libpango1.0-dev \
                      libgdk-pixbuf2.0-dev libfftw3-dev \
                      git flex bison libzmq3-dev libmsgpack-dev \
                      libeigen3-dev libccp4-dev \
                      meson ninja-build
$ cd /home/user/downloads/crystfel
$ meson build
$ ninja -C build
$ sudo ninja -C build install
$ sudo ldconfig
```


Installation problems and solutions
-----------------------------------

* **Problem**: After installation, CrystFEL programs fail to start with an error about a missing shared object file:
    ```
    $ crystfel
    crystfel: error while loading shared libraries: libxgandalf.so: cannot open shared object file: No such file directory
    ```

    The error message might mention a different shared library, such as `libcrystfel.so`.

  **Explanation**: The default installation location for CrystFEL is `/usr/local`.
  Some Linux distributions (including Fedora) don't include this location in
  the default search path for libraries.  This problem will affect any program
  you install under `/usr/local`, not just CrystFEL.

  **Solution**: The best solution is to create (as root) a new file called
  `/etc/ld.so.conf.d/local.conf`, with the following two lines:
  ```
  /usr/local/lib
  /usr/local/lib64
  ```
  Then run `sudo ldconfig`.

  You can also simply run `sudo ldconfig /usr/local/lib /usr/local/lib64`, but
  the effect will not be permanent, e.g. across system updates.

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

  **Solution**: Two options:

  Option 1: Do not use your distribution's HDF5 package.  Instead, compile
  HDF5 yourself.  If you compile HDF5 with autotools (`./configure` et al.),
  take care to use `H5_CFLAGS` instead of `CFLAGS` to add any compiler flags.
  Ensure that the `h5cc` corresponding to the correct version of HDF5 is found,
  by making sure that its location comes first in `PATH` when compiling
  CrystFEL.

  Option 2: Compile HDF5 using CMake instead of autotools.  In this case,
  HDF5 will install a `pkg-config` file which Meson will use in preference to
  `h5cc`.


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

