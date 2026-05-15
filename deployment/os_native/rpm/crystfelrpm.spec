Name: crystfel
Version: %{version_from_tag}
Release: 1%{?dist}
Summary: CrystFEL - Data processing for serial crystallography

License: GPLv3.0
URL: https://gitlab.desy.de/thomas.white/%{name}
%dnl Source0: https://gitlab.desy.de/thomas.white/%{name}.git

Prefix: /usr/local
%global _prefix /usr/local
%global debug_package %{nil}
%global __requires_exclude ^perl\\(
%global CFPREFIX /software/crystfel/devel
%global PKG_CONFIG_PATH %{_prefix}/%{_libdir}/%{name}/pkgconfig:%{_prefix}/lib/%{name}/pkgconfig
%global LD_LIBRARY_PATH %{_prefix}/%{_libdir}/%{name}:%{_prefix}/lib%{name}
%global LD_RUN_PATH %{_prefix}/%{_libdir}/%{name}:%{_prefix}/lib/%{name}
%global PATH %{_prefix}/bin/:$PATH
%global QA_RPATHS 0x0002

BuildRequires: which
BuildRequires: wget
BuildRequires: diffutils 
BuildRequires: sed
BuildRequires: flex
BuildRequires: bison
BuildRequires: cmake
BuildRequires: meson
BuildRequires: ninja-build
BuildRequires: gcc
BuildRequires: gcc-c++
BuildRequires: gcc-gfortran
BuildRequires: glibc-devel
BuildRequires: cppcheck
BuildRequires: cppcheck-htmlreport
BuildRequires: gtk3-devel
BuildRequires: cairo-devel
BuildRequires: pango-devel
BuildRequires: gdk-pixbuf2-devel
BuildRequires: fftw-devel
BuildRequires: libpng-devel
BuildRequires: zeromq-devel
BuildRequires: python3-devel
BuildRequires: lz4-devel
BuildRequires: bzip2-devel
BuildRequires: libcurl-devel
BuildRequires: lz4
BuildRequires: bzip2-libs
BuildRequires: pandoc


Requires: gtk3
Requires: cairo
Requires: pango
Requires: gdk-pixbuf2
Requires: libcurl-minimal
Requires: bzip2-libs
Requires: lz4
Requires: hdf5
Requires: python3
Requires: perl-libs


%description
CrystFEL is a suite of programs for processing data from serial crystallography experiments,
performed at synchrotron and X-ray free-electron laser facilities, as well as in your home lab using an electron microscope.

%prep
%dnl %auto_setup

%build
%dnl %configure
%dnl %make_build

meson setup build --prefix %{CFPREFIX}
ninja -v -C build

%check
ninja -v -C build test

%install
ninja -C build install
     
# Tweak CrystFEL GUI so that it can find syminfo.lib
mv %{CFPREFIX}/bin/crystfel %{CFPREFIX}/bin/crystfel.real
echo '#!/bin/sh' > %{CFPREFIX}/bin/crystfel
echo "export SYMINFO=%{_prefix}/share/ccp4/syminfo.lib" >> %{CFPREFIX}/bin/crystfel
echo "%{_prefix}/bin/crystfel.real \"\$@\"" >> %{CFPREFIX}/bin/crystfel
chmod +x %{CFPREFIX}/bin/crystfel
    
# Tweak get_hkl in the same way
mv %{CFPREFIX}/bin/get_hkl %{CFPREFIX}/bin/get_hkl.real
echo '#!/bin/sh' > %{CFPREFIX}/bin/get_hkl
echo "export SYMINFO=%{_prefix}/share/ccp4/syminfo.lib" >> %{CFPREFIX}/bin/get_hkl
echo "%{_prefix}/bin/get_hkl.real \"\$@\"" >> %{CFPREFIX}/bin/get_hkl
chmod +x %{CFPREFIX}/bin/get_hkl
    
# Mosflm (tweaked to find syminfo.lib and not load environ.def/default.def)
wget -nv https://www.mrc-lmb.cam.ac.uk/harry/imosflm/ver740/downloads/imosflm-7.4.0-linux-64.zip
unzip imosflm-7.4.0-linux-64.zip imosflm/bin/mosflm
mv imosflm/bin/mosflm %{CFPREFIX}/bin/mosflm.real
echo '#!/bin/sh' > %{CFPREFIX}/bin/mosflm
echo "export SYMINFO=%{_prefix}/share/ccp4/syminfo.lib" >> %{CFPREFIX}/bin/mosflm
echo "%{_prefix}/bin/mosflm.real -n \"\$@\"" >> %{CFPREFIX}/bin/mosflm
chmod +x %{CFPREFIX}/bin/mosflm

sed -i 's/python/python3/'   %{CFPREFIX}/bin/fg-graph
grep "python" %{CFPREFIX}/bin/fg-graph
sed -i 's/python/python3/'   %{CFPREFIX}/bin/peakogram-stream
grep "python" %{CFPREFIX}/bin/peakogram-stream

export QA_RPATHS=%{QA_RPATHS}
%dnl rm -rf $RPM_BUILD_ROOT
rm -rf %{buildroot}

cd %{CFPREFIX}

if [ -d bin ] ;
then
    mkdir -p %{buildroot}/%{_bindir} ;
    cp -rp bin/*  -t %{buildroot}/%{_bindir} ;
fi ;

if [ -d lib ] ;
then
    mkdir -p %{buildroot}/usr/local/lib/%{name} ;
    cp -rp lib/*  -t %{buildroot}/usr/local/lib/%{name} ;
fi ;

if [ -d lib64 ] ;
then
    mkdir -p %{buildroot}/%{_libdir}/%{name}/pkgconfig ;
    cp -rp lib64/*  -t %{buildroot}/%{_libdir}/%{name} ;
    cp -rp lib64/pkgconfig/*  -t %{buildroot}/%{_libdir}/%{name}/pkgconfig ;
fi ;

if [ -d include ] ;
then
    mkdir -p %{buildroot}/%{_includedir}/%{name} ;
    cp -rp include/*  -t %{buildroot}/%{_includedir}/%{name} ;
fi ;

if [ -d share  ] ;
then
    mkdir -p %{buildroot}/%{_datarootdir}/%{name} ;
    cp -rp share/*  -t %{buildroot}/%{_datarootdir}/%{name} ;
fi ;

%dnl install -m 0755 bin/%{name} %{buildroot}/%{_bindir}/%{name}

%files
# Package all subdirectories and files under /usr/local/*
/usr/local/*

# OR
# Package individual directories.

%dnl %license <path-to-file>
%dnl %doc <path-to-file>
%dnl %{_bindir}/*
%dnl /usr/local/lib/%{name}/*
%dnl %{_libdir}/%{name}/*
%dnl %{_includedir}/%{name}/*
%dnl %{_datarootdir}/%{name}/*
%dnl %{_libdir}/%{name}/pkgconfig/*

%changelog
* Wed May 8 2026 tirumala
- Initial version
