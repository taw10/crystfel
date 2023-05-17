# frozen_string_literal: true

# Suite for processing Serial Femtosecond Crystallography (SFX) data
class Crystfel < Formula
  desc 'Suite for processing Serial Femtosecond Crystallography (SFX) data'
  homepage 'https://www.desy.de/~twhite/crystfel/index.html'

  stable do
    url 'https://www.desy.de/~twhite/crystfel/crystfel-0.10.2.tar.gz'
    sha256 '9c23bd9dd0ca4b9e1b54df48845062095373aa4ec1b029635c5ace9a5c7eb0fe'
  end

  head do
    url 'https://gitlab.desy.de/thomas.white/crystfel.git'
  end

  depends_on 'bison' => :build
  depends_on 'meson' => :build
  depends_on 'ninja' => :build
  depends_on 'pkg-config' => :build
  depends_on 'argp-standalone' => :build
  depends_on 'cairo'
  depends_on 'fftw'
  depends_on 'gdk-pixbuf'
  depends_on 'gsl'
  depends_on 'gtk+3'
  depends_on 'hdf5'
  depends_on 'libpng'
  depends_on 'libtiff'
  depends_on 'pango'

  def install
      args = std_meson_args
      if OS.mac? 
        args << "-Dc_args=-I/opt/argp-standalone/include"
        args << "-Dc_link_args=-L/usr/local/opt/argp-standalone/lib -largp"
      end
      system 'meson', 'setup', *args, 'build'
      system 'ninja','-C','build'
      system 'ninja','-C','build','test'
      system 'ninja','-C','build','install'
  end

  test do
    # All functionality tests are done in the above section with make test
    system bin / 'ambigator', '--help'
    system bin / 'cell_explorer', '--help'
    system bin / 'cell_tool', '--help'
    system bin / 'check_hkl', '--help'
    system bin / 'compare_hkl', '--help'
    system bin / 'crystfel', '--help'
    system bin / 'detector-shift', '--help'
    system bin / 'get_hkl', '--help'
    system bin / 'indexamajig', '--help'
    system bin / 'list_events', '--help'
    system bin / 'make_pixelmap', '--help'
    system bin / 'partialator', '--help'
    system bin / 'peakogram-stream', '--help'
    system bin / 'process_hkl', '--help'
    system bin / 'render_hkl', '--help'
    system bin / 'whirligig', '--help'
  end
end