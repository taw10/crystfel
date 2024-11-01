CrystFEL - Data processing for serial crystallography
=====================================================

Overview
--------

CrystFEL is a suite of programs for processing data from [serial
crystallography experiments](https://en.wikipedia.org/wiki/Serial_Femtosecond_Crystallography),
performed at synchrotron and X-ray free-electron laser facilities, as well as
in your home lab using an electron microscope.


Getting started
---------------

See [INSTALL.md](INSTALL.md) for installation instructions, including
our container registry, installation via package manager and details of
pre-existing installations at X-ray facilities around the world.

CrystFEL can be used from the command line or via a graphical user interface.
To start the graphical user interface, run ```crystfel```.

There is a [video tutorial](https://vimeo.com/585412404), as well as a [text
tutorial](doc/articles/tutorial.rst) to get you started with processing via the
GUI.

For command-line use, standard ```man``` pages are available.  Start with
```man crystfel```.  The manual pages are also
[available on the web](https://www.desy.de/~twhite/crystfel/manual.html).


Documentation
-------------

* [Basic Tutorial](doc/articles/tutorial.rst)
* [How to choose the right point group for merging](doc/articles/pointgroup.rst)
* [How to increase data processing speed](doc/articles/speed.rst)
* [Real-time data processing](doc/articles/online.rst)
* [Processing electron diffraction data](doc/articles/electrons.rst)
* [Symmetry classification for serial crystallography](doc/twin-calculator.pdf)
* [Matrix conventions used in CrystFEL code](doc/matrix-notation.pdf) - for
  developers, written mostly for my own benefit.
* [Hit rate graph](doc/hitrate.png)
* [Examples folder](doc/examples) - contains some template input files.
* [Contributing to CrystFEL](CONTRIBUTING.md) - including how to cite CrystFEL
  and how to find good first issues to work on.
* [Citation list](https://www.desy.de/~twhite/crystfel/citations.html) - please
  send us details of your paper, if it's missing!
* [Scripts folder](scripts) - a miscellany of smaller programs to help at
  various stages of data processing.


Journal articles and book chapters
----------------------------------

* [Processing serial crystallography data with CrystFEL: a step-by-step
  guide](https://doi.org/10.1107/S205979831801238X) - covers command-line
  processing only (pre-dates the GUI).
* [Recent developments in CrystFEL](http://dx.doi.org/10.1107/S1600576716004751) -
  now somewhat out of date, but contains some useful information about the
  algorithms used.
* [Crystallography and Molecular Imaging using X-ray
  Lasers](https://doi.org/10.23730/CYRSP-2018-001.605) - an introduction to the
  biological aspects and possibilities, written for physicists (in contrast to
  most other articles, which introduce the physical aspects for biologists!).
* [Original paper about CrystFEL](http://dx.doi.org/10.1107/S0021889812002312)
  from 2012.  Not open access, but a "reprint" is available
  [here](https://www.desy.de/~twhite/crystfel/db5097-reprint.pdf).
* [Climbing the Data Mountain: Processing of SFX
  Data](https://link.springer.com/chapter/10.1007/978-3-030-00551-1_7) -
  emphasizes data volume issues for XFELs.  Unfortunately not open access.
* [Processing of XFEL
  Data](https://link.springer.com/protocol/10.1007/978-1-4939-7000-1_13) -
  describes the entire processing pipeline.  Unfortunately not open access.


Awards
------

In 2017, the development of CrystFEL was recognised with the [Max von Laue
Prize](https://www.desy.de/news/news_search/index_eng.html?openDirectAnchor=1202)
from the [German Society for Crystallography (DGK)](https://dgk-home.de/en/).


Citing CrystFEL
---------------

Please see [CONTRIBUTING.md](CONTRIBUTING.md) for citation instructions.


Related software
----------------

[OnDA Monitor](https://www.ondamonitor.com/) for real-time monitoring of data
quality.  Read [the paper](https://doi.org/10.1107/s1600576716007469).

[Pixel Anomaly Detection Tool](https://github.com/gihankaushyal/PixelAnomalyDetectorTool)
uses machine learning techniques to find misbehaving detector pixels.
Read [the paper](https://doi.org/10.1107/s1600576724000116).

[Careless](https://github.com/rs-station/careless) scales and merges data
using variational inference.  Use its `stream2mtz` script to import data from
CrystFEL.  Read [the paper](http://dx.doi.org/10.1038/s41467-022-35280-8).

[DatView](https://github.com/nstander/DatView) helps with multivariate analysis
of large datasets.  Read [the paper](https://doi.org/10.1107/s1600576719012044).
[More information here](https://sites.google.com/view/zatsepinlab/resources/datview)
including a tutorial and manual.


Funding acknowledgements
------------------------

Development of CrystFEL is primarily funded by the
[Helmholtz Association](https://www.helmholtz.de/) via
[DESY](https://www.desy.de/).

Partial funding for CrystFEL has been provided by:

* The consortium DAPHNE4NFDI in the context of the work of the NFDI e.V. The
  consortium is funded by the Deutsche Forschungsgemeinschaft (DFG, German
  Research Foundation) - project number 460248799 (2022-).

* European Union’s Horizon 2020 research and innovation programme under grant
  agreement No 857641 ([ExPaNDS](https://expands.eu/)) (2019-2023).

* [X-Probe](http://x-probe.org/), a project of the European Union's 2020
  Research and Innovation Program Under the Marie Skłodowska-Curie grant
  agreement 637295 (2015-2018).

* The [BMBF](https://www.bmbf.de/) German-Russian Cooperation
  [SyncFELMed](http://www.syncfelmed.org/), grant 05K14CHA (2014-2017).

* [BioStruct-X](https://www.biostruct-x.eu/), a project funded by the Seventh
  Framework Programme (FP7) of the European Commission (2011-2016).


Licence
-------

Copyright © 2012-2024 Deutsches Elektronen-Synchrotron DESY, a research centre
of the Helmholtz Association.

See [AUTHORS](AUTHORS) as well as individual source code files for full details
of contributors.

CrystFEL is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CrystFEL is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CrystFEL.  If not, see <http://www.gnu.org/licenses/>.


