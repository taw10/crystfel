CrystFEL - Data processing for serial crystallography
=====================================================

Copyright © 2012-2023 Deutsches Elektronen-Synchrotron DESY,
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

See INSTALL.md for installation instructions.

The best way to get started, after installation, is to run command ```crystfel```
to start the graphical user interface.


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

- European Union’s Horizon 2020 research and innovation programme under grant
  agreement No 857641 (ExPaNDS) (2019-2023).

- "X-Probe", a project of the European Union's 2020 Research and Innovation
  Program Under the Marie Skłodowska-Curie grant agreement 637295 (2015-2018).

- The BMBF German-Russian Cooperation "SyncFELMed", grant 05K14CHA (2014-2017).

- BioStruct-X, a project funded by the Seventh Framework Programme (FP7) of the
  European Commission (2011-2016).

