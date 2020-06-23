#!/usr/bin/env python

import h5py
import numpy

blank = numpy.zeros((1,1), dtype=float)

with h5py.File('ev_enum1.h5', 'w') as fh:
    fh.create_dataset('/data/panelA/ev_1/panel_data1t/dataABCset/array', data=blank)
    fh.create_dataset('/data/panelA/ev_2/panel_data1t/dataDEFset/array', data=blank)
    fh.create_dataset('/data/panelA/ev_3/panel_data1t/dataGHIset/array', data=blank)
    fh.create_dataset('/data/panelA/ev_4/panel_data1t/dataKLMset/nomatch', data=blank)
    fh.create_dataset('/data/panelA/ev_5/panel_data1t/dataNOPset/array', data=blank)
    fh.create_dataset('/data/panelB/ev_1/panel_data1t/dataABCset/array', data=blank)
    fh.create_dataset('/data/panelB/ev_2/panel_data1t/dataDEFset/array', data=blank)
    fh.create_dataset('/data/panelB/ev_3/panel_data1t/dataGHIset/array', data=blank)
    fh.create_dataset('/data/panelB/ev_4/panel_data1t/dataKLMset/array', data=blank)
    fh.create_dataset('/data/panelB/ev_5/panel_data1t/dataNOPset/array', data=blank)
    fh.create_dataset('/data/panelB/ev_1/dataABCset/array', data=blank)
    fh.create_dataset('/data/panelB/ev_2/dataDEFset/array', data=blank)
    fh.create_dataset('/data/panelB/ev_3/dataGHIset/array', data=blank)
    fh.create_dataset('/data/panelB/ev_4/dataKLMset/array', data=blank)
    fh.create_dataset('/data/panelB/ev_5/dataNOPset/array', data=blank)
