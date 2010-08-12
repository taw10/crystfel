n_panels = 2

; Upper panel (found to the "right" in the HDF5, nearest the beam)
0/min_x = 512
0/max_x = 1023
0/min_y = 0
0/max_y = 1023
0/cx = 459.0
0/cy = 511.0
0/clen = 64.6e-3
0/res = 13333.3  ; 75 micron pixel size
0/badrow_direction = x

; Lower panel (found to the "left" in the HDF5, furthest from the beam)
1/min_x = 0
1/max_x = 511
1/min_y = 0
1/max_y = 1023
1/cx = 901.0
1/cy = 519.0
1/clen = 67.4e-3
1/res = 13333.3  ; 75 micron pixel size
1/badrow_direction = x
