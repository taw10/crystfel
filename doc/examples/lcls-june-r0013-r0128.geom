mask = /processing/hitfinder/masks
mask_good = 0x27
mask_bad = 0x00

; These default values will be used unless overridden by the per-panel values
peak_sep = 50.0
integr_radius = 10.0

; Upper panel (nearest the beam)
0/min_fs = 0
0/max_fs = 1023
0/min_ss = 512
0/max_ss = 1023
0/corner_x = -512.90
0/corner_y = 53.00
0/fs = x
0/ss = y
0/clen = 64.78e-3
0/res = 13333.3  ; 75 micron pixel size
0/badrow_direction = y
0/rigid_group = 0

; Lower panel (furthest from the beam)
1/min_fs = 0
1/max_fs = 1023
1/min_ss = 0
1/max_ss = 511
1/corner_x = -519.00
1/corner_y = -901.00
1/fs = x
1/ss = y
1/clen = 67.73e-3
1/res = 13333.3  ; 75 micron pixel size
1/badrow_direction = y
1/rigid_group = 1

bad_jetneartobeam/min_x = -15.0
bad_jetneartobeam/max_x = +15.0
bad_jetneartobeam/min_y = 71.3
bad_jetneartobeam/max_y = 159.3

bad_jetfarfrombeam/min_x = -25.0
bad_jetfarfrombeam/max_x = +25.0
bad_jetfarfrombeam/min_y = 159.3
bad_jetfarfrombeam/max_y = 600.0

bad_rectangle/min_x = -520.0
bad_rectangle/max_x = +520.0
bad_rectangle/min_y = -469.0
bad_rectangle/max_y = 0.0
