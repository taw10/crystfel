mask = /processing/hitfinder/masks
mask_good = 0x07
mask_bad = 0x00

; These default values will be used unless overridden by the per-panel values
peak_sep = 50.0
integr_radius = 10.0

upper/min_fs = 0
upper/max_fs = 1023
upper/min_ss = 512
upper/max_ss = 1023
upper/corner_x = -491.90
upper/corner_y = 71.30
upper/fs = x
upper/ss = y
upper/clen = 67.8e-3
upper/res = 13333.3  ; 75 micron pixel size
upper/badrow_direction = y

lower/min_fs = 0
lower/max_fs = 1023
lower/min_ss = 0
lower/max_ss = 511
lower/corner_x = -492.00
lower/corner_y = -779.70
lower/fs = x
lower/ss = y
lower/clen = 70.8e-3
lower/res = 13333.3  ; 75 micron pixel size
lower/badrow_direction = y

bad_jetneartobeam/min_x = -15.0
bad_jetneartobeam/max_x = +15.0
bad_jetneartobeam/min_y = 71.3
bad_jetneartobeam/max_y = 159.3

bad_jetfarfrombeam/min_x = -25.0
bad_jetfarfrombeam/max_x = +25.0
bad_jetfarfrombeam/min_y = 159.3
bad_jetfarfrombeam/max_y = 600.0
