; Example geometry file for Eiger 16M detector, using its native file format
; and binning 2.

; Camera length (in m) and photon energy (eV)
clen = 0.1
photon_energy = 22000

; adu_per_photon needs a relatively recent CrystFEL version.  If your version is
; older, change it to adu_per_eV and set it to one over the photon energy in eV
adu_per_photon = 1
res = 13333.3   ; 75 micron pixel size

; These lines describe the data layout for the Eiger native multi-event files
dim0 = %
dim1 = ss
dim2 = fs
data = /entry/data/data

; Mask out strips between panels
bad_v0/min_fs = 1030
bad_v0/min_ss = 0
bad_v0/max_fs = 1039
bad_v0/max_ss = 2166

bad_h0/min_fs = 0
bad_h0/min_ss = 514
bad_h0/max_fs = 2069
bad_h0/max_ss = 550

bad_h1/min_fs = 0
bad_h1/min_ss = 1065
bad_h1/max_fs = 2069
bad_h1/max_ss = 1101

bad_h2/min_fs = 0
bad_h2/min_ss = 1616
bad_h2/max_fs = 2069
bad_h2/max_ss = 1652

; Uncomment these lines if you have a separate bad pixel map (recommended!)
;mask_file = eiger-badmap.h5
;mask = /data/data
;mask_good = 0x0
;mask_bad = 0x1

; corner_{x,y} set the position of the corner of the detector (in pixels)
; relative to the beam
panel0/min_fs = 0
panel0/min_ss = 0
panel0/max_fs = 2069
panel0/max_ss = 2166
panel0/corner_x = -1000.0
panel0/corner_y = -1000.0
panel0/fs = x
panel0/ss = y
