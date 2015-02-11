; Example of a geometry file for CSPAD data output by Cheetah in CXI format

photon_energy = /LCLS/photon_energy_eV
clen = /LCLS/detector_1/EncoderValue
coffset = 0.573224
adu_per_eV = 0.00338
res = 9097.53


; The following lines define how to interpret the three-dimensional data array
; containing the image data.  The first dimension is the event number, the
; second and third are spatial dimensions.

dim0 = %
dim1 = ss
dim2 = fs
data = /entry_1/data_1/data


; These following lines define where to find a "bad pixel mask" for each event,
; and how to interpret its contents.

mask = /entry_1/data_1/mask
mask_good = 0x0000
mask_bad = 0xffff


; The following lines define "rigid groups" which express the physical
; construction of the detector.  This is used when refining the detector
; geometry.

rigid_group_q0 = q0a0,q0a1,q0a2,q0a3,q0a4,q0a5,q0a6,q0a7,q0a8,q0a9,q0a10,q0a11,q0a12,q0a13,q0a14,q0a15
rigid_group_q1 = q1a0,q1a1,q1a2,q1a3,q1a4,q1a5,q1a6,q1a7,q1a8,q1a9,q1a10,q1a11,q1a12,q1a13,q1a14,q1a15
rigid_group_q2 = q2a0,q2a1,q2a2,q2a3,q2a4,q2a5,q2a6,q2a7,q2a8,q2a9,q2a10,q2a11,q2a12,q2a13,q2a14,q2a15
rigid_group_q3 = q3a0,q3a1,q3a2,q3a3,q3a4,q3a5,q3a6,q3a7,q3a8,q3a9,q3a10,q3a11,q3a12,q3a13,q3a14,q3a15

rigid_group_a0 = q0a0,q0a1
rigid_group_a1 = q0a2,q0a3
rigid_group_a2 = q0a4,q0a5
rigid_group_a3 = q0a6,q0a7
rigid_group_a4 = q0a8,q0a9
rigid_group_a5 = q0a10,q0a11
rigid_group_a6 = q0a12,q0a13
rigid_group_a7 = q0a14,q0a15
rigid_group_a8 = q1a0,q1a1
rigid_group_a9 = q1a2,q1a3
rigid_group_a10 = q1a4,q1a5
rigid_group_a11 = q1a6,q1a7
rigid_group_a12 = q1a8,q1a9
rigid_group_a13 = q1a10,q1a11
rigid_group_a14 = q1a12,q1a13
rigid_group_a15 = q1a14,q1a15
rigid_group_a16 = q2a0,q2a1
rigid_group_a17 = q2a2,q2a3
rigid_group_a18 = q2a4,q2a5
rigid_group_a19 = q2a6,q2a7
rigid_group_a20 = q2a8,q2a9
rigid_group_a21 = q2a10,q2a11
rigid_group_a22 = q2a12,q2a13
rigid_group_a23 = q2a14,q2a15
rigid_group_a24 = q3a0,q3a1
rigid_group_a25 = q3a2,q3a3
rigid_group_a26 = q3a4,q3a5
rigid_group_a27 = q3a6,q3a7
rigid_group_a28 = q3a8,q3a9
rigid_group_a29 = q3a10,q3a11
rigid_group_a30 = q3a12,q3a13
rigid_group_a31 = q3a14,q3a15

rigid_group_collection_quadrants = q0,q1,q2,q3
rigid_group_collection_asics = a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31


; The geometrical parameters for the panels follow.
; The above parameters can also be set for each panel individually, but usually
; it makes more sense to define them once for all panels.

q0a0/min_fs = 0
q0a0/min_ss = 0
q0a0/max_fs = 193
q0a0/max_ss = 184
q0a0/fs = +0.004806x +0.999989y
q0a0/ss = -0.999989x +0.004806y
q0a0/corner_x = 443.819
q0a0/corner_y = -49.8719

q0a1/min_fs = 194
q0a1/min_ss = 0
q0a1/max_fs = 387
q0a1/max_ss = 184
q0a1/fs = +0.004806x +0.999989y
q0a1/ss = -0.999989x +0.004806y
q0a1/corner_x = 444.766
q0a1/corner_y = 147.126

q0a2/min_fs = 0
q0a2/min_ss = 185
q0a2/max_fs = 193
q0a2/max_ss = 369
q0a2/fs = +0.003265x +0.999995y
q0a2/ss = -0.999995x +0.003265y
q0a2/corner_x = 239.8
q0a2/corner_y = -49.3504

q0a3/min_fs = 194
q0a3/min_ss = 185
q0a3/max_fs = 387
q0a3/max_ss = 369
q0a3/fs = +0.003265x +0.999995y
q0a3/ss = -0.999995x +0.003265y
q0a3/corner_x = 240.444
q0a3/corner_y = 147.649

q0a4/min_fs = 0
q0a4/min_ss = 370
q0a4/max_fs = 193
q0a4/max_ss = 554
q0a4/fs = -0.999997x +0.002424y
q0a4/ss = -0.002424x -0.999997y
q0a4/corner_x = 872.219
q0a4/corner_y = 342.054

q0a5/min_fs = 194
q0a5/min_ss = 370
q0a5/max_fs = 387
q0a5/max_ss = 554
q0a5/fs = -0.999997x +0.002424y
q0a5/ss = -0.002424x -0.999997y
q0a5/corner_x = 675.22
q0a5/corner_y = 342.532

q0a6/min_fs = 0
q0a6/min_ss = 555
q0a6/max_fs = 193
q0a6/max_ss = 739
q0a6/fs = -0.999997x +0.002685y
q0a6/ss = -0.002685x -0.999997y
q0a6/corner_x = 871.381
q0a6/corner_y = 135.836

q0a7/min_fs = 194
q0a7/min_ss = 555
q0a7/max_fs = 387
q0a7/max_ss = 739
q0a7/fs = -0.999997x +0.002685y
q0a7/ss = -0.002685x -0.999997y
q0a7/corner_x = 674.382
q0a7/corner_y = 136.365

q0a8/min_fs = 0
q0a8/min_ss = 740
q0a8/max_fs = 193
q0a8/max_ss = 924
q0a8/fs = -0.000078x -0.999999y
q0a8/ss = +0.999999x -0.000078y
q0a8/corner_x = 480.758
q0a8/corner_y = 769.64

q0a9/min_fs = 194
q0a9/min_ss = 740
q0a9/max_fs = 387
q0a9/max_ss = 924
q0a9/fs = -0.000078x -0.999999y
q0a9/ss = +0.999999x -0.000078y
q0a9/corner_x = 480.743
q0a9/corner_y = 572.64

q0a10/min_fs = 0
q0a10/min_ss = 925
q0a10/max_fs = 193
q0a10/max_ss = 1109
q0a10/fs = +0.001551x -0.999999y
q0a10/ss = +0.999999x +0.001551y
q0a10/corner_x = 689.447
q0a10/corner_y = 770.295

q0a11/min_fs = 194
q0a11/min_ss = 925
q0a11/max_fs = 387
q0a11/max_ss = 1109
q0a11/fs = +0.001551x -0.999999y
q0a11/ss = +0.999999x +0.001551y
q0a11/corner_x = 689.752
q0a11/corner_y = 573.296

q0a12/min_fs = 0
q0a12/min_ss = 1110
q0a12/max_fs = 193
q0a12/max_ss = 1294
q0a12/fs = -0.999998x -0.002161y
q0a12/ss = +0.002161x -0.999998y
q0a12/corner_x = 445.672
q0a12/corner_y = 751.701

q0a13/min_fs = 194
q0a13/min_ss = 1110
q0a13/max_fs = 387
q0a13/max_ss = 1294
q0a13/fs = -0.999998x -0.002161y
q0a13/ss = +0.002161x -0.999998y
q0a13/corner_x = 248.672
q0a13/corner_y = 751.276

q0a14/min_fs = 0
q0a14/min_ss = 1295
q0a14/max_fs = 193
q0a14/max_ss = 1479
q0a14/fs = -0.999999x -0.000074y
q0a14/ss = +0.000074x -0.999999y
q0a14/corner_x = 445.151
q0a14/corner_y = 541.081

q0a15/min_fs = 194
q0a15/min_ss = 1295
q0a15/max_fs = 387
q0a15/max_ss = 1479
q0a15/fs = -0.999999x -0.000074y
q0a15/ss = +0.000074x -0.999999y
q0a15/corner_x = 248.151
q0a15/corner_y = 541.066

q1a0/min_fs = 388
q1a0/min_ss = 0
q1a0/max_fs = 581
q1a0/max_ss = 184
q1a0/fs = -0.999990x -0.004167y
q1a0/ss = +0.004167x -0.999990y
q1a0/corner_x = 28.4776
q1a0/corner_y = 436.83

q1a1/min_fs = 582
q1a1/min_ss = 0
q1a1/max_fs = 775
q1a1/max_ss = 184
q1a1/fs = -0.999990x -0.004167y
q1a1/ss = +0.004167x -0.999990y
q1a1/corner_x = -168.52
q1a1/corner_y = 436.009

q1a2/min_fs = 388
q1a2/min_ss = 185
q1a2/max_fs = 581
q1a2/max_ss = 369
q1a2/fs = -1.000001x +0.000385y
q1a2/ss = -0.000385x -1.000001y
q1a2/corner_x = 29.3559
q1a2/corner_y = 226.978

q1a3/min_fs = 582
q1a3/min_ss = 185
q1a3/max_fs = 775
q1a3/max_ss = 369
q1a3/fs = -1.000001x +0.000385y
q1a3/ss = -0.000385x -1.000001y
q1a3/corner_x = -167.644
q1a3/corner_y = 227.054

q1a4/min_fs = 388
q1a4/min_ss = 370
q1a4/max_fs = 581
q1a4/max_ss = 554
q1a4/fs = +0.000539x -1.000000y
q1a4/ss = +1.000000x +0.000539y
q1a4/corner_x = -364.144
q1a4/corner_y = 859.163

q1a5/min_fs = 582
q1a5/min_ss = 370
q1a5/max_fs = 775
q1a5/max_ss = 554
q1a5/fs = +0.000539x -1.000000y
q1a5/ss = +1.000000x +0.000539y
q1a5/corner_x = -364.038
q1a5/corner_y = 662.163

q1a6/min_fs = 388
q1a6/min_ss = 555
q1a6/max_fs = 581
q1a6/max_ss = 739
q1a6/fs = -0.000337x -1.000000y
q1a6/ss = +1.000000x -0.000337y
q1a6/corner_x = -156.511
q1a6/corner_y = 857.902

q1a7/min_fs = 582
q1a7/min_ss = 555
q1a7/max_fs = 775
q1a7/max_ss = 739
q1a7/fs = -0.000337x -1.000000y
q1a7/ss = +1.000000x -0.000337y
q1a7/corner_x = -156.577
q1a7/corner_y = 660.902

q1a8/min_fs = 388
q1a8/min_ss = 740
q1a8/max_fs = 581
q1a8/max_ss = 924
q1a8/fs = +0.999996x +0.002303y
q1a8/ss = -0.002303x +0.999996y
q1a8/corner_x = -786.718
q1a8/corner_y = 463.506

q1a9/min_fs = 582
q1a9/min_ss = 740
q1a9/max_fs = 775
q1a9/max_ss = 924
q1a9/fs = +0.999996x +0.002303y
q1a9/ss = -0.002303x +0.999996y
q1a9/corner_x = -589.719
q1a9/corner_y = 463.959

q1a10/min_fs = 388
q1a10/min_ss = 925
q1a10/max_fs = 581
q1a10/max_ss = 1109
q1a10/fs = +0.999997x +0.001741y
q1a10/ss = -0.001741x +0.999997y
q1a10/corner_x = -787.022
q1a10/corner_y = 668.135

q1a11/min_fs = 582
q1a11/min_ss = 925
q1a11/max_fs = 775
q1a11/max_ss = 1109
q1a11/fs = +0.999997x +0.001741y
q1a11/ss = -0.001741x +0.999997y
q1a11/corner_x = -590.022
q1a11/corner_y = 668.478

q1a12/min_fs = 388
q1a12/min_ss = 1110
q1a12/max_fs = 581
q1a12/max_ss = 1294
q1a12/fs = -0.000201x -0.999999y
q1a12/ss = +0.999999x -0.000201y
q1a12/corner_x = -761.085
q1a12/corner_y = 428.541

q1a13/min_fs = 582
q1a13/min_ss = 1110
q1a13/max_fs = 775
q1a13/max_ss = 1294
q1a13/fs = -0.000201x -0.999999y
q1a13/ss = +0.999999x -0.000201y
q1a13/corner_x = -761.125
q1a13/corner_y = 231.541

q1a14/min_fs = 388
q1a14/min_ss = 1295
q1a14/max_fs = 581
q1a14/max_ss = 1479
q1a14/fs = +0.003097x -0.999995y
q1a14/ss = +0.999995x +0.003097y
q1a14/corner_x = -559.624
q1a14/corner_y = 428.347

q1a15/min_fs = 582
q1a15/min_ss = 1295
q1a15/max_fs = 775
q1a15/max_ss = 1479
q1a15/fs = +0.003097x -0.999995y
q1a15/ss = +0.999995x +0.003097y
q1a15/corner_x = -559.014
q1a15/corner_y = 231.348

q2a0/min_fs = 776
q2a0/min_ss = 0
q2a0/max_fs = 969
q2a0/max_ss = 184
q2a0/fs = -0.004086x -0.999991y
q2a0/ss = +0.999991x -0.004086y
q2a0/corner_x = -442.346
q2a0/corner_y = 20.3382

q2a1/min_fs = 970
q2a1/min_ss = 0
q2a1/max_fs = 1163
q2a1/max_ss = 184
q2a1/fs = -0.004086x -0.999991y
q2a1/ss = +0.999991x -0.004086y
q2a1/corner_x = -443.151
q2a1/corner_y = -176.66

q2a2/min_fs = 776
q2a2/min_ss = 185
q2a2/max_fs = 969
q2a2/max_ss = 369
q2a2/fs = +0.000302x -1.000000y
q2a2/ss = +1.000000x +0.000302y
q2a2/corner_x = -235.519
q2a2/corner_y = 19.2312

q2a3/min_fs = 970
q2a3/min_ss = 185
q2a3/max_fs = 1163
q2a3/max_ss = 369
q2a3/fs = +0.000302x -1.000000y
q2a3/ss = +1.000000x +0.000302y
q2a3/corner_x = -235.459
q2a3/corner_y = -177.769

q2a4/min_fs = 776
q2a4/min_ss = 370
q2a4/max_fs = 969
q2a4/max_ss = 554
q2a4/fs = +0.999997x -0.002037y
q2a4/ss = +0.002037x +0.999997y
q2a4/corner_x = -863.817
q2a4/corner_y = -370.344

q2a5/min_fs = 970
q2a5/min_ss = 370
q2a5/max_fs = 1163
q2a5/max_ss = 554
q2a5/fs = +0.999997x -0.002037y
q2a5/ss = +0.002037x +0.999997y
q2a5/corner_x = -666.817
q2a5/corner_y = -370.746

q2a6/min_fs = 776
q2a6/min_ss = 555
q2a6/max_fs = 969
q2a6/max_ss = 739
q2a6/fs = +1.000000x -0.001155y
q2a6/ss = +0.001155x +1.000000y
q2a6/corner_x = -863.549
q2a6/corner_y = -165.126

q2a7/min_fs = 970
q2a7/min_ss = 555
q2a7/max_fs = 1163
q2a7/max_ss = 739
q2a7/fs = +1.000000x -0.001155y
q2a7/ss = +0.001155x +1.000000y
q2a7/corner_x = -666.549
q2a7/corner_y = -165.353

q2a8/min_fs = 776
q2a8/min_ss = 740
q2a8/max_fs = 969
q2a8/max_ss = 924
q2a8/fs = +0.002076x +0.999998y
q2a8/ss = -0.999998x +0.002076y
q2a8/corner_x = -473.62
q2a8/corner_y = -793.473

q2a9/min_fs = 970
q2a9/min_ss = 740
q2a9/max_fs = 1163
q2a9/max_ss = 924
q2a9/fs = +0.002076x +0.999998y
q2a9/ss = -0.999998x +0.002076y
q2a9/corner_x = -473.211
q2a9/corner_y = -596.474

q2a10/min_fs = 776
q2a10/min_ss = 925
q2a10/max_fs = 969
q2a10/max_ss = 1109
q2a10/fs = +0.004134x +0.999991y
q2a10/ss = -0.999991x +0.004134y
q2a10/corner_x = -676.809
q2a10/corner_y = -792.653

q2a11/min_fs = 970
q2a11/min_ss = 925
q2a11/max_fs = 1163
q2a11/max_ss = 1109
q2a11/fs = +0.004134x +0.999991y
q2a11/ss = -0.999991x +0.004134y
q2a11/corner_x = -675.995
q2a11/corner_y = -595.655

q2a12/min_fs = 776
q2a12/min_ss = 1110
q2a12/max_fs = 969
q2a12/max_ss = 1294
q2a12/fs = +0.999981x -0.006417y
q2a12/ss = +0.006417x +0.999981y
q2a12/corner_x = -442.034
q2a12/corner_y = -769.447

q2a13/min_fs = 970
q2a13/min_ss = 1110
q2a13/max_fs = 1163
q2a13/max_ss = 1294
q2a13/fs = +0.999981x -0.006417y
q2a13/ss = +0.006417x +0.999981y
q2a13/corner_x = -245.038
q2a13/corner_y = -770.711

q2a14/min_fs = 776
q2a14/min_ss = 1295
q2a14/max_fs = 969
q2a14/max_ss = 1479
q2a14/fs = +0.999996x -0.002727y
q2a14/ss = +0.002727x +0.999996y
q2a14/corner_x = -441.283
q2a14/corner_y = -566.627

q2a15/min_fs = 970
q2a15/min_ss = 1295
q2a15/max_fs = 1163
q2a15/max_ss = 1479
q2a15/fs = +0.999996x -0.002727y
q2a15/ss = +0.002727x +0.999996y
q2a15/corner_x = -244.283
q2a15/corner_y = -567.164

q3a0/min_fs = 1164
q3a0/min_ss = 0
q3a0/max_fs = 1357
q3a0/max_ss = 184
q3a0/fs = +0.999988x -0.004965y
q3a0/ss = +0.004965x +0.999988y
q3a0/corner_x = -33.3507
q3a0/corner_y = -458.693

q3a1/min_fs = 1358
q3a1/min_ss = 0
q3a1/max_fs = 1551
q3a1/max_ss = 184
q3a1/fs = +0.999988x -0.004965y
q3a1/ss = +0.004965x +0.999988y
q3a1/corner_x = 163.647
q3a1/corner_y = -459.671

q3a2/min_fs = 1164
q3a2/min_ss = 185
q3a2/max_fs = 1357
q3a2/max_ss = 369
q3a2/fs = +0.999998x -0.002316y
q3a2/ss = +0.002316x +0.999998y
q3a2/corner_x = -31.8316
q3a2/corner_y = -254.931

q3a3/min_fs = 1358
q3a3/min_ss = 185
q3a3/max_fs = 1551
q3a3/max_ss = 369
q3a3/fs = +0.999998x -0.002316y
q3a3/ss = +0.002316x +0.999998y
q3a3/corner_x = 165.168
q3a3/corner_y = -255.388

q3a4/min_fs = 1164
q3a4/min_ss = 370
q3a4/max_fs = 1357
q3a4/max_ss = 554
q3a4/fs = +0.002474x +0.999997y
q3a4/ss = -0.999997x +0.002474y
q3a4/corner_x = 359.553
q3a4/corner_y = -886.512

q3a5/min_fs = 1358
q3a5/min_ss = 370
q3a5/max_fs = 1551
q3a5/max_ss = 554
q3a5/fs = +0.002474x +0.999997y
q3a5/ss = -0.999997x +0.002474y
q3a5/corner_x = 360.04
q3a5/corner_y = -689.512

q3a6/min_fs = 1164
q3a6/min_ss = 555
q3a6/max_fs = 1357
q3a6/max_ss = 739
q3a6/fs = +0.000059x +1.000000y
q3a6/ss = -1.000000x +0.000059y
q3a6/corner_x = 154.142
q3a6/corner_y = -884.763

q3a7/min_fs = 1358
q3a7/min_ss = 555
q3a7/max_fs = 1551
q3a7/max_ss = 739
q3a7/fs = +0.000059x +1.000000y
q3a7/ss = -1.000000x +0.000059y
q3a7/corner_x = 154.154
q3a7/corner_y = -687.763

q3a8/min_fs = 1164
q3a8/min_ss = 740
q3a8/max_fs = 1357
q3a8/max_ss = 924
q3a8/fs = -0.999993x +0.004040y
q3a8/ss = -0.004040x -0.999993y
q3a8/corner_x = 784.877
q3a8/corner_y = -492.935

q3a9/min_fs = 1358
q3a9/min_ss = 740
q3a9/max_fs = 1551
q3a9/max_ss = 924
q3a9/fs = -0.999993x +0.004040y
q3a9/ss = -0.004040x -0.999993y
q3a9/corner_x = 587.878
q3a9/corner_y = -492.139

q3a10/min_fs = 1164
q3a10/min_ss = 925
q3a10/max_fs = 1357
q3a10/max_ss = 1109
q3a10/fs = -0.999971x +0.007529y
q3a10/ss = -0.007529x -0.999971y
q3a10/corner_x = 784.254
q3a10/corner_y = -699.59

q3a11/min_fs = 1358
q3a11/min_ss = 925
q3a11/max_fs = 1551
q3a11/max_ss = 1109
q3a11/fs = -0.999971x +0.007529y
q3a11/ss = -0.007529x -0.999971y
q3a11/corner_x = 587.26
q3a11/corner_y = -698.107

q3a12/min_fs = 1164
q3a12/min_ss = 1110
q3a12/max_fs = 1357
q3a12/max_ss = 1294
q3a12/fs = +0.004516x +0.999990y
q3a12/ss = -0.999990x +0.004516y
q3a12/corner_x = 769.176
q3a12/corner_y = -460.51

q3a13/min_fs = 1358
q3a13/min_ss = 1110
q3a13/max_fs = 1551
q3a13/max_ss = 1294
q3a13/fs = +0.004516x +0.999990y
q3a13/ss = -0.999990x +0.004516y
q3a13/corner_x = 770.066
q3a13/corner_y = -263.512

q3a14/min_fs = 1164
q3a14/min_ss = 1295
q3a14/max_fs = 1357
q3a14/max_ss = 1479
q3a14/fs = +0.004918x +0.999989y
q3a14/ss = -0.999989x +0.004918y
q3a14/corner_x = 554.764
q3a14/corner_y = -460.25

q3a15/min_fs = 1358
q3a15/min_ss = 1295
q3a15/max_fs = 1551
q3a15/max_ss = 1479
q3a15/fs = +0.004918x +0.999989y
q3a15/ss = -0.999989x +0.004918y
q3a15/corner_x = 555.732
q3a15/corner_y = -263.253
