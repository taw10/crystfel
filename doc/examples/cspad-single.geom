; Example of a CrystFEL geometry file for CSPAD data in single-event HDF5
; format, with the panel data arranged as by Cheetah

adu_per_eV = 0.00105                   ; correct value for CSPAD 1.0 (i.e. old versions only)
                                       ; for newer versions (since about 2013), use 0.00338 instead
res = 9090.91                          ; pixels per metre
clen = /LCLS/detectorPosition          ; camera length from HDF5 file
coffset = 0.0                          ; no adjustment to camera length from HDF5 file
photon_energy = /LCLS/photon_energy_eV ; photon energy (in eV) from HDF5 file
data = /data/rawdata                   ; where to find the image data in the HDF5 file


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

q0a0/min_fs =    0
q0a0/min_ss =    0
q0a0/max_fs =  193
q0a0/max_ss =  184
q0a0/fs = -0.0057225772x +0.9999836087y
q0a0/ss = -0.9999836087x -0.0057225772y
q0a0/corner_x =       424.384
q0a0/corner_y =      -10.6473

q0a1/min_fs =  194
q0a1/min_ss =    0
q0a1/max_fs =  387
q0a1/max_ss =  184
q0a1/fs = -0.0057225772x +0.9999836087y
q0a1/ss = -0.9999836087x -0.0057225772y
q0a1/corner_x =       423.257
q0a1/corner_y =       186.349

q1a0/min_fs =  388
q1a0/min_ss =    0
q1a0/max_fs =  581
q1a0/max_ss =  184
q1a0/fs = -0.9999815226x +0.0060744919y
q1a0/ss = -0.0060744919x -0.9999815226y
q1a0/corner_x =       13.1196
q1a0/corner_y =       422.727

q1a1/min_fs =  582
q1a1/min_ss =    0
q1a1/max_fs =  775
q1a1/max_ss =  184
q1a1/fs = -0.9999815226x +0.0060744919y
q1a1/ss = -0.0060744919x -0.9999815226y
q1a1/corner_x =      -183.877
q1a1/corner_y =       423.924

q2a0/min_fs =  776
q2a0/min_ss =    0
q2a0/max_fs =  969
q2a0/max_ss =  184
q2a0/fs = +0.0048825117x -0.9999880791y
q2a0/ss = +0.9999880791x +0.0048825117y
q2a0/corner_x =      -428.519
q2a0/corner_y =       8.29208

q2a1/min_fs =  970
q2a1/min_ss =    0
q2a1/max_fs = 1163
q2a1/max_ss =  184
q2a1/fs = +0.0048825117x -0.9999880791y
q2a1/ss = +0.9999880791x +0.0048825117y
q2a1/corner_x =      -427.557
q2a1/corner_y =      -188.706

q3a0/min_fs = 1164
q3a0/min_ss =    0
q3a0/max_fs = 1357
q3a0/max_ss =  184
q3a0/fs = +0.9999964833x +0.0026521713y
q3a0/ss = -0.0026521713x +0.9999964833y
q3a0/corner_x =      -12.1295
q3a0/corner_y =      -424.046

q3a1/min_fs = 1358
q3a1/min_ss =    0
q3a1/max_fs = 1551
q3a1/max_ss =  184
q3a1/fs = +0.9999964833x +0.0026521713y
q3a1/ss = -0.0026521713x +0.9999964833y
q3a1/corner_x =       184.870
q3a1/corner_y =      -423.523

q0a2/min_fs =    0
q0a2/min_ss =  185
q0a2/max_fs =  193
q0a2/max_ss =  369
q0a2/fs = -0.0003575958x +0.9999999404y
q0a2/ss = -0.9999999404x -0.0003575958y
q0a2/corner_x =       211.039
q0a2/corner_y =      -10.9889

q0a3/min_fs =  194
q0a3/min_ss =  185
q0a3/max_fs =  387
q0a3/max_ss =  369
q0a3/fs = -0.0003575958x +0.9999999404y
q0a3/ss = -0.9999999404x -0.0003575958y
q0a3/corner_x =       210.969
q0a3/corner_y =       186.011

q1a2/min_fs =  388
q1a2/min_ss =  185
q1a2/max_fs =  581
q1a2/max_ss =  369
q1a2/fs = -0.9999999404x +0.0003117446y
q1a2/ss = -0.0003117446x -0.9999999404y
q1a2/corner_x =       13.3322
q1a2/corner_y =       208.973

q1a3/min_fs =  582
q1a3/min_ss =  185
q1a3/max_fs =  775
q1a3/max_ss =  369
q1a3/fs = -0.9999999404x +0.0003117446y
q1a3/ss = -0.0003117446x -0.9999999404y
q1a3/corner_x =      -183.668
q1a3/corner_y =       209.034

q2a2/min_fs =  776
q2a2/min_ss =  185
q2a2/max_fs =  969
q2a2/max_ss =  369
q2a2/fs = -0.0004957085x -0.9999998808y
q2a2/ss = +0.9999998808x -0.0004957085y
q2a2/corner_x =      -215.108
q2a2/corner_y =       8.41637

q2a3/min_fs =  970
q2a3/min_ss =  185
q2a3/max_fs = 1163
q2a3/max_ss =  369
q2a3/fs = -0.0004957085x -0.9999998808y
q2a3/ss = +0.9999998808x -0.0004957085y
q2a3/corner_x =      -215.206
q2a3/corner_y =      -188.584

q3a2/min_fs = 1164
q3a2/min_ss =  185
q3a2/max_fs = 1357
q3a2/max_ss =  369
q3a2/fs = +1.0000000000x -0.0001447762y
q3a2/ss = +0.0001447762x +1.0000000000y
q3a2/corner_x =      -12.1720
q3a2/corner_y =      -210.466

q3a3/min_fs = 1358
q3a3/min_ss =  185
q3a3/max_fs = 1551
q3a3/max_ss =  369
q3a3/fs = +1.0000000000x -0.0001447762y
q3a3/ss = +0.0001447762x +1.0000000000y
q3a3/corner_x =       184.828
q3a3/corner_y =      -210.495

q0a4/min_fs =    0
q0a4/min_ss =  370
q0a4/max_fs =  193
q0a4/max_ss =  554
q0a4/fs = -0.9999362230x -0.0112929661y
q0a4/ss = +0.0112929661x -0.9999362230y
q0a4/corner_x =       840.468
q0a4/corner_y =       392.435

q0a5/min_fs =  194
q0a5/min_ss =  370
q0a5/max_fs =  387
q0a5/max_ss =  554
q0a5/fs = -0.9999362230x -0.0112929661y
q0a5/ss = +0.0112929661x -0.9999362230y
q0a5/corner_x =       643.481
q0a5/corner_y =       390.210

q1a4/min_fs =  388
q1a4/min_ss =  370
q1a4/max_fs =  581
q1a4/max_ss =  554
q1a4/fs = -0.0003113084x -0.9999999404y
q1a4/ss = +0.9999999404x -0.0003113084y
q1a4/corner_x =      -387.060
q1a4/corner_y =       844.837

q1a5/min_fs =  582
q1a5/min_ss =  370
q1a5/max_fs =  775
q1a5/max_ss =  554
q1a5/fs = -0.0003113084x -0.9999999404y
q1a5/ss = +0.9999999404x -0.0003113084y
q1a5/corner_x =      -387.121
q1a5/corner_y =       647.837

q2a4/min_fs =  776
q2a4/min_ss =  370
q2a4/max_fs =  969
q2a4/max_ss =  554
q2a4/fs = +0.9999818802x +0.0060196919y
q2a4/ss = -0.0060196919x +0.9999818802y
q2a4/corner_x =      -849.683
q2a4/corner_y =      -392.889

q2a5/min_fs =  970
q2a5/min_ss =  370
q2a5/max_fs = 1163
q2a5/max_ss =  554
q2a5/fs = +0.9999818802x +0.0060196919y
q2a5/ss = -0.0060196919x +0.9999818802y
q2a5/corner_x =      -652.687
q2a5/corner_y =      -391.703

q3a4/min_fs = 1164
q3a4/min_ss =  370
q3a4/max_fs = 1357
q3a4/max_ss =  554
q3a4/fs = -0.0023373319x +0.9999972582y
q3a4/ss = -0.9999972582x -0.0023373319y
q3a4/corner_x =       388.452
q3a4/corner_y =      -844.271

q3a5/min_fs = 1358
q3a5/min_ss =  370
q3a5/max_fs = 1551
q3a5/max_ss =  554
q3a5/fs = -0.0023373319x +0.9999972582y
q3a5/ss = -0.9999972582x -0.0023373319y
q3a5/corner_x =       387.991
q3a5/corner_y =      -647.271

q0a6/min_fs =    0
q0a6/min_ss =  555
q0a6/max_fs =  193
q0a6/max_ss =  739
q0a6/fs = -0.9999533892x -0.0096549187y
q0a6/ss = +0.0096549187x -0.9999533892y
q0a6/corner_x =       842.954
q0a6/corner_y =       182.508

q0a7/min_fs =  194
q0a7/min_ss =  555
q0a7/max_fs =  387
q0a7/max_ss =  739
q0a7/fs = -0.9999533892x -0.0096549187y
q0a7/ss = +0.0096549187x -0.9999533892y
q0a7/corner_x =       645.963
q0a7/corner_y =       180.606

q1a6/min_fs =  388
q1a6/min_ss =  555
q1a6/max_fs =  581
q1a6/max_ss =  739
q1a6/fs = +0.0007214849x -0.9999997616y
q1a6/ss = +0.9999997616x +0.0007214849y
q1a6/corner_x =      -174.672
q1a6/corner_y =       844.923

q1a7/min_fs =  582
q1a7/min_ss =  555
q1a7/max_fs =  775
q1a7/max_ss =  739
q1a7/fs = +0.0007214849x -0.9999997616y
q1a7/ss = +0.9999997616x +0.0007214849y
q1a7/corner_x =      -174.530
q1a7/corner_y =       647.923

q2a6/min_fs =  776
q2a6/min_ss =  555
q2a6/max_fs =  969
q2a6/max_ss =  739
q2a6/fs = +0.9999990463x +0.0013817062y
q2a6/ss = -0.0013817062x +0.9999990463y
q2a6/corner_x =      -850.035
q2a6/corner_y =      -179.939

q2a7/min_fs =  970
q2a7/min_ss =  555
q2a7/max_fs = 1163
q2a7/max_ss =  739
q2a7/fs = +0.9999990463x +0.0013817062y
q2a7/ss = -0.0013817062x +0.9999990463y
q2a7/corner_x =      -653.035
q2a7/corner_y =      -179.667

q3a6/min_fs = 1164
q3a6/min_ss =  555
q3a6/max_fs = 1357
q3a6/max_ss =  739
q3a6/fs = -0.0021083837x +0.9999977946y
q3a6/ss = -0.9999977946x -0.0021083837y
q3a6/corner_x =       175.688
q3a6/corner_y =      -844.681

q3a7/min_fs = 1358
q3a7/min_ss =  555
q3a7/max_fs = 1551
q3a7/max_ss =  739
q3a7/fs = -0.0021083837x +0.9999977946y
q3a7/ss = -0.9999977946x -0.0021083837y
q3a7/corner_x =       175.273
q3a7/corner_y =      -647.681

q0a8/min_fs =    0
q0a8/min_ss =  740
q0a8/max_fs =  193
q0a8/max_ss =  924
q0a8/fs = +0.0078456579x -0.9999692440y
q0a8/ss = +0.9999692440x +0.0078456579y
q0a8/corner_x =       440.576
q0a8/corner_y =       811.780

q0a9/min_fs =  194
q0a9/min_ss =  740
q0a9/max_fs =  387
q0a9/max_ss =  924
q0a9/fs = +0.0078456579x -0.9999692440y
q0a9/ss = +0.9999692440x +0.0078456579y
q0a9/corner_x =       442.122
q0a9/corner_y =       614.786

q1a8/min_fs =  388
q1a8/min_ss =  740
q1a8/max_fs =  581
q1a8/max_ss =  924
q1a8/fs = +0.9999817014x +0.0060460833y
q1a8/ss = -0.0060460833x +0.9999817014y
q1a8/corner_x =      -808.912
q1a8/corner_y =       453.082

q1a9/min_fs =  582
q1a9/min_ss =  740
q1a9/max_fs =  775
q1a9/max_ss =  924
q1a9/fs = +0.9999817014x +0.0060460833y
q1a9/ss = -0.0060460833x +0.9999817014y
q1a9/corner_x =      -611.915
q1a9/corner_y =       454.273

q2a8/min_fs =  776
q2a8/min_ss =  740
q2a8/max_fs =  969
q2a8/max_ss =  924
q2a8/fs = -0.0027005512x +0.9999963641y
q2a8/ss = -0.9999963641x -0.0027005512y
q2a8/corner_x =      -451.888
q2a8/corner_y =      -813.482

q2a9/min_fs =  970
q2a9/min_ss =  740
q2a9/max_fs = 1163
q2a9/max_ss =  924
q2a9/fs = -0.0027005512x +0.9999963641y
q2a9/ss = -0.9999963641x -0.0027005512y
q2a9/corner_x =      -452.420
q2a9/corner_y =      -616.482

q3a8/min_fs = 1164
q3a8/min_ss =  740
q3a8/max_fs = 1357
q3a8/max_ss =  924
q3a8/fs = -0.9999983311x -0.0018213630y
q3a8/ss = +0.0018213630x -0.9999983311y
q3a8/corner_x =       806.687
q3a8/corner_y =      -442.882

q3a9/min_fs = 1358
q3a9/min_ss =  740
q3a9/max_fs = 1551
q3a9/max_ss =  924
q3a9/fs = -0.9999983311x -0.0018213630y
q3a9/ss = +0.0018213630x -0.9999983311y
q3a9/corner_x =       609.687
q3a9/corner_y =      -443.241

q0a10/min_fs =    0
q0a10/min_ss =  925
q0a10/max_fs =  193
q0a10/max_ss = 1109
q0a10/fs = +0.0069182217x -0.9999760985y
q0a10/ss = +0.9999760985x +0.0069182217y
q0a10/corner_x =       653.083
q0a10/corner_y =       813.670

q0a11/min_fs =  194
q0a11/min_ss =  925
q0a11/max_fs =  387
q0a11/max_ss = 1109
q0a11/fs = +0.0069182217x -0.9999760985y
q0a11/ss = +0.9999760985x +0.0069182217y
q0a11/corner_x =       654.445
q0a11/corner_y =       616.674

q1a10/min_fs =  388
q1a10/min_ss =  925
q1a10/max_fs =  581
q1a10/max_ss = 1109
q1a10/fs = +0.9999817014x +0.0060461056y
q1a10/ss = -0.0060461056x +0.9999817014y
q1a10/corner_x =      -808.912
q1a10/corner_y =       658.082

q1a11/min_fs =  582
q1a11/min_ss =  925
q1a11/max_fs =  775
q1a11/max_ss = 1109
q1a11/fs = +0.9999817014x +0.0060461056y
q1a11/ss = -0.0060461056x +0.9999817014y
q1a11/corner_x =      -611.915
q1a11/corner_y =       659.273

q2a10/min_fs =  776
q2a10/min_ss =  925
q2a10/max_fs =  969
q2a10/max_ss = 1109
q2a10/fs = -0.0027005058x +0.9999963641y
q2a10/ss = -0.9999963641x -0.0027005058y
q2a10/corner_x =      -656.888
q2a10/corner_y =      -813.482

q2a11/min_fs =  970
q2a11/min_ss =  925
q2a11/max_fs = 1163
q2a11/max_ss = 1109
q2a11/fs = -0.0027005058x +0.9999963641y
q2a11/ss = -0.9999963641x -0.0027005058y
q2a11/corner_x =      -657.420
q2a11/corner_y =      -616.482

q3a10/min_fs = 1164
q3a10/min_ss =  925
q3a10/max_fs = 1357
q3a10/max_ss = 1109
q3a10/fs = -0.9999980927x -0.0019562324y
q3a10/ss = +0.0019562324x -0.9999980927y
q3a10/corner_x =       806.906
q3a10/corner_y =      -655.522

q3a11/min_fs = 1358
q3a11/min_ss =  925
q3a11/max_fs = 1551
q3a11/max_ss = 1109
q3a11/fs = -0.9999980927x -0.0019562324y
q3a11/ss = +0.0019562324x -0.9999980927y
q3a11/corner_x =       609.906
q3a11/corner_y =      -655.907

q0a12/min_fs =    0
q0a12/min_ss = 1110
q0a12/max_fs =  193
q0a12/max_ss = 1294
q0a12/fs = -0.9999986887x -0.0016366633y
q0a12/ss = +0.0016366633x -0.9999986887y
q0a12/corner_x =       416.434
q0a12/corner_y =       791.587

q0a13/min_fs =  194
q0a13/min_ss = 1110
q0a13/max_fs =  387
q0a13/max_ss = 1294
q0a13/fs = -0.9999986887x -0.0016366633y
q0a13/ss = +0.0016366633x -0.9999986887y
q0a13/corner_x =       219.435
q0a13/corner_y =       791.265

q1a12/min_fs =  388
q1a12/min_ss = 1110
q1a12/max_fs =  581
q1a12/max_ss = 1294
q1a12/fs = +0.0016421415x -0.9999986291y
q1a12/ss = +0.9999986291x +0.0016421415y
q1a12/corner_x =      -781.411
q1a12/corner_y =       419.856

q1a13/min_fs =  582
q1a13/min_ss = 1110
q1a13/max_fs =  775
q1a13/max_ss = 1294
q1a13/fs = +0.0016421415x -0.9999986291y
q1a13/ss = +0.9999986291x +0.0016421415y
q1a13/corner_x =      -781.088
q1a13/corner_y =       222.856

q2a12/min_fs =  776
q2a12/min_ss = 1110
q2a12/max_fs =  969
q2a12/max_ss = 1294
q2a12/fs = +0.9999987483x -0.0015812991y
q2a12/ss = +0.0015812991x +0.9999987483y
q2a12/corner_x =      -423.530
q2a12/corner_y =      -793.420

q2a13/min_fs =  970
q2a13/min_ss = 1110
q2a13/max_fs = 1163
q2a13/max_ss = 1294
q2a13/fs = +0.9999987483x -0.0015812991y
q2a13/ss = +0.0015812991x +0.9999987483y
q2a13/corner_x =      -226.531
q2a13/corner_y =      -793.732

q3a12/min_fs = 1164
q3a12/min_ss = 1110
q3a12/max_fs = 1357
q3a12/max_ss = 1294
q3a12/fs = +0.0000857039x +1.0000000000y
q3a12/ss = -1.0000000000x +0.0000857039y
q3a12/corner_x =       787.680
q3a12/corner_y =      -418.886

q3a13/min_fs = 1358
q3a13/min_ss = 1110
q3a13/max_fs = 1551
q3a13/max_ss = 1294
q3a13/fs = +0.0000857039x +1.0000000000y
q3a13/ss = -1.0000000000x +0.0000857039y
q3a13/corner_x =       787.697
q3a13/corner_y =      -221.886

q0a14/min_fs =    0
q0a14/min_ss = 1295
q0a14/max_fs =  193
q0a14/max_ss = 1479
q0a14/fs = -0.9999981523x -0.0019349456y
q0a14/ss = +0.0019349456x -0.9999981523y
q0a14/corner_x =       416.424
q0a14/corner_y =       579.070

q0a15/min_fs =  194
q0a15/min_ss = 1295
q0a15/max_fs =  387
q0a15/max_ss = 1479
q0a15/fs = -0.9999981523x -0.0019349456y
q0a15/ss = +0.0019349456x -0.9999981523y
q0a15/corner_x =       219.424
q0a15/corner_y =       578.689

q1a14/min_fs =  388
q1a14/min_ss = 1295
q1a14/max_fs =  581
q1a14/max_ss = 1479
q1a14/fs = +0.0016421643x -0.9999986291y
q1a14/ss = +0.9999986291x +0.0016421643y
q1a14/corner_x =      -576.411
q1a14/corner_y =       419.856

q1a15/min_fs =  582
q1a15/min_ss = 1295
q1a15/max_fs =  775
q1a15/max_ss = 1479
q1a15/fs = +0.0016421643x -0.9999986291y
q1a15/ss = +0.9999986291x +0.0016421643y
q1a15/corner_x =      -576.088
q1a15/corner_y =       222.856

q2a14/min_fs =  776
q2a14/min_ss = 1295
q2a14/max_fs =  969
q2a14/max_ss = 1479
q2a14/fs = +0.9999768734x -0.0068050129y
q2a14/ss = +0.0068050129x +0.9999768734y
q2a14/corner_x =      -423.430
q2a14/corner_y =      -582.404

q2a15/min_fs =  970
q2a15/min_ss = 1295
q2a15/max_fs = 1163
q2a15/max_ss = 1479
q2a15/fs = +0.9999768734x -0.0068050129y
q2a15/ss = +0.0068050129x +0.9999768734y
q2a15/corner_x =      -226.435
q2a15/corner_y =      -583.745

q3a14/min_fs = 1164
q3a14/min_ss = 1295
q3a14/max_fs = 1357
q3a14/max_ss = 1479
q3a14/fs = -0.0010058292x +0.9999995232y
q3a14/ss = -0.9999995232x -0.0010058292y
q3a14/corner_x =       575.673
q3a14/corner_y =      -418.865

q3a15/min_fs = 1358
q3a15/min_ss = 1295
q3a15/max_fs = 1551
q3a15/max_ss = 1479
q3a15/fs = -0.0010058292x +0.9999995232y
q3a15/ss = -0.9999995232x -0.0010058292y
q3a15/corner_x =       575.475
q3a15/corner_y =      -221.866
