adu_per_photon = 1
res = 10000
clen = 100.0 mm
photon_energy = 9000 eV

dim0 = %
dim1 = ss
dim2 = fs
data = /data/data

q0/dim0 = 0
q0/min_fs = 0
q0/min_ss = 0
q0/max_fs = 1024
q0/max_ss = 1024
q0/fs = x
q0/ss = y
q0/corner_x = -1056
q0/corner_y = 30

q1/dim0 = 1
q1/min_fs = 0
q1/min_ss = 0
q1/max_fs = 1024
q1/max_ss = 1024
q1/fs = x
q1/ss = y
q1/corner_x = 30
q1/corner_y = 30

q2/dim0 = 2
q2/min_fs = 0
q2/min_ss = 0
q2/max_fs = 1024
q2/max_ss = 1024
q2/fs = x
q2/ss = y
q2/corner_x = 30
q2/corner_y = -1054

q3/dim0 = 3
q3/min_fs = 0
q3/min_ss = 0
q3/max_fs = 1024
q3/max_ss = 1024
q3/fs = x
q3/ss = y
q3/corner_x = -1054
q3/corner_y = -1054

group_all = q0,q1,q2,q3
