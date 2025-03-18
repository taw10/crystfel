using HDF5

h5open("header_read.h5", "w") do f

    write(f, "bigint_single", Int64(1234567890123456789))
    write(f, "data_single", zeros(Int8, 1024, 1024))

    v = Vector{Int64}(undef, 3)
    v[1] = 1234567890123456785
    v[2] = 1234567890123456789
    v[3] = 1234567890123456782
    write(f, "bigint_arr", v)
    write(f, "data_arr", zeros(Int8, 1024, 1024, 3))

end
