using CrystFEL
using Plots
using LinearAlgebra

dtempl = loaddatatemplate("electrons.geom")

w = 256
h = 256
xcenters = Matrix{Float32}(undef, w, h)
ycenters = Matrix{Float32}(undef, w, h)
brightfield = Matrix{Float32}(undef, w, h)
darkfield = Matrix{Float32}(undef, w, h)

function centroid(a::AbstractArray)
    total = 0.0
    moments = zeros(Float64, ndims(a))
    for p in eachindex(a)
        moments = moments .+ collect(Tuple(p))*a[p]
        total += a[p]
    end
    return moments ./ total
end

for fnstr in eachline("files.lst")

    filename = split(fnstr)[1]
    ev = split(fnstr)[2]
    image = try
        Image(dtempl, filename, ev)
    catch e
        println("Failed to load image")
        println(((n-1) ÷ w)+1, "   ", ((n-1) % w)+1)
        continue
    end

    @GC.preserve image begin
        region = view(image.dp[1], 235:290, 200:255)
        cen = centroid(region)
        x = Int(image.scan_coords[1])+1
        y = Int(image.scan_coords[2])+1
        xcenters[x,y],ycenters[x,y] = Tuple(cen).+(235,200)
        brightfield[x,y] = sum(region)
        darkfield[x,y] = sum(image.dp[1]) - brightfield[x,y]
    end
end

function fitplane(plane,flags)

    n = 1
    A = Matrix{Float32}(undef, 3, length(plane))
    v = Vector{Float32}(undef, 0)

    for x in 1:size(plane)[1]
        for y in 1:size(plane)[2]
            if flags[x,y] > 400000
                A[1,n] = x
                A[2,n] = y
                A[3,n] = 1
                push!(v, plane[x,y])
                n += 1
            end
        end
    end

    B = Matrix{Float32}(undef, 3, n-1)
    for y in 1:n-1
        B[:,y] = A[:,y]
    end

    return transpose(pinv(B))*v
end


x = fitplane(xcenters, darkfield)
y = fitplane(ycenters, darkfield)

println("panel/corner_x = ",-x[3])
println("panel/corner_y = ",-y[3])
println("scanv_beamxy = ", -x[1], ",", -x[2], ",", -y[1], ",", -y[2])

p1 = heatmap(xcenters, title="Beam center x")
p2 = heatmap(ycenters, title="Beam center y")
p3 = heatmap(brightfield, title="Bright field image")
p4 = heatmap(darkfield, title="Dark field image")
plot(p1,p2,p3,p4, layout=(2,2))
