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

n = 1
A = Matrix{Float32}(undef, 3, w*h)
xv = Vector{Float32}(undef, 0)
yv = Vector{Float32}(undef, 0)

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

        if darkfield[x,y] > 400000
            A[1,n] = image.scan_coords[1]
            A[2,n] = image.scan_coords[2]
            A[3,n] = 1
            push!(xv, xcenters[x,y])
            push!(yv, ycenters[x,y])
            global n += 1
        end
    end
end

B = Matrix{Float32}(undef, 3, n-1)
for y in 1:n-1
    B[:,y] = A[:,y]
end

xans = transpose(pinv(B))*xv
yans = transpose(pinv(B))*yv

println("panel/corner_x = ",-xans[3]+0.5)
println("panel/corner_y = ",-yans[3]+0.5)
println("scanv_beamxy = ", -xans[1], ",", -xans[2], ",", -yans[1], ",", -yans[2])

p1 = heatmap(xcenters, title="Beam center x")
p2 = heatmap(ycenters, title="Beam center y")
p3 = heatmap(brightfield, title="Bright field image")
p4 = heatmap(darkfield, title="Dark field image")
plot(p1,p2,p3,p4, layout=(2,2))
