module Images
export Image

guardian = []

function protect(guardian, obj)
    push!(guardian, obj)
    obj
end

function unprotect(guardian, obj)
    let pos = findfirst(==(obj), guardian)
        if pos !== nothing
            deleteat!(guardian, pos)
        end
    end
end


mutable struct Image
end


"""
    Image()

Create a CrystFEL image structure
"""
function Image(panels)


end

end  # of module
