using CrystFEL

let st = Stream("input.stream", "r")
    n = 0
    for (cr,reflections) in allcrystals(st)
        n += 1
    end
    println("Read ", n, " chunks")
end
