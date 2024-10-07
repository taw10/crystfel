using CrystFEL

let st = Stream("/home/twhite/experiments/cxidb-193/short.stream", "r"),
    merged = mergereflections(allcrystals(st), SymOpList("mmm"))
    savereflist!(merged, "merged.hkl")
end
