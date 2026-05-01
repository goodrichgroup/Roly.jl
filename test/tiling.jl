@testset "tiling" begin

    
end

begin
    # sys = AssemblySystem([1 1 2 3; 2 4 3 2; 3 4 1 2], UnitSquare)
    sys = AssemblySystem([1 3 1 1; 1 2 1 4], UnitSquare)
    strs = polygen(sys; maxsize=3)
    ms = Vector{Int}[]
end

function tile_add(s, args...)
    m = composition(s)
    push!(ms, m)
    tbonds = tile_bonds(s)
    if !isnothing(tbonds)
        np = nspecies(assemblysystem(s))
        tile_m = copy(m)
        for b in tbonds
            tile_m[np+b] += 1
        end
    end
end

begin
    metasys = AssemblySystem([strs[8]])
    ps = Roly.species(metasys, 1)
    metastrs = polygen(metasys; maxsize=3)
end

