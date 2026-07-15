@testset "tiling" begin
    # Polymino square: all 4 sites same color; opposite faces bond (1↔3, 2↔4).
    sq_sys = BindingRules([1 1 1 3; 1 2 1 4], PolygonParticleSpecies(4; labels=[1, 1, 1, 1]))
    sq = Polyform(sq_sys, 1)

    # Polyiamond triangle: all 3 sites same color; all faces bond.
    tri_sys = BindingRules([1 1 1 1; 1 2 1 2; 1 3 1 3], PolygonParticleSpecies(3; labels=[1, 1, 1]))
    tri = Polyform(tri_sys, 1)

    # Square with only 1 active site → n_opensites = 1, cannot tile.
    sys_1bond = BindingRules([1 1 1 1], UnitSquare)
    sq_1 = Polyform(sys_1bond, 1)

    # LatticeIter
    iter = LatticeIter(sq)
    @test Base.IteratorSize(typeof(iter)) == Base.SizeUnknown()
    @test Base.eltype(typeof(iter)) ==
            Tuple{Vector{SVector{2,Float64}},Vector{Vector{NTuple{2,UnitRange{Int}}}}}

    io = IOBuffer()
    show(io, iter)
    @test contains(String(take!(io)), "LatticeIter")

    @test isnothing(iterate(LatticeIter(sq_1)))

    tilings = collect(iter)
    @test !isempty(tilings)
    for (lvs, contacts) in tilings
        @test length(lvs) == length(contacts)
    end

    # Square tiling
    result = iterate(LatticeIter(sq))
    @test !isnothing(result)
    (latvecs, contacts), _ = result

    # Square tiles 2D: 2 orthogonal unit vectors.
    @test length(latvecs) == 2
    @test abs(latvecs[1]' * latvecs[2]) < 1e-10
    @test all(v -> isapprox(norm(v), 1.0; atol=1e-10), latvecs)

    # One contact group per lattice direction; 2 inter-cell bonds total.
    @test length(contacts) == 2
    @test sum(length, contacts) == 2
    for dir_contacts in contacts, c in dir_contacts
        @test c isa NTuple{2,UnitRange{Int}}
    end


    # Triangle tiling
    result = iterate(LatticeIter(tri))
    @test !isnothing(result)
    (latvecs, _), _ = result
    @test length(latvecs) == 2
    @test all(v -> isapprox(norm(v), 1.0; atol=1e-10), latvecs)

    # # istranslationtile
    @test istranslationtile(sq)
    @test istranslationtile(tri)
    @test !istranslationtile(sq_1)

    # # cantile
    @test cantile(sq; maxtilesize=1)
    @test cantile(tri; maxtilesize=1)
end
