@testset "tiling" begin

    
end

begin
    sys = AssemblySystem([1 1 2 3; 2 4 3 2; 3 4 1 2], UnitSquare)
    strs = polygen(sys; maxsize=3)
end
begin
    metasys = AssemblySystem([strs[8]])
    ps = Roly.species(metasys, 1)
    metastrs = polygen(metasys; maxsize=3)
end

