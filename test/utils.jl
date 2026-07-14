using Roly: is_cutset, blockdiag!

@testset "utils" begin
    # is_cutset
    g = NautyGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)

    @test is_cutset(g, [2]) == true
    @test is_cutset(g, [1]) == false

    add_edge!(g, 1, 3)

    @test is_cutset(g, [2]) == false
    @test is_cutset(g, [1, 2]) == false
    @test is_cutset(g, [1, 2, 3]) == false

    # blockdiag!
    h = NautyGraph(; vertex_labels=[1,2,3,4])
    add_edge!(h, 1, 3)
    add_edge!(h, 2, 4)

    g_old = copy(g)
    b = blockdiag!(g, h)

    @test b[1:3] == g_old
    @test b[4:7] == h
    @test labels(b) == vcat(labels(g_old), labels(h))
    for i in 1:3, j in 4:7
        @test has_edge(g, i, j) == has_edge(g, j, i) == false
    end
end