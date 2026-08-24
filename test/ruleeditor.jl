using Roly: EditorState, SquareLat, TriangleLat, HexLat, Key,
            inferred_bonds, bonds_matrix, _handle_key!, _visible_species_count,
            _ensure_species!

# Convenience: build a fresh EditorState with just the fields tests care about.
function _mkstate(lat, sp; cells=Dict{Tuple{Int,Int},Tuple{Int,Int}}(),
                  cursor=(1, 1), active_species=1, active_rot=0,
                  grid_rows=8, grid_cols=8, species=[copy(sp)])
    return EditorState(lat, sp, species, cells, cursor, active_species, active_rot,
                       grid_rows, grid_cols, :basic, "",
                       Dict{Tuple{Int,Int},Tuple{Int,Int}}(), (0, 0), 0, -1, 0, (0, 0), true)
end

@testset "ruleeditor" begin
   # bond inference
    s = _mkstate(SquareLat(), UnitSquare;
                    cells=Dict((1,1) => (1, 0), (1,2) => (1, 0)))
    @test inferred_bonds(s) == [(2, 4)]

    s = _mkstate(SquareLat(), UnitSquare;
                    cells=Dict((1,1) => (1, 0), (2,1) => (1, 0)))
    @test inferred_bonds(s) == [(1, 3)]

    # Empty grid -> no bonds
    s = _mkstate(SquareLat(), UnitSquare)
    @test isempty(inferred_bonds(s))
    @test size(bonds_matrix(s)) == (0, 4)

    s = _mkstate(TriangleLat(), UnitTriangle;
                    cells=Dict((1,1) => (1, 0), (1,2) => (1, 0)))
    @test inferred_bonds(s) == [(3, 3)]

    s = _mkstate(TriangleLat(), UnitTriangle;
                    cells=Dict((1,1) => (1, 0), (2,1) => (1, 0)))
    @test inferred_bonds(s) == [(1, 1)]

    s = _mkstate(HexLat(), UnitHexagon; cells=Dict((1,1) => (1, 0), (1,2) => (1, 0)))
    @test inferred_bonds(s) == [(3, 6)]

    s = _mkstate(HexLat(), UnitHexagon;
                    cells=Dict((1,1) => (1, 0), (2,1) => (1, 0)))
    @test inferred_bonds(s) == [(2, 5)]

    # end-to-end:
    cells = Dict((1,1) => (1, 0), (2,1) => (1, 0), (2,2) => (1, 0))
    s = _mkstate(SquareLat(), UnitSquare; cells=cells)
    bonds = bonds_matrix(s)
    rules = BindingRules(bonds, UnitSquare)
    r = polyenum(rules; maxsize=3, maxstrs=100)
    @test r.nstructures > 0

     #_visible_species_count"
    s = _mkstate(SquareLat(), UnitSquare)
    @test _visible_species_count(s) == 1

    s.active_species = 5
    @test _visible_species_count(s) == 5

    # After switching back with no cells, contracts.
    s.active_species = 1
    @test _visible_species_count(s) == 1

    # If a cell uses species 3, count must include it.
    s.cells[(1,1)] = (3, 0)
    @test _visible_species_count(s) == 3
    s.active_species = 2
    @test _visible_species_count(s) == 3  # max(active, used_max)

    "cursor/keys"
    s = _mkstate(SquareLat(), UnitSquare; cursor=(3, 3))
    _handle_key!(s, Key(:arrow, :up));    @test s.cursor == (2, 3)
    _handle_key!(s, Key(:arrow, :down));  @test s.cursor == (3, 3)
    _handle_key!(s, Key(:arrow, :left));  @test s.cursor == (3, 2)
    _handle_key!(s, Key(:arrow, :right)); @test s.cursor == (3, 3)

    # Clamp at boundaries.
    s.cursor = (1, 1)
    _handle_key!(s, Key(:arrow, :up));    @test s.cursor == (1, 1)
    _handle_key!(s, Key(:arrow, :left));  @test s.cursor == (1, 1)
    s.cursor = (s.grid_rows, s.grid_cols)
    _handle_key!(s, Key(:arrow, :down));  @test s.cursor == (s.grid_rows, s.grid_cols)
    _handle_key!(s, Key(:arrow, :right)); @test s.cursor == (s.grid_rows, s.grid_cols)

    s = _mkstate(SquareLat(), UnitSquare; cursor=(2, 3), active_species=1, active_rot=2)
    _handle_key!(s, Key(:enter, nothing))
    @test s.cells[(2, 3)] == (1, 2)

    # Placing on an occupied cell is a no-op.
    s.active_species = 2
    _ensure_species!(s, 2)
    _handle_key!(s, Key(:enter, nothing))
    @test s.cells[(2, 3)] == (1, 2)  # unchanged

    # Space erases.
    _handle_key!(s, Key(:char, ' '))
    @test !haskey(s.cells, (2, 3))

    # Then Enter places (now that cell is empty).
    _handle_key!(s, Key(:enter, nothing))
    @test s.cells[(2, 3)] == (2, 2)  # active_species=2 now

    # Cursor on empty cell rotates active_rot only.
    s = _mkstate(SquareLat(), UnitSquare; cursor=(1, 1), active_rot=0)
    _handle_key!(s, Key(:char, 'r'))
    @test s.active_rot == 1
    _handle_key!(s, Key(:char, 'R'))
    @test s.active_rot == 0

    # Cursor on a placed tile rotates the tile AND syncs active_rot.
    s.cells[(1, 1)] = (1, 2)
    _handle_key!(s, Key(:char, 'r'))
    @test s.cells[(1, 1)] == (1, 3)
    @test s.active_rot == 3
    _handle_key!(s, Key(:char, 'R'))
    @test s.cells[(1, 1)] == (1, 2)
    @test s.active_rot == 2

    # Wraps around n_rotations (4 for square).
    s.cells[(1, 1)] = (1, 3)
    _handle_key!(s, Key(:char, 'r'))
    @test s.cells[(1, 1)] == (1, 0)

    s = _mkstate(SquareLat(), UnitSquare)
    @test length(s.species) == 1
    @test s.active_species == 1

    _handle_key!(s, Key(:char, '3'))
    @test s.active_species == 3
    @test length(s.species) == 3  # extended by _ensure_species!

    _handle_key!(s, Key(:char, '1'))
    @test s.active_species == 1
    @test length(s.species) == 3  # doesn't shrink; visibility handled separately

    # Digits outside 1-9 are ignored.
    _handle_key!(s, Key(:char, '0'))
    @test s.active_species == 1

    s = _mkstate(SquareLat(), UnitSquare;
                    cells=Dict((1,1) => (1, 0), (2,2) => (1, 0)))
    _handle_key!(s, Key(:char, 'c'))
    @test isempty(s.cells)

    # 'q' returns true (quit signal).
    @test _handle_key!(s, Key(:char, 'q')) == true

    # Other keys return false.
    @test _handle_key!(s, Key(:arrow, :up)) == false
end
