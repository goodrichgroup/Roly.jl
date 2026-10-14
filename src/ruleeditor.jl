"""
    RuleEditor

A terminal-based geometric editor for constructing [`BindingRules`](@ref) by placing particles
on a lattice.

Everything here is internal to the editor except its entry point [`ruleeditor`](@ref), which
`Roly` re-exports, so `Roly.ruleeditor(species)` is the way to reach it.
This submodule was written by Claude Opus 5.
"""
module RuleEditor

using REPL
using REPL.Terminals: TTYTerminal, raw!
using StaticArrays: SVector
using Rotations: Angle2d
using ..Roly: BindingRules, ParticleSpecies, Pose, bindingsites, color, interactionmatrix, isaligned,
              istouching, nsites

export ruleeditor

const CSI = "\e["
const CLR_EOL = "\e[K"

alt_screen_on()   = print(CSI, "?1049h")
alt_screen_off()  = print(CSI, "?1049l")
hide_cursor()     = print(CSI, "?25l")
show_cursor()     = print(CSI, "?25h")
clear_screen()    = print(CSI, "2J")

# --- Palette / color depth ------------------------------------------------------------

const _WONG_MAKIE_MID = [
    (0x1A, 0x78, 0xC6), (0xE7, 0x54, 0x51), (0xFF, 0xB1, 0x2F), (0xB5, 0x71, 0xC4),
    (0x43, 0xCF, 0xE2), (0x19, 0xB6, 0x84), (0xFF, 0x88, 0x35), (0x7C, 0x85, 0xFF),
    (0xB9, 0xD6, 0x3E),
]
const _WONG_MAKIE_256 = [33, 203, 214, 141, 87, 79, 208, 105, 149]
const _WONG_BASIC     = [34, 31, 33, 35, 36, 32, 33, 34, 33]

function _color_depth()
    ct = get(ENV, "COLORTERM", "")
    (ct == "truecolor" || ct == "24bit") && return :truecolor
    occursin("256", get(ENV, "TERM", "")) && return :color256
    return :basic
end

function species_style(i::Int, depth::Symbol)
    if depth == :truecolor
        r, g, b = _WONG_MAKIE_MID[mod1(i, length(_WONG_MAKIE_MID))]
        return string(CSI, "38;2;", r, ";", g, ";", b, "m")
    elseif depth == :color256
        return string(CSI, "38;5;", _WONG_MAKIE_256[mod1(i, length(_WONG_MAKIE_256))], "m")
    else
        return string(CSI, _WONG_BASIC[mod1(i, length(_WONG_BASIC))], "m")
    end
end

# Count visible codepoints in a string containing CSI escapes (which we strip).
_visible_length(s::AbstractString) = length(replace(s, r"\e\[[0-9;]*[a-zA-Z]" => ""))
_pad_to(s::AbstractString, width::Int) = string(s, ' '^max(0, width - _visible_length(s)))

# --- Keyboard --------------------------------------------------------------------------

struct Key
    kind::Symbol
    value::Any
end

# POSIX poll+read on stdin (fd 0). This bypasses Julia's libuv layer so there is no
# lingering background read to clean up on exit -- after the editor returns, Julia's REPL
# reclaims stdin normally.
mutable struct _PollFd
    fd::Cint
    events::Cshort
    revents::Cshort
end
const _POLLIN = Cshort(1)

function _poll_stdin(timeout_ms::Integer)
    pfd = _PollFd(Cint(0), _POLLIN, Cshort(0))
    ret = ccall(:poll, Cint, (Ref{_PollFd}, Culong, Cint), pfd, Culong(1), Cint(timeout_ms))
    return ret > 0 && (pfd.revents & _POLLIN) != 0
end

function _read_byte()
    buf = Ref{UInt8}(0)
    n = ccall(:read, Cssize_t, (Cint, Ref{UInt8}, Csize_t), Cint(0), buf, Csize_t(1))
    n > 0 || return nothing
    return Char(buf[])
end

function readkey()
    c = _read_byte()
    c === nothing && return Key(:esc, nothing)
    if c == '\e'
        # Look for an escape sequence with a short timeout.
        if _poll_stdin(2)
            next = _read_byte()
            if next == '['
                if _poll_stdin(2)
                    code = _read_byte()
                    code == 'A' && return Key(:arrow, :up)
                    code == 'B' && return Key(:arrow, :down)
                    code == 'C' && return Key(:arrow, :right)
                    code == 'D' && return Key(:arrow, :left)
                    return Key(:other, something(code, '\0'))
                end
            end
            return Key(:other, something(next, '\0'))
        end
        return Key(:esc, nothing)
    end
    (c == '\r' || c == '\n') && return Key(:enter, nothing)
    c == '\x7f' && return Key(:backspace, nothing)
    return Key(:char, c)
end

# --- Braille rendering ----------------------------------------------------------------
# Each cell of a braille sprite is 2 dot-columns wide x 4 dot-rows tall. This gives
# sub-cell resolution for drawing outlines of shapes that don't align to a char grid.

# Bit index for each (dot-row, dot-col) in a braille character:
#   Layout:  row 1 col 1 = 0, col 2 = 3
#            row 2 col 1 = 1, col 2 = 4
#            row 3 col 1 = 2, col 2 = 5
#            row 4 col 1 = 6, col 2 = 7
const _BRAILLE_BITS = ((0, 3), (1, 4), (2, 5), (6, 7))

function _braille_char(dots::AbstractMatrix{Bool})
    v = 0
    for r in 1:4
        dots[r, 1] && (v |= 1 << _BRAILLE_BITS[r][1])
        dots[r, 2] && (v |= 1 << _BRAILLE_BITS[r][2])
    end
    return Char(0x2800 + v)
end

function _braille_line!(dots, r0, c0, r1, c1)
    # Bresenham. Coordinates are 0-based; we shift to Julia 1-based on write.
    dr = abs(r1 - r0); dc = abs(c1 - c0)
    sr = r0 < r1 ? 1 : -1
    sc = c0 < c1 ? 1 : -1
    err = dr - dc
    r, c = r0, c0
    while true
        rr, cc = r + 1, c + 1
        1 <= rr <= size(dots, 1) && 1 <= cc <= size(dots, 2) && (dots[rr, cc] = true)
        r == r1 && c == c1 && break
        e2 = 2 * err
        if e2 > -dc; err -= dc; r += sr; end
        if e2 <  dr; err += dr; c += sc; end
    end
end

function _braille_sprite(w::Int, h::Int, segments)
    dots = falses(h * 4, w * 2)
    for ((r0, c0), (r1, c1)) in segments
        _braille_line!(dots, r0, c0, r1, c1)
    end
    rows = String[]
    for row in 1:h
        chars = Char[]
        for col in 1:w
            sub = @view dots[(row-1)*4+1:row*4, (col-1)*2+1:col*2]
            push!(chars, _braille_char(sub))
        end
        push!(rows, String(chars))
    end
    return rows
end

function _overlay!(sprite::Vector{String}, r::Int, c::Int, chr::Char)
    chars = collect(sprite[r])
    if 1 <= c <= length(chars)
        chars[c] = chr
    end
    sprite[r] = String(chars)
    return sprite
end

# --- Lattice trait ---------------------------------------------------------------------
# Each lattice defines: cell size in screen chars/rows, per-cell world pose,
# and per-cell sprite (a Vector{String}, one row per line of the tile).

abstract type Lattice end
struct SquareLat   <: Lattice end
struct TriangleLat <: Lattice end
struct HexLat      <: Lattice end

lattice_for(n::Integer) = n == 3 ? TriangleLat() :
                          n == 4 ? SquareLat()   :
                          n == 6 ? HexLat()      :
                          throw(ArgumentError("editor only supports 3-, 4-, 6-gons"))

# Cell footprint on the terminal: (cols, rows). Heights are odd so single-arrow overlays
# sit exactly at the vertical center.
cell_size(::SquareLat)   = (7, 5)
cell_size(::TriangleLat) = (7, 4)
cell_size(::HexLat)      = (9, 5)

n_rotations(::SquareLat)   = 4
n_rotations(::TriangleLat) = 3
n_rotations(::HexLat)      = 6

# Top-left screen position (in canvas-local coords, 1-based).
function screen_pos(lat::SquareLat, r, c)
    cw, ch = cell_size(lat)
    return ((r - 1) * ch + 1, (c - 1) * cw + 1)
end
function screen_pos(lat::TriangleLat, r, c)
    cw, ch = cell_size(lat)
    # Triangles already alternate up/down within a row via (r+c) parity, so no
    # row offset is needed to tessellate; adjacent rows have flipped parity.
    return ((r - 1) * ch + 1, (c - 1) * cw + 1)
end
function screen_pos(lat::HexLat, r, c)
    cw, ch = cell_size(lat)
    row_off = isodd(r) ? cw ÷ 2 : 0
    return ((r - 1) * ch + 1, (c - 1) * cw + 1 + row_off)
end

# World-space pose (for bond inference; independent of screen layout).
function world_pose(::SquareLat, r, c, k, F=Float64)
    x = SVector{2,F}(F(c), F(-r))
    # Clockwise rotation: for k = 0..3 site 1 lands at south, west, north, east.
    ψ = Angle2d{F}(-F(π / 2) * k)
    return Pose(x, ψ)
end
function world_pose(::TriangleLat, r, c, k, F=Float64)
    h = sqrt(F(3)) / 2
    up = iseven(r + c)
    y_off = up ? -h / 6 : h / 6
    x = SVector{2,F}(F(c) / 2, -F(r) * h + y_off)
    base = up ? F(0) : F(π)
    ψ = Angle2d{F}(base - F(2π / 3) * k)
    return Pose(x, ψ)
end
function world_pose(::HexLat, r, c, k, F=Float64)
    h = sqrt(F(3))
    x = SVector{2,F}(F(c) * h + (isodd(r) ? h / 2 : F(0)), -F(r) * F(3 / 2))
    ψ = Angle2d{F}(F(π / 6) - F(π / 3) * k)
    return Pose(x, ψ)
end

# --- Sprites (Vector{String}, one entry per screen row) ---------------------------------

const _SQ_ARROWS      = ('↑', '→', '↓', '←')
const _TR_UP_ARROWS   = ('↑', '↘', '↙')
const _TR_DOWN_ARROWS = ('↓', '↖', '↗')
const _HX_ARROWS      = ('↑', '↗', '↘', '↓', '↙', '↖')

function tile_sprite(::SquareLat, k, _r, _c)
    arrow = _SQ_ARROWS[mod1(k + 1, 4)]
    return ["┌─────┐", "│     │", string("│  ", arrow, "  │"), "│     │", "└─────┘"]
end
empty_sprite(::SquareLat, _r, _c) = ["┌─────┐", "│     │", "│     │", "│     │", "└─────┘"]

# --- Triangle sprites (braille, 7×4) ---
# 14 dots wide × 16 dots tall. Aspect ratio 7:8 in real units, near-perfect
# equilateral in a 1:2 char aspect terminal.
const _TRI_UP_SEGMENTS = (
    ((0, 7), (15, 0)),   # left slant
    ((0, 7), (15, 13)),  # right slant
    ((15, 0), (15, 13)), # base
)
const _TRI_DOWN_SEGMENTS = (
    ((0, 0), (0, 13)),   # top edge
    ((0, 0), (15, 6)),   # left slant
    ((0, 13), (15, 7)),  # right slant
)

function _tri_up_sprite()
    return _braille_sprite(7, 4, _TRI_UP_SEGMENTS)
end
function _tri_down_sprite()
    return _braille_sprite(7, 4, _TRI_DOWN_SEGMENTS)
end

function tile_sprite(::TriangleLat, k, r, c)
    down = isodd(r + c)
    arrow = down ? _TR_DOWN_ARROWS[mod1(k + 1, 3)] : _TR_UP_ARROWS[mod1(k + 1, 3)]
    sprite = down ? _tri_down_sprite() : _tri_up_sprite()
    _overlay!(sprite, down ? 2 : 3, 4, arrow)
    return sprite
end
function empty_sprite(::TriangleLat, r, c)
    down = isodd(r + c)
    return down ? _tri_down_sprite() : _tri_up_sprite()
end

# --- Hexagon sprites (braille, 9×5) ---
# 6 vertices of a pointy-top hex in an 18×20 dot grid.
const _HEX_SEGMENTS = let N=(0,9), NE=(5,17), SE=(14,17), S=(19,9), SW=(14,0), NW=(5,0)
    (
        (N,  NE),   # top-right slant
        (NE, SE),   # right vertical
        (SE, S),    # bottom-right slant
        (S,  SW),   # bottom-left slant
        (SW, NW),   # left vertical
        (NW, N),    # top-left slant
    )
end

function _hex_sprite()
    return _braille_sprite(9, 5, _HEX_SEGMENTS)
end

function tile_sprite(::HexLat, k, _r, _c)
    arrow = _HX_ARROWS[mod1(k + 1, 6)]
    sprite = _hex_sprite()
    _overlay!(sprite, 3, 5, arrow)  # center of hex
    return sprite
end
empty_sprite(::HexLat, _r, _c) = _hex_sprite()

# --- Editor state ---------------------------------------------------------------------

mutable struct EditorState{L<:Lattice, S}
    lat::L
    base_species::S
    species::Vector{S}                                # extended lazily by digit keys
    cells::Dict{Tuple{Int,Int}, Tuple{Int,Int}}       # (r,c) => (species_idx, rot_step)
    cursor::Tuple{Int,Int}
    active_species::Int
    active_rot::Int
    grid_rows::Int
    grid_cols::Int
    depth::Symbol
    message::String
    # Rendering snapshot for diff redraws:
    prev_cells::Dict{Tuple{Int,Int}, Tuple{Int,Int}}
    prev_cursor::Tuple{Int,Int}
    prev_active_species::Int
    prev_active_rot::Int
    prev_num_species::Int
    prev_term_size::Tuple{Int,Int}
    force_full::Bool
end

function _ensure_species!(state::EditorState, n::Int)
    while length(state.species) < n
        push!(state.species, copy(state.base_species))
    end
end

# --- Rendering ------------------------------------------------------------------------

# Layout: sidebar and canvas are each drawn inside a bounding box.
# The sidebar has a fixed inner width; the canvas takes what fits, capped so no more
# than MAX_GRID_COLS particles are shown horizontally.
const SIDEBAR_INNER_WIDTH = 18
const SIDEBAR_TOTAL_WIDTH = SIDEBAR_INNER_WIDTH + 2  # + 2 box borders
const CANVAS_GAP = 1
const MAX_GRID_COLS = 20
const MAX_GRID_ROWS = 40

# Dim + dark gray; renders the empty-cell "lattice" as faint outlines.
const FAINT_STYLE = "\e[2;38;5;238m"

function _draw_box(io, r0, c0, w, h; title="")
    inner = w - 2
    top = if isempty(title)
        "┌" * "─"^inner * "┐"
    else
        pad = max(1, inner - length(title) - 2)
        "┌ " * title * " " * "─"^pad * "┐"
    end
    print(io, CSI, r0, ";", c0, "H", top)
    for r in (r0+1):(r0+h-2)
        print(io, CSI, r, ";", c0, "H", "│")
        print(io, CSI, r, ";", c0 + w - 1, "H", "│")
    end
    print(io, CSI, r0+h-1, ";", c0, "H", "└", "─"^inner, "┘")
end

function _visible_species_count(state)
    used_max = 0
    for (spidx, _) in values(state.cells)
        spidx > used_max && (used_max = spidx)
    end
    return max(state.active_species, used_max)
end

function _sidebar_lines(state)
    depth = state.depth
    lines = String[]
    push!(lines, "── Species ──")
    for i in 1:_visible_species_count(state)
        marker = i == state.active_species ? "▶" : " "
        style = species_style(i, depth)
        push!(lines, string(marker, " ", style, "■", "\e[0m", "  species ", i))
    end
    push!(lines, "")
    push!(lines, "── Keys ─────")
    push!(lines, "arrows  move")
    push!(lines, "enter   place")
    push!(lines, "space   erase")
    push!(lines, "r / R   rotate")
    push!(lines, "1-9     species")
    push!(lines, "c       clear")
    push!(lines, "q       accept")
    return lines
end

# Layout accessors given the current terminal size.
_sidebar_box(term_rows) = (r0=1, c0=1, w=SIDEBAR_TOTAL_WIDTH, h=term_rows)
function _canvas_box(state, term_rows, term_cols)
    cw, _ = cell_size(state.lat)
    reserve = state.lat isa HexLat ? cw ÷ 2 : 0
    inner_cols = state.grid_cols * cw + reserve
    c0 = SIDEBAR_TOTAL_WIDTH + CANVAS_GAP + 1
    w = min(inner_cols + 2, term_cols - c0 + 1)
    return (r0=1, c0=c0, w=w, h=term_rows)
end

function _draw_sidebar(io, state, term_rows)
    box = _sidebar_box(term_rows)
    _draw_box(io, box.r0, box.c0, box.w, box.h; title="Editor")
    lines = _sidebar_lines(state)
    inner_r0 = box.r0 + 1
    inner_c0 = box.c0 + 1
    for i in 1:(box.h - 2)
        line = i <= length(lines) ? lines[i] : ""
        padded = _pad_to(line, SIDEBAR_INNER_WIDTH)
        print(io, CSI, inner_r0 + i - 1, ";", inner_c0, "H", padded)
    end
end


function _draw_cell(io, state, r, c, canvas_inner_r, canvas_inner_c)
    row0, col0 = screen_pos(state.lat, r, c)
    is_cursor = (r, c) == state.cursor

    sprite = if haskey(state.cells, (r, c))
        sp_idx, rot = state.cells[(r, c)]
        rows = tile_sprite(state.lat, rot, r, c)
        style = species_style(sp_idx, state.depth)
        prefix = is_cursor ? "\e[7m" * style : style
        [prefix * ln * "\e[0m" for ln in rows]
    elseif is_cursor
        rows = tile_sprite(state.lat, state.active_rot, r, c)
        style = species_style(state.active_species, state.depth)
        [string("\e[2;7m", style, ln, "\e[0m") for ln in rows]
    else
        [FAINT_STYLE * ln * "\e[0m" for ln in empty_sprite(state.lat, r, c)]
    end

    for (i, line) in enumerate(sprite)
        print(io, CSI, canvas_inner_r + row0 + i - 2, ";", canvas_inner_c + col0 - 1, "H", line)
    end
end

function _draw_full_canvas(io, state, term_rows, term_cols)
    box = _canvas_box(state, term_rows, term_cols)
    _draw_box(io, box.r0, box.c0, box.w, box.h; title="Particles")
    inner_r = box.r0 + 1
    inner_c = box.c0 + 1
    # Clear canvas interior
    inner_w = box.w - 2
    for r in inner_r:(box.r0 + box.h - 2)
        print(io, CSI, r, ";", inner_c, "H", " "^inner_w)
    end
    for r in 1:state.grid_rows, c in 1:state.grid_cols
        _draw_cell(io, state, r, c, inner_r, inner_c)
    end
end

function _draw_diff_canvas(io, state, term_rows, term_cols)
    box = _canvas_box(state, term_rows, term_cols)
    inner_r = box.r0 + 1
    inner_c = box.c0 + 1
    dirty = Set{Tuple{Int,Int}}()
    push!(dirty, state.cursor, state.prev_cursor)
    for k in keys(state.cells)
        get(state.prev_cells, k, nothing) == state.cells[k] || push!(dirty, k)
    end
    for k in keys(state.prev_cells)
        haskey(state.cells, k) || push!(dirty, k)
    end
    for (r, c) in dirty
        1 <= r <= state.grid_rows && 1 <= c <= state.grid_cols || continue
        _draw_cell(io, state, r, c, inner_r, inner_c)
    end
end

function _redraw(state)
    term_rows, term_cols = displaysize(stdout)
    resized = (term_rows, term_cols) != state.prev_term_size
    _recompute_grid!(state, term_rows, term_cols)
    io = IOBuffer(sizehint=8192)

    full = state.force_full || resized
    sidebar_dirty = full ||
                    state.active_species != state.prev_active_species ||
                    state.active_rot != state.prev_active_rot ||
                    _visible_species_count(state) != state.prev_num_species

    full && print(io, CSI, "2J")
    sidebar_dirty && _draw_sidebar(io, state, term_rows)

    if full
        _draw_full_canvas(io, state, term_rows, term_cols)
    else
        _draw_diff_canvas(io, state, term_rows, term_cols)
    end
    print(io, CSI, term_rows, ";1H")  # park cursor
    write(stdout, take!(io))

    state.prev_cells = copy(state.cells)
    state.prev_cursor = state.cursor
    state.prev_active_species = state.active_species
    state.prev_active_rot = state.active_rot
    state.prev_num_species = _visible_species_count(state)
    state.prev_term_size = (term_rows, term_cols)
    state.force_full = false
end

function _recompute_grid!(state, term_rows, term_cols)
    # Canvas inner region sits inside its bounding box, to the right of the sidebar box.
    canvas_c0 = SIDEBAR_TOTAL_WIDTH + CANVAS_GAP + 1
    canvas_inner_cols = term_cols - canvas_c0 - 1  # -1 for right border
    canvas_inner_rows = term_rows - 2  # -2 top/bottom borders
    cw, ch = cell_size(state.lat)
    reserve = state.lat isa HexLat ? cw ÷ 2 : 0
    fit_rows = max(1, canvas_inner_rows ÷ ch)
    fit_cols = max(1, (canvas_inner_cols - reserve) ÷ cw)
    new_rows = min(fit_rows, MAX_GRID_ROWS)
    new_cols = min(fit_cols, MAX_GRID_COLS)
    if (new_rows, new_cols) != (state.grid_rows, state.grid_cols)
        state.grid_rows = new_rows
        state.grid_cols = new_cols
        state.force_full = true
    end
    r, c = state.cursor
    state.cursor = (clamp(r, 1, state.grid_rows), clamp(c, 1, state.grid_cols))
end

# --- Bond inference ------------------------------------------------------------------

function _collect_sites(state)
    lat = state.lat
    F = Float64
    tiles = [(rc, spidx, rot) for (rc, (spidx, rot)) in state.cells]
    sites = map(tiles) do (rc, spidx, rot)
        r, c = rc
        pose = world_pose(lat, r, c, rot, F)
        sp = state.species[spidx]
        [pose * s for s in bindingsites(sp)]
    end
    return tiles, sites
end

function inferred_bonds(state::EditorState)
    tiles, sites = _collect_sites(state)
    bondpairs = Set{Tuple{Int,Int}}()
    for i in eachindex(tiles), j in eachindex(tiles)
        j > i || continue
        for s1 in sites[i], s2 in sites[j]
            if istouching(s1, s2) && isaligned(s1, s2)
                a, b = minmax(color(s1), color(s2))
                push!(bondpairs, (a, b))
            end
        end
    end
    return sort(collect(bondpairs))
end

function bonds_matrix(state::EditorState)
    tiles, sites = _collect_sites(state)
    seen = Set{NTuple{4,Int}}()
    rows = NTuple{4,Int}[]
    for i in eachindex(tiles), j in eachindex(tiles)
        j > i || continue
        (_, sp_i, _) = tiles[i]
        (_, sp_j, _) = tiles[j]
        for (idx1, s1) in enumerate(sites[i]), (idx2, s2) in enumerate(sites[j])
            if istouching(s1, s2) && isaligned(s1, s2)
                a = (sp_i, idx1, sp_j, idx2)
                b = (sp_j, idx2, sp_i, idx1)
                if a ∉ seen && b ∉ seen
                    push!(rows, a); push!(seen, a); push!(seen, b)
                end
            end
        end
    end
    return isempty(rows) ? zeros(Int, 0, 4) : reduce(vcat, [collect(r)' for r in rows])
end

# --- Entry point ---------------------------------------------------------------------

"""
    ruleeditor(species::ParticleSpecies; output=:rules)

Open a terminal-based geometric editor for constructing binding rules by placing particles
on a lattice. `species` must be a regular polygon (triangle, square, or hexagon); the
lattice is inferred from it.

The editor starts with one "instance" of this species. Press digits `1`-`9` to add more
instances (of the same shape but distinct species colors) and switch between them; e.g.
pressing `3` for the first time creates species 3 and makes it active. Any species with
at least one placed particle contributes to the returned output.

Controls: arrow keys move the cursor, Enter places a particle (only on empty cells),
space erases, `r`/`R` rotate (the particle under the cursor if any, else the active
rotation), digits `1`-`9` pick the active species, `c` clears, `q` accepts. Bonds are
inferred from touching site pairs across all placed particles.

`output` controls the return value:
- `:rules` (default): a `BindingRules` object, ready for `polyenum`/`polygen`.
- `:bonds`: the `nx4` integer bonds matrix `[species1 site1 species2 site2; ...]`,
  suitable for pasting back into code as `BindingRules(bonds, species)`.
- `:matrix`: the `Symmetric{Bool}` color-indexed interaction matrix, suitable for
  pasting back as `BindingRules(intmat, species)`.

If nothing is placed, returns `nothing` regardless of `output`.
"""
function ruleeditor(species::ParticleSpecies; output::Symbol=:rules)
    output in (:rules, :bonds, :matrix) ||
        throw(ArgumentError("output must be :rules, :bonds, or :matrix"))
    lat = lattice_for(nsites(species))
    depth = _color_depth()
    state = EditorState(lat, species, [copy(species)],
                        Dict{Tuple{Int,Int},Tuple{Int,Int}}(),
                        (1, 1), 1, 0, 1, 1, depth, "",
                        Dict{Tuple{Int,Int},Tuple{Int,Int}}(),
                        (0, 0), 0, -1, 0, (0, 0), true)

    term = TTYTerminal(get(ENV, "TERM", "dumb"), stdin, stdout, stderr)
    raw!(term, true); hide_cursor(); alt_screen_on()

    try
        _redraw(state)
        last_size = displaysize(stdout)
        while true
            # Wait for input, but wake every 30ms to check for a terminal resize.
            while !_poll_stdin(30)
                current = displaysize(stdout)
                if current != last_size
                    last_size = current
                    state.force_full = true
                    _redraw(state)
                end
            end
            # Drain all buffered keys so held keys collapse into one redraw.
            should_quit = false
            while true
                if _handle_key!(state, readkey())
                    should_quit = true
                    break
                end
                _poll_stdin(0) || break
            end
            should_quit && break
            _redraw(state)
        end
    finally
        print(CSI, "0m")
        alt_screen_off(); show_cursor(); raw!(term, false)
    end

    # If nothing was placed, there's nothing to infer, return nothing.
    isempty(state.cells) && return nothing

    used = Set(sp_idx for (sp_idx, _) in values(state.cells))
    sorted_used = sort!(collect(used))
    remap = Dict(sp => i for (i, sp) in enumerate(sorted_used))
    bonds = bonds_matrix(state)
    for row in axes(bonds, 1)
        bonds[row, 1] = remap[bonds[row, 1]]
        bonds[row, 3] = remap[bonds[row, 3]]
    end
    kept_species = [state.species[sp] for sp in sorted_used]
    output == :bonds && return bonds
    rules = length(kept_species) == 1 ?
            BindingRules(bonds, kept_species[1]) :
            BindingRules(bonds, kept_species)
    output == :matrix && return interactionmatrix(rules)
    return rules
end

function _handle_key!(state::EditorState, key::Key)
    r, c = state.cursor
    state.message = ""
    if key.kind == :char && key.value == 'q'
        return true
    elseif key.kind == :arrow
        key.value == :up    && (state.cursor = (max(1, r - 1), c))
        key.value == :down  && (state.cursor = (min(state.grid_rows, r + 1), c))
        key.value == :left  && (state.cursor = (r, max(1, c - 1)))
        key.value == :right && (state.cursor = (r, min(state.grid_cols, c + 1)))
    elseif key.kind == :enter
        haskey(state.cells, (r, c)) || (state.cells[(r, c)] = (state.active_species, state.active_rot))
    elseif key.kind == :char && key.value == ' '
        delete!(state.cells, (r, c))
    elseif key.kind == :char && (key.value == 'r' || key.value == 'R')
        step = key.value == 'r' ? 1 : -1
        n = n_rotations(state.lat)
        if haskey(state.cells, (r, c))
            sp_idx, rot = state.cells[(r, c)]
            new_rot = mod(rot + step, n)
            state.cells[(r, c)] = (sp_idx, new_rot)
            state.active_rot = new_rot
        else
            state.active_rot = mod(state.active_rot + step, n)
        end
    elseif key.kind == :char && key.value == 'c'
        empty!(state.cells)
    elseif key.kind == :char && isdigit(key.value)
        d = parse(Int, key.value)
        if 1 <= d <= 9
            _ensure_species!(state, d)
            state.active_species = d
        end
    end
    return false
end

end # module
