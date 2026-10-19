module Roly

using LinearAlgebra, StaticArrays, SparseArrays, Rotations, Statistics, Random
using Base.Iterators, DataStructures
using Graphs, NautyGraphs
using ReverseSearch

# Geometry primitives
export Pose, dimension, translate, translate!
public numtype, posetype, sitetype, particletype
export Rotation, Angle2d, RotXYZ, RotMatrix3, rotation_angle, rotation_axis, SVector

# Binding sites
export BindingSite, color

# Bodies and their symmetry groups
export Polyhedron
export Tetrahedron, Cube, Octahedron, Dodecahedron, Icosahedron, Pyramid, Prism, Antiprism
export RotationGroup, Cyclic, Dihedral, Tetrahedral, Octahedral, Icosahedral, grouporder

# Particle species
export ParticleSpecies
public SpeciesAndPose
export nsites, bindingsite, bindingsites, graphrep, symmetrynumber

# Assembly system
export BindingRules, SpeciesSiteLoc, interactionmatrix
export ncolors, nspecies, nbonds, isinert, species
public bonded_colors, bonded_sites, bonded_species, sitesofcolor, speciesofcolor

# Polyforms
export AbstractPolyform, Polyform, ParticleSiteLoc, nparticles, bindingrules, composition, bonds
export exposedsites, opensites, exposedsitelocs, opensitelocs

# The graph a polyform carries, and the vertices it is written in, are an implementation detail:
# every one of these takes or yields bare graph vertices, whose numbering (canonical or original)
# a caller has to get right. Nothing exported does, so that hazard stops at the boundary.
public bondindex, interior_edges, exterior_edges, subpolyform
public canonbindingsite, canonbindingsites, rotationcenter

export bondtypes

# Enumeration
export ACCEPT, REJECT, BREAK
export EnumerationStatus, Complete, MaxSizeReached, MaxPolysReached, Interrupted
export polyenum, polygen, polycount, PolyformCount, isbounded

# Species
export PolygonParticleSpecies, UnitNgon, UnitTriangle, UnitSquare, UnitHexagon
export SymmetricUnitTriangle, SymmetricUnitSquare, SymmetricUnitHexagon
export PolyhedronParticleSpecies
export UnitTetrahedron, UnitCube, UnitOctahedron, UnitDodecahedron, UnitIcosahedron
export SymmetricUnitTetrahedron, SymmetricUnitCube, SymmetricUnitOctahedron
export SymmetricUnitDodecahedron, SymmetricUnitIcosahedron
export UnitPyramid, UnitPrism, UnitAntiprism
export PatchyParticleSpecies, PatchyDisk, PatchySphere
export MetaParticleSpecies, polyform, inducedrules, recast

# Public, but not exported: reach for these as `Roly.faces(p)`, or import them by name.
# They are stable API, but specific enough to a body, an encoding or a species that putting
# them in every user's namespace is not worth it.

# Polyhedron geometry
public corners, faces, facevertices, nfaces, nedges
public facecentroid, facecentroids, facenormal, facenormals, edgemidpoint
public inradius, minedgelength

# Graph encodings and the symmetry they record
public dartencoding, cycleencoding
public rotationgroup, faceorbits, facesym, siteorbits, stabilizerorders, sitelabel
public permutationgroup, check_encoding

# What a bond fixes about relative orientation
public contact_pairing, standard_twist, twistfreedom, twist, Contact

# The `ParticleSpecies` interface, implemented rather than called
public isconvex, bounding_radius, could_contact, overlap
public sat_overlap, edgenormals, polyhedron

include("utils.jl")
include("pose.jl")
include("bindingsite.jl")
include("particlespecies.jl")
include("encoding.jl")
include("bindingrules.jl")
include("particle.jl")
include("polyform.jl")
include("enumeration.jl")

include("species/polygonparticlespecies.jl")
include("species/polyhedronparticlespecies.jl")
include("species/patchyparticlespecies.jl")
include("species/metaparticlespecies.jl")

export editrules, editpolyform
export render, polyformplot, polyformplot!

"""
    editrules(species::ParticleSpecies; output=:rules)
    editrules(catalog::AbstractVector; output=:rules)
    editrules(rules::BindingRules; output=:rules)

Open a terminal editor that builds binding rules geometrically, and return them once you
accept.

Given a vector of species, the editor may set any of its species to any of those shapes, which
is how a system of differently shaped particles is designed. Given a `BindingRules`, it opens on
that system's species with its bonds already set, so an existing one can be worked on further. Provided by the Tachikoma extension, so it needs `import Tachikoma`, which is preferable
to `using` because Tachikoma exports `render` and `Rect` as well.

Rules can be built two ways, and the editor keeps both on screen at once. Placing particles reads
bonds off the geometry: each step attaches a copy of a species to one of the structure's free
binding sites, choosing which free site to attach to, which site of the incoming particle meets
it, and in 3D which twist of the bond to use. The rules then come from every pair of placed
particles, so contacts a particle makes with neighbors it was not attached to count too, which
is how an arrangement reveals rules you did not set out to specify. Editing the interaction
matrix sets bonds directly, and an edit holds against whatever the geometry goes on to produce,
so the two can be used together.

The editor opens on two panes, the rules and the polyforms they enumerate, with a third for
building structures that `b` brings in. `tab` moves the focus between the panes on screen and the
arrow keys and enter act on whichever has it. The rules are always shown, since they are what the
editor produces.

The editor starts with one instance of `species`. In the rules pane `a` adds another, of the same
shape but distinct colors, and `d` drops the last; in the construction pane, which shows them as
a strip marking the active one, digits `1`-`9` pick which to place. Every instance gets a color
range of its own in the rules whether or not it is ever placed, and one carrying a bond or a
placement is not dropped. There is no limit on how many there may be: a table that then outgrows
the window says so rather than the editor refusing to go on.

`f` turns on fully addressable mode, in which every placement takes a species of its own, the
editor adding one as each particle goes down and dropping it again as each comes back up. No two
particles then share a binding site color, so each bond in the rules names one particular pair of
particles rather than a pair of species, which is what a fully addressable structure means.
Turning the mode on over a structure already built renumbers that structure rather than asking
for it to be rebuilt. Where the editor was given several shapes, `p` and `P` set the active
species to another of them.

With the rules focused, the arrow keys move a cursor over the interaction matrix, rows with up
and down and columns with left and right, and enter turns that pair of colors on or off. The two
colors are picked out in the species drawings above the matrix, so moving the cursor is how two
sites are selected and enter bonds them. The matrix and the list of bonding pairs are both shown,
the matrix saying what is possible and the list what is set; the cursor drives the matrix, being
the editable one, and the list marks whichever pair it is on. A matrix too big to sit beside the
list moves to a full-height pane of its own, so long as it leaves the other panes their minimum
width. A matrix grows as the square of the colors, so a fully addressable structure asks for most
of the screen; past that the list of bonding pairs takes over and the rules pane is titled with
the color count.

A site is drawn in a shade of its species' own hue, so its color says which particle it belongs
to, and a matrix cell is split between the colors of the two sites it joins, the row's above and
the column's below. A color is written as one character while one is enough to tell them all
apart, `1`-`9` then `a`-`z` then `A`-`Z`, and as two decimal digits past those 61; the editor
stops adding species at 99 colors rather than naming two different sites the same way. Each particle carries a small filled triangle at its center pointing away
from its first binding site, showing which way it is turned, left off wherever the particle is
drawn too small for it to read as a direction.

With the construction focused, each arrow key takes the free site lying next in that direction,
out of every site ahead rather than out of the two nearest round the perimeter: the one taken
minimizes its distance along the direction plus a weight times its offset to the side, so a
nearly aligned site wins over a nearer one off to one side. The distance along the direction
strictly increases at every step, so holding one arrow down marches across the structure and then
stops, with nothing ahead left to take, rather than oscillating between a pair each of which was
the better match from the other. `,` and `.` walk the perimeter in order, which always moves.
`r`/`R` turn the pending particle by changing which of its sites meets the cursor, `t`/`T` pick
the twist in 3D, enter attaches and backspace undoes.

Structures are built in construction windows, and `n` opens another while `w` and `W` cycle
through them. The rules are read off every window at once, so an arrangement can be kept on the
books while the next is built beside it. `c` clears the window in front, taking with it the bonds
that only that window produced and leaving the ones another window or a hand edit also accounts
for; clearing an already empty window closes it, and the last one is always kept.

Seating a particle against the chosen site can put it against several at once, which is how a
ring closes. Every contact the pending particle would make is labelled before it is committed
and the sidebar counts the ones beyond the site aimed at, so a closure shows in advance. A
placement that would overlap an existing particle is refused, and is drawn crossed out in gray
beforehand rather than looking legal until enter does nothing.

The enumeration pane draws the polyforms the rules allow. It runs on `e` rather than
automatically, since a permissive rules set can take long enough to be felt between keystrokes;
changing the rules marks the last run out of date rather than repeating it. `s`/`S` set the
largest structure to look for and `x`/`X` the number to stop after, which is what bounds the
work. Neither has a ceiling, what is worth waiting for being the user's to judge.

Each run first asks `isbounded` whether the rules admit finitely many structures at all and
`polycount` how many there are, exactly while the count fits its own budget and by subsampling
the search tree beyond it. The pane title then carries the found count, how many there are to
find when that is more, and a tag for each bound in force: `16, complete, bounded` for a listing
that holds every structure the rules admit, `15 of 16, size ≤ 4, bounded` for one the size bound
stops short, and `200 of ≈688k, size ≤ 12, unbounded` under rules that admit arbitrarily large
structures. Counts above ten thousand are abbreviated with `k`, `M` or `B`. The sidebar carries the same verdict on the rules themselves,
`bounded, 16 total` against `unbounded`. Reporting only how many were found leaves no way to tell
a finished design from a partial view of one. Structures are numbered by their position in the run; with the pane focused the arrow keys
pick one, and the box beside the grid redraws it as large as the pane allows, titled with its
number, naming every site in its own color and writing both colors where two meet, with its
composition vector underneath.

`-` and `=` zoom the construction view and `0` refits it, those being the only pane with a camera.
`q` accepts from anywhere. A line along the foot reports what the editor last did, or why it
declined to, and counts the species, bonds and placed particles.

A 3D species is drawn in isometric projection, as the edges of the faces that turn toward the
camera, each face's edges in the shade its orientation gives it so that the drawing reads as a lit
solid. Hidden edges go in two steps: within a particle, dropping the faces turned away removes
exactly the edges on its far side, the body being convex; between particles, the drawing runs
farthest first and each particle erases what its silhouette covers of everything already drawn
behind it. A species with no polyhedron behind it is drawn as its bounding sphere's silhouette,
the way a 2D species with no corners is drawn as a circle.

`[` and `]` step the camera through eight viewpoints, one per octant, shared by every pane. The
camera never moves on its own; what follows it is the cursor, the sites offered to attach to
being the ones the camera can see, so the arrow keys walk the near side of the structure and
turning the camera offers a different set. A face of a convex particle is visible from four of
the eight, so every site is reachable. A small set of labelled axes in the corner of the
construction pane and of the inspector says which way the camera is looking. Only the sites on
the near side are named, and a bond between two placed particles is left unlabelled, its two
faces meeting inside the solid where nothing of them shows.

The construction view holds its zoom while you build, so attaching a particle leaves the rest of
the structure where it was, and it zooms out only when the structure would leave the pane.
Zooming by hand holds the cursor still rather than the middle of the structure, and switches the
automatic fitting off until `0` asks for it again.

Binding sites are drawn with the color the returned rules give them, so a label names the same
row the interaction matrix does. Where two sites meet, both colors are written side by side, the
pending bond included, and the pending particle's remaining sites are named where they sit so
that turning it turns the labels with it.

`output` controls the return value:
- `:rules` (default): a `BindingRules` object, ready for `polyenum`/`polygen`.
- `:bonds`: the `nx4` integer bonds matrix `[species1 site1 species2 site2; ...]`,
  suitable for pasting back into code as `BindingRules(bonds, species)`.
- `:matrix`: the `Symmetric{Bool}` color-indexed interaction matrix, suitable for
  pasting back as `BindingRules(intmat, species)`.

If no bond was formed, returns `nothing` regardless of `output`.
"""
function editrules end

"""
    editpolyform(species::ParticleSpecies; addressable=false)
    editpolyform(catalog::AbstractVector; addressable=false)
    editpolyform(poly::Polyform; addressable=false)

Open the construction pane of [`editrules`](@ref) on its own and return the structure built in
the window in front, as a [`Polyform`](@ref) carrying the binding rules its own contacts imply.
Provided by the Tachikoma extension, so it needs `import Tachikoma`.

Given a `Polyform`, the editor opens on that structure and the system it was built under, so one
can be carried on from or taken apart. A window holds one connected structure, so a particle the
rest hangs off is not removed.

The keys are the construction keys of [`editrules`](@ref), the panes it has no use for being
left off. `addressable=true` opens in fully addressable mode, which gives every particle a
species of its own, so that the rules the structure implies name one particular pair of particles
per bond rather than a pair of species.

Returns `nothing` if nothing was built, or if the window does not hold one connected structure.
The polyform is grown through Roly's own attachment, in an order read off the contacts, so it
arrives with the graph, bonds and canonical particle ordering that every other polyform has. That
canonical ordering is also why its particles are not in the order they were placed, and why a
structure that goes round the editor comes back equal by its canonical graph rather than by
`==`, which asks for the same rules object.

The rules come from the structure, so they are exactly as permissive as its contacts require and
generally admit other structures too. Use [`polycount`](@ref) or [`polyenum`](@ref) on
`bindingrules(poly)` to see what else they allow.
"""
function editpolyform end

"""
    render(p; hidedecorations=true, kwargs...)

Draw a [`Polyform`](@ref) or [`ParticleSpecies`](@ref) into a fresh figure, picking a 2D or 3D
axis to match its dimension. Provided by the Makie extension, so it needs a backend loaded.

3D output needs a backend with a depth buffer, GLMakie or WGLMakie. CairoMakie sorts
primitives instead of depth-testing them, so 3D polyforms show artifacts where faces meet.

`bindingrules=rules` draws sites no bond can use as inert, which works for a bare species too.
"""
function render end

"""
    polyformplot(p; kwargs...)

Makie recipe behind [`render`](@ref), drawing a polyform or species into a new axis.
Provided by the Makie extension.
"""
function polyformplot end

"""
    polyformplot!(ax, p; kwargs...)

Draw a polyform or species onto an axis you already have, the mutating form of
[`polyformplot`](@ref). Provided by the Makie extension.
"""
function polyformplot! end

end
