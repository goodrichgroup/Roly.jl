# Orientation and phases

A bond in Roly fixes more than which two sites touch. It fixes how the two particles are turned
relative to one another, and on a symmetric particle it may fix it in more than one way. This
page explains the model: what a binding site's frame means, how many distinct attachments a
bond admits, and why the rules that align frames across a solid's faces are what they are.

None of it has to be read to use the package. It matters when a species does not enumerate what
you expected, because the failure mode is silent — wrong counts, not errors.

## What a binding site's pose says

A site's pose carries three things, with different status:

- **position** — physical.
- **outward normal**, the pose's local x — physical.
- **twist reference**, the freedom left in the frame. For a polyhedron face this is "local z
  points at the midpoint of the face's first edge", so it is really a choice of which corner the
  face's corner list starts at.

Three integers attach to a site. Only the first is visible in the graph:

| quantity | meaning | for a square side face of a triangular prism |
|---|---|---|
| `degree` | number of darts, i.e. graph vertices under [`dartencoding`](@ref) | 4 |
| `gauge` | order of the *face's own* rotational symmetry about its normal | 4 |
| `stab` | order of the stabiliser of that site in the *particle's* symmetry group | 2 |

`stab ∣ gauge ∣ degree`. They collapse together for regular faces on highly symmetric bodies,
which is why cubes hide every distinction below: a cube's face has all three equal to 4.

[`facegauge`](@ref) computes the middle one from a face's geometry, [`sitestabilisers`](@ref)
the last one from the particle's.

## Phases

A bond fixes the partner's position and normal. What remains is a turn about the bond axis, and
the question is how many values of it give genuinely different structures.

It is tempting to say: the twist reference is only defined up to `C_gauge`, so a bond must admit
all `lcm(gauge₁, gauge₂)` turns. That is wrong. A face turn that is not a symmetry of the whole
particle carries the frame to a different edge of the *body*, and so describes a different
relative orientation to everything else on that particle. Nothing about it is redundant.

Only a rotation carrying the particle onto itself produces two frames that describe it equally
well. So, with `q` the site's [`twistfreedom`](@ref):

```
L = lcm(q₁, q₂)        q = stab  (locking)  or  gauge  (rotation-free)
R = L / stab₂
```

`L` is [`bondperiod`](@ref) and `R` is [`nphases`](@ref), the number of distinct
attachments, with `stab₂` belonging to the particle being attached. Phases come from
**symmetry**, and an unsymmetric particle has exactly one per bond — which is the ordinary
reading of an oriented binding site. The quotient by `stab₂` is an optimisation rather than a
correctness requirement; those duplicates would be merged by canonical form anyway.

`locking` is the per-site modelling choice, a keyword of
[`PolyhedronParticleSpecies`](@ref) and [`PatchySphere`](@ref). A locking site — the default —
holds its partner in the orientation its frame names. A rotation-free site declines to fix the
twist and admits every orientation its face geometrically permits. Setting a triangular prism's
square side faces rotation-free is how to get both the coplanar and the out-of-plane lattice;
leaving them locking gives the coplanar one only.

Two conditions this has to satisfy:

- `L ∣ lcm(degree₁, degree₂)`, so the phase can be recorded in the graph. This holds
  because `q ∣ gauge ∣ degree`.
- [`cycleencoding`](@ref), one vertex per site, can only express `L = 1`, so it is valid exactly
  when every twist freedom is 1. Distinct labels already force every `stab` to 1, so for locking
  sites this comes free — which is why most species get the cheap encoding.

### Choosing a phase rather than counting them

`twists` turns a face's site about its own normal, by an angle in radians. This is the bond
*phase*: which relative orientation the bond means, as opposed to how many it admits. Whole
dart steps rotate the face's corner list, so the frame and the graph vertices move together;
whatever is left over turns the frame alone.

Note that turning both faces of a bond by the same amount does not cancel — it turns the partner
by twice that, the offset appearing on both sides of the face-to-face flip, since
`Δ·Rx(-θ) = Rx(θ)·Δ`.

Turning one face of a symmetry orbit differently from its fellows splits the orbit and lowers
the symmetry number. That is what it is for, and the graph records the break: the twist is
folded into the key [`siteorbits`](@ref) groups by. A twist shared across an orbit leaves the
symmetry intact, since turns about a site's own normal commute with its stabiliser.

Selecting a single phase on a *symmetric* particle is not possible, and that is correct
rather than missing: the two attachments are related by a symmetry of the particle, so a model
that keeps the symmetry cannot prefer one. Take both with `locking=false`, or break the symmetry
with colors.

## The two face rules

Two functions rotate a solid's corner lists, and they enforce different — and unequal —
conditions.

[`Roly._canonical_faces`](@ref) reads one face's own shape: rotate to the lexicographically
least cyclic word of `(edge length, interior angle)`. That word determines a planar polygon up
to congruence and its cyclic symmetries are exactly the face's rotations, so the minimising
positions form exactly one gauge orbit, and congruent faces pick corresponding positions.

This aligns references **up to `gauge`**, which is all [`site_symmetry`](@ref) needs — it
compares frames up to a face turn, because a quarter turn about a cube face is a real symmetry
that shifts darts. It is *not* enough for bonds.

[`Roly._propagate_faces`](@ref) aligns references **up to `stab`**, which is what bonds need.
Where `gauge > stab` the leftover turns are real: two frames off by a gauge turn that is not a
stab turn describe different bonds.

The case that exposes the difference is a triangular prism with `h = a`, all three square sides
alike. Under the face rule alone the three references come out at `[0.87,-0.5,0]`, `[0,0,-1]`,
`[0,1,0]` — one of them a quarter turn from the others. Legal by gauge (4), illegal by stab (2).
Bonding face 2 to face 3 then turns the neighbour out of the plane where face 2 to face 4 does
not, and the prisms enumerate as `[1,3,6,22,73,357]` instead of the polyiamonds
`[1,2,3,6,10,22]`. With the same solid at `h = 2` the sides are rectangles, the face rule has
something to bite on, all three references come out vertical, and the count is right without
propagation.

## Why propagation is correct

Let `G` be the rotations of the solid that preserve the labelling. This is a group: composing
two label-preserving symmetries preserves labels.

**It must be `G`, not the solid's full rotation group.** Propagating under the full group relates
two faces by a rotation that need not lie in `G`, leaving their frames off by an element of the
*solid's* stabiliser rather than the species'. A cube with four sticky sides and two caps is the
counterexample: the solid's group is all 24 rotations with every face in one orbit, but the
species' group is the 8 preserving sides-vs-caps. Propagating under the 24 leaves adjacent side
faces a quarter turn apart — inside gauge, outside stab — and polyominoes come out as
`[1,3,9,39,208,1402]` instead of `[1,2,4,9,21,56]`. Since labels come from the coloring, this is
species-level information, and the propagation belongs there rather than on
[`Polyhedron`](@ref).

**Construction.** Process faces in index order; each inherits its reference from the first
earlier face in its `G`-orbit, carried over by some element of `G`.

**Lemma.** Every face `k` ends with `frameₖ = hₖ · frame_r`, where `r` is the least-indexed face
of its orbit and `hₖ ∈ G` satisfies `hₖ(r) = k`.

*Proof.* Induction on index. The orbit minimum `r` is earlier than every other member, so each
face has someone to inherit from. If `j` inherits from `i` via `g ∈ G`, then by the induction
hypothesis `frameⱼ = g · frameᵢ = g hᵢ · frame_r`, and `g hᵢ ∈ G` carries `r` to `j`. ∎

**Theorem.** For *every* `g ∈ G` — not only the ones used in the construction — and every face
`i` with `j = g(i)`, the frames `frameⱼ` and `g · frameᵢ` differ by an element of `Stab_G(j)`, a
turn about face `j`'s normal of order `stab_j`.

*Proof.* Put `t = frameⱼ · (g · frameᵢ)⁻¹`. By the lemma,

```
t = hⱼ · frame_r · frame_r⁻¹ · hᵢ⁻¹ · g⁻¹ = hⱼ hᵢ⁻¹ g⁻¹ ∈ G
```

and `t(j) = hⱼ(hᵢ⁻¹(g⁻¹(j))) = hⱼ(hᵢ⁻¹(i)) = hⱼ(r) = j`, so `t ∈ Stab_G(j)`. A proper rotation
fixing a face fixes its plane and its centroid, hence is a turn about its normal, and
`Stab_G(j)` is cyclic of order `stab_j`. ∎

This is also the best achievable: no assignment of frames can be pinned finer than the
stabiliser, since the stabiliser's elements permute the darts of a face that the particle cannot
tell apart.

**Two side conditions.**

1. The target dart always exists: `g` is a congruence of the solid, so it carries face `i`'s
   edge midpoints onto face `j`'s.
2. Propagation moves dart 1 only *within its gauge orbit*, so it changes frames only by gauge
   turns. `_canonical_faces` has already put dart 1 in the class minimising a
   congruence-invariant word, and `g` is a congruence, so the image lands in `j`'s minimising
   class. This is what makes it sound to compute [`facegauge`](@ref) and [`siteorbits`](@ref) —
   both of which compare frames up to gauge — *before* propagating, which the ordering requires,
   since the labels are needed to know `G`.

## What checks this, and what cannot

[`site_symmetry`](@ref) compares frames up to `gauge`, so it is blind to precisely the choices
propagation makes. It is still a complete backstop for its own failure mode, and the asymmetry
is worth relying on: a badly chosen twist reference can only make `site_symmetry` too **small**,
never too large. A rotation matching every site's position and normal maps every face's
supporting plane to another, hence maps the intersection of the half-spaces — the solid, *by
convexity* — onto itself, so it is a genuine symmetry, and no spurious ones can be invented.
[`Roly._check_encoding`](@ref) compares it against the graph's symmetry number and therefore
catches every under-count.

Convexity carries real weight there, and elsewhere: [`PolyhedronParticleSpecies`](@ref) reports
`isconvex` unconditionally and tests overlap by separating axes. `Polyhedron` checks it on
construction.
