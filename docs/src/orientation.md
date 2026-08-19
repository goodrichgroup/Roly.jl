# Orientation and phases

A bond fixes which two sites touch, and how the two particles are turned relative to one another.
On a symmetric particle it may fix that turn in more than one way.
Each way is a *phase*.

Read this page if a species enumerates the wrong number of structures, since that failure is silent.

## What a site's pose says

A binding site's pose carries its position, its outward normal (the pose's local x axis), and its *twist reference*, the freedom left in the frame.
For a polyhedron face the twist reference is "local z points at the midpoint of the face's first edge", which is a choice of where the face's corner list starts.

## Three numbers per site

| quantity | meaning | square side face of a triangular prism |
|---|---|---|
| `degree` | graph vertices the site owns under [`dartencoding`](@ref) | 4 |
| `gauge` | order of the face's own rotational symmetry about its normal | 4 |
| `stab` | order of that site's stabiliser in the particle's symmetry group | 2 |

Each divides the next.
All three are equal for a cube face, which is why cubes hide the distinction.
[`facegauge`](@ref) computes `gauge` from the face, [`sitestabilisers`](@ref) computes `stab` from the particle.

## How many phases a bond has

A bond fixes the partner's position and normal, leaving a turn about the bond axis.
Only a rotation carrying the particle onto itself gives two frames that describe it equally well.
A face turn that is not such a rotation points the frame at a different edge of the body, so it means a different orientation.

Writing `q` for a site's [`twistfreedom`](@ref):

```
L = lcm(q₁, q₂)          q = stab if the site is locking, gauge if not
R = L / stab₂
```

`L` is [`bondperiod`](@ref), the turns the bond admits.
`R` is [`nphases`](@ref), those giving different structures, with `stab₂` the stabiliser of the site being attached.
An unsymmetric particle has one phase per bond.

## Choosing how a bond turns

[`PolyhedronParticleSpecies`](@ref) and [`PatchySphere`](@ref) take two keywords for this, each one value for the species or one per face.

`locking` decides whether a site pins its partner.
The default `true` holds it in the orientation the site's frame names; `false` admits every orientation the face allows.
Making a triangular prism's square sides rotation-free gives both the flat and the out-of-plane lattice, where locking gives the flat one only.

`twists` turns a site about its own normal by an angle in radians, picking which orientation the bond means.
Turning both faces of a bond by the same angle does not cancel; it turns the partner by twice that angle.
Turning one face of a symmetry orbit differently from the others splits the orbit and lowers the symmetry number, while turning a whole orbit leaves the symmetry alone.

You cannot select a single phase on a symmetric particle.
The two attachments are related by a symmetry, so a model that keeps the symmetry cannot prefer one.
Take both with `locking=false`, or break the symmetry with colors.

## Aligning frames across faces

Faces related by a symmetry of the particle must start their corner lists at corresponding corners.
Roly does this in two steps: first each face is rotated to a starting corner chosen from its own shape, which aligns frames up to `gauge`, then frames are aligned up to `stab`, which is what bonds need.
Where `gauge` exceeds `stab` the leftover turns are real, and two frames off by one describe different bonds.

A triangular prism with `h = a` shows why the second step is needed.
Its three square sides are alike, and the first step alone leaves one reference a quarter turn from the other two: legal by `gauge` (4), not by `stab` (2).
Bonding face 2 to face 3 then turns the neighbour out of the plane where face 2 to face 4 does not, and the prisms enumerate as `[1,3,6,22,73,357]` instead of the polyiamonds `[1,2,3,6,10,22]`.
With `h = 2` the sides are rectangles, the first step has something to work with, and the count is right without the second.

Alignment uses the rotations preserving the *coloring*, not all rotations of the solid.
A cube with four sticky sides and two caps has 24 rotations as a solid but only the 8 that keep sides apart from caps.
Using all 24 leaves adjacent sides a quarter turn apart and gives `[1,3,9,39,208,1402]` instead of the polyominoes `[1,2,4,9,21,56]`.

## What catches a mistake

[`site_symmetry`](@ref) compares frames up to `gauge`, so it cannot see these choices.
It is still useful, because a bad twist reference can only make it too small, never too large.
[`check_encoding`](@ref) compares it against the graph's symmetry number, so it catches every undercount.
