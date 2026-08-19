# Orientation and phases

A bond fixes which two sites touch, and how the two particles are turned relative to one another.
On a symmetric particle it can fix that turn in more than one way, and each way is a *phase*.

This page matters for 3D species whose faces have their own symmetry.
In 2D every bond has exactly one phase.

## What a site's pose says

A binding site's pose carries its position, its outward normal (the pose's local x axis), and its *twist reference*, the remaining freedom in the frame.
For a polyhedron face the twist reference is "local z points at the midpoint of the face's first edge", so it is a choice of where the face's corner list starts.

## Three numbers per site

| quantity | meaning | square side face of a triangular prism |
|---|---|---|
| `degree` | graph vertices the site owns under [`dartencoding`](@ref) | 4 |
| `gauge` | order of the face's own rotational symmetry about its normal | 4 |
| `stab` | order of that site's stabilizer in the particle's symmetry group | 2 |

Each divides the next, and all three are equal for a cube face.
[`facegauge`](@ref) computes `gauge` from the face, [`sitestabilizers`](@ref) computes `stab` from the particle.

## How many phases a bond has

A bond fixes the partner's position and normal, leaving a turn about the bond axis.
Two turns describe the same structure when a rotation of the particle relates them, so the count comes from the particle's symmetry.

Writing `q` for a site's [`twistfreedom`](@ref):

```
L = lcm(q₁, q₂)          q = stab if the site is locking, gauge if not
R = L / stab₂
```

`L` is [`bondperiod`](@ref), the turns the bond admits.
`R` is [`nphases`](@ref), those giving different structures, with `stab₂` the stabilizer of the site being attached.
A particle with no symmetry has one phase per bond.

## Choosing how a bond turns

[`PolyhedronParticleSpecies`](@ref) and [`PatchySphere`](@ref) take two keywords for this, each one value for the species or one per face.

`locking` decides whether a site pins its partner.
The default `true` holds it in the orientation the site's frame names, while `false` admits every orientation the face allows.
Making a triangular prism's square sides rotation-free gives both the flat and the out-of-plane lattice; locking gives the flat one only.

`twists` turns a site about its own normal by an angle in radians, picking which orientation the bond means.
Turning both faces of a bond by the same angle turns the partner by twice that angle.
Turning one face of a symmetry orbit differently from the others splits the orbit and lowers the symmetry number, while turning a whole orbit leaves the symmetry alone.

On a symmetric particle the phases come as a set, since a symmetry relates them.
Take them all with `locking=false`, or give the faces different colors to single one out.

## Frames on symmetric faces

Faces that a symmetry of the particle relates start their corner lists at corresponding corners, so a bond means the same thing at each of them.
Roly arranges this when the species is built, using the rotations that preserve the coloring.

## Checking a species

[`site_symmetry`](@ref) counts the rotations that carry every site onto a site of the same color.
[`check_encoding`](@ref) compares that against the graph's symmetry number and throws if they disagree.
Every built-in species runs it in its constructor; call it in your own if you write graph labels by hand.
