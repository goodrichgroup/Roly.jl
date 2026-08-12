# Helpers shared across the test suite.

# Which graph encoding a species gets is not a modelling choice, so the public constructors pick
# it -- the sparse one whenever it provably carries everything the dart encoding would -- and do
# not expose the choice. These reach the internal builders to force one, which is how the suite
# checks that the two agree wherever both are valid, and that the sparse one is rejected where
# it is not.
dartspecies(p; colors=1:Roly.nfaces(p), locking=true, twists=0) =
    Roly._polyhedronspecies(p, colors, locking, twists, false)
cyclespecies(p; colors=1:Roly.nfaces(p), locking=true, twists=0) =
    Roly._polyhedronspecies(p, colors, locking, twists, true)

dartsphere(p, r=1; colors=1:Roly.nfaces(p), locking=true, twists=0) =
    Roly._patchysphere(p, r, colors, locking, twists, false)
cyclesphere(p, r=1; colors=1:Roly.nfaces(p), locking=true, twists=0) =
    Roly._patchysphere(p, r, colors, locking, twists, true)

# The same rules with the on-lattice shortcut switched off, so a system that takes it can be
# enumerated both ways and the two answers compared. `_onlattice` is derived in the constructor
# and there is deliberately no way to set it, hence the field surgery.
function withoutlattice(sys::BindingRules)
    vals = Any[getfield(sys, f) for f in fieldnames(typeof(sys))]
    vals[findfirst(==(:_onlattice), collect(fieldnames(typeof(sys))))] = false
    return typeof(sys)(vals...)
end
