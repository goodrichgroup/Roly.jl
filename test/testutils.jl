# Helpers shared across the test suite.

# helper functions to override encoding choice
dartspecies(p; colors=1:Roly.nfaces(p), locking=true, twists=0) =
    Roly._polyhedronspecies(p, colors, locking, twists, false)
cyclespecies(p; colors=1:Roly.nfaces(p), locking=true, twists=0) =
    Roly._polyhedronspecies(p, colors, locking, twists, true)

dartsphere(p, r=1; colors=1:Roly.nfaces(p), locking=true, twists=0) =
    Roly._patchysphere(p, r, colors, locking, twists, false)
cyclesphere(p, r=1; colors=1:Roly.nfaces(p), locking=true, twists=0) =
    Roly._patchysphere(p, r, colors, locking, twists, true)

# helper function to disable the onlattice check
function withoutlattice(sys::BindingRules)
    vals = Any[getfield(sys, f) for f in fieldnames(typeof(sys))]
    vals[findfirst(==(:_onlattice), collect(fieldnames(typeof(sys))))] = false
    return typeof(sys)(vals...)
end
