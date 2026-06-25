"""
    Pose{D,F,R<:Rotation{D,F}}

A `Pose` describes the configuration of a `D`-dimensional rigid body.
It consists of translational coordinates `pose.x` and orientational
coordinates `pose.psi`.
"""
struct Pose{D,F,R<:Rotation{D,F}}
    x::SVector{D,F}
    psi::R
end

"""
    Pose(x, psi)

Create a pose from a position vector `x` and orientation `psi`.
"""
Pose(x::AbstractVector, psi::R) where {R<:Rotation} = Pose{length(x),eltype(x),R}(SVector{length(x),eltype(x)}(x), psi)
Pose{D,F}(x::AbstractVector, psi::R) where {D,F,R<:Rotation{D,F}} = Pose{D,F,R}(SVector{D,F}(x), psi)

"""
    one(::Type{Pose{D,F,R}})
    one(p::Pose)

Return the identity pose (origin, identity rotation) of the given type.
"""
Base.one(::Type{Pose{D,F,R}}) where {D,F,R} = Pose{D,F,R}(zero(SVector{D,F}), one(R))
Base.one(p::Pose) = one(typeof(p))

"""
    Pose{2,F}(; orientationtype=Angle2d)
    Pose{3,F}(; orientationtype=RotXYZ)
    Pose{D}(; ...)

Return the identity pose of eltype `F` in `D` dimensions, using the default rotation type
(`Angle2d` in 2D, `RotXYZ` in 3D). The rotation type can be overridden with `orientationtype`.
"""
Pose{2,F}(; orientationtype=Angle2d) where {F} = one(Pose{2,F,orientationtype{F}})
Pose{3,F}(; orientationtype=RotXYZ) where {F} = one(Pose{3,F,orientationtype{F}})
Pose{D,F}(; kwargs...) where {D,F} = throw(ArgumentError("no default pose for D=$D, use D=2 or D=3"))
Pose{D}(; kwargs...) where {D} = Pose{D,Float64}(; kwargs...)

"""
    Pose{3,F}(p::Pose{2,F})

Lift a 2D pose into 3D by appending a zero z-coordinate and wrapping the 2D rotation angle
into the z-axis rotation of the 3D rotation type (default `RotXYZ`).
"""
function Pose{3,F}(p::Pose{2,F}; orientationtype=RotXYZ) where {F}
    return Pose(SVector{3,F}(p.x..., F(0)), orientationtype(F(0), F(0), rotation_angle(p.psi)))
end
function Pose{3}(p::Pose{2,F}; orientationtype=RotXYZ) where {F}
    return Pose{3,F}(p; orientationtype)
end

Base.:(==)(p1::Pose, p2::Pose) = p1.x == p2.x && p1.psi == p2.psi
Base.isapprox(p1::Pose, p2::Pose; kwargs...) = isapprox(p1.x, p2.x; kwargs...) && isapprox(p1.psi, p2.psi; kwargs...)

Base.:*(p1::Pose, p2::Pose) = Pose(p1.psi * p2.x + p1.x, p1.psi * p2.psi)
Base.:*(p::Pose, v::AbstractVector) = p.psi * v + p.x
Base.:*(R::Rotation, p::Pose) = typeof(p)(R * p.x, R * p.psi)
Base.:*(p::Pose, R::Rotation) = typeof(p)(p.x, p.psi * R)
Base.:+(v, p::Pose) = typeof(p)(p.x + v, p.psi)
Base.:+(p::Pose, v) = v + p

Base.inv(p::Pose) =
    let ψinv = inv(p.psi)
        typeof(p)(-ψinv * p.x, ψinv)
    end
Base.:/(p1::Pose, p2) = p1 * inv(p2)
Base.:\(p1, p2::Pose) = inv(p1) * p2

Base.Matrix(p::Pose{D,F}) where {D,F} = [p.psi p.x; zeros(F, D)' F(1)]

Base.eltype(::Type{<:Pose{D,F}}) where {D,F} = F

dimension(::Pose{D}) where {D} = D
dimension(::Type{<:Pose{D}}) where {D} = D

function Base.show(io::Core.IO, p::Pose{D,F,R}) where {D,F,R}
    if D == 2
        print(io, "[x: $(p.x); ψ: $(rotation_angle(p.psi) / π)π]")
    elseif D == 3
        print(io, "[x: $(p.x); ψ: $(rotation_angle(p.psi) / π)π, $(rotation_axis(p.psi))]")
    end
end
function Base.show(io::Core.IO, ::MIME"text/plain", p::Pose{D,F,R}) where {D,F,R}
    print(io, "$D-dimensional Pose{$D,$F,$(string(nameof(R)))}:\n")
    print(io, " - x: $(p.x)\n")
    if D == 2
        print(io, " - psi: $(rotation_angle(p.psi) / π)π")
    elseif D == 3
        print(io, " - psi: $(rotation_angle(p.psi) / π)π, $(rotation_axis(p.psi))")
    end
end

function cart2pol(x::F, y::Real) where {F}
    y = convert(F, y)
    return SVector(sqrt(x^2 + y^2), atan(y, x))
end
function pol2cart(r::F, psi::Real) where {F}
    psi = convert(F, psi)
    return SVector(r * cos(psi), r * sin(psi))
end
