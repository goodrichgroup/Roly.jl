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
    Pose{D,F}([; orientationtype)) where {D,F<:Real}

Create a `D`-dimensional pose of eltype `F` at the origin and standard orientation.
The type of the orientational coordinates defaults to `Angle2d` in 2d
and `RotXYZ` in 3d, but can optionally be overridden.
"""
function Pose{D,F}() where {D,F} 
    throw(ArgumentError("Cannot create D=$D-dimensional pose, use D=2 or D=3."))
end
function Pose{2,F}(; orientationtype=Angle2d) where {F} 
    return Pose{2,F,orientationtype{F}}(zeros(2), orientationtype(0.0))
end
function Pose{3,F}(; orientationtype=RotXYZ) where {F} 
    return Pose{3,F,orientationtype{F}}(zeros(3), orientationtype(0.0, 0.0, 0.0))
end

"""
    Pose{D}([; orientationtype)) where {D}

Create a `D`-dimensional pose of eltype `Float64` at the origin and standard orientation.
The type of the orientational coordinates defaults to `Angle2d` in 2d
and `RotXYZ` in 3d, but can optionally be overridden.
"""
Pose{D}(; kwargs...) where {D} = Pose{D,Float64}(; kwargs...)

"""
    Pose(x, psi)

Create a pose from a position vector `x` and orientation `psi`.
"""
function Pose{D,F}(x, psi) where {D, F} 
    return Pose{D,eltype(x),typeof(psi)}(x, psi)
end
function Pose(x, psi)
    return Pose{length(x),eltype(x),typeof(psi)}(x, psi)
end


function Pose{3,F}(p::Pose{2,F}; orientationtype=RotXYZ) where F
    return Pose(SVector{3,F}(p.x..., F(0)), orientationtype(F(0), F(0), rotation_angle(p.psi)))
end
function Pose{3}(p::Pose{2,F}; orientationtype=RotXYZ) where F
    return Pose{3,F}(p; orientationtype)
end

promote_rule(::Type{<:Pose{D,F1}}, ::Type{<:Pose{D,F2}}) where {D, F1, F2} = promote_rule(F1, F2)

Base.:(==)(p1::Pose, p2::Pose) = p1.x == p2.x && p1.psi == p2.psi
Base.isapprox(p1::Pose, p2::Pose; kwargs...) = isapprox(p1.x, p2.x; kwargs...) && isapprox(p1.psi, p2.psi; kwargs...)

Base.:*(p1::Pose, p2::Pose) = Pose(p1.psi * p2.x + p1.x, p1.psi * p2.psi)
Base.:*(p::Pose, v::AbstractVector) = p.psi * v + p.x
Base.:*(R::Rotation, p::Pose) = typeof(p)(R * p.x, R * p.psi)
Base.:*(p::Pose, R::Rotation) = typeof(p)(p.x, p.psi * R)
Base.:+(v, p::Pose) = typeof(p)(p.x + v, p.psi)
Base.:+(p::Pose, v) = v + p

Base.inv(p::Pose) = let ψinv = inv(p.psi); typeof(p)(-ψinv * p.x, ψinv) end
Base.:/(p1::Pose, p2) = p1 * inv(p2)
Base.:\(p1, p2::Pose) = inv(p1) * p2

Base.Matrix(p::Pose{D,F}) where {D,F} = [p.psi p.x; zeros(F, D)' F(1)]

Base.eltype(::Type{<:Pose{D,F}}) where {D,F} = F

dimension(::Pose{D}) where D = D
dimension(::Type{<:Pose{D}}) where D = D

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
