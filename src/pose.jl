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
    return  Pose{3,F,orientationtype{F}}(zeros(3), orientationtype(0.0, 0.0, 0.0))
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
function Pose(x, psi)
    return Pose{length(x),eltype(x),typeof(psi)}(x, psi)
end

Base.:(==)(p1::Pose, p2::Pose) = p1.x == p2.x && p1.psi == p2.psi
Base.isapprox(p1::Pose, p2::Pose; kwargs...) = isapprox(p1.x, p2.x; kwargs...) && isapprox(p1.psi, p2.psi; kwargs...)
Base.:*(p1::Pose, p2::Pose) = typeof(p1)(p2.psi * p1.x + p2.x, p1.psi * p2.psi)
Base.:*(R::Rotation, p::Pose) = typeof(p)(R * p.x, R * p.psi)
Base.:*(p::Pose, R::Rotation) = typeof(p)(p.x, p.psi * R)
Base.inv(p::Pose) = let ψinv = inv(p.psi); typeof(p)(-ψinv * p.x, ψinv) end
Base.:/(p1::Pose, p2) = p1 * inv(p2)
Base.:\(p1, p2::Pose) = inv(p1) * p2
Base.eltype(::Type{<:Pose{D,F}}) where {D,F} = F
dimension(::Pose{D}) where D = D
dimension(::Type{<:Pose{D}}) where D = D

"""
    tomatrix(p::Pose{D,F}) where {D,F}

Convert a pose `p` into a matrix of the form `[R(p.psi) p.x; 0 1]`.
"""
tomatrix(p::Pose{2,F}) where F = [p.psi p.x; zeros(F, 2)' F(1)]
tomatrix(p::Pose{3,F}) where F = [p.psi p.x; zeros(F, 3)' F(1)]

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
    return [sqrt(x^2 + y^2), atan(y, x)]
end
function pol2cart(r::F, psi::Real) where {F}
    psi = convert(F, psi)
    return [r * cos(psi),  r * sin(psi)]
end
