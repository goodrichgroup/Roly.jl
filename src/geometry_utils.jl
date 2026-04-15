struct Pose{D,F,R<:Rotation{D,F}}
    x::SVector{D,F}
    ψ::R
end
function Pose{D}() where {D} 
    if D == 3
        return  Pose{D,Float64,RotXYZ{Float64}}(zeros(3), RotXYZ(0.0, 0.0, 0.0))
    elseif D == 2 
        return Pose{D,Float64,Angle2d{Float64}}(zeros(2), Angle2d(0.0))
    else
        throw(ArgumentError("Pose must be 2d or 3d"))
    end
end
function Pose(x, psi)
    return Pose{length(x),eltype(x),typeof(psi)}(x, psi)
end

Base.:(==)(p1::Pose, p2::Pose) = p1.x == p2.x && p1.ψ == p2.ψ
Base.isapprox(p1::Pose, p2::Pose; kwargs...) = isapprox(p1.x, p2.x; kwargs...) && isapprox(p1.ψ, p2.ψ; kwargs...)
Base.:*(p1::Pose, p2::Pose) = typeof(p1)(p2.ψ * p1.x + p2.x, p1.ψ * p2.ψ)
Base.:*(R::Rotation, p::Pose) = typeof(p)(R * p.x, R * p.ψ)
Base.:*(p::Pose, R::Rotation) = typeof(p)(p.x, p.ψ * R)

Base.inv(p::Pose) = let ψinv = inv(p.ψ); typeof(p)(-ψinv * p.x, ψinv) end
Base.:/(p1::Pose, p2::Pose) = p1 * inv(p2)
Base.:\(p1::Pose, p2::Pose) = inv(p1) * p2

tomatrix(p::Pose{2,F}) where F = [p.ψ p.x; zeros(F, 2)' F(1)]
tomatrix(p::Pose{3,F}) where F = [p.ψ p.x; zeros(F, 3)' F(1)]

function Base.show(io::Core.IO, p::Pose{D,F,R}) where {D,F,R}
    print(io, "$D-dimensional Pose{$D,$F,$(string(nameof(R)))}:\n")
    print(io, " - x: $(p.x)\n")
    if D == 2
        print(io, " - ψ: $(rotation_angle(p.ψ) / π)π")
    elseif D == 3
        print(io, " - ψ: $(rotation_angle(p.ψ) / π)π, $(rotation_axis(p.ψ))")
    end
end
dimension(::Pose{D}) where D = D
dimension(::Type{<:Pose{D}}) where D = D

function cart2pol(x::F, y::Real) where {F}
    y = convert(F, y)
    return [sqrt(x^2 + y^2), atan(y, x)]
end
function pol2cart(r::F, ψ::Real) where {F}
    ψ = convert(F, ψ)
    return [r * cos(ψ),  r * sin(ψ)]
end