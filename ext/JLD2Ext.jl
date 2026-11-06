module JLD2Ext

using Roly, JLD2

# The contents of a file are written as HDF5 groups and plain arrays, so that any HDF5 implementation
# can read them: a dictionary becomes a group, a list of equal-length rows a matrix with one column
# per row, and any other list a group with members "1", "2", .... The polyforms of a file are written
# as flat columns together with their lengths, since a file may hold very many of them.

function tojld2(file::AbstractString, d::AbstractDict)
    jldopen(file, "w") do f
        for (k, v) in pairs(d)
            k == "polyforms" ? _writepolyforms(f, k, v) : _write(f, k, v)
        end
    end
    return nothing
end

function fromjld2(file::AbstractString)
    return jldopen(file, "r") do f
        Dict{String,Any}(k => k == "polyforms" ? _readpolyforms(f[k]) : _read(f[k]) for k in keys(f))
    end
end

_write(f, path, x) = (f[path] = x)

function _write(f, path, x::AbstractDict)
    for (k, v) in pairs(x)
        _write(f, "$path/$k", v)
    end
    return nothing
end

function _write(f, path, x::AbstractVector)
    if isempty(x)
        f[path] = Int[]
    elseif all(v -> v isa Union{Number,AbstractString}, x)
        f[path] = collect(x)
    elseif all(v -> v isa AbstractVector{<:Number} && length(v) == length(first(x)) > 0, x)
        f[path] = reduce(hcat, x)
    else
        for (i, v) in enumerate(x)
            _write(f, "$path/$i", v)
        end
    end
    return nothing
end

function _writepolyforms(f, path, ps)
    rows = [r for p in ps for r in p["bonds"]]
    f["$path/species"] = reduce(vcat, (p["species"] for p in ps); init=Int[])
    f["$path/nparticles"] = [length(p["species"]) for p in ps]
    f["$path/bonds"] = isempty(rows) ? zeros(Int, 5, 0) : reduce(hcat, rows)
    f["$path/nbonds"] = [length(p["bonds"]) for p in ps]
    return nothing
end

_read(x) = x
_read(x::AbstractMatrix) = [collect(c) for c in eachcol(x)]

function _read(g::JLD2.Group)
    ks = keys(g)
    all(k -> all(isdigit, k), ks) && return [_read(g[string(i)]) for i in 1:length(ks)]
    return Dict{String,Any}(k => _read(g[k]) for k in ks)
end

function _readpolyforms(g::JLD2.Group)
    species, nparticles, bonds, nbonds = g["species"], g["nparticles"], g["bonds"], g["nbonds"]
    ps = Vector{Dict{String,Any}}(undef, length(nparticles))
    a = b = 0
    for i in eachindex(nparticles)
        ps[i] = Dict{String,Any}(
            "species" => species[(a + 1):(a + nparticles[i])],
            "bonds" => [bonds[:, j] for j in (b + 1):(b + nbonds[i])],
        )
        a += nparticles[i]
        b += nbonds[i]
    end
    return ps
end

end
