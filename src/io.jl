### Writing binding rules and polyforms to files, and reading them back.
#
# A file holds a set of binding rules, a polyform, or a vector of polyforms, and always the rules:
# a polyform means nothing without them. Everything is written through `todict` and read through
# `fromdict`, which reduce rules and polyforms to named tuples of strings, numbers and arrays, so a
# file format only has to store those. TOML is built in, and JLD2 is provided by an extension once
# JLD2.jl is loaded.

const _FILEVERSION = 1

"""
    todict(rules::BindingRules)
    todict(p::Polyform; rules=false)

Serialize the object into a `Dict{String,Any}` of strings, numbers and arrays.

A polyform refers to the species and sites of its binding rules by number, and its rules are
included only with `rules=true`. [`fromdict`](@ref) reads the result back.
"""
function todict end

"""
    fromdict(::Type{BindingRules}, d)
    fromdict(::Type{Polyform}, d; rules=nothing)
    fromdict(::Type{ParticleSpecies}, d)

Deserialize the object from a dictionary written by [`todict`](@ref), or read from a file.

A polyform is read within `rules`, or within the rules its dictionary includes.
"""
function fromdict end

const _FileContent = Union{BindingRules,Polyform,AbstractVector{<:Polyform}}

"""
    tofile(file, rules::BindingRules; format=nothing)
    tofile(file, p::Polyform; format=nothing)
    tofile(file, ps::AbstractVector{<:Polyform}; format=nothing)

Write binding rules, a polyform, or a vector of polyforms to `file`, together with the binding
rules of the polyforms, and return `file`. All polyforms of a vector have to share one set of
binding rules.

`format` is `:toml` or `:jld2`. When it is not given, files ending in `.rlyb` or `.jld2` are
written as JLD2 and all others as TOML. JLD2 needs JLD2.jl to be loaded. JLD2 files are HDF5 files,
which other HDF5 implementations can read.
"""
function tofile(file::AbstractString, x::_FileContent; format=nothing)
    d = _filecontent(x)
    if _fileformat(file, format) == :toml
        open(io -> TOML.print(io, d; sorted=true, by=_keyorder), file, "w")
    else
        _jld2ext().tojld2(file, d)
    end
    return file
end

const _KEYPRIORITY = Dict("contenttype" => 1, "formatversion" => 2, "rolyversion" => 3, "rules" => 4, "type" => 5)
_keyorder(k) = (get(_KEYPRIORITY, k, 100), k)


"""
    fromfile(file; format=nothing)

Read the binding rules, polyform, or vector of polyforms that [`tofile`](@ref) wrote to `file`,
and return it.

`format` is `:toml` or `:jld2`, and is taken from the extension of `file` when not given, as for
[`tofile`](@ref).
"""
function fromfile(file::AbstractString; format=nothing)
    d = _fileformat(file, format) == :toml ? TOML.parsefile(file) : _jld2ext().fromjld2(file)

    contenttype = get(d, "contenttype", nothing)
    (contenttype isa AbstractString && startswith(contenttype, "roly/")) ||
        throw(ArgumentError("$file is not a Roly file."))
    d["formatversion"] <= _FILEVERSION || throw(
        ArgumentError(
            "Format too new: $file was written in file format version $(d["formatversion"]). This version of Roly can only " *
            "read up to version $_FILEVERSION.",
        ),
    )

    rules = fromdict(BindingRules, d["rules"])
    content = Symbol(chopprefix(contenttype, "roly/"))
    content == :rulesonly && return rules
    content == :polysingle && return fromdict(Polyform, d["polyform"]; rules)
    content == :polylist && return [fromdict(Polyform, p; rules) for p in d["polyforms"]]
    throw(ArgumentError("$file holds content of type $contenttype, which this version of Roly cannot read."))
end

function _fileheader(content::Symbol)
    return Dict{String,Any}(
        "contenttype" => "roly/$content",
        "formatversion" => _FILEVERSION,
        "rolyversion" => string(pkgversion(@__MODULE__)),
    )
end

_filecontent(rules::BindingRules) = merge!(_fileheader(:rulesonly), Dict("rules" => todict(rules)))

function _filecontent(p::Polyform)
    return merge!(_fileheader(:polysingle), Dict("rules" => todict(bindingrules(p)), "polyform" => todict(p; rules=false)))
end

function _filecontent(ps::AbstractVector{<:Polyform})
    isempty(ps) && throw(ArgumentError("An empty vector of polyforms has no binding rules to write."))
    rules = bindingrules(first(ps))
    all(p -> bindingrules(p) === rules, ps) ||
        throw(ArgumentError("All polyforms written to one file have to share one set of binding rules."))
    return merge!(_fileheader(:polylist),
                  Dict("rules" => todict(rules), "polyforms" => [todict(p; rules=false) for p in ps]))
end


function _fileformat(file::AbstractString, format)
    if !isnothing(format)
        format in (:toml, :jld2) || throw(ArgumentError("The file format has to be :toml or :jld2, not $(repr(format))."))
        return format
    end
    return lowercase(last(splitext(file))) in (".rlyb", ".jld2") ? :jld2 : :toml
end

function _jld2ext()
    ext = Base.get_extension(@__MODULE__, :JLD2Ext)
    isnothing(ext) && throw(ArgumentError("Reading and writing JLD2 files requires JLD2.jl to be loaded."))
    return ext
end
