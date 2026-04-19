function are_cutvertices(g::AbstractNautyGraph, vs::AbstractVector{<:Integer})
    n = nv(g)
    return are_cutvertices!(g, vs, zeros(Bool, n), zeros(Bool, n), zeros(Bool, n), zeros(Cint, n))
end

function are_cutvertices!(g::AbstractNautyGraph, vs::AbstractVector{<:Integer},
    is_target::AbstractVector{Bool}, forbidden::AbstractVector{Bool},
    explored::AbstractVector{Bool}, queue::AbstractVector)

    fill!(is_target, false)
    fill!(forbidden, false)
    forbidden[vs] .= true
    for v in vs
        is_target .|= NautyGraphs.adjrow(g, v)
    end
    is_target[vs] .= false

    n_targets = count(is_target)
    n_targets < 2 && return false

    v0 = findfirst(is_target)
    fill!(explored, false)
    explored[v0] = true
    n_found = 1

    queue[1] = v0
    q_start = 1
    q_end = 2

    while q_end > q_start
        v = queue[q_start]
        q_start += 1

        outneighs = NautyGraphs.adjrow(g, v)
        for neigh in eachindex(outneighs)
            outneighs[neigh] || continue
            forbidden[neigh] && continue
            if !explored[neigh]
                explored[neigh] = true
                if is_target[neigh]
                    n_found += 1
                    n_found == n_targets && return false
                end
                queue[q_end] = neigh
                q_end += 1
            end
        end
    end

    return true
end


function blockdiag!(g::AbstractNautyGraph, h::AbstractNautyGraph)
    ng, nh = nv(g), nv(h)
    add_vertices!(g, nh; vertex_labels=labels(h))
    for e in edges(h)
        add_edge!(g, e.src+ng, e.dst+ng)
    end
    return g
end