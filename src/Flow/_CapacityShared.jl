function _solve_with_algorithm(
    algorithm::Symbol,
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    tol::Float64=1e-10,
    validate::Bool=true
)
    if algorithm === :edmonds_karp
        return solve_max_flow_edmonds_karp(
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            source_nodes,
            sink_nodes;
            tol=tol,
            validate=validate
        )
    elseif algorithm === :dinic
        return solve_max_flow_dinic(
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            source_nodes,
            sink_nodes;
            tol=tol,
            validate=validate
        )
    elseif algorithm === :push_relabel
        return solve_max_flow_push_relabel(
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            source_nodes,
            sink_nodes;
            tol=tol,
            validate=validate
        )
    else
        throw(ArgumentError("Unsupported algorithm: $algorithm. Expected :edmonds_karp, :dinic, or :push_relabel."))
    end
end

function _copy_capacities(capacities::Dict{Tuple{Int64,Int64},Float64})
    return Dict{Tuple{Int64,Int64},Float64}(capacities)
end

function _backward_reachable_residual(
    sink::Int64,
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    flow::Dict{Tuple{Int64,Int64},Float64},
    tol::Float64;
    finite_caps_only::Bool=false
)::Set{Int64}
    reachable = Set{Int64}([sink])
    queue = Int64[sink]
    head = 1

    while head <= length(queue)
        v = queue[head]
        head += 1

        for u in get(incoming_index, v, Set{Int64}())
            residual = capacities[(u, v)] - get(flow, (u, v), 0.0)
            if residual > tol && !(u in reachable)
                push!(reachable, u)
                push!(queue, u)
            end
        end

        for w in get(outgoing_index, v, Set{Int64}())
            if finite_caps_only && !isfinite(get(capacities, (v, w), Inf))
                continue
            end
            residual = get(flow, (v, w), 0.0)
            if residual > tol && !(w in reachable)
                push!(reachable, w)
                push!(queue, w)
            end
        end
    end

    return reachable
end

function _edges_in_some_mincut(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    flow_result::FlowSolveResult;
    tol::Float64=1e-10
)::Vector{Tuple{Int64,Int64}}
    aug_out = flow_result.augmented_outgoing
    aug_in = flow_result.augmented_incoming
    aug_caps = flow_result.augmented_capacities
    super_sink = flow_result.super_sink
    all_nodes = union(Set(first.(edgelist)), Set(last.(edgelist)))
    all_aug_nodes = union(Set(keys(aug_out)), Set(keys(aug_in)), all_nodes, Set([flow_result.super_source, super_sink]))
    aug_flow = flow_result.augmented_flow
    can_reach_sink = _backward_reachable_residual(
        super_sink,
        aug_out,
        aug_in,
        aug_caps,
        aug_flow,
        tol;
        finite_caps_only=true
    )
    S_star = flow_result.mincut_S
    S_double_star = setdiff(all_aug_nodes, can_reach_sink)

    candidates = Tuple{Int64,Int64}[]
    for e in edgelist
        u, v = e
        saturated = abs(get(flow_result.flow, e, 0.0) - capacities[e]) <= tol
        if saturated && (u in S_double_star) && !(v in S_star)
            push!(candidates, e)
        end
    end

    sort!(candidates)
    return candidates
end

# ── Shared utility: bounded baseline guard ────────────────
# Used identically in 8 modules. Defined once here.
function _require_bounded_baseline(
    flow_result::FlowSolveResult,
    context::String="Analysis"
)::Nothing
    flow_result.is_unbounded && throw(ArgumentError(
        "$context is undefined for an unbounded flow result."
    ))
    nothing
end

# ── Shared utility: graph node extraction ─────────────────
# Canonical single implementation used by all modules.
function _graph_nodes_set(
    edgelist::Vector{Tuple{Int64,Int64}}
)::Set{Int64}
    return union(Set(first.(edgelist)), Set(last.(edgelist)))
end