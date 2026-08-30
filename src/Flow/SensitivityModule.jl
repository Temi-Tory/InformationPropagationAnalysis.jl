module SensitivityModule

if isdefined(parentmodule(@__MODULE__), :FlowModule)
    const FlowModule = parentmodule(@__MODULE__).FlowModule
else
    include("FlowModule.jl")
end
using .FlowModule

include("_CapacityShared.jl")
if isdefined(parentmodule(@__MODULE__), :CapacityTypes)
    const CapacityTypes = parentmodule(@__MODULE__).CapacityTypes
else
    include("CapacityTypes.jl")
end
using .CapacityTypes

export SensitivityResult,
       critical_edge_ranking,
       marginal_capacity_values,
       birnbaum_importance,
       analyze_sensitivity

struct SensitivityResult
    critical_edges::Vector{CriticalEdgeRecord}
    marginal_capacity::Dict{Tuple{Int64,Int64},Float64}
    birnbaum::Dict{Tuple{Int64,Int64},Float64}
end

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

        # Reverse-residual traversal from sink:
        # from current node v we can step to u iff residual edge u->v exists.
        # Case 1 (forward residual): original edge (u,v) has c(u,v)-f(u,v)>0.
        for u in get(incoming_index, v, Set{Int64}())
            residual = capacities[(u, v)] - get(flow, (u, v), 0.0)
            if residual > tol && !(u in reachable)
                push!(reachable, u)
                push!(queue, u)
            end
        end

        # Case 2 (backward residual): original edge (v,w) has f(v,w)>0,
        # yielding residual edge w->v, so reverse traversal steps v->w.
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
    super_source = flow_result.super_source
    super_sink = flow_result.super_sink
    all_nodes = union(Set(first.(edgelist)), Set(last.(edgelist)))
    all_aug_nodes = union(Set(keys(aug_out)), Set(keys(aug_in)), all_nodes, Set([super_source, super_sink]))
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

    return candidates
end

function _zero_capacity_reruns(
    edges::Vector{Tuple{Int64,Int64}},
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10
)::Dict{Tuple{Int64,Int64},FlowSolveResult}
    cache = Dict{Tuple{Int64,Int64},FlowSolveResult}()
    for edge in edges
        modified = _copy_capacities(capacities)
        modified[edge] = 0.0
        cache[edge] = _solve_with_algorithm(
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            modified,
            source_nodes,
            sink_nodes;
            tol=tol,
            validate=true
        )
    end
    return cache
end

function _critical_edge_ranking_from_cache(
    edges::Vector{Tuple{Int64,Int64}},
    zero_cache::Dict{Tuple{Int64,Int64},FlowSolveResult},
    baseline::Float64
)::Vector{CriticalEdgeRecord}
    rankings = CriticalEdgeRecord[]
    for edge in edges
        rerun = zero_cache[edge]
        drop = baseline - rerun.max_flow
        push!(rankings, CriticalEdgeRecord(edge, baseline, rerun.max_flow, drop))
    end
    sort!(rankings; by=x -> (-x.drop, x.edge))
    return rankings
end

function _birnbaum_from_cache(
    mincut_edges::Vector{Tuple{Int64,Int64}},
    zero_cache::Dict{Tuple{Int64,Int64},FlowSolveResult},
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10
)::Dict{Tuple{Int64,Int64},Float64}
    importance = Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in edgelist)

    for edge in mincut_edges
        haskey(zero_cache, edge) || throw(ArgumentError("Missing zero-capacity cache entry for edge $edge"))
        cap_inf = _copy_capacities(capacities)
        cap_inf[edge] = Inf

        rerun_inf = _solve_with_algorithm(
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            cap_inf,
            source_nodes,
            sink_nodes;
            tol=tol,
            validate=true
        )
        rerun_zero = zero_cache[edge]

        if rerun_inf.is_unbounded && rerun_zero.is_unbounded
            throw(ArgumentError("Birnbaum importance is undefined for edge $edge because both c_e=Inf and c_e=0 scenarios are unbounded."))
        elseif rerun_inf.is_unbounded
            importance[edge] = Inf
        elseif rerun_zero.is_unbounded
            throw(ArgumentError("Birnbaum importance monotonicity violated for edge $edge: c_e=0 is unbounded while c_e=Inf is bounded."))
        else
            importance[edge] = rerun_inf.max_flow - rerun_zero.max_flow
        end
    end

    return importance
end

"""
    critical_edge_ranking(...)

Rank saturated candidate edges by the exact drop in max flow when the edge capacity
is set to 0. Only saturated edges from the solved flow are tested.
"""
function critical_edge_ranking(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10
)::Vector{CriticalEdgeRecord}
    _require_bounded_baseline(flow_result)
    baseline = flow_result.max_flow
    zero_cache = _zero_capacity_reruns(
        flow_result.saturated_edges,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes;
        algorithm=algorithm,
        tol=tol
    )
    return _critical_edge_ranking_from_cache(flow_result.saturated_edges, zero_cache, baseline)
end

"""
    marginal_capacity_values(...)

Compute exact marginal capacity values by increasing each saturated edge capacity by 1
and recomputing max flow. Unsaturated edges are assigned 0 by theorem-guided pruning.
"""
function marginal_capacity_values(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    delta::Float64=1.0
)::Dict{Tuple{Int64,Int64},Float64}
    _require_bounded_baseline(flow_result)
    delta <= 0.0 && throw(ArgumentError("delta must be positive."))
    values = Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in edgelist)
    baseline = flow_result.max_flow

    for edge in flow_result.saturated_edges
        modified = _copy_capacities(capacities)
        modified[edge] = capacities[edge] + delta
        rerun = _solve_with_algorithm(
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            modified,
            source_nodes,
            sink_nodes;
            tol=tol,
            validate=true
        )
        values[edge] = (rerun.max_flow - baseline) / delta
    end

    return values
end

"""
    birnbaum_importance(...)

Compute exact Birnbaum importance:
    I(e) = max_flow(c_e = Inf) - max_flow(c_e = 0)

Only edges in the solved min-cut partition crossing from S to T are evaluated; all
other edges are assigned 0 by pruning.
"""
function birnbaum_importance(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10
)::Dict{Tuple{Int64,Int64},Float64}
    _require_bounded_baseline(flow_result)
    mincut_edges = _edges_in_some_mincut(
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        flow_result;
        tol=tol
    )
    zero_cache = _zero_capacity_reruns(
        mincut_edges,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes;
        algorithm=algorithm,
        tol=tol
    )
    return _birnbaum_from_cache(
        mincut_edges,
        zero_cache,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes;
        algorithm=algorithm,
        tol=tol
    )
end

"""
    analyze_sensitivity(...)

Run the exact sensitivity module on top of a solved flow result.
"""
function analyze_sensitivity(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    delta::Float64=1.0
)::SensitivityResult
    _require_bounded_baseline(flow_result)
    some_mincut_edges = _edges_in_some_mincut(
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        flow_result;
        tol=tol
    )
    zero_cache = _zero_capacity_reruns(
        flow_result.saturated_edges,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes;
        algorithm=algorithm,
        tol=tol
    )
    critical = _critical_edge_ranking_from_cache(flow_result.saturated_edges, zero_cache, flow_result.max_flow)
    # Marginal (+delta) perturbations are distinct from zero-capacity perturbations,
    # so they require separate reruns and cannot reuse zero_cache.
    marginal = marginal_capacity_values(
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes,
        flow_result;
        algorithm=algorithm,
        tol=tol,
        delta=delta
    )
    birnbaum = _birnbaum_from_cache(
        some_mincut_edges,
        zero_cache,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes;
        algorithm=algorithm,
        tol=tol
    )

    return SensitivityResult(critical, marginal, birnbaum)
end

end
