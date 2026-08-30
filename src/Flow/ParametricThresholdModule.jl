module ParametricThresholdModule

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

export DegradationThreshold,
       UpgradeThreshold,
       ParametricThresholdResult,
       find_degradation_threshold,
       find_upgrade_threshold,
       find_all_degradation_thresholds,
       analyze_parametric_thresholds

struct DegradationThreshold
    target_edge::Tuple{Int64,Int64}
    original_capacity::Float64
    threshold_capacity::Float64
    degradation_margin::Float64
    target_flow::Float64
    baseline_flow::Float64
    target_achievable::Bool
    target_reachable_at_zero::Bool
    solver_calls::Int64
end

struct UpgradeThreshold
    target_edge::Tuple{Int64,Int64}
    original_capacity::Float64
    required_capacity::Float64
    required_increase::Float64
    target_flow::Float64
    baseline_flow::Float64
    already_sufficient::Bool
    upgrade_ineffective::Bool
    solver_calls::Int64
end

struct ParametricThresholdResult
    degradation_thresholds::Vector{DegradationThreshold}
    target_flow::Float64
    baseline_flow::Float64
end

function _validate_target_edge(
    edgelist::Vector{Tuple{Int64,Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    target_edge::Tuple{Int64,Int64}
)::Float64
    target_edge in Set(edgelist) || throw(ArgumentError("target_edge $target_edge is not present in edgelist."))
    haskey(capacities, target_edge) || throw(ArgumentError("Missing capacity for target_edge $target_edge."))
    c = capacities[target_edge]
    (isnan(c) || c < 0.0) && throw(ArgumentError("Invalid capacity for target_edge $target_edge: $c"))
    isfinite(c) || throw(ArgumentError("Parametric threshold analysis currently requires a finite original capacity for target_edge $target_edge."))
    return c
end

function _solve_and_count(
    call_count::Vector{Int64},
    algorithm::Symbol,
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    tol::Float64=1e-10
)::FlowSolveResult
    call_count[1] += 1
    return _solve_with_algorithm(
        algorithm,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes;
        tol=tol,
        validate=true
    )
end

function _solve_edge_capacity_and_count(
    call_count::Vector{Int64},
    algorithm::Symbol,
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    target_edge::Tuple{Int64,Int64},
    edge_capacity::Float64;
    tol::Float64=1e-10
)::FlowSolveResult
    modified = _copy_capacities(capacities)
    modified[target_edge] = edge_capacity
    return _solve_and_count(
        call_count,
        algorithm,
        edgelist,
        outgoing_index,
        incoming_index,
        modified,
        source_nodes,
        sink_nodes;
        tol=tol
    )
end

_partition_equal(a::Set{Int64}, b::Set{Int64}) = a == b

function _closed_form_threshold(
    c_lo::Float64,
    c_hi::Float64,
    flow_lo::Float64,
    flow_hi::Float64,
    target_flow::Float64,
    tol::Float64
)::Tuple{Float64,Bool}
    c_hi < c_lo && throw(ArgumentError("Invalid threshold interval: c_lo=$c_lo, c_hi=$c_hi."))
    abs(c_hi - c_lo) <= tol && return c_hi, false

    slope = (flow_hi - flow_lo) / (c_hi - c_lo)
    if slope <= tol
        return c_hi, true
    end

    threshold = c_lo + (target_flow - flow_lo) / slope
    return clamp(threshold, c_lo, c_hi), false
end

function _degradation_recurse(
    call_count::Vector{Int64},
    algorithm::Symbol,
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    target_edge::Tuple{Int64,Int64},
    target_flow::Float64,
    c_lo::Float64,
    c_hi::Float64,
    flow_lo::Float64,
    flow_hi::Float64,
    part_lo::Set{Int64},
    part_hi::Set{Int64},
    depth::Int;
    tol::Float64=1e-10,
    max_depth::Int=64
)::Float64
    if _partition_equal(part_lo, part_hi)
        threshold, flat = _closed_form_threshold(c_lo, c_hi, flow_lo, flow_hi, target_flow, tol)
        return flat ? c_hi : threshold
    end

    depth > max_depth && throw(ArgumentError("Parametric threshold recursion exceeded max_depth=$max_depth at interval [c_lo=$c_lo, c_hi=$c_hi]. This may indicate degenerate capacity structure or nearly-equal partition boundaries. Increase max_depth or inspect the graph."))

    c_mid = (c_lo + c_hi) / 2.0
    mid = _solve_edge_capacity_and_count(
        call_count,
        algorithm,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes,
        target_edge,
        c_mid;
        tol=tol
    )
    mid.is_unbounded && throw(ArgumentError("Parametric threshold solve became unbounded at c_e=$c_mid for edge $target_edge."))

    if mid.max_flow >= target_flow
        return _degradation_recurse(
            call_count,
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            source_nodes,
            sink_nodes,
            target_edge,
            target_flow,
            c_lo,
            c_mid,
            flow_lo,
            mid.max_flow,
            part_lo,
            mid.mincut_S,
            depth + 1;
            tol=tol,
            max_depth=max_depth
        )
    else
        return _degradation_recurse(
            call_count,
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            source_nodes,
            sink_nodes,
            target_edge,
            target_flow,
            c_mid,
            c_hi,
            mid.max_flow,
            flow_hi,
            mid.mincut_S,
            part_hi,
            depth + 1;
            tol=tol,
            max_depth=max_depth
        )
    end
end

function _upgrade_recurse(
    call_count::Vector{Int64},
    algorithm::Symbol,
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    target_edge::Tuple{Int64,Int64},
    target_flow::Float64,
    c_lo::Float64,
    c_hi::Float64,
    flow_lo::Float64,
    flow_hi::Float64,
    part_lo::Set{Int64},
    part_hi::Set{Int64},
    depth::Int;
    tol::Float64=1e-10,
    max_depth::Int=64
)::Tuple{Float64,Bool}
    if _partition_equal(part_lo, part_hi)
        threshold, flat = _closed_form_threshold(c_lo, c_hi, flow_lo, flow_hi, target_flow, tol)
        return flat ? (Inf, true) : (threshold, false)
    end

    depth > max_depth && throw(ArgumentError("Parametric threshold recursion exceeded max_depth=$max_depth at interval [c_lo=$c_lo, c_hi=$c_hi]. This may indicate degenerate capacity structure or nearly-equal partition boundaries. Increase max_depth or inspect the graph."))

    c_mid = (c_lo + c_hi) / 2.0
    mid = _solve_edge_capacity_and_count(
        call_count,
        algorithm,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes,
        target_edge,
        c_mid;
        tol=tol
    )
    mid.is_unbounded && throw(ArgumentError("Parametric threshold solve became unbounded at c_e=$c_mid for edge $target_edge."))

    if mid.max_flow >= target_flow
        return _upgrade_recurse(
            call_count,
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            source_nodes,
            sink_nodes,
            target_edge,
            target_flow,
            c_lo,
            c_mid,
            flow_lo,
            mid.max_flow,
            part_lo,
            mid.mincut_S,
            depth + 1;
            tol=tol,
            max_depth=max_depth
        )
    else
        return _upgrade_recurse(
            call_count,
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            source_nodes,
            sink_nodes,
            target_edge,
            target_flow,
            c_mid,
            c_hi,
            mid.max_flow,
            flow_hi,
            mid.mincut_S,
            part_hi,
            depth + 1;
            tol=tol,
            max_depth=max_depth
        )
    end
end

"""
    find_degradation_threshold(edgelist, outgoing_index, incoming_index,
                               capacities, source_nodes, sink_nodes,
                               target_edge, target_flow;
                               algorithm=:dinic, tol=1e-10, max_depth=64)

Find the exact minimum capacity of `target_edge` such that max flow remains at least
`target_flow`. The threshold is exact: partition changes are located recursively, and
on any fixed-partition segment the threshold is solved analytically in closed form.
"""
function find_degradation_threshold(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    target_edge::Tuple{Int64,Int64},
    target_flow::Float64;
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    max_depth::Int=64
)::DegradationThreshold
    target_flow > 0.0 || throw(ArgumentError("target_flow must be positive."))
    max_depth >= 0 || throw(ArgumentError("max_depth must be nonnegative."))

    original_capacity = _validate_target_edge(edgelist, capacities, target_edge)
    call_count = Int64[0]

    zero_result = _solve_edge_capacity_and_count(
        call_count,
        algorithm,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes,
        target_edge,
        0.0;
        tol=tol
    )
    zero_result.is_unbounded && throw(ArgumentError("Degradation threshold is undefined because the c_e=0 solve is unbounded for edge $target_edge."))

    baseline_result = _solve_edge_capacity_and_count(
        call_count,
        algorithm,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes,
        target_edge,
        original_capacity;
        tol=tol
    )
    baseline_result.is_unbounded && throw(ArgumentError("Degradation threshold is undefined for an unbounded baseline solve on edge $target_edge."))

    if zero_result.max_flow >= target_flow
        return DegradationThreshold(
            target_edge,
            original_capacity,
            0.0,
            original_capacity,
            target_flow,
            baseline_result.max_flow,
            true,
            true,
            call_count[1]
        )
    end

    # threshold_capacity=Inf is a sentinel value indicating the
    # target is unachievable even at full capacity. Callers must
    # check target_achievable=false before using threshold_capacity
    # or degradation_margin arithmetically.
    if baseline_result.max_flow < target_flow
        return DegradationThreshold(
            target_edge,
            original_capacity,
            Inf,
            0.0,
            target_flow,
            baseline_result.max_flow,
            false,
            false,
            call_count[1]
        )
    end

    threshold = _degradation_recurse(
        call_count,
        algorithm,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes,
        target_edge,
        target_flow,
        0.0,
        original_capacity,
        zero_result.max_flow,
        baseline_result.max_flow,
        zero_result.mincut_S,
        baseline_result.mincut_S,
        0;
        tol=tol,
        max_depth=max_depth
    )

    return DegradationThreshold(
        target_edge,
        original_capacity,
        threshold,
        original_capacity - threshold,
        target_flow,
        baseline_result.max_flow,
        true,
        false,
        call_count[1]
    )
end

"""
    find_upgrade_threshold(edgelist, outgoing_index, incoming_index,
                           capacities, source_nodes, sink_nodes,
                           target_edge, target_flow;
                           algorithm=:dinic, tol=1e-10, max_depth=64)

Find the exact minimum capacity `target_edge` must have for max flow to first reach at
least `target_flow`. A doubling search finds an exact enclosing interval, then the
threshold is solved recursively using exact partition-change recursion.
"""
function find_upgrade_threshold(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    target_edge::Tuple{Int64,Int64},
    target_flow::Float64;
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    max_depth::Int=64
)::UpgradeThreshold
    target_flow > 0.0 || throw(ArgumentError("target_flow must be positive."))
    max_depth >= 0 || throw(ArgumentError("max_depth must be nonnegative."))

    original_capacity = _validate_target_edge(edgelist, capacities, target_edge)
    call_count = Int64[0]

    baseline_result = _solve_edge_capacity_and_count(
        call_count,
        algorithm,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes,
        target_edge,
        original_capacity;
        tol=tol
    )
    baseline_result.is_unbounded && throw(ArgumentError("Upgrade threshold is undefined for an unbounded baseline solve on edge $target_edge."))

    if baseline_result.max_flow >= target_flow
        return UpgradeThreshold(
            target_edge,
            original_capacity,
            original_capacity,
            0.0,
            target_flow,
            baseline_result.max_flow,
            true,
            false,
            call_count[1]
        )
    end

    base_scale = max(1.0, original_capacity)
    max_scale = (2.0 ^ 32) * base_scale
    search_cap = base_scale
    prev_flow = baseline_result.max_flow
    c_hi = original_capacity + search_cap

    while true
        search_cap > max_scale && throw(ArgumentError("Upgrade doubling search exceeded 2^32 * max(1.0, original_capacity) for edge $target_edge without reaching target_flow=$target_flow."))

        hi_result = _solve_edge_capacity_and_count(
            call_count,
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            source_nodes,
            sink_nodes,
            target_edge,
            c_hi;
            tol=tol
        )
        hi_result.is_unbounded && throw(ArgumentError("Upgrade threshold solve became unbounded at c_e=$c_hi for edge $target_edge."))

        if hi_result.max_flow >= target_flow
            required_capacity, ineffective = _upgrade_recurse(
                call_count,
                algorithm,
                edgelist,
                outgoing_index,
                incoming_index,
                capacities,
                source_nodes,
                sink_nodes,
                target_edge,
                target_flow,
                original_capacity,
                c_hi,
                baseline_result.max_flow,
                hi_result.max_flow,
                baseline_result.mincut_S,
                hi_result.mincut_S,
                0;
                tol=tol,
                max_depth=max_depth
            )
            return UpgradeThreshold(
                target_edge,
                original_capacity,
                required_capacity,
                ineffective ? Inf : (required_capacity - original_capacity),
                target_flow,
                baseline_result.max_flow,
                false,
                ineffective,
                call_count[1]
            )
        end

        if abs(hi_result.max_flow - prev_flow) <= tol
            # required_capacity=Inf and required_increase=Inf are sentinel
            # values indicating this edge alone cannot achieve target_flow
            # regardless of capacity. Callers must check upgrade_ineffective
            # before using required_capacity or required_increase arithmetically.
            return UpgradeThreshold(
                target_edge,
                original_capacity,
                Inf,
                Inf,
                target_flow,
                baseline_result.max_flow,
                false,
                true,
                call_count[1]
            )
        end

        prev_flow = hi_result.max_flow
        search_cap *= 2.0
        c_hi = original_capacity + search_cap
    end
end

"""
    find_all_degradation_thresholds(edgelist, outgoing_index, incoming_index,
                                    capacities, source_nodes, sink_nodes,
                                    target_flow;
                                    candidate_edges=nothing,
                                    algorithm=:dinic, tol=1e-10, max_depth=64)

Run exact degradation-threshold analysis for each candidate edge. Results are sorted by
ascending `degradation_margin` (smallest margin = most vulnerable edge first).
"""
function find_all_degradation_thresholds(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    target_flow::Float64;
    candidate_edges=nothing,
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    max_depth::Int=64
)::Vector{DegradationThreshold}
    edge_set = Set(edgelist)
    candidates = if candidate_edges === nothing
        [e for e in edgelist if isfinite(capacities[e])]
    else
        converted = Vector{Tuple{Int64,Int64}}(candidate_edges)
        for edge in converted
            edge in edge_set || throw(ArgumentError("candidate edge $edge is not in edgelist."))
        end
        converted
    end

    isempty(candidates) && throw(ArgumentError("No finite-capacity candidate edges available for degradation threshold analysis."))

    results = DegradationThreshold[]
    for edge in candidates
        push!(results, find_degradation_threshold(
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            source_nodes,
            sink_nodes,
            edge,
            target_flow;
            algorithm=algorithm,
            tol=tol,
            max_depth=max_depth
        ))
    end

    sort!(results; by=x -> (x.degradation_margin, x.target_edge))
    return results
end

"""
    analyze_parametric_thresholds(edgelist, outgoing_index, incoming_index,
                                  capacities, source_nodes, sink_nodes,
                                  flow_result;
                                  target_flow=nothing,
                                  candidate_edges=nothing,
                                  algorithm=:dinic, tol=1e-10, max_depth=64)

Run aggregate exact degradation-threshold analysis. If `target_flow` is `nothing`, the
module uses `0.9 * flow_result.max_flow` by default (10% degradation tolerance).
"""
function analyze_parametric_thresholds(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    target_flow=nothing,
    candidate_edges=nothing,
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    max_depth::Int=64
)::ParametricThresholdResult
    _require_bounded_baseline(flow_result)

    selected_target_flow = target_flow === nothing ? 0.9 * flow_result.max_flow : Float64(target_flow)
    selected_target_flow > 0.0 || throw(ArgumentError("target_flow must be positive."))
    selected_target_flow <= flow_result.max_flow + tol || throw(ArgumentError("target_flow=$selected_target_flow exceeds baseline max flow=$(flow_result.max_flow)."))

    thresholds = find_all_degradation_thresholds(
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes,
        selected_target_flow;
        candidate_edges=candidate_edges,
        algorithm=algorithm,
        tol=tol,
        max_depth=max_depth
    )

    return ParametricThresholdResult(thresholds, selected_target_flow, flow_result.max_flow)
end

end
