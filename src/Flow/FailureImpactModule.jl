module FailureImpactModule

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

export FailureImpactResult,
       extract_min_cut_sets,
       analyze_single_edge_failures,
       analyze_k_edge_failures,
       analyze_capacity_degradation,
       analyze_failure_impact

struct FailureImpactResult
    min_cut_edges::Vector{Tuple{Int64,Int64}}
    single_edge_failures::Vector{SingleEdgeFailureRecord}
    k_edge_failures::Vector{KEdgeFailureRecord}
    degradation_results::Vector{DegradationScenarioRecord}
end

function _drop_from_baseline(
    baseline_flow::Float64,
    perturbed_result::FlowSolveResult
)::Tuple{Float64,Bool}
    if perturbed_result.is_unbounded
        return -Inf, true
    end
    return baseline_flow - perturbed_result.max_flow, false
end

function _build_scenario_capacities(
    capacities::Dict{Tuple{Int64,Int64},Float64},
    scenario
)::Dict{Tuple{Int64,Int64},Float64}
    modified = _copy_capacities(capacities)

    if scenario isa AbstractDict
        for (edge, raw_value) in scenario
            edge isa Tuple{Int64,Int64} || throw(ArgumentError("Scenario override edge key must be Tuple{Int64,Int64}; got $(typeof(edge))."))
            haskey(modified, edge) || throw(ArgumentError("Scenario override references unknown edge $edge."))
            value = Float64(raw_value)
            (isnan(value) || value < 0.0) && throw(ArgumentError("Scenario override for edge $edge must be nonnegative and not NaN."))
            modified[edge] = value
        end
        return modified
    elseif scenario isa Real
        factor = Float64(scenario)
        (isnan(factor) || factor < 0.0) && throw(ArgumentError("Scenario scale factor must be nonnegative and not NaN."))
        for edge in keys(modified)
            cap = modified[edge]
            if isfinite(cap)
                modified[edge] = factor * cap
            end
        end
        return modified
    else
        throw(ArgumentError("Invalid scenario type $(typeof(scenario)). Expected an edge-capacity dictionary or a numeric scale factor."))
    end
end

function _index_combinations(n::Int, k::Int)::Vector{Vector{Int}}
    if k == 0
        return [Int[]]
    end
    result = Vector{Vector{Int}}()
    current = Vector{Int}(undef, k)

    function backtrack(start::Int, depth::Int)
        if depth > k
            push!(result, copy(current))
            return
        end

        max_i = n - (k - depth)
        for i in start:max_i
            current[depth] = i
            backtrack(i + 1, depth + 1)
        end
    end

    backtrack(1, 1)
    return result
end

"""
    extract_min_cut_sets(edgelist, flow_result)

Return the representative minimum-cut edge set induced by the solved partition
`(mincut_S, mincut_T)` from `flow_result`, i.e. all original edges `(u,v)` with
`u ∈ mincut_S` and `v ∈ mincut_T`.

This returns one valid solved minimum cut and does not enumerate all minimum cuts.
"""
function extract_min_cut_sets(
    edgelist::Vector{Tuple{Int64,Int64}},
    flow_result::FlowSolveResult
)::Vector{Tuple{Int64,Int64}}
    _require_bounded_baseline(flow_result)
    edges = Tuple{Int64,Int64}[]
    for edge in edgelist
        u, v = edge
        if (u in flow_result.mincut_S) && (v in flow_result.mincut_T)
            push!(edges, edge)
        end
    end
    sort!(edges)
    return edges
end

function extract_min_cut_sets(
    flow_result::FlowSolveResult
)::Vector{Tuple{Int64,Int64}}
    _require_bounded_baseline(flow_result)
    deterministic_edgelist = collect(keys(flow_result.flow))
    sort!(deterministic_edgelist)
    return extract_min_cut_sets(deterministic_edgelist, flow_result)
end

"""
    analyze_single_edge_failures(...; algorithm=:dinic, tol=1e-10)

For each edge in the exact candidate set of edges that appear in some minimum cut,
set its capacity to zero, rerun exact max flow with validation, and report the
flow impact sorted deterministically by descending drop then edge lexicographic order.
"""
function analyze_single_edge_failures(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10
)::Vector{SingleEdgeFailureRecord}
    _require_bounded_baseline(flow_result)
    baseline_flow = flow_result.max_flow
    candidates = _edges_in_some_mincut(
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        flow_result;
        tol=tol
    )

    results = SingleEdgeFailureRecord[]
    for edge in candidates
        modified = _copy_capacities(capacities)
        modified[edge] = 0.0
        perturbed = _solve_with_algorithm(
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
        drop, is_unbounded = _drop_from_baseline(baseline_flow, perturbed)
        push!(results, SingleEdgeFailureRecord(
            edge,
            baseline_flow,
            perturbed.max_flow,
            drop,
            drop > 0.0,
            is_unbounded
        ))
    end

    sort!(results; by=x -> (-x.drop, x.edge))
    return results
end

"""
    analyze_k_edge_failures(...; k=2, algorithm=:dinic, tol=1e-10, combination_limit=10_000)

Enumerate all exact `k`-edge removals from the theorem-safe candidate set of edges
that appear in some minimum cut. Each combination sets all selected capacities to
zero simultaneously and reruns exact max flow with validation.
"""
function analyze_k_edge_failures(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    k::Int=2,
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    combination_limit::Int=10_000
)::Vector{KEdgeFailureRecord}
    _require_bounded_baseline(flow_result)
    k < 0 && throw(ArgumentError("k must be nonnegative."))
    if k == 0
        return KEdgeFailureRecord[]
    end

    candidates = _edges_in_some_mincut(
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        flow_result;
        tol=tol
    )
    n = length(candidates)
    if k > n
        return KEdgeFailureRecord[]
    end

    total_combinations = binomial(n, k)
    if total_combinations > combination_limit
        throw(ArgumentError("k-edge failure analysis would evaluate $total_combinations combinations (C($n,$k)), exceeding combination_limit=$combination_limit. Reduce k or use a smaller candidate set."))
    end

    baseline_flow = flow_result.max_flow
    staged = KEdgeFailureRecord[]

    for idxs in _index_combinations(n, k)
        combo_edges = [candidates[i] for i in idxs]
        combo_key = Tuple(combo_edges)
        modified = _copy_capacities(capacities)
        for edge in combo_edges
            modified[edge] = 0.0
        end

        perturbed = _solve_with_algorithm(
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
        drop, is_unbounded = _drop_from_baseline(baseline_flow, perturbed)
        push!(staged, KEdgeFailureRecord(
            combo_key,
            baseline_flow,
            perturbed.max_flow,
            drop,
            is_unbounded
        ))
    end

    sort!(staged; by=x -> (-x.drop, x.edges))
    return staged
end

"""
    analyze_capacity_degradation(...; scenarios, algorithm=:dinic, tol=1e-10)

Run exact max-flow analyses for each degradation scenario.
Scenarios may be either explicit edge-capacity override dictionaries or Float64
uniform scale factors applied to finite capacities.
"""
function analyze_capacity_degradation(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    scenarios,
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    baseline_result::Union{Nothing,FlowSolveResult}=nothing
)::Vector{DegradationScenarioRecord}
    baseline = baseline_result === nothing ? _solve_with_algorithm(
        algorithm,
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes;
        tol=tol,
        validate=true
    ) : baseline_result
    _require_bounded_baseline(baseline)
    baseline_flow = baseline.max_flow

    results = DegradationScenarioRecord[]
    for (idx, scenario) in enumerate(scenarios)
        scenario_capacities = _build_scenario_capacities(capacities, scenario)
        rerun = _solve_with_algorithm(
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            scenario_capacities,
            source_nodes,
            sink_nodes;
            tol=tol,
            validate=true
        )
        drop, is_unbounded = _drop_from_baseline(baseline_flow, rerun)
        push!(results, DegradationScenarioRecord(
            idx,
            scenario_capacities,
            rerun.max_flow,
            rerun.sink_flow,
            rerun.saturated_edges,
            drop,
            is_unbounded
        ))
    end

    sort!(results; by=x -> x.scenario_id)
    return results
end

"""
    analyze_failure_impact(...; k=2, scenarios=nothing, algorithm=:dinic, tol=1e-10)

Run the complete exact failure-impact analysis pipeline and return a typed aggregate
result with representative min-cut edges, single-edge failures, k-edge failures,
and optional capacity degradation analysis.
"""
function analyze_failure_impact(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    k::Int=2,
    scenarios=nothing,
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    combination_limit::Int=10_000
)::FailureImpactResult
    min_cut_edges = extract_min_cut_sets(edgelist, flow_result)
    single_edge_failures = analyze_single_edge_failures(
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes,
        flow_result;
        algorithm=algorithm,
        tol=tol
    )

    k_edge_failures = if k == 0
        KEdgeFailureRecord[]
    else
        analyze_k_edge_failures(
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            source_nodes,
            sink_nodes,
            flow_result;
            k=k,
            algorithm=algorithm,
            tol=tol,
            combination_limit=combination_limit
        )
    end

    degradation_results = if scenarios === nothing
        DegradationScenarioRecord[]
    else
        analyze_capacity_degradation(
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            source_nodes,
            sink_nodes;
            scenarios=scenarios,
            algorithm=algorithm,
            tol=tol,
            baseline_result=flow_result
        )
    end

    return FailureImpactResult(min_cut_edges, single_edge_failures, k_edge_failures, degradation_results)
end

end
