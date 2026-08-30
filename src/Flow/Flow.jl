module Flow

include("FlowModule.jl")
using .FlowModule

include("_CapacityShared.jl")
include("CapacityTypes.jl")
using .CapacityTypes

include("SensitivityModule.jl")
using .SensitivityModule

include("FailureImpactModule.jl")
using .FailureImpactModule

include("StructuralModule.jl")
using .StructuralModule

include("FlowDecompositionModule.jl")
using .FlowDecompositionModule

include("ParametricThresholdModule.jl")
using .ParametricThresholdModule

include("NodeCapacitatedFlowModule.jl")
using .NodeCapacitatedFlowModule

include("MinCutUtilitiesModule.jl")
using .MinCutUtilitiesModule

include("GlobalConnectivityModule.jl")
using .GlobalConnectivityModule

# ─────────────────────────────────────────────────────────────────────────────
# Public API
#
# `analyze_all` runs every analysis and returns a `FlowCapacityResult`. The
# per-concern `analyze_*` entry points and the max-flow solvers are the other
# primary calls; their result types and the core min-cut / connectivity queries
# round out the export set. Secondary helpers and low-level record types are
# `public` (importable, documented) but not pulled in by `using`.
# ─────────────────────────────────────────────────────────────────────────────
export analyze_all, FlowCapacityResult,
       # max-flow solvers
       solve_max_flow_dinic, solve_max_flow_edmonds_karp, solve_max_flow_push_relabel,
       FlowSolveResult,
       # per-concern entry points + their result types
       analyze_sensitivity, SensitivityResult,
       analyze_failure_impact, FailureImpactResult,
       analyze_structure, StructuralResult,
       decompose_flow, FlowDecomposition,
       analyze_parametric_thresholds, ParametricThresholdResult,
       analyze_node_capacitated_flow, NodeCapacitatedAnalysisResult,
       analyze_min_cuts, MinCutAnalysis,
       analyze_global_connectivity, GlobalConnectivityResult,
       # core min-cut / connectivity queries
       edges_in_every_mincut, edges_in_some_mincut, mincut_partition,
       minimum_st_cut_edges, minimum_st_cut_capacity,
       edge_connectivity, node_connectivity, global_min_cut

public sink_flows, node_inflow, node_outflow,
       critical_edge_ranking, marginal_capacity_values, birnbaum_importance,
       extract_min_cut_sets, analyze_single_edge_failures, analyze_k_edge_failures,
       analyze_capacity_degradation,
       identify_spof_edges, identify_spof_nodes, enumerate_paths,
       path_flow_contributions, bottleneck_ranking, node_topological_positions,
       edge_redundancy_scores,
       find_degradation_threshold, find_upgrade_threshold, find_all_degradation_thresholds,
       DegradationThreshold, UpgradeThreshold,
       build_node_split_graph, solve_node_capacitated_flow, node_capacitated_spof_nodes,
       NodeCapacitatedFlowResult, NodeSplitGraph,
       enumerate_min_cuts, MinCut, MinCutEnumeration,
       EdgeConnectivityResult, NodeConnectivityResult, GlobalMinCutResult,
       FlowPathComponent
# Internal (reach as `Flow.x`): validate_capacity_constraints, validate_flow_conservation,
# validate_maxflow_mincut, validate_exactness, validate_decomposition, CriticalEdgeRecord,
# SingleEdgeFailureRecord, KEdgeFailureRecord, DegradationScenarioRecord,
# PathFlowContribution, BottleneckRecord

# ─────────────────────────────────────────────────────────────────────────────
# Result struct
# ─────────────────────────────────────────────────────────────────────────────

"""
    FlowCapacityResult

Unified result returned by a single `analyze_all` call.
Every field is produced by a specific sub-module; no algorithm logic lives here.

Fields (in pipeline order):

- `flow`                  (`FlowModule`)               Baseline max-flow solve result.
                                                        All downstream analyses share this single solve.
- `sensitivity`           (`SensitivityModule`)         Critical-edge ranking, marginal capacity values,
                                                        and Birnbaum importance scores.
- `failure_impact`        (`FailureImpactModule`)        Min-cut sets, single-edge failure impacts,
                                                        k-edge failure combos, and capacity-degradation
                                                        scenario results.
- `structure`             (`StructuralModule`)           SPOF edges/nodes, enumerated source-to-sink paths,
                                                        bottleneck ranking, node topological positions,
                                                        and edge redundancy scores.
- `flow_decomposition`    (`FlowDecompositionModule`)   Path-flow decomposition of the baseline flow into
                                                        source-to-sink components.
- `parametric_thresholds` (`ParametricThresholdModule`) Per-edge degradation and upgrade thresholds
                                                        relative to a target flow level.
- `node_capacitated`      (`NodeCapacitatedFlowModule`) Node-capacity-constrained max-flow result.
                                                        `nothing` if `node_capacities` was not provided
                                                        to `analyze_all`.
- `min_cut_analysis`      (`MinCutUtilitiesModule`)     Lattice-based full enumeration of minimum cuts,
                                                        plus edges in some / every minimum cut.
- `global_connectivity`   (`GlobalConnectivityModule`)  Exact edge connectivity, node connectivity,
                                                        and global minimum cut.

Metadata fields:
- `algorithm`             Symbol identifying the max-flow algorithm used throughout (e.g. `:dinic`).
- `tol`                   Numerical tolerance used for all flow comparisons.
- `baseline_max_flow`     `flow.max_flow` stored at top level for quick access without field chaining.
"""
struct FlowCapacityResult
    # Step 1
    flow                     :: FlowSolveResult

    # Step 2
    sensitivity              :: SensitivityResult

    # Step 3
    failure_impact           :: FailureImpactResult

    # Step 4
    structure                :: StructuralResult

    # Step 5
    flow_decomposition       :: FlowDecomposition

    # Step 6
    parametric_thresholds    :: ParametricThresholdResult

    # Step 7 — nothing if node_capacities was not provided to analyze_all
    node_capacitated         :: Union{NodeCapacitatedAnalysisResult, Nothing}

    # Step 8
    min_cut_analysis         :: MinCutAnalysis

    # Step 9
    global_connectivity      :: GlobalConnectivityResult

    # Metadata
    algorithm                :: Symbol
    tol                      :: Float64
    baseline_max_flow        :: Float64  # = flow.max_flow, stored at top level for quick access
end

# ─────────────────────────────────────────────────────────────────────────────
# analyze_all
# ─────────────────────────────────────────────────────────────────────────────

"""
    analyze_all(
        edgelist, outgoing_index, incoming_index,
        capacities, source_nodes, sink_nodes;
        node_capacities        = nothing,
        target_flow            = nothing,
        k_failure              = 2,
        degradation_scenarios  = nothing,
        cut_limit              = 1000,
        path_limit             = 10_000,
        redundancy_candidates  = nothing,
        combination_limit      = 10_000,
        algorithm              = :dinic,
        tol                    = 1e-10,
        max_depth              = 64
    ) -> FlowCapacityResult

Run the complete capacity-analysis pipeline in a single call.

The baseline max-flow solve (Step 1) executes exactly once. Every downstream
analysis step receives the same `FlowSolveResult` and does not re-solve the
baseline independently. Steps that require additional solves for their own
algorithmic purposes (e.g. ParametricThresholdModule performing per-edge binary
search, SensitivityModule probing zero-capacity perturbations) do so internally.
All exactness guarantees from individual modules are preserved — no approximations
are introduced by this aggregator.

# Positional arguments
- `edgelist`        : `Vector{Tuple{Int64,Int64}}` — directed edge list.
- `outgoing_index`  : `Dict{Int64,Set{Int64}}` — adjacency map (u => neighbours).
- `incoming_index`  : `Dict{Int64,Set{Int64}}` — reverse adjacency map (v => predecessors).
- `capacities`      : `Dict{Tuple{Int64,Int64},Float64}` — per-edge capacity.
- `source_nodes`    : `Vector{Int64}` — source nodes.
- `sink_nodes`      : `Vector{Int64}` — sink nodes.

# Keyword arguments
- `node_capacities`         : `Union{Dict{Int64,Float64}, Nothing}` (default `nothing`).
  Per-node capacity constraints. If provided, Step 7 (node-capacitated flow analysis)
  is executed; `result.node_capacitated` will be a `NodeCapacitatedAnalysisResult`.
  If `nothing`, Step 7 is skipped and `result.node_capacitated` is `nothing`.
- `target_flow`             : `Union{Float64, Nothing}` (default `nothing`).
  Target flow level for Step 6 (parametric threshold analysis).
  If `nothing`, ParametricThresholdModule defaults to `0.9 * max_flow`.
- `k_failure`               : `Int` (default `2`).
  Maximum combination size for k-edge failure analysis in Step 3.
- `degradation_scenarios`   : `Union{Vector, Nothing}` (default `nothing`).
  Capacity-degradation scenarios for Step 3. If `nothing`,
  `result.failure_impact.degradation_results` will be an empty vector.
- `cut_limit`               : `Int` (default `1000`).
  Maximum number of minimum cuts to enumerate in Step 8.
- `path_limit`              : `Int` (default `10_000`).
  Maximum number of source-to-sink paths to enumerate in Step 4.
- `redundancy_candidates`   : (default `nothing`).
  Candidate edges for redundancy scoring in Step 4. If `nothing`, all edges
  are considered candidates.
- `combination_limit`       : `Int` (default `10_000`).
  Passed to FailureImpactModule (k-edge combos) and StructuralModule.
- `algorithm`               : `Symbol` (default `:dinic`).
  Max-flow algorithm for all solver calls. Accepted: `:dinic`, `:edmonds_karp`,
  `:push_relabel`. Propagated to every sub-module that accepts it.
- `tol`                     : `Float64` (default `1e-10`).
  Numerical tolerance for flow comparisons, propagated throughout.
- `max_depth`               : `Int` (default `64`).
  Maximum binary-search depth in Step 6 (ParametricThresholdModule).

# Pipeline steps
1. **FlowModule**               — baseline max-flow solve.
2. **SensitivityModule**        — `analyze_sensitivity`.
3. **FailureImpactModule**      — `analyze_failure_impact`.
4. **StructuralModule**         — `analyze_structure`.
5. **FlowDecompositionModule**  — `decompose_flow`.
6. **ParametricThresholdModule**— `analyze_parametric_thresholds`.
7. **NodeCapacitatedFlowModule**— `analyze_node_capacitated_flow` (conditional on `node_capacities`).
8. **MinCutUtilitiesModule**    — `analyze_min_cuts`.
9. **GlobalConnectivityModule** — `analyze_global_connectivity`.

# Errors
- Throws `ArgumentError` if `edgelist`, `source_nodes`, `sink_nodes`, or `capacities` are empty.
- Throws `ArgumentError` with the message
  "analyze_all cannot proceed: baseline max-flow is unbounded. Check edge capacities for infinite paths."
  if the baseline solve is unbounded. No downstream steps are executed in that case.
- Any error thrown by a sub-module propagates unchanged with its original message.
"""
function analyze_all(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    node_capacities::Union{Dict{Int64,Float64}, Nothing}=nothing,
    target_flow::Union{Float64, Nothing}=nothing,
    k_failure::Int=2,
    degradation_scenarios=nothing,
    cut_limit::Int=1000,
    path_limit::Int=10_000,
    redundancy_candidates=nothing,
    combination_limit::Int=10_000,
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    max_depth::Int=64
)::FlowCapacityResult
    # ── Input validation ─────────────────────────────────────────────────────
    isempty(edgelist)      && throw(ArgumentError("edgelist must be non-empty."))
    isempty(source_nodes)  && throw(ArgumentError("source_nodes must be non-empty."))
    isempty(sink_nodes)    && throw(ArgumentError("sink_nodes must be non-empty."))
    isempty(capacities)    && throw(ArgumentError("capacities must be non-empty."))

    # ── Step 1: Baseline max-flow solve ──────────────────────────────────────
    flow = _solve_with_algorithm(
        algorithm,
        edgelist, outgoing_index, incoming_index, capacities,
        source_nodes, sink_nodes;
        tol=tol,
        validate=true
    )

    flow.is_unbounded && throw(ArgumentError(
        "analyze_all cannot proceed: baseline max-flow is unbounded. " *
        "Check edge capacities for infinite paths."
    ))

    # ── Step 2: Sensitivity analysis ─────────────────────────────────────────
    sensitivity = analyze_sensitivity(
        edgelist, outgoing_index, incoming_index, capacities,
      source_nodes, sink_nodes, flow;
        algorithm=algorithm,
        tol=tol
    )

    # ── Step 3: Failure impact analysis ──────────────────────────────────────
    failure_impact = analyze_failure_impact(
        edgelist, outgoing_index, incoming_index, capacities,
        source_nodes, sink_nodes, flow;
        k=k_failure,
        scenarios=degradation_scenarios,
        algorithm=algorithm,
        tol=tol,
        combination_limit=combination_limit
    )

    # ── Step 4: Structural analysis ───────────────────────────────────────────
    structure = analyze_structure(
        edgelist, outgoing_index, incoming_index, capacities,
        source_nodes, sink_nodes, flow;
        algorithm=algorithm,
        tol=tol,
        path_limit=path_limit,
        redundancy_candidates=redundancy_candidates,
        combination_limit=combination_limit
    )

    # ── Step 5: Flow decomposition ────────────────────────────────────────────
    flow_decomposition = decompose_flow(
        edgelist, source_nodes, sink_nodes, flow;
        tol=tol
    )

    # ── Step 6: Parametric threshold analysis ─────────────────────────────────
    parametric_thresholds = analyze_parametric_thresholds(
        edgelist, outgoing_index, incoming_index, capacities,
        source_nodes, sink_nodes, flow;
        target_flow=target_flow,
        algorithm=algorithm,
        tol=tol,
        max_depth=max_depth
    )

    # ── Step 7: Node-capacitated flow (conditional) ───────────────────────────
    node_capacitated = if node_capacities !== nothing
        analyze_node_capacitated_flow(
            edgelist, outgoing_index, incoming_index, capacities,
            source_nodes, sink_nodes, node_capacities;
            algorithm=algorithm,
            tol=tol
        )
    else
        nothing
    end

    # ── Step 8: Min-cut analysis ──────────────────────────────────────────────
    min_cut_analysis = analyze_min_cuts(
        edgelist, outgoing_index, incoming_index, capacities,
        source_nodes, sink_nodes, flow;
        cut_limit=cut_limit,
        tol=tol
    )

    # ── Step 9: Global connectivity ───────────────────────────────────────────
    global_connectivity = analyze_global_connectivity(
        edgelist, outgoing_index, incoming_index, capacities,
        source_nodes, sink_nodes;
        algorithm=algorithm,
        tol=tol
    )

    return FlowCapacityResult(
        flow,
        sensitivity,
        failure_impact,
        structure,
        flow_decomposition,
        parametric_thresholds,
        node_capacitated,
        min_cut_analysis,
        global_connectivity,
        algorithm,
        tol,
        flow.max_flow
    )
end

end # module Flow
