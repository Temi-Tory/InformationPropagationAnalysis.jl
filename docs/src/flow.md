# [Flow](@id flow)

`Flow` answers the capacity questions about a DAG-structured flow network: how
much can get through, where the bottleneck is, how much a component failure
costs, and how sensitive the throughput is to each edge.

## `analyze_all`

One call runs every analysis and returns a `FlowCapacityResult`:

```
result = analyze_all(edgelist, outgoing_index, incoming_index, capacities,
                     source_nodes, sink_nodes;
                     node_capacities = nothing, algorithm = :dinic,
                     k_failure = 2, path_limit = 10_000, ...)
```

`source_nodes` and `sink_nodes` are `Vector{Int64}`. `capacities` is
`Dict{Tuple{Int64,Int64},Float64}`; pass `node_capacities` to add per-node
throughput limits (solved by node-splitting). The keyword limits
(`path_limit`, `cut_limit`, `combination_limit`) bound the enumerations on dense
graphs.

```@example flow
using InformationPropagationAnalysis
const IPA = InformationPropagationAnalysis

edgelist = Tuple{Int64,Int64}[(1,2), (1,3), (2,4), (3,4)]
outgoing = Dict{Int64,Set{Int64}}(1 => Set([2,3]), 2 => Set([4]), 3 => Set([4]))
incoming = Dict{Int64,Set{Int64}}(2 => Set([1]), 3 => Set([1]), 4 => Set([2,3]))
capacities = Dict{Tuple{Int64,Int64},Float64}((1,2) => 10.0, (2,4) => 5.0,
                                              (1,3) => 7.0,  (3,4) => 8.0)

mf = IPA.Flow.solve_max_flow_dinic(edgelist, outgoing, incoming, capacities,
                                   Int64[1], Int64[4])
(max_flow = mf.max_flow, mincut_capacity = mf.mincut_capacity, cut = mf.mincut_S)
```

The `(2,4)` edge caps the top path at 5 and `(1,3)` caps the bottom at 7, so the
s–t max-flow is 12 — equal to the min-cut capacity, as it must be.

## Individual analyses

Each is also callable on its own, taking the edge list, indices, capacities,
source/sink vectors and a solved `FlowSolveResult`:

| Entry point | Result type | What it gives |
|---|---|---|
| `solve_max_flow_dinic` / `_edmonds_karp` / `_push_relabel` | `FlowSolveResult` | max-flow value, flow dict, the min cut |
| `analyze_structure` | `StructuralResult` | single points of failure, path contributions, edge redundancy |
| `analyze_sensitivity` | `SensitivityResult` | Birnbaum importance, marginal capacity values, critical-edge ranking |
| `analyze_failure_impact` | `FailureImpactResult` | single- and *k*-edge failure throughput loss, degradation scenarios |
| `analyze_min_cuts` | `MinCutAnalysis` | all edges in *some* / *every* minimum cut, cut enumeration |
| `analyze_global_connectivity` | `GlobalConnectivityResult` | edge- and node-connectivity, global min cut |
| `analyze_parametric_thresholds` | `ParametricThresholdResult` | the degradation / upgrade thresholds |
| `analyze_node_capacitated_flow` | `NodeCapacitatedAnalysisResult` | max-flow and SPOF nodes under node capacities |
| `decompose_flow` | `FlowDecomposition` | the max-flow solution split into path-flows |

Granular helpers (`birnbaum_importance`, `identify_spof_nodes`, `enumerate_paths`,
the threshold finders, …) are `public` but not exported — reach them as
`Flow.birnbaum_importance` and friends.

## Docstrings

```@docs
InformationPropagationAnalysis.Flow.analyze_all
InformationPropagationAnalysis.Flow.FlowCapacityResult
```
