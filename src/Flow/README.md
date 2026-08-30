# CapacityAnalysisKit — Exact Network Flow Analysis

You have a network. You need to know: how much can flow through it, where will it fail, how vulnerable is it, and what happens when components degrade. CapacityAnalysisKit answers those questions exactly.

Designed for: DAG-structured flow networks (directed acyclic graphs with clear source-to-sink directionality), which is the dominant modeling pattern for infrastructure reliability systems such as supply chains, distribution networks, pipeline systems, communication backbones, and transportation corridors.

---

## Section 1: What Questions Does This Framework Answer?

This section is organized by engineering concern, not by module boundaries.

### 1.1 Flow and Capacity

**Q1. What is the maximum throughput of this network?**
- Use: `analyze_all(...)`
- Read: `result.flow.max_flow` or `result.baseline_max_flow`

**Q2. How much flow reaches each sink?**
- Use: `analyze_all(...)`
- Read: `result.flow.sink_flow`

**Q3. Which paths actually carry flow and how much?**
- Use: `analyze_all(...)` or `decompose_flow(...)`
- Read: `result.flow_decomposition.components`

### 1.2 Bottlenecks and Constraints

**Q1. Where is flow most constrained?**
- Use: `analyze_all(...)`
- Read: `result.flow.mincut_S`, `result.flow.mincut_T`, `result.structure.bottleneck_ranking`

**Q2. Which edges are fully saturated?**
- Use: `analyze_all(...)`
- Read: `result.flow.saturated_edges`

**Q3. What is the weakest point in the network?**
- Use: `analyze_all(...)`
- Read: `result.min_cut_analysis.representative_cut`, `result.structure.spof_edges`, `result.structure.spof_nodes`

### 1.3 Single Points of Failure

**Q1. Which components, if removed, stop all flow?**
- Use: `analyze_all(...)`
- Read: `result.structure.spof_edges`, `result.structure.spof_nodes`

**Q2. Which edges appear in every possible failure mode (minimum-cut mode)?**
- Use: `analyze_all(...)` or `analyze_min_cuts(...)`
- Read: `result.min_cut_analysis.edges_in_every_cut`

**Q3. Which nodes are structural chokepoints?**
- Use: `analyze_all(...; node_capacities=...)`
- Read: `result.node_capacitated.spof_nodes`

### 1.4 Resilience and Degradation

**Q1. How much can edge `e` degrade before throughput drops below target `T`?**
- Use: `find_degradation_threshold(...)` or `analyze_parametric_thresholds(...)`
- Read: `DegradationThreshold.degradation_margin`, `DegradationThreshold.threshold_capacity`

**Q2. What is the resilience margin of each component?**
- Use: `analyze_all(...)`
- Read: `result.parametric_thresholds.degradation_thresholds`

**Q3. What happens if multiple components fail together?**
- Use: `analyze_failure_impact(...; k=...)`, `analyze_capacity_degradation(...; scenarios=...)`
- Read: `result.failure_impact.k_edge_failures`, `result.failure_impact.degradation_results`

### 1.5 Upgrade and Investment Planning

**Q1. What capacity increase on edge `e` achieves target throughput `T`?**
- Use: `find_upgrade_threshold(...)`
- Read: `UpgradeThreshold.required_increase`, `UpgradeThreshold.required_capacity`

**Q2. Which edges have highest marginal value per unit capacity?**
- Use: `analyze_sensitivity(...)`
- Read: `result.sensitivity.marginal_capacity`

**Q3. Which upgrade gives most flow per unit cost?**
- Status: not in current scope (requires explicit cost model / min-cost flow layer).
- Planned direction: Phase 4 (min-cost max-flow integration).

### 1.6 Failure Mode Enumeration

**Q1. What are all ways this system can fail at maximum capacity (minimum cuts)?**
- Use: `enumerate_min_cuts(...)` or `analyze_all(...)`
- Read: `result.min_cut_analysis.enumeration.cuts`

**Q2. How many distinct failure modes exist?**
- Use: `enumerate_min_cuts(...)`
- Read: `MinCutEnumeration.total_cuts`, `MinCutEnumeration.is_complete`

**Practical note on examples**
- Multiple minimum cuts require a topology that actually induces a nonempty free zone $F = S^{**} \setminus S^*$.
- Rich real-looking networks often still have a unique minimum cut, even when they show strong sensitivity and nontrivial multi-edge failures.
- In the example data:
    - `network.edges` is the small end-to-end sanity example,
    - `network_complex.edges` is the richer analysis example,
    - `network_lattice.edges` is the dedicated minimum-cut lattice enumeration example.

**Q3. What is the smallest edge set whose failure disables max throughput?**
- Use: `minimum_st_cut_edges(...)`, `minimum_st_cut_capacity(...)`
- Read: cut edges and cut capacity outputs

### 1.7 Network Redundancy

**Q1. How many independent paths exist between source and sink regions?**
- Use: `edge_redundancy_scores(...)`, `edge_connectivity(...)`
- Read: `result.structure.edge_redundancy`, `result.global_connectivity.edge_connectivity.lambda`

**Q2. How many edges must fail to disconnect the network?**
- Use: `edge_connectivity(...)`
- Read: `EdgeConnectivityResult.lambda`

**Q3. How many nodes must fail to disconnect the network?**
- Use: `node_connectivity(...)`
- Read: `NodeConnectivityResult.kappa`

**Q4. Why do `lambda` and `kappa` sometimes return 0 on DAG models? Is that a bug?**
- Usually not a bug.
- In a directed DAG, many node pairs are not mutually reachable (the graph is not strongly connected by design).
- Global directed connectivity scans all relevant node-pair cuts; if any directed pair has no path, the corresponding cut value is 0, so global `lambda`/`kappa` can be 0.
- Interpretation: this is mathematically correct for global directed connectivity on DAGs.
- For DAG reliability studies, the most meaningful resilience signal is usually source-to-sink throughput and cut structure (max-flow/min-cut, SPOFs, failure impact), not strong-connectivity-style global metrics.

### 1.8 Component Importance

**Q1. Which components matter most to performance?**
- Use: `analyze_sensitivity(...)`
- Read: `result.sensitivity.critical_edges`, `result.sensitivity.birnbaum`

**Q2. What is Birnbaum importance of each edge?**
- Use: `birnbaum_importance(...)` or `analyze_sensitivity(...)`
- Read: `result.sensitivity.birnbaum`

**Q3. Which one-unit edge upgrade improves flow most?**
- Use: `marginal_capacity_values(...)`
- Read: `result.sensitivity.marginal_capacity`

---

## Section 2: What This Framework Does NOT Answer

Clear scope boundaries are essential for reliability-engineering credibility.

**Q: What is the probability throughput exceeds `T`?**
- Not answered here.
- Requires probabilistic reliability modeling.
- Recommended workflow: use deterministic outputs from this framework as inputs/features to a separate probabilistic model.

**Q: What is optimal routing under cost constraints?**
- Not answered in current version.
- Requires min-cost max-flow / cost-aware optimization.
- Planned direction: Phase 4.

**Q: What happens as demand/capacity changes over time?**
- Not answered here.
- Requires dynamic or time-expanded flow modeling.

**Q: What if capacities are uncertain intervals/distributions?**
- Not directly modeled.
- Recommended deterministic envelope workflow:
  - run best-case capacities,
  - run worst-case capacities,
  - compare exact outcomes for each scenario.

**Q: My network has directed cycles. Can I still use this?**
- Core max-flow solving remains valid on general directed graphs.
- This framework is intentionally designed for DAG-structured reliability models.
- Path enumeration and flow decomposition modules assume acyclicity; with cycles, those analyses are not guaranteed to behave as intended.
- For cyclic systems, either (1) transform into an acyclic stage-expanded model, or (2) use only modules that do not rely on DAG structure.

---

## Section 3: Quick Start

Minimal complete usage example (build inputs, run full pipeline, read key outputs).

```julia
include("CapacityAnalysisKit.jl")
using .CapacityAnalysisKit

edgelist = Tuple{Int64,Int64}[(1,2),(1,3),(2,4),(3,4)]
capacities = Dict{Tuple{Int64,Int64},Float64}((1,2)=>5.0,(1,3)=>3.0,(2,4)=>4.0,(3,4)=>4.0)

outgoing_index = Dict{Int64,Set{Int64}}()
incoming_index = Dict{Int64,Set{Int64}}()
for (u,v) in edgelist
    push!(get!(outgoing_index, u, Set{Int64}()), v)
    push!(get!(incoming_index, v, Set{Int64}()), u)
end

result = analyze_all(edgelist, outgoing_index, incoming_index, capacities, Int64[1], Int64[4])
println("max flow = ", result.baseline_max_flow)
println("bottleneck edges = ", result.flow.saturated_edges)
println("edge connectivity λ = ", result.global_connectivity.edge_connectivity.lambda)
```

---

## Section 4: Mathematical Foundations

For each theorem/result: formal statement, guarantee in framework, using modules.

### 4.1 Max-Flow Min-Cut Theorem
**Formal statement**
$$
\max\;\text{flow}(s\to t)=\min\;\text{cut\_capacity}(s,t)
$$

**Framework guarantee**
- Global optimality of every flow solve.
- Correctness of min-cut-derived outputs.

**Using module(s)**
- `FlowModule`, `FailureImpactModule`, `MinCutUtilitiesModule`, `GlobalConnectivityModule`, `CapacityAnalysisKit`.

### 4.2 Flow Conservation
**Formal statement**
For each non-terminal node $v$:
$$
\sum_u f(u,v)=\sum_w f(v,w)
$$

**Framework guarantee**
- Feasibility of returned flows.
- Correctness basis for decomposition and residual reasoning.

**Using module(s)**
- `FlowModule`, `FlowDecompositionModule`, `NodeCapacitatedFlowModule`.

### 4.3 Integrality Theorem
**Formal statement**
If capacities are integers, an optimal max-flow exists with integer edge flows.

**Framework guarantee**
- Integer `lambda`/`kappa` in unit/integer capacity connectivity constructions.
- Integer component values in decomposition under integer-flow instances.

**Using module(s)**
- `FlowModule`, `GlobalConnectivityModule`, `FlowDecompositionModule`.

### 4.4 Menger's Theorem (edge version)
**Formal statement**
Maximum number of edge-disjoint $s\to t$ paths equals minimum $s$-$t$ edge-cut size.
With unit capacities:
$$
\text{max-flow}(s,t)=\#\text{edge-disjoint paths}(s,t)
$$

**Framework guarantee**
- Correctness of edge redundancy and edge connectivity metrics.

**Using module(s)**
- `StructuralModule` (`edge_redundancy_scores`), `GlobalConnectivityModule` (`edge_connectivity`).

### 4.5 Menger's Theorem (node version)
**Formal statement**
Maximum number of internally node-disjoint $s\to t$ paths equals minimum internal node-cut size.

**Framework guarantee**
- Correctness of node connectivity.
- Correctness basis for node SPOF analysis under node constraints.

**Using module(s)**
- `GlobalConnectivityModule` (`node_connectivity`), `NodeCapacitatedFlowModule` (`node_capacitated_spof_nodes`).

### 4.6 Min-Cut Lattice Structure
**Formal statement**
Let:
- $S^*$ = forward residual BFS source-reachable set,
- $S^{**}$ = complement of backward residual BFS sink-reachable set.

Every minimum-cut source side $S$ satisfies:
$$
S^*\subseteq S\subseteq S^{**}
$$

Characterizations:
- Edge in **some** minimum cut iff saturated and $u\in S^{**}$ and $v\notin S^*$.
- Edge in **every** minimum cut iff saturated and $u\in S^*$ and $v\in T^{**}$.

Free zone:
$$
F=S^{**}\setminus S^*
$$

**Framework guarantee**
- Correctness of `edges_in_some_mincut`, `edges_in_every_mincut`, `enumerate_min_cuts`.
- Correctness basis for cut-critical structural edge identification.

**Using module(s)**
- `MinCutUtilitiesModule`, `StructuralModule` (`identify_spof_edges`).

### 4.7 Node-Splitting Bijection
**Formal statement**
Node capacity $c(v)$ is enforced by replacing node $v$ with edge $v_{in}\to v_{out}$ capacity $c(v)$, preserving feasible-flow correspondence and objective value.

**Framework guarantee**
- Exact equivalence between node-capacitated and split edge-capacitated formulations.

**Using module(s)**
- `NodeCapacitatedFlowModule`.

### 4.8 Flow Decomposition Theorem
**Formal statement**
Any feasible conserved flow decomposes into positive source-to-sink path components (not unique), summing exactly to total flow.

**Framework guarantee**
- Correctness of `decompose_flow` components and totals.

**Using module(s)**
- `FlowDecompositionModule`.

### 4.9 Piecewise-Linear Parametric Structure
**Formal statement**
For single-edge capacity parameter $c_e$, $F(c_e)=\text{max-flow}(c_e)$ is monotone non-decreasing, concave, piecewise linear with finite breakpoints.

**Framework guarantee**
- Exact threshold results in degradation/upgrade analyses.

**Using module(s)**
- `ParametricThresholdModule` (`find_degradation_threshold`, `find_upgrade_threshold`).

---

## Section 5: Input Contract

### 5.1 `edgelist :: Vector{Tuple{Int64,Int64}}`
- **Type:** `Vector{Tuple{Int64,Int64}}`
- **Must contain:** directed edges `(u,v)`.
- **Valid values:** non-empty list of node-ID pairs.
- **Constraints:** duplicate edges unsupported; self-loops unsupported.
- **Design target:** DAG-structured networks (directed acyclic graphs).
- **Why this is intentional:** acyclicity guarantees path enumeration termination, decomposition termination (no flow-cycles to unwind), well-defined topological ordering, and tractable structural analysis behavior.
- **Cyclic-graph note:** max-flow solve remains mathematically valid, but DAG-dependent analyses (especially path enumeration/decomposition) are outside intended scope.

### 5.2 `outgoing_index :: Dict{Int64, Set{Int64}}`
- **Type:** `Dict{Int64,Set{Int64}}`
- **Must contain:** direct successors per node.
- **Coverage:** should include all graph nodes.
- **Missing keys:** treated as empty outgoing set by `get(..., Set())` patterns.

Build loop:
```julia
outgoing = Dict{Int64,Set{Int64}}()
for (u,v) in edgelist
    push!(get!(outgoing, u, Set{Int64}()), v)
end
```

### 5.3 `incoming_index :: Dict{Int64, Set{Int64}}`
- **Type:** `Dict{Int64,Set{Int64}}`
- **Must contain:** direct predecessors per node.
- **Coverage:** should include all graph nodes.

Build loop:
```julia
incoming = Dict{Int64,Set{Int64}}()
for (u,v) in edgelist
    push!(get!(incoming, v, Set{Int64}()), u)
end
```

### 5.4 `capacities :: Dict{Tuple{Int64,Int64}, Float64}`
- **Type:** `Dict{Tuple{Int64,Int64},Float64}`
- **Must contain:** capacity for every edge in `edgelist`.
- **Valid values:** finite `Float64 >= 0` or `Inf`.
- **Invalid values:** `NaN` or negative (throws `ArgumentError`).
- **Meaning of `Inf`:** unconstrained edge capacity.
- **Min-cut implication:** `Inf` edges cannot belong to finite minimum cuts.

### 5.5 `source_nodes :: Vector{Int64}`
- **Type:** `Vector{Int64}`
- **Must contain:** source IDs present in graph.
- **Constraint:** non-empty.
- **Handling:** multi-source via internal super-source.

### 5.6 `sink_nodes :: Vector{Int64}`
- **Type:** `Vector{Int64}`
- **Must contain:** sink IDs present in graph.
- **Constraint:** non-empty.
- **Handling:** multi-sink via internal super-sink.

### 5.7 `node_capacities :: Dict{Int64, Float64}` (optional)
- **Type:** `Dict{Int64,Float64}`
- **Must contain:** constrained-node capacities.
- **Constraints:** keys in graph, values finite/nonnegative/non-NaN.
- **Absent nodes:** treated as unconstrained.
- **In `analyze_all`:** if `nothing`, node-capacitated stage is skipped.

### Complete minimal build example
```julia
edgelist = Tuple{Int64,Int64}[(1,2),(1,3),(2,4),(3,4)]

outgoing_index = Dict{Int64,Set{Int64}}()
incoming_index = Dict{Int64,Set{Int64}}()
for (u,v) in edgelist
    push!(get!(outgoing_index, u, Set{Int64}()), v)
    push!(get!(incoming_index, v, Set{Int64}()), u)
end

capacities = Dict{Tuple{Int64,Int64},Float64}(
    (1,2)=>5.0,
    (1,3)=>3.0,
    (2,4)=>4.0,
    (3,4)=>4.0,
)

source_nodes = Int64[1]
sink_nodes   = Int64[4]
```

---

## Section 6: Module Reference

### 6.1 FlowModule

| Field | Content |
|---|---|
| Purpose | Solve exact max-flow/min-cut and validate flow/cut correctness. |
| Algorithm | Edmonds-Karp, Dinic, Push-Relabel. |
| Theorem | Max-Flow Min-Cut, Flow Conservation, Integrality. |
| Solver calls | 1 solve per `solve_max_flow_*` call. |
| Key inputs | Standard contract only. |
| Key outputs | `FlowSolveResult` (value, cut, residual, saturation, boundedness). |
| Example call | `r = solve_max_flow_dinic(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes)` |

### 6.2 SensitivityModule

| Field | Content |
|---|---|
| Purpose | Rank critical edges and quantify capacity-value sensitivity. |
| Algorithm | Exact perturb-and-resolve against baseline flow. |
| Theorem | Max-Flow Min-Cut, lattice saturation logic. |
| Solver calls | Up to `|saturated_edges|` zero-capacity reruns + up to `|min_cut_candidates|` Birnbaum reruns + up to `|saturated_edges|` marginal-delta reruns; with `|saturated_edges|<=|E|`, `|min_cut_candidates|<=|E|`, overall `O(E)`. |
| Key inputs | Baseline `flow_result` from `FlowModule`. |
| Key outputs | `SensitivityResult` (`critical_edges`, `marginal_capacity`, `birnbaum`). |
| Example call | `s = analyze_sensitivity(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes, flow_result)` |

### 6.3 FailureImpactModule

| Field | Content |
|---|---|
| Purpose | Quantify throughput impact of failures/degradation scenarios. |
| Algorithm | Exact reruns on perturbed capacity/failure configurations. |
| Theorem | Max-Flow Min-Cut and solver-level feasibility guarantees. |
| Solver calls | Single-edge: up to `|min_cut_candidates|` with `|min_cut_candidates|<=|E|`; k-edge: combinatorial in `k` (bounded by `combination_limit`); degradation: one rerun per scenario entry. |
| Key inputs | Baseline `flow_result`; optional `k`, `scenarios`, `combination_limit`. |
| Key outputs | `FailureImpactResult` (single-edge, k-edge, degradation outcomes). |
| Example call | `f = analyze_failure_impact(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes, flow_result; k=2)` |

### 6.4 StructuralModule

| Field | Content |
|---|---|
| Purpose | Identify bottlenecks, SPOFs, path/topology diagnostics, redundancy. |
| Algorithm | Traversal + exact reruns for redundancy scoring. |
| Theorem | Menger (edge), min-cut lattice logic, max-flow optimality. |
| Solver calls | SPOF/path/bottleneck parts are traversal-based; redundancy stage up to `O(|candidate_edges|)` reruns. |
| Key inputs | Baseline `flow_result`; optional `path_limit`, `redundancy_candidates`. |
| Key outputs | `StructuralResult` (`spof_edges`, `spof_nodes`, `paths`, `bottleneck_ranking`, `edge_redundancy`, etc.). |
| Example call | `st = analyze_structure(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes, flow_result; path_limit=10_000)` |

### 6.5 FlowDecompositionModule

| Field | Content |
|---|---|
| Purpose | Decompose solved flow into explicit positive path components. |
| Algorithm | Deterministic path extraction/subtraction on positive-flow edges. |
| Theorem | Flow Decomposition Theorem + Flow Conservation. |
| Solver calls | 0 (pure graph/flow traversal). |
| Key inputs | Baseline `flow_result` and graph path structure. |
| Key outputs | `FlowDecomposition` (`components`, `total_flow`, `is_unique`). |
| Example call | `d = decompose_flow(edgelist, source_nodes, sink_nodes, flow_result)` |

### 6.6 ParametricThresholdModule

| Field | Content |
|---|---|
| Purpose | Compute exact degradation/upgrade thresholds for target throughput. |
| Algorithm | Piecewise-linear parametric logic + bounded exact reruns (`max_depth`). |
| Theorem | Piecewise-linear parametric max-flow structure + Max-Flow Min-Cut. |
| Solver calls | Candidate-dependent, bounded by roughly `O(|candidate_edges| * max_depth)`. |
| Key inputs | Baseline `flow_result`; optional `target_flow`, `candidate_edges`, `max_depth`. |
| Key outputs | `ParametricThresholdResult` + per-edge threshold structs. |
| Example call | `p = analyze_parametric_thresholds(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes, flow_result; target_flow=nothing)` |

### 6.7 NodeCapacitatedFlowModule

| Field | Content |
|---|---|
| Purpose | Exact node-capacitated max-flow and node-SPOF analysis via splitting. |
| Algorithm | Node-splitting transform + delegated edge-cap max-flow solve. |
| Theorem | Node-Splitting Bijection; Menger node version. |
| Solver calls | Baseline node-cap solve `+1`; SPOF reruns up to constrained-node set size. |
| Key inputs | `node_capacities::Dict{Int64,Float64}`. |
| Key outputs | `NodeCapacitatedAnalysisResult` (`flow_result`, `spof_nodes`). |
| Example call | `nc = analyze_node_capacitated_flow(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes, node_capacities)` |

### 6.8 MinCutUtilitiesModule

| Field | Content |
|---|---|
| Purpose | Min-cut extraction, edge characterizations, bounded min-cut enumeration. |
| Algorithm | Residual BFS lattice characterization + free-zone enumeration. |
| Theorem | Min-cut lattice (`S*`, `S**`) + Max-Flow Min-Cut. |
| Solver calls | 0 additional solves (consumes baseline `flow_result`). |
| Key inputs | Baseline `flow_result`; optional `cut_limit`. |
| Key outputs | `MinCutAnalysis`, `MinCutEnumeration`, cut-edge families. |
| Example call | `mc = analyze_min_cuts(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes, flow_result; cut_limit=1000)` |

### 6.9 GlobalConnectivityModule

| Field | Content |
|---|---|
| Purpose | Exact directed edge/node connectivity and global min-cut. |
| Algorithm | Super-sink O(V) formulations + two-pass directed global-min-cut solves. |
| Theorem | Menger (edge/node), Max-Flow Min-Cut, Integrality. |
| Solver calls | `edge_connectivity`: `O(V)`; `node_connectivity`: `O(V)`; `global_min_cut`: `2(V-1)`. |
| Key inputs | Standard contract (+ capacities for weighted cuts/connectivity). |
| Key outputs | `GlobalConnectivityResult` (`lambda`, `kappa`, global cut witness). |
| Example call | `g = analyze_global_connectivity(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes)` |

### 6.10 CapacityAnalysisKit (`analyze_all`)

| Field | Content |
|---|---|
| Purpose | Single-call full analysis pipeline returning unified typed output. |
| Algorithm | Pure orchestration; no additional optimization logic. |
| Theorem | Preserves exactness guarantees of delegated modules. |
| Solver calls | 1 baseline solve plus module-required additional reruns. |
| Key inputs | Standard contract + optional knobs (`node_capacities`, `target_flow`, `k_failure`, limits, `algorithm`, `tol`). |
| Key outputs | `CapacityAnalysisKitResult`. |
| Example call | `kit = analyze_all(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes; algorithm=:dinic)` |

---

## Section 7: Output Structure Reference

### FlowSolveResult
Exact result of a max-flow solve with augmented super-terminals.

| Field | Type | Description |
|---|---|---|
| max_flow | Float64 | Optimal total flow value. |
| flow | Dict{Tuple{Int64,Int64},Float64} | Flow on original edges. |
| augmented_flow | Dict{Tuple{Int64,Int64},Float64} | Flow on augmented graph. |
| augmented_outgoing | Dict{Int64,Set{Int64}} | Augmented outgoing adjacency. |
| augmented_incoming | Dict{Int64,Set{Int64}} | Augmented incoming adjacency. |
| augmented_capacities | Dict{Tuple{Int64,Int64},Float64} | Augmented capacities. |
| residual_capacity | Dict{Tuple{Int64,Int64},Float64} | Post-solve residual capacities. |
| node_flow | Dict{Int64,Float64} | Node throughput summary. |
| sources | Vector{Int64} | Source IDs in solve. |
| sinks | Vector{Int64} | Sink IDs in solve. |
| super_source | Int64 | Internal super-source ID. |
| super_sink | Int64 | Internal super-sink ID. |
| mincut_S | Set{Int64} | One source-side min-cut set. |
| mincut_T | Set{Int64} | One sink-side min-cut set. |
| mincut_capacity | Float64 | Capacity of reported minimum cut. |
| saturated_edges | Vector{Tuple{Int64,Int64}} | Edges saturated within tolerance. |
| sink_flow | Dict{Int64,Float64} | Per-sink inflow. |
| is_unbounded | Bool | True if max-flow is unbounded. |

### SensitivityResult

| Field | Type | Description |
|---|---|---|
| critical_edges | Vector{NamedTuple} | Criticality ranking tuples. |
| marginal_capacity | Dict{Tuple{Int64,Int64},Float64} | Marginal value per unit edge capacity. |
| birnbaum | Dict{Tuple{Int64,Int64},Float64} | Birnbaum importance scores. |

### FailureImpactResult

| Field | Type | Description |
|---|---|---|
| min_cut_edges | Vector{Tuple{Int64,Int64}} | Baseline min-cut candidate edge set. |
| single_edge_failures | Vector{NamedTuple} | Single-edge failure impacts. |
| k_edge_failures | Vector{NamedTuple} | k-edge failure combination impacts. |
| degradation_results | Vector{NamedTuple} | Scenario degradation impacts. |

### StructuralResult

| Field | Type | Description |
|---|---|---|
| spof_edges | Vector{Tuple{Int64,Int64}} | Structural edge SPOFs. |
| spof_nodes | Vector{Int64} | Structural node SPOFs. |
| paths | Vector{Vector{Int64}} | Enumerated source-sink paths. |
| path_flow_contributions | Vector{NamedTuple} | Contribution tuples per path. |
| bottleneck_ranking | Vector{NamedTuple} | Bottleneck ranking tuples. |
| node_positions | Dict{Int64,Symbol} | Node position labels. |
| edge_redundancy | Dict{Tuple{Int64,Int64},Int64} | Edge redundancy map. |

### FlowPathComponent

| Field | Type | Description |
|---|---|---|
| path | Vector{Int64} | Ordered path nodes. |
| flow_value | Float64 | Path component flow. |
| bottleneck_edge | Tuple{Int64,Int64} | Bottleneck edge for this component. |

### FlowDecomposition

| Field | Type | Description |
|---|---|---|
| components | Vector{FlowPathComponent} | Decomposition components. |
| total_flow | Float64 | Sum of component flows. |
| is_unique | Bool | Module-reported uniqueness flag. |

### DegradationThreshold

| Field | Type | Description |
|---|---|---|
| target_edge | Tuple{Int64,Int64} | Edge under degradation analysis. |
| original_capacity | Float64 | Baseline capacity. |
| threshold_capacity | Float64 | Exact threshold capacity for target flow. |
| degradation_margin | Float64 | Allowed degradation from baseline. |
| target_flow | Float64 | Throughput target. |
| baseline_flow | Float64 | Baseline max-flow value. |
| target_achievable | Bool | Whether target is achievable in tested range. |
| target_reachable_at_zero | Bool | Whether target holds at zero edge capacity. |
| solver_calls | Int64 | Solver calls used for this edge threshold. |

### UpgradeThreshold

| Field | Type | Description |
|---|---|---|
| target_edge | Tuple{Int64,Int64} | Edge under upgrade analysis. |
| original_capacity | Float64 | Baseline capacity. |
| required_capacity | Float64 | Required capacity to hit target flow. |
| required_increase | Float64 | Additional required capacity. |
| target_flow | Float64 | Throughput target. |
| baseline_flow | Float64 | Baseline max-flow value. |
| already_sufficient | Bool | Whether baseline already satisfies target. |
| upgrade_ineffective | Bool | Whether upgrade cannot achieve target under constraints. |
| solver_calls | Int64 | Solver calls used for this edge threshold. |

### ParametricThresholdResult

| Field | Type | Description |
|---|---|---|
| degradation_thresholds | Vector{DegradationThreshold} | Per-edge degradation thresholds. |
| target_flow | Float64 | Effective target flow used. |
| baseline_flow | Float64 | Baseline max-flow value. |

### NodeSplitGraph

| Field | Type | Description |
|---|---|---|
| split_edgelist | Vector{Tuple{Int64,Int64}} | Split-graph edges. |
| split_outgoing | Dict{Int64,Set{Int64}} | Split-graph outgoing adjacency. |
| split_incoming | Dict{Int64,Set{Int64}} | Split-graph incoming adjacency. |
| split_capacities | Dict{Tuple{Int64,Int64},Float64} | Split-graph capacities. |
| split_sources | Vector{Int64} | Split-graph sources. |
| split_sinks | Vector{Int64} | Split-graph sinks. |
| node_to_in | Dict{Int64,Int64} | Original node -> split in-node ID. |
| node_to_out | Dict{Int64,Int64} | Original node -> split out-node ID. |
| split_to_original | Dict{Int64,Int64} | Split node -> original node. |
| original_edgelist | Vector{Tuple{Int64,Int64}} | Stored original edge list. |
| split_edge_to_original | Dict{Tuple{Int64,Int64},Tuple{Int64,Int64}} | Split edge -> original edge map. |
| node_internal_edges | Dict{Int64,Tuple{Int64,Int64}} | Original node -> internal split edge. |

### NodeCapacitatedFlowResult

| Field | Type | Description |
|---|---|---|
| max_flow | Float64 | Node-capacitated optimal flow. |
| flow | Dict{Tuple{Int64,Int64},Float64} | Original-edge flow map. |
| node_flow | Dict{Int64,Float64} | Throughput by original node. |
| sources | Vector{Int64} | Source nodes. |
| sinks | Vector{Int64} | Sink nodes. |
| sink_flow | Dict{Int64,Float64} | Per-sink flow. |
| saturated_edges | Vector{Tuple{Int64,Int64}} | Saturated original edges. |
| saturated_nodes | Vector{Int64} | Saturated constrained nodes. |
| mincut_S | Set{Int64} | Source-side cut nodes. |
| mincut_T | Set{Int64} | Sink-side cut nodes. |
| mincut_capacity | Float64 | Node-cap min-cut capacity. |
| is_unbounded | Bool | True if unbounded. |
| node_split_graph | NodeSplitGraph | Stored split graph/mapping. |

### NodeCapacitatedAnalysisResult

| Field | Type | Description |
|---|---|---|
| flow_result | NodeCapacitatedFlowResult | Baseline node-cap flow solve result. |
| spof_nodes | Vector{Int64} | Node SPOFs under node-cap analysis. |

### MinCut

| Field | Type | Description |
|---|---|---|
| S | Set{Int64} | Source-side cut set. |
| T | Set{Int64} | Sink-side cut set. |
| crossing_edges | Vector{Tuple{Int64,Int64}} | Edges crossing S->T. |
| capacity | Float64 | Cut capacity. |

### MinCutEnumeration

| Field | Type | Description |
|---|---|---|
| cuts | Vector{MinCut} | Enumerated min cuts. |
| total_cuts | Int64 | Total cuts if complete; truncated count otherwise. |
| is_complete | Bool | True if full family enumerated. |
| free_zone_size | Int64 | Size of free-zone set `S** \ S*`. |

### MinCutAnalysis

| Field | Type | Description |
|---|---|---|
| representative_cut | MinCut | Representative minimum cut. |
| edges_in_some_cut | Vector{Tuple{Int64,Int64}} | Edges in at least one minimum cut. |
| edges_in_every_cut | Vector{Tuple{Int64,Int64}} | Edges in every minimum cut. |
| enumeration | MinCutEnumeration | Enumeration summary/results. |
| max_flow | Float64 | Baseline max flow used by analysis. |
| min_cut_capacity | Float64 | Baseline min-cut capacity. |

### EdgeConnectivityResult

| Field | Type | Description |
|---|---|---|
| lambda | Int64 | Directed edge connectivity value. |
| achieving_source | Int64 | Source achieving minimum. |
| achieving_sink | Int64 | Sink achieving minimum. |
| min_cut_edges | Vector{Tuple{Int64,Int64}} | Witness edge cut. |
| solver_calls | Int64 | Solver calls consumed. |

### NodeConnectivityResult

| Field | Type | Description |
|---|---|---|
| kappa | Int64 | Directed node connectivity value. |
| achieving_source | Int64 | Source achieving minimum. |
| achieving_sink | Int64 | Sink achieving minimum. |
| min_cut_nodes | Vector{Int64} | Witness node cut. |
| solver_calls | Int64 | Solver calls consumed. |

### GlobalMinCutResult

| Field | Type | Description |
|---|---|---|
| min_cut_capacity | Float64 | Global directed minimum-cut capacity. |
| achieving_source | Int64 | Source of witness pair. |
| achieving_sink | Int64 | Sink of witness pair. |
| min_cut_edges | Vector{Tuple{Int64,Int64}} | Witness crossing edges. |
| cut_S | Set{Int64} | Source-side witness cut set. |
| cut_T | Set{Int64} | Sink-side witness cut set. |
| solver_calls | Int64 | Solver calls consumed. |

### GlobalConnectivityResult

| Field | Type | Description |
|---|---|---|
| edge_connectivity | EdgeConnectivityResult | Edge connectivity result. |
| node_connectivity | NodeConnectivityResult | Node connectivity result. |
| global_min_cut | GlobalMinCutResult | Global directed min-cut result. |

### CapacityAnalysisKitResult

| Field | Type | Description |
|---|---|---|
| flow | FlowSolveResult | Baseline flow solve (Step 1). |
| sensitivity | SensitivityResult | Sensitivity analysis (Step 2). |
| failure_impact | FailureImpactResult | Failure impact analysis (Step 3). |
| structure | StructuralResult | Structural analysis (Step 4). |
| flow_decomposition | FlowDecomposition | Flow decomposition (Step 5). |
| parametric_thresholds | ParametricThresholdResult | Parametric thresholds (Step 6). |
| node_capacitated | Union{NodeCapacitatedAnalysisResult, Nothing} | Node-cap analysis (Step 7) or `nothing`. |
| min_cut_analysis | MinCutAnalysis | Min-cut analysis (Step 8). |
| global_connectivity | GlobalConnectivityResult | Global connectivity (Step 9). |
| algorithm | Symbol | Chosen flow algorithm. |
| tol | Float64 | Floating-point comparison tolerance. |
| baseline_max_flow | Float64 | Alias of `flow.max_flow`. |

---

## Section 8: Algorithm Selection Guide

| Algorithm | Complexity | Best for | Worst for |
|---|---|---|---|
| `:dinic` | `O(V²E)` | General, DAG-heavy and mixed networks | Very dense graphs |
| `:edmonds_karp` | `O(VE²)` | Small/sparse graphs and augmenting-path tracing | Large graphs |
| `:push_relabel` | `O(V²√E)` | Dense infrastructure graphs | Some unit-capacity instances where Dinic can be faster |

**Recommendation**
- Default `:dinic` for most workloads.
- Use `:push_relabel` for very dense graphs.
- Use `:edmonds_karp` when path-trace interpretability matters.

---

## Section 9: Exactness Statement

### What “exact” means here
Every returned value is one of:
1. A direct exact max-flow solve output, validated.
2. A closed-form derivation from proved structural results.
3. An exact graph traversal result (BFS/DFS) by construction.

No approximation, sampling, or heuristic convergence mechanism is introduced in the computation chain.

### What `tol` means
- `tol` is used only for floating-point comparison checks (saturation/equality conditions).
- `tol` is **not** an approximation tolerance on objective values.
- Default `tol=1e-10` is the standard setting.
- Practical tuning: `tol=1e-12` for high-precision integer-capacity cases; `tol=1e-8` for large floating-point-capacity networks when additional numerical robustness is needed.

### Validation suite checks
For each `FlowSolveResult`, validation includes:
- `validate_capacity_constraints`: verifies `0 <= f(e) <= c(e)`.
- `validate_flow_conservation`: verifies node-wise conservation at non-terminals.
- `validate_maxflow_mincut`: verifies max-flow value equals min-cut capacity within comparison tolerance.

### Explicitly out of scope
- Probabilistic component reliability.
- Interval/uncertain capacities (use deterministic scenario envelopes).
- Time-varying/dynamic network flow.
- Min-cost max-flow.
- Approximate/heuristic solvers.

---

## Appendix A: `analyze_all` Step-to-Output and Solver Cost Map

Let $V$ be node count and $E$ edge count.

| Step | `analyze_all` call | Result field written | Incremental solver calls (beyond prior steps) |
|---|---|---|---|
| 1 | `_solve_with_algorithm(...)` | `result.flow`, `result.baseline_max_flow` | `+1` |
| 2 | `analyze_sensitivity(...)` | `result.sensitivity` | up to `|saturated_edges|` zero-capacity + up to `|min_cut_candidates|` Birnbaum + up to `|saturated_edges|` marginal reruns (overall `O(E)`) |
| 3 | `analyze_failure_impact(...)` | `result.failure_impact` | single-edge up to `|min_cut_candidates|`; k-edge combinational bounded by `combination_limit`; scenario reruns one per scenario |
| 4 | `analyze_structure(...)` | `result.structure` | traversal-heavy + candidate-dependent redundancy reruns |
| 5 | `decompose_flow(...)` | `result.flow_decomposition` | `+0` |
| 6 | `analyze_parametric_thresholds(...)` | `result.parametric_thresholds` | roughly bounded by `O(|candidate_edges|*max_depth)` |
| 7 | `analyze_node_capacitated_flow(...)` (conditional) | `result.node_capacitated` | if no node capacities: `+0`; else baseline node-cap solve `+1` plus SPOF reruns |
| 8 | `analyze_min_cuts(...)` | `result.min_cut_analysis` | `+0` |
| 9 | `analyze_global_connectivity(...)` | `result.global_connectivity` | `+V` (edge connectivity) + `+V` (node connectivity) + `+2(V-1)` (global min-cut) = `+(4V-2)` |

### Notes
- Baseline solve is executed once in Step 1.
- Downstream modules reuse the same baseline result when appropriate.
- Additional solves are required by module logic, not by aggregator duplication.

---

## Appendix B: Question-to-Function Quick Reference

| Engineering question | Primary function(s) | Result field(s) to read |
|---|---|---|
| Maximum throughput? | `analyze_all`, `solve_max_flow_*` | `result.baseline_max_flow`, `result.flow.max_flow` |
| Flow to each sink? | `analyze_all` | `result.flow.sink_flow` |
| Active routing paths? | `decompose_flow`, `analyze_all` | `result.flow_decomposition.components` |
| Bottleneck cut location? | `analyze_all` | `result.flow.mincut_S`, `result.flow.mincut_T`, `result.structure.bottleneck_ranking` |
| Fully saturated edges? | `analyze_all` | `result.flow.saturated_edges` |
| Single points of failure? | `analyze_structure`, `analyze_node_capacitated_flow` | `result.structure.spof_edges`, `result.structure.spof_nodes`, `result.node_capacitated.spof_nodes` |
| Degradation margin per edge? | `analyze_parametric_thresholds` | `result.parametric_thresholds.degradation_thresholds` |
| Required upgrade per edge? | `find_upgrade_threshold` | `UpgradeThreshold.required_increase`, `UpgradeThreshold.required_capacity` |
| Multi-failure impact? | `analyze_failure_impact` | `result.failure_impact.k_edge_failures`, `result.failure_impact.degradation_results` |
| All minimum-cut failure modes? | `enumerate_min_cuts`, `analyze_min_cuts` | `MinCutEnumeration.cuts`, `MinCutEnumeration.total_cuts` |
| Edge/node disconnect robustness? | `edge_connectivity`, `node_connectivity`, `analyze_global_connectivity` | `lambda`, `kappa`, `result.global_connectivity.*` |
| Most important components? | `analyze_sensitivity` | `result.sensitivity.birnbaum`, `result.sensitivity.critical_edges`, `result.sensitivity.marginal_capacity` |
