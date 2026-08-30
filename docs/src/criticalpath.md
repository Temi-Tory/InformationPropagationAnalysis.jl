# [CriticalPath](@id criticalpath)

`CriticalPath` is a mode-based longest / shortest / max-scaling path analysis.
Every mode carries **forward and backward** semantics, so one call gives the
project value, the per-node margins (slack), and — for the additive modes — the
early-start / late-finish / late-start schedule.

## Modes

| Mode | Combine | Compound | Margin | Additive schedule |
|---|---|---|---|---|
| `LONGEST_PATH` | `maximum` | `+` | slack | ✓ |
| `SHORTEST_PATH` | `minimum` | `+` | margin | ✓ |
| `MAX_SCALING` | `maximum` | `*` | ratio-slack | — |

## `analyze` / `critical_path`

```
result = CriticalPath.analyze(iteration_sets, outgoing_index, incoming_index,
                              source_nodes, node_values, edge_values;
                              mode = CriticalPath.LONGEST_PATH)
```

`critical_path` (top-level) is an alias for `CriticalPath.analyze`. The result is
a `PathResult{T}` with fields `forward`, `reverse_completion`, `project_value`,
`through`, `margin`, `critical` (the zero-margin nodes), and the schedule dicts
`early_start` / `late_finish` / `late_start`.

```@example cpm
using InformationPropagationAnalysis
const IPA = InformationPropagationAnalysis

# 1 → 2 → 4 and 1 → 3 → 4 ; durations on the nodes
edgelist = Tuple{Int64,Int64}[(1,2), (1,3), (2,4), (3,4)]
outgoing = Dict{Int64,Set{Int64}}(1 => Set([2,3]), 2 => Set([4]), 3 => Set([4]))
incoming = Dict{Int64,Set{Int64}}(2 => Set([1]), 3 => Set([1]), 4 => Set([2,3]))
sources  = Set{Int64}([1])
itersets, _, _ = IPA.Input.find_iteration_sets(edgelist, outgoing, incoming)

durations = Dict{Int64,Float64}(1 => 0.0, 2 => 5.0, 3 => 2.0, 4 => 3.0)
delays    = Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in edgelist)

r = critical_path(itersets, outgoing, incoming, sources, durations, delays;
                  mode = CriticalPath.LONGEST_PATH)

(makespan = r.project_value, critical = r.critical, slack_of_3 = r.margin[3])
```

The 1→2→4 path (0 + 5 + 3 = 8) is longer than 1→3→4 (0 + 2 + 3 = 5), so node 3
carries 3 units of slack and the critical nodes are `[1, 2, 4]`.

## Interval durations

`interval_analyze` runs the two corner problems (all-low, all-high). For
`LONGEST_PATH`, `interval_analyze_split` does better: it subdivides the
uncertain nodes to get tighter, still-rigorous bounds on each node's margin,
declining (`SplitDeclined`) rather than exceeding `max_runs`.

Interval **edge** delays are not supported by the split — subdivide such an edge
into a duration node first.

## Docstrings

```@docs
InformationPropagationAnalysis.CriticalPath.analyze
InformationPropagationAnalysis.CriticalPath.accumulation_analysis
InformationPropagationAnalysis.CriticalPath.interval_analyze
InformationPropagationAnalysis.CriticalPath.interval_analyze_exact
InformationPropagationAnalysis.CriticalPath.interval_analyze_split
InformationPropagationAnalysis.CriticalPath.PathResult
InformationPropagationAnalysis.CriticalPath.AccumulationResult
InformationPropagationAnalysis.CriticalPath.SplitDeclined
```
