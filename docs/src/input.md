# [Input](@id input)

`Input` reads networks and their numeric data from files. It also holds the
generic uncertainty-arithmetic layer (`add_values`, `multiply_values`, …) that
lets the other toolkits stay agnostic over `Float64` / `Interval` / `pbox`.

## Graph structure

`read_graph_to_dict(path)` takes an `.EDGES` file — a CSV with a
`source,destination` header — and returns

```
(edgelist, outgoing_index, incoming_index, source_nodes)
```

`outgoing_index` and `incoming_index` are `Dict{Int64,Set{Int64}}` adjacency
maps; `source_nodes` is the `Set` of nodes with no incoming edge.

`find_iteration_sets(edgelist, outgoing_index, incoming_index)` returns
`(iteration_sets, ancestors, descendants)` — the topological layering plus
transitive ancestor / descendant sets that every downstream toolkit needs.
`identify_fork_and_join_nodes(outgoing_index, incoming_index)` returns
`(fork_nodes, join_nodes)` — the nodes with out-degree ≥ 2 and in-degree ≥ 2.

```@example input
using InformationPropagationAnalysis
const IPA = InformationPropagationAnalysis

edgelist = Tuple{Int64,Int64}[(1,2), (1,3), (2,4), (3,4)]
outgoing = Dict{Int64,Set{Int64}}(1 => Set([2,3]), 2 => Set([4]), 3 => Set([4]))
incoming = Dict{Int64,Set{Int64}}(2 => Set([1]), 3 => Set([1]), 4 => Set([2,3]))

itersets, anc, desc = IPA.Input.find_iteration_sets(edgelist, outgoing, incoming)
forks, joins        = IPA.Input.identify_fork_and_join_nodes(outgoing, incoming)
(itersets, forks, joins)
```

## Node priors and edge probabilities

Each quantity has a **generic** reader that auto-detects the value type from the
file's `data_type` field:

- `read_node_priors_from_json(path)` → `Dict{Int64,T}`
- `read_edge_probabilities_from_json(path)` → `Dict{Tuple{Int64,Int64},T}`

and **type-specific** readers with a guaranteed return type — the suffix picks
the type:

| suffix | node priors | edge probabilities |
|---|---|---|
| `_float64` | `read_node_priors_from_json_float64` | `read_edge_probabilities_from_json_float64` |
| `_interval` | `read_node_priors_from_json_interval` | `read_edge_probabilities_from_json_interval` |
| `_pbox` | `read_node_priors_from_json_pbox` | `read_edge_probabilities_from_json_pbox` |

`read_complete_network(edges_file, node_priors_file, edge_probs_file)` is the
one-call convenience that wires all three together.

## Capacities

- `read_edge_capacities_from_json(path)` → `Dict{Tuple{Int64,Int64},Float64}`,
  reading `{"edges": [{"source": …, "destination": …, "capacity": …}]}`
  (`"Inf"` is accepted for an uncapacitated edge).
- `read_node_capacities_from_json(path)` → `Dict{Int64,Float64}`.
- `read_capacities_input(path)` reads the combined node + edge capacity schema.

## Value types and p-boxes

`Interval(lo, hi)` is the package's own interval type (a validated
`lo ≤ hi` pair). `pbox` and `PBA` are re-exports of
`ProbabilityBoundsAnalysis.pbox` and the `ProbabilityBoundsAnalysis` module —
use `PBA.normal(μ, σ)`, `PBA.uniform(a, b)` and friends to build parametric
p-boxes.

## Docstrings

```@autodocs
Modules = [InformationPropagationAnalysis.Input]
```
