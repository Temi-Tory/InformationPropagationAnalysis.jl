# InformationPropagationAnalysis.jl

Exact analysis of how information, probability, time and flow propagate through a
**directed acyclic network**, with every quantity carried as a point value
(`Float64`), an interval, or a probability box (`pbox`) through the *same*
algorithm.

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://Temi-Tory.github.io/InformationPropagationAnalysis.jl/)

## Toolkits

The package is five toolkits, each a submodule you can `using` on its own:

| Toolkit | Question it answers |
|---|---|
| `Input` | Read a network, node priors, edge probabilities and capacities from files. |
| `Diamonds` | Where do paths reconverge, and what must be conditioned on to decorrelate them? |
| `Probability` | The exact probability each node is operational **and** reachable from a source. |
| `CriticalPath` | Longest / shortest / max-scaling path, forward and backward, with interval bounds. |
| `Flow` | Max-flow, min-cut, bottlenecks, failure impact, sensitivity. |

## Installation

```julia
using Pkg
Pkg.add("InformationPropagationAnalysis")
```

Requires Julia 1.12 or newer. `ProbabilityBoundsAnalysis` does not precompile on
Julia 1.12 (upstream); the package still loads and runs, just interpreted — see
the [Reproducibility](https://Temi-Tory.github.io/InformationPropagationAnalysis.jl/reproducibility/)
docs.

## Quick start — exact reachability reliability

```julia
using InformationPropagationAnalysis
const IPA = InformationPropagationAnalysis

# a four-node diamond: 1 forks to {2,3}, which both feed the join at 4
edgelist = Tuple{Int64,Int64}[(1,2), (1,3), (2,4), (3,4)]
outgoing = Dict{Int64,Set{Int64}}(1 => Set([2,3]), 2 => Set([4]), 3 => Set([4]))
incoming = Dict{Int64,Set{Int64}}(2 => Set([1]), 3 => Set([1]), 4 => Set([2,3]))
sources  = Set{Int64}([1])

node_priors = Dict{Int64,Float64}(n => 0.9  for n in 1:4)
edge_probs  = Dict{Tuple{Int64,Int64},Float64}(e => 0.85 for e in edgelist)

itersets, anc, desc = IPA.Input.find_iteration_sets(edgelist, outgoing, incoming)
forks, joins        = IPA.Input.identify_fork_and_join_nodes(outgoing, incoming)

structure, lookup = new_identify(edgelist, node_priors, edge_probs, sources,
                                 forks, joins, anc, desc, itersets)

belief = update_beliefs_iterative(edgelist, itersets, outgoing, incoming, sources,
                                  node_priors, edge_probs, desc, anc,
                                  structure, joins, forks, lookup)

belief[4]   # exact P(node 4 up and reachable from 1)
```

`using InformationPropagationAnalysis` brings in the five toolkit names, the
value types `Interval` and `pbox`, and one blessed entry point per toolkit —
`new_identify`, `update_beliefs_iterative`, `critical_path`, `analyze_all`.
Everything else is reached through its toolkit (`Flow.analyze_structure`,
`CriticalPath.LONGEST_PATH`, `Input.read_graph_to_dict`, …).

## File formats

Node priors and edge probabilities are JSON keyed by node id / `"(u,v)"` edge
string, with a `data_type` field:

```json
{ "data_type": "Float64",  "nodes": {"1": 0.8, "2": 0.7} }
{ "data_type": "Interval", "nodes": {"1": {"type": "interval", "lower": 0.7, "upper": 0.9}} }
{ "data_type": "pbox",     "nodes": {"1": {"type": "pbox", "construction": "parametric",
                                           "shape": "normal", "params": [0.8, 0.1]}} }
```

Networks are `.EDGES` files (CSV, `source,destination` header) or 0/1 adjacency
matrices; capacities are `{"edges": [{"source": …, "destination": …, "capacity": …}]}`.

## Documentation

Full documentation, one page per toolkit, at
<https://Temi-Tory.github.io/InformationPropagationAnalysis.jl/>.

## Applications

[**information-propagation-no-code**](https://github.com/Temi-Tory/information-propagation-no-code)
— a no-code workbench (HTTP server + Angular UI) built on this package: upload a
network, run every analysis from the browser, no Julia required.

## Citation

```bibtex
@software{InformationPropagationAnalysis,
  title  = {InformationPropagationAnalysis.jl: Exact Reachability Analysis with Diamond Structure Optimization},
  author = {T. Ohiani and E. Patelli},
  year   = {2024},
  url    = {https://github.com/Temi-Tory/InformationPropagationAnalysis.jl}
}
```

The theoretical foundation:

```bibtex
@inproceedings{ohiani2023information,
  title        = {The Information Propagation Method for Efficient Network Reliability Analysis},
  author       = {T. Ohiani and E. Patelli},
  booktitle    = {2023 7th International Conference on System Reliability and Safety (ICSRS)},
  pages        = {580--584},
  year         = {2023},
  organization = {IEEE},
  address      = {Bologna, Italy},
  doi          = {10.1109/ICSRS59833.2023.10381157}
}
```

## Acknowledgments

T. Ohiani is supported by the UK National Decommissioning Authority Bursary
Studentship "Developing a resilience framework for decommissioning plan of a
nuclear facility".

## License

MIT. See [LICENSE](LICENSE).
