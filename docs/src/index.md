# InformationPropagationAnalysis.jl

Exact analysis of how information, probability, time and flow propagate through a
**directed acyclic network**, with every quantity carried as a point value
(`Float64`), an interval, or a probability box (`pbox`) through the *same*
algorithm.

The package is five toolkits, each a submodule you can `using` on its own:

| Toolkit | Question it answers |
|---|---|
| [`Input`](@ref input) | Read a network, node priors, edge probabilities and capacities from files. |
| [`Diamonds`](@ref diamonds) | Where do paths reconverge, and what must be conditioned on to decorrelate them? |
| [`Probability`](@ref probability) | What is the exact probability each node is operational **and** reachable from a source? |
| [`CriticalPath`](@ref criticalpath) | Longest / shortest / max-scaling path, forward and backward, with interval bounds. |
| [`Flow`](@ref flow) | Max-flow, min-cut, bottlenecks, failure impact, sensitivity. |

## Installation

```julia
julia> using Pkg; Pkg.add("InformationPropagationAnalysis")
```

Requires Julia 1.12 or newer. See [Reproducibility](@ref) for a note on the
`ProbabilityBoundsAnalysis` precompilation caveat on 1.12.

## The API at a glance

`using InformationPropagationAnalysis` brings in the five toolkit names, the two
universal value types (`Interval`, `pbox`), and one blessed entry point per
toolkit: `new_identify`, `update_beliefs_iterative`, `critical_path`,
`analyze_all`. Everything else is reached through its toolkit, e.g.
`Flow.analyze_structure`, `CriticalPath.LONGEST_PATH`.

```jldoctest
julia> using InformationPropagationAnalysis

julia> critical_path === CriticalPath.analyze
true

julia> sort(String.(filter(n -> n != :InformationPropagationAnalysis, names(InformationPropagationAnalysis))))
11-element Vector{String}:
 "CriticalPath"
 "Diamonds"
 "Flow"
 "Input"
 "Interval"
 "Probability"
 "analyze_all"
 "critical_path"
 "new_identify"
 "pbox"
 "update_beliefs_iterative"
```

## A first analysis — exact reachability reliability

Take a four-node diamond: node 1 forks to 2 and 3, which both feed the join at
node 4. Every node and edge is up with some independent probability; the belief
of a node is the probability it is up *and* reachable from node 1.

```@example diamond
using InformationPropagationAnalysis
const IPA = InformationPropagationAnalysis

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

round(belief[4]; digits = 6)
```

The join at node 4 has two parents that share the fork ancestor 1, so their
reachability is correlated — `Diamonds` finds that structure and
`Probability` conditions on it to get the exact value (inclusion–exclusion on
the raw parents would over- or under-count the shared path).
