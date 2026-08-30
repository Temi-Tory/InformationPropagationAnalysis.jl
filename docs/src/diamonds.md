# [Diamonds](@id diamonds)

At a join node, two parents can share an upstream fork — their reachability
events are then **correlated**, and combining them with independent
inclusion–exclusion is wrong. That correlated substructure is a **diamond**.
`Diamonds` finds every diamond in a DAG and computes, for each, a
**conditioning set**: a set of nodes such that, conditioned on their states, the
join's parents become mutually independent.

## `new_identify`

`new_identify` is the sole diamond producer — a factorized,
correct-by-construction identification that emits the root diamonds and the full
nested (`unique`) diamond set together, ready for
[`update_beliefs_iterative`](@ref probability).

```
structure, lookup = new_identify(edgelist, node_priors, link_probs, source_nodes,
                                 fork_nodes, join_nodes, ancestors, descendants,
                                 iteration_sets)
```

`link_probs` is accepted for signature symmetry with propagation but is not used
by identification itself. The call is generic over
`T ∈ {Float64, Interval, pbox}` (the priors' element type selects the path).

```@example diamonds
using InformationPropagationAnalysis
const IPA = InformationPropagationAnalysis

# node 4 joins {2,3}, both descended from the fork at node 1
edgelist = Tuple{Int64,Int64}[(1,2), (1,3), (2,4), (3,4)]
outgoing = Dict{Int64,Set{Int64}}(1 => Set([2,3]), 2 => Set([4]), 3 => Set([4]))
incoming = Dict{Int64,Set{Int64}}(2 => Set([1]), 3 => Set([1]), 4 => Set([2,3]))
sources  = Set{Int64}([1])
priors   = Dict{Int64,Float64}(n => 0.9 for n in 1:4)
probs    = Dict{Tuple{Int64,Int64},Float64}(e => 0.9 for e in edgelist)

itersets, anc, desc = IPA.Input.find_iteration_sets(edgelist, outgoing, incoming)
forks, joins        = IPA.Input.identify_fork_and_join_nodes(outgoing, incoming)

structure, lookup = new_identify(edgelist, priors, probs, sources,
                                 forks, joins, anc, desc, itersets)

collect(keys(structure))   # join nodes that carry a diamond
```

## Result types

- `Diamond` — a single diamond: its relevant nodes, conditioning nodes and edge list.
- `DiamondsAtNode` — the diamonds attached to one join, plus its non-diamond parents.
- `DiamondComputationData` — the per-diamond payload `update_beliefs_iterative` consumes.

`Diamonds.create_diamond_hash_key` (a `public`, non-exported helper) is the
context key used to cache diamond computations across the nested recursion.

## Docstrings

```@docs
InformationPropagationAnalysis.Diamonds.new_identify
InformationPropagationAnalysis.Diamonds.Diamond
InformationPropagationAnalysis.Diamonds.DiamondsAtNode
InformationPropagationAnalysis.Diamonds.DiamondComputationData
InformationPropagationAnalysis.Diamonds.create_diamond_hash_key
```
