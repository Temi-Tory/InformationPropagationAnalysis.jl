# [Probability](@id probability)

`Probability` computes, for **every** node `v`,

```
belief(v) = P(v is operational AND reachable from some source)
```

with every node prior and edge probability an independent random variable. Paths
overlap, so this is not a simple product — computing it exactly is #P-hard in
general. The package gets the exact value by **diamond conditioning**: wherever
[`Diamonds`](@ref diamonds) finds correlated parents at a join, propagation sums
over the states of that diamond's conditioning set, within which the parents *are*
independent, and recurses for nested diamonds.

## `update_beliefs_iterative`

```
belief = update_beliefs_iterative(edgelist, iteration_sets, outgoing_index,
                                  incoming_index, source_nodes, node_priors,
                                  link_probability, descendants, ancestors,
                                  diamond_structures, join_nodes, fork_nodes,
                                  computation_lookup)
```

`diamond_structures` and `computation_lookup` are the two objects returned by
[`new_identify`](@ref diamonds). The result is `Dict{Int64,T}` with the same
element type `T ∈ {Float64, Interval, pbox}` as the priors — one propagation,
any value form.

`validate_network_data` runs the same topology / probability contract checks that
`update_beliefs_iterative` performs internally, so you can fail fast before a long
run.

## Point, interval and p-box in one call

```@example prob
using InformationPropagationAnalysis
const IPA = InformationPropagationAnalysis

edgelist = Tuple{Int64,Int64}[(1,2), (1,3), (2,4), (3,4)]
outgoing = Dict{Int64,Set{Int64}}(1 => Set([2,3]), 2 => Set([4]), 3 => Set([4]))
incoming = Dict{Int64,Set{Int64}}(2 => Set([1]), 3 => Set([1]), 4 => Set([2,3]))
sources  = Set{Int64}([1])
itersets, anc, desc = IPA.Input.find_iteration_sets(edgelist, outgoing, incoming)
forks, joins        = IPA.Input.identify_fork_and_join_nodes(outgoing, incoming)

function belief4(priors, probs)
    s, l = new_identify(edgelist, priors, probs, sources, forks, joins, anc, desc, itersets)
    update_beliefs_iterative(edgelist, itersets, outgoing, incoming, sources,
                             priors, probs, desc, anc, s, joins, forks, l)[4]
end

# point
p  = belief4(Dict{Int64,Float64}(n => 0.9 for n in 1:4),
             Dict{Tuple{Int64,Int64},Float64}(e => 0.9 for e in edgelist))
# interval — priors known only to ±0.02
iv = belief4(Dict{Int64,Interval}(n => Interval(0.88, 0.92) for n in 1:4),
             Dict{Tuple{Int64,Int64},Interval}(e => Interval(0.88, 0.92) for e in edgelist))

(point = round(p; digits = 5), interval = (round(iv.lower; digits = 5), round(iv.upper; digits = 5)))
```

The interval result brackets the point result: same algorithm, the uncertainty
in the inputs carried through to the output.

## Docstrings

```@docs
InformationPropagationAnalysis.Probability.update_beliefs_iterative
InformationPropagationAnalysis.Probability.calculate_regular_belief
InformationPropagationAnalysis.Probability.inclusion_exclusion
```
