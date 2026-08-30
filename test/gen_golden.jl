# Regenerates the golden-master belief vectors used as regression guards by runtests.jl.
# Run from the package root:  julia --project=. test/gen_golden.jl
#
# The values are the output of the current build. They are only meaningful as a golden
# master because the underlying algorithm is independently validated elsewhere (see
# test/fixtures/README.md for the oracle provenance of each fixture).

using InformationPropagationAnalysis
const IPA = InformationPropagationAnalysis

const FIX = joinpath(@__DIR__, "fixtures")

function belief_float(dir, stem)
    el, out, inc, src = IPA.Input.read_graph_to_dict(joinpath(dir, "$stem.EDGES"))
    itersets, anc, desc = IPA.Input.find_iteration_sets(el, out, inc)
    fk, jn = IPA.Input.identify_fork_and_join_nodes(out, inc)
    np = IPA.Input.read_node_priors_from_json(joinpath(dir, "float", "$stem-nodepriors.json"))
    ep = IPA.Input.read_edge_probabilities_from_json(joinpath(dir, "float", "$stem-linkprobabilities.json"))
    r, u = IPA.new_identify(el, np, ep, src, fk, jn, anc, desc, itersets)
    IPA.update_beliefs_iterative(el, itersets, out, inc, src, np, ep, desc, anc, r, jn, fk, u)
end

let dir = joinpath(FIX, "power-network")
    bel = belief_float(dir, "power-network")
    open(joinpath(dir, "expected-float-beliefs.csv"), "w") do io
        println(io, "node,belief")
        for n in sort(collect(keys(bel)))
            println(io, "$n,$(bel[n])")
        end
    end
    println("wrote power-network/expected-float-beliefs.csv ($(length(bel)) nodes)")
end
