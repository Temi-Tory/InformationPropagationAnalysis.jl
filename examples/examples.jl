# Runnable examples over the networks in examples/data/.
#
#   julia --project=. examples/examples.jl
#
# (activate an environment that has InformationPropagationAnalysis added, or run
# from the package root with --project=.)

using InformationPropagationAnalysis
const IPA = InformationPropagationAnalysis

const DATA = joinpath(@__DIR__, "data")

# ── helpers ──────────────────────────────────────────────────────────────────

function load_network(name)
    dir = joinpath(DATA, name)
    el, out, inc, src = IPA.Input.read_graph_to_dict(joinpath(dir, "$name.EDGES"))
    itersets, anc, desc = IPA.Input.find_iteration_sets(el, out, inc)
    forks, joins = IPA.Input.identify_fork_and_join_nodes(out, inc)
    (; dir, el, out, inc, src, itersets, anc, desc, forks, joins)
end

function reachability(n, priors, probs)
    s, l = new_identify(n.el, priors, probs, n.src, n.forks, n.joins, n.anc, n.desc, n.itersets)
    update_beliefs_iterative(n.el, n.itersets, n.out, n.inc, n.src,
                             priors, probs, n.desc, n.anc, s, n.joins, n.forks, l)
end

# ── 1. exact reachability reliability, three value forms ──────────────────────

let name = "power-network", n = load_network("power-network")
    println("── $name : $(length(n.src)) sources, $(length(n.el)) edges ──")

    np = IPA.Input.read_node_priors_from_json(joinpath(n.dir, "float", "$name-nodepriors.json"))
    ep = IPA.Input.read_edge_probabilities_from_json(joinpath(n.dir, "float", "$name-linkprobabilities.json"))
    bel = reachability(n, np, ep)
    println("  float  belief, worst node: ", round(minimum(values(bel)); digits = 4))

    npi = IPA.Input.read_node_priors_from_json_interval(joinpath(n.dir, "interval", "$name-nodepriors.json"))
    epi = IPA.Input.read_edge_probabilities_from_json_interval(joinpath(n.dir, "interval", "$name-linkprobabilities.json"))
    beli = reachability(n, npi, epi)
    w = argmin(v -> bel[v], collect(keys(bel)))
    println("  interval belief at node $w: [", round(beli[w].lower; digits = 4), ", ",
            round(beli[w].upper; digits = 4), "]  (float ", round(bel[w]; digits = 4), ")")
end

# ── 2. diamond structure ─────────────────────────────────────────────────────

let name = "munin-dag", n = load_network("munin-dag")
    np = IPA.Input.read_node_priors_from_json(joinpath(n.dir, "float", "$name-nodepriors.json"))
    ep = IPA.Input.read_edge_probabilities_from_json(joinpath(n.dir, "float", "$name-linkprobabilities.json"))
    structure, lookup = new_identify(n.el, np, ep, n.src, n.forks, n.joins, n.anc, n.desc, n.itersets)
    println("── $name : $(length(n.el)) edges, $(length(structure)) joins carry a diamond, ",
            "$(length(lookup)) unique diamonds ──")
end

# ── 3. flow capacity ─────────────────────────────────────────────────────────

let name = "power-network", n = load_network("power-network")
    # synthesise capacities: every edge 10.0
    caps = Dict{Tuple{Int64,Int64},Float64}(e => 10.0 for e in n.el)
    sinks = Int64[v for v in union(keys(n.out), keys(n.inc)) if !haskey(n.out, v) || isempty(n.out[v])]
    mf = IPA.Flow.solve_max_flow_dinic(n.el, n.out, n.inc, caps, collect(n.src), sinks)
    println("── $name flow : max-flow = ", mf.max_flow, ", min-cut capacity = ", mf.mincut_capacity, " ──")
end
