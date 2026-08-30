# Test-only oracles and helpers. Pure Julia.
#
# `brute_reliability` is exact reachability reliability by state enumeration over the
# 2^(|nodes|+|edges|) component-up/down space — lifted from validation/oracles.jl. Only
# tractable for tiny graphs (|nodes|+|edges| ≲ 22), which is what the synthetic diamond
# fixtures are for. Larger fixtures assert against recorded values instead.

module Oracle

using JSON

export brute_reliability, build_indices, parse_cpm_inputs, parse_cpm_inputs_interval,
       single_source, read_json_dict

"Parse a JSON object file into a `Dict`."
read_json_dict(path::AbstractString) = JSON.parsefile(path)

"Adjacency (outgoing, incoming) from an edge list."
function build_indices(edgelist::Vector{Tuple{Int,Int}})
    out = Dict{Int,Set{Int}}()
    inc = Dict{Int,Set{Int}}()
    for (u, v) in edgelist
        push!(get!(out, u, Set{Int}()), v)
        push!(get!(inc, v, Set{Int}()), u)
    end
    return out, inc
end

"The unique node with no incoming edges (errors if not unique)."
function single_source(edgelist::Vector{Tuple{Int,Int}})
    _, inc = build_indices(edgelist)
    nodes = union(Set(first.(edgelist)), Set(last.(edgelist)))
    srcs = [n for n in nodes if !haskey(inc, n) || isempty(inc[n])]
    length(srcs) == 1 || error("expected exactly one source, got $srcs")
    return srcs[1]
end

"""
    brute_reliability(edgelist, node_priors, edge_probs, sources) -> Dict(node => belief)

Exact P(node operational AND reachable from some source), every component independent.
"""
function brute_reliability(edgelist::Vector{Tuple{Int,Int}},
                           node_priors::Dict{Int,Float64},
                           edge_probs::Dict{Tuple{Int,Int},Float64},
                           sources)
    nodes = sort(collect(union(Set(first.(edgelist)), Set(last.(edgelist)))))
    edges = collect(keys(edge_probs))
    nid = Dict(n => i for (i, n) in enumerate(nodes))
    eid = Dict(e => length(nodes) + k for (k, e) in enumerate(edges))
    V = length(nodes) + length(edges)
    V <= 24 || error("brute_reliability: $V bits is too large; use a recorded value")
    srcset = Set(sources)

    bel = Dict(n => 0.0 for n in nodes)
    for s in 0:(UInt64(1) << V - 1)
        up_node(n) = (s >> (nid[n] - 1)) & 1 == 1
        up_edge(e) = (s >> (eid[e] - 1)) & 1 == 1
        w = 1.0
        for n in nodes
            w *= up_node(n) ? node_priors[n] : (1 - node_priors[n])
        end
        for e in edges
            w *= up_edge(e) ? edge_probs[e] : (1 - edge_probs[e])
        end
        w == 0.0 && continue
        reach = Set{Int}(n for n in srcset if up_node(n))
        changed = true
        while changed
            changed = false
            for (u, v) in edgelist
                if u in reach && !(v in reach) && up_node(v) && up_edge((u, v))
                    push!(reach, v); changed = true
                end
            end
        end
        for n in reach
            bel[n] += w
        end
    end
    return bel
end

_edgekey(k) = (m = match(r"\(\s*(\d+)\s*,\s*(\d+)\s*\)", k); (parse(Int, m[1]), parse(Int, m[2])))

"""
    parse_cpm_inputs(path) -> (node_durations::Dict{Int,Float64}, edge_delays::Dict{Tuple{Int,Int},Float64})

Reads the `time_analysis.node_durations` / `edge_delays` schema of the `*-cpm-inputs.json`
fixtures (scalar Float64 variant).
"""
function parse_cpm_inputs(path::String)
    ta = JSON.parsefile(path)["time_analysis"]
    nd = Dict{Int,Float64}(parse(Int, k) => Float64(v) for (k, v) in ta["node_durations"])
    ed = Dict{Tuple{Int,Int},Float64}()
    if haskey(ta, "edge_delays")
        for (k, v) in ta["edge_delays"]
            ed[_edgekey(k)] = Float64(v)
        end
    end
    return nd, ed
end

"""
    parse_cpm_inputs_interval(path) -> (node_durations::Dict{Int,Tuple{Float64,Float64}}, edge_delays)

Interval variant: node durations are `{"type":"interval","lower":lo,"upper":hi}`.
"""
function parse_cpm_inputs_interval(path::String)
    ta = JSON.parsefile(path)["time_analysis"]
    nd = Dict{Int,Tuple{Float64,Float64}}()
    for (k, v) in ta["node_durations"]
        nd[parse(Int, k)] = v isa Number ? (Float64(v), Float64(v)) :
                            (Float64(v["lower"]), Float64(v["upper"]))
    end
    ed = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    if haskey(ta, "edge_delays")
        for (k, v) in ta["edge_delays"]
            ed[_edgekey(k)] = v isa Number ? (Float64(v), Float64(v)) :
                              (Float64(v["lower"]), Float64(v["upper"]))
        end
    end
    return nd, ed
end

end # module Oracle
