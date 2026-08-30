#
# DOMINATION SPLIT — exact interval floats for the longest-path mode.
#
# Theory and proofs: validation/cpm_v2/DOMINATION_SPLIT_THEORY.md. Per node v, coordinates
# split into incomparable (float nondecreasing), dominated (float nonincreasing) and the
# bypass set H_v; the monotone coordinates are pinned at their extremal ends and only H_v
# corners are enumerated. Never worse than exhaustive: if the split's run count exceeds the
# shared exhaustive sweep's, the exhaustive driver is used instead.
#
# Current scope (matching what is proven): LONGEST_PATH only; edge delays must be crisp
# (interval delays reduce to node durations by edge subdivision — not yet implemented).
#

"""
The split (or its internal exhaustive delegation) declining because the run
budget is exceeded — a legitimate, expected outcome, not a bug. Kept distinct
from every other exception `interval_analyze_split` can throw (a genuine
usage error, or an implementation bug) so a caller's fallback can catch
*this* specifically and report the real reason, instead of a generic
`ArgumentError` catch silently mislabelling an unrelated failure as "the
budget was exceeded" too.
"""
struct SplitDeclined <: Exception
    needed::Int
    max_runs::Int
    exhaustive_k::Int
end

function Base.showerror(io::IO, e::SplitDeclined)
    print(io, "domination split declined: needs $(e.needed) runs against a limit of ",
          "$(e.max_runs) (exhaustive would need 2^$(e.exhaustive_k))")
end

function _closure_sets(iteration_sets, incoming_index, outgoing_index)
    anc = Dict{Int64,Set{Int64}}()
    for layer in iteration_sets, n in layer
        s = Set{Int64}()
        for p in get(incoming_index, n, Set{Int64}())
            push!(s, p); union!(s, anc[p])
        end
        anc[n] = s
    end
    desc = Dict{Int64,Set{Int64}}()
    for i in length(iteration_sets):-1:1, n in iteration_sets[i]
        s = Set{Int64}()
        for c in get(outgoing_index, n, Set{Int64}())
            push!(s, c); union!(s, desc[c])
        end
        desc[n] = s
    end
    (anc, desc)
end

function _bypass_sets(v, nodes, outgoing_index, incoming_index)
    sinks = [n for n in nodes if n != v && isempty(get(outgoing_index, n, Set{Int64}()))]
    S = Set{Int64}(sinks); stack = copy(sinks)
    while !isempty(stack)
        n = pop!(stack)
        for p in get(incoming_index, n, Set{Int64}())
            (p == v || p in S) && continue
            push!(S, p); push!(stack, p)
        end
    end
    srcs = [n for n in nodes if n != v && isempty(get(incoming_index, n, Set{Int64}()))]
    T = Set{Int64}(srcs); stack = copy(srcs)
    while !isempty(stack)
        n = pop!(stack)
        for s in get(outgoing_index, n, Set{Int64}())
            (s == v || s in T) && continue
            push!(T, s); push!(stack, s)
        end
    end
    (S, T)
end

"""
Exact interval floats by the domination split (longest-path mode). Falls back to the shared
exhaustive sweep when that is cheaper; throws if both exceed `max_runs`. Forward, through and
project values are the exact corner-pair bounds; margins and the criticality classification
are exact.
"""
function interval_analyze_split(
    iteration_sets::Vector{Set{Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    node_values::Dict{Int64,ValueInterval},
    edge_values::Dict{Tuple{Int64,Int64},ValueInterval};
    mode::AnalysisMode = LONGEST_PATH,
    atol::Float64 = 1e-9,
    max_runs::Int = 4_000_000,
    force_split::Bool = false
)::IntervalPathResult
    mode.name == :longest_path ||
        throw(ArgumentError("domination split is proven for LONGEST_PATH only; requested $(mode.name)"))
    any(!is_degenerate(v) for v in values(edge_values)) &&
        throw(ArgumentError("interval edge delays not supported yet; subdivide edges into duration nodes"))

    nodes = sort!(collect(Set(n for layer in iteration_sets for n in layer)))
    varn = sort!([k for (k, v) in node_values if !is_degenerate(v)])
    k = length(varn)
    anc, desc = _closure_sets(iteration_sets, incoming_index, outgoing_index)

    plans = Dict{Int64,Tuple{Vector{Int64},Dict{Int64,Symbol}}}()
    total = 0
    for v in nodes
        S, T = _bypass_sets(v, nodes, outgoing_index, incoming_index)
        H = Int64[]; rule = Dict{Int64,Symbol}()
        for u in varn
            if u == v || (u in anc[v] && !(u in S)) || (u in desc[v] && !(u in T))
                rule[u] = :dominated
            elseif !(u in anc[v]) && !(u in desc[v])
                rule[u] = :incomparable
            else
                push!(H, u)
            end
        end
        plans[v] = (H, rule)
        total += 2 * (1 << length(H))
    end

    if !force_split && k <= 62 && 2.0^k <= total && 2.0^k <= max_runs
        return interval_analyze_exact(iteration_sets, outgoing_index, incoming_index, source_nodes,
                                      node_values, edge_values; mode = mode, atol = atol,
                                      max_corners = max_runs)
    end
    total <= max_runs || throw(SplitDeclined(total, max_runs, k))

    base = interval_analyze(iteration_sets, outgoing_index, incoming_index, source_nodes,
                            node_values, edge_values; mode = mode, atol = atol)
    w0 = Dict{Tuple{Int64,Int64},Float64}(e => v.lo for (e, v) in edge_values)
    d = Dict{Int64,Float64}(kk => vv.lo for (kk, vv) in node_values)
    fl = Dict{Int64,Float64}(); fu = Dict{Int64,Float64}()
    for v in nodes
        H, rule = plans[v]
        fplus = -Inf; fminus = Inf
        for (target_hi, store) in ((true, :plus), (false, :minus))
            for u in varn
                r = get(rule, u, :H)
                if r == :incomparable
                    d[u] = target_hi ? node_values[u].hi : node_values[u].lo
                elseif r == :dominated
                    d[u] = target_hi ? node_values[u].lo : node_values[u].hi
                end
            end
            for mask in 0:(1 << length(H)) - 1
                for (i, u) in enumerate(H)
                    d[u] = (mask >> (i - 1)) & 1 == 1 ? node_values[u].hi : node_values[u].lo
                end
                r = analyze(iteration_sets, outgoing_index, incoming_index, source_nodes, d, w0;
                            mode = mode, atol = atol)
                if store == :plus
                    fplus = max(fplus, r.margin[v])
                else
                    fminus = min(fminus, r.margin[v])
                end
            end
        end
        # fminus and fplus come from two INDEPENDENT corner sweeps (the
        # dominated/incomparable coordinates are pinned oppositely in each),
        # so — unlike a single min/max pass over one stream of values —
        # nothing forces fminus <= fplus by construction. A node whose true
        # margin is exactly 0 (every source node, trivially) can see the two
        # sweeps sum a different arrangement of the same Float64 inputs and
        # land a few ULPs apart, tripping `ValueInterval`'s strict lo<=hi
        # check on a value that is analytically zero-width. Snap to the same
        # `atol` the criticality classification already uses, rather than
        # loosening `ValueInterval` itself — the constructor stays strict,
        # which is where strictness belongs; the fix is that the value
        # crossing the tolerance band is now honestly zero before it gets
        # there, not that the check quietly stops looking.
        abs(fminus) <= atol && (fminus = 0.0)
        abs(fplus) <= atol && (fplus = 0.0)
        (fplus < fminus && fminus - fplus < atol) && (fplus = fminus)
        fl[v] = fminus; fu[v] = fplus
    end

    margin = Dict{Int64,ValueInterval}(v => ValueInterval(fl[v], fu[v]) for v in nodes)
    necessarily = sort!([v for v in nodes if fu[v] <= atol])
    possibly = sort!([v for v in nodes if fl[v] <= atol])

    return IntervalPathResult(mode.name, :exact_domination_split, base.forward, base.through,
                              base.project_value, margin, mode.margin_name,
                              necessarily, possibly, total)
end
