# new_identify — correct-by-construction diamond IDENTIFICATION (replaces the buggy completeness loop +
# hybrid reuse of the old Pipeline). Recursive conditioning with independent-diamond FACTORIZATION:
# at a correlated join, partition the parents into INDEPENDENT groups (disjoint un-conditioned ancestry);
# each correlated group becomes a diamond conditioning on ONE new shared fork; independent groups combine
# by inclusion-exclusion at the join. Recurse per group in context E∪{f} until parents are independent.
#
# CONTEXT-AWARE IDENTITY: a diamond's identity is create_diamond_hash_key = (edgelist, conditioning). When
# an upstream fork is a graph SOURCE, cutting its incoming edges is a no-op, so the SAME (edgelist,{f})
# diamond arises whether that fork is FREE or already FIXED upstream, yet the correct nesting differs. So
# conditioning_nodes = {f} ∪ (E ∩ relevant(edges)) — the new fork PLUS already-conditioned upstream forks
# present in this diamond's edges — which makes the identity encode the context (distinct hkeys per
# context). The extra (already-fixed) conditioning nodes are dropped for free by updateDiamondJoin's
# zero-weight conditioning skip, so each level enumerates only its one free fork (no 2^depth blow-up).
#
# Produces the two objects the propagation consumes:
#   root_diamonds   :: Dict{Int64, Vector{DiamondsAtNode}}   (one entry per independent group per join)
#   unique_diamonds :: Dict{UInt64, DiamondComputationData{T}}
# Reference: validation/rc_core.jl (recursion) + rc_core_factored.jl (factorization). Validated exact vs
# CUDD on 129 random+mutant DAGs + structured networks; fanin-k op-count 2^k -> 2k+1. See ROADMAP.md.

# --- build the sub_* fields of DiamondComputationData from a diamond edgelist ---
function _subgraph_structure(edges::Vector{Tuple{Int64,Int64}}, join_node::Int64,
                             conditioning::Set{Int64}, node_priors::Dict{Int64,T},
                             full_anc::Dict{Int64,Set{Int64}}, full_desc::Dict{Int64,Set{Int64}},
                             iteration_sets::Vector{Set{Int64}}) where {T}
    sub_out = Dict{Int64,Set{Int64}}(); sub_in = Dict{Int64,Set{Int64}}()
    for (i,j) in edges
        push!(get!(sub_out,i,Set{Int64}()), j); push!(get!(sub_in,j,Set{Int64}()), i)
    end
    relevant = Set{Int64}()
    for (i,j) in edges; push!(relevant,i); push!(relevant,j); end
    sub_sources = Set(n for n in relevant if !haskey(sub_in,n) || isempty(sub_in[n]))
    sub_forks = Set(n for (n,t) in sub_out if length(t) > 1)
    sub_joins = Set(n for (n,s) in sub_in if length(s) > 1)
    sub_anc = Dict{Int64,Set{Int64}}(); sub_desc = Dict{Int64,Set{Int64}}()
    for n in relevant
        sub_anc[n]  = intersect(get(full_anc, n, Set{Int64}()), relevant)
        sub_desc[n] = intersect(get(full_desc, n, Set{Int64}()), relevant)
    end
    sub_iter = Vector{Set{Int64}}()
    for s in iteration_sets
        fs = intersect(s, relevant); isempty(fs) || push!(sub_iter, fs)
    end
    # sub_node_priors is REBUILT by updateDiamondJoin (Case 1/2/3); stored copy is ignored.
    sub_np = Dict{Int64,T}(n => node_priors[n] for n in relevant)
    (sub_out, sub_in, sub_sources, sub_forks, sub_joins, sub_anc, sub_desc, sub_iter, sub_np)
end

"""
    new_identify(edgelist, node_priors, link_probs, source_nodes, fork_nodes, join_nodes,
                 ancestors, descendants, iteration_sets)
        -> (root_diamonds::Dict{Int64,Vector{DiamondsAtNode}},
            unique_diamonds::Dict{UInt64,DiamondComputationData{T}})

Correct diamond identification with independent-diamond factorization. Drop-in producer for the two
objects `update_beliefs_iterative` consumes (replaces `identify_and_group_diamonds` +
`build_unique_diamond_storage_depth_first_parallel`). Type-generic over T ∈ {Float64, Interval, pbox}.
`link_probs` is accepted for signature symmetry but not needed by identification.
"""
function new_identify(edgelist::Vector{Tuple{Int64,Int64}}, node_priors::Dict{Int64,T},
                      link_probs::Dict{Tuple{Int64,Int64},T}, source_nodes::Set{Int64},
                      fork_nodes::Set{Int64}, join_nodes::Set{Int64},
                      ancestors::Dict{Int64,Set{Int64}}, descendants::Dict{Int64,Set{Int64}},
                      iteration_sets::Vector{Set{Int64}};
                      is_det_override::Union{Nothing,Function} = nothing) where {T}

    incoming = Dict{Int64,Set{Int64}}()
    for (u,v) in edgelist; push!(get!(incoming,v,Set{Int64}()), u); end

    topi = Dict{Int64,Int}()
    let k = 0
        for s in iteration_sets, n in sort(collect(s)); k += 1; topi[n] = k; end
    end
    # A node may be excluded from conditioning only if its REACHABILITY is deterministic: a dead node
    # (prior 0, never reachable) or a certain SOURCE (prior 1, always reachable). A prior-1 NON-source
    # fork still has UNCERTAIN reachability (via its edges) and MUST remain conditionable — excluding it
    # (the old prior∈{0,1} test) drops all diamonds when every node prior is 1.0 (e.g. the grid benchmark)
    # and yields wrong reliabilities.
    #
    # Type-generic over Float64/Interval/pbox (2026-08-17 fix; was Float64-only, so Interval/pbox NEVER
    # got this exclusion even for genuinely certain sources -- confirmed to inflate diamond count/maxcond
    # substantially on real inputs, e.g. drone concentrated-minimal K=8: 687 diamonds/maxcond=10 measured
    # with pbox priors vs 284/maxcond=9 with Float64-midpoint priors on the IDENTICAL network. Safe by
    # construction, not just empirically: excluding an already-degenerate node from conditioning cannot
    # change the result (Lemma 1, conditional invariance -- P(R_A=c) is 1 for exactly one state and 0 for
    # the rest when A is degenerate, so the total-probability sum collapses to a single term whether or
    # not A is enumerated; conditioning on it is redundant work, never a different answer). Exactness-gate
    # this before trusting it in the manuscript: rerun the corpus and confirm belief VALUES are unchanged
    # (only diamond count/maxcond/timing should move) -- see MASTER_FINDINGS.md corrections ledger row 11.
    _is_zero_val(v::Float64) = v == 0.0
    _is_zero_val(v::Interval) = v.lower == 0.0 && v.upper == 0.0
    _is_zero_val(v::pbox) = v.ml == 0.0 && v.mh == 0.0
    _is_one_val(v::Float64) = v == 1.0
    _is_one_val(v::Interval) = v.lower == 1.0 && v.upper == 1.0
    _is_one_val(v::pbox) = v.ml == 1.0 && v.mh == 1.0
    default_is_det(n) = begin
        v = get(node_priors, n, nothing)
        v === nothing && return false
        _is_zero_val(v) && return true
        (_is_one_val(v) && n in source_nodes) && return true
        false
    end
    # is_det_override swaps the settledness predicate (the module's single value-consulting
    # rule) for consumers with different value semantics, e.g. interval CPM excluding
    # zero-width durations. The recursion is untouched; with no override, behaviour is
    # identical to before this keyword existed.
    is_det(n) = is_det_override === nothing ? default_is_det(n) : Bool(is_det_override(n))

    # un-conditioned influence set of a parent p given E (p + its ancestors, minus conditioned nodes)
    infl(p, E) = begin
        s = Set{Int64}(); (p in E) || push!(s, p)
        for a in get(ancestors, p, Set{Int64}()); (a in E) || push!(s, a); end
        s
    end
    # partition parent list P into INDEPENDENT groups (disjoint un-conditioned ancestry) via union-find
    function components(P, E)
        pv = collect(P); n = length(pv); par = collect(1:n)
        find(i) = (while par[i] != i; par[i] = par[par[i]]; i = par[i]; end; i)
        infls = [infl(p, E) for p in pv]
        for i in 1:n, j in i+1:n
            isempty(intersect(infls[i], infls[j])) || (par[find(i)] = find(j))
        end
        groups = Dict{Int,Vector{Int64}}(); for i in 1:n; push!(get!(groups, find(i), Int64[]), pv[i]); end
        collect(values(groups))
    end
    # highest un-conditioned fork shared by >=2 parents in G (or an asymmetric parent), else nothing
    function shared_fork(G, E)
        length(G) < 2 && return nothing
        cover = Dict{Int64,Int}()
        for p in G, a in get(ancestors, p, Set{Int64}())
            (a in fork_nodes) && !(a in E) && !is_det(a) && (cover[a] = get(cover,a,0)+1)
        end
        for p in G, q in G
            p == q && continue
            (p in get(ancestors,q,Set{Int64}())) && !(p in E) && !is_det(p) && (cover[p] = get(cover,p,0)+1)
        end
        best = nothing; bestk = typemax(Int)
        for (a,c) in cover; c < 2 && continue; get(topi,a,typemax(Int)) < bestk && (bestk = topi[a]; best = a); end
        best
    end

    unique_diamonds = Dict{UInt64, DiamondComputationData{T}}()

    # one diamond for join v restricted to correlated parent-group G, conditioning on fork f, context E
    function group_diamond(v, E, G, f)
        B = union(E, Set([f])); Gset = Set(G)
        RN = Set{Int64}([v]); for p in G; push!(RN, p); union!(RN, get(ancestors, p, Set{Int64}())); end
        # induced edges on RN, cutting into B, and cutting v's incoming from parents OUTSIDE this group
        # (so this diamond computes P(v reached via group G) alone; other groups combine by IE at v).
        edges = Tuple{Int64,Int64}[]
        for (u,w) in edgelist
            (u in RN && w in RN && !(w in B) && !(w == v && !(u in Gset))) && push!(edges, (u,w))
        end
        relevant = Set{Int64}(); for (u,w) in edges; push!(relevant,u); push!(relevant,w); end
        cond = union(Set([f]), intersect(E, relevant))   # context-aware identity (see header)
        diamond = Diamond(RN, cond, edges)
        hkey = create_diamond_hash_key(diamond)
        haskey(unique_diamonds, hkey) && return DiamondsAtNode(diamond, Set{Int64}(), v)
        sub_in = Dict{Int64,Set{Int64}}(); for (u,w) in edges; push!(get!(sub_in,w,Set{Int64}()), u); end
        inner_joins = Set(n for (n,s) in sub_in if length(s) > 1)
        sub_structs = Dict{Int64,Vector{DiamondsAtNode}}()
        so,si,ss,sf,sj,sa,sd,siter,snp = _subgraph_structure(edges, v, cond, node_priors, ancestors, descendants, iteration_sets)
        unique_diamonds[hkey] = DiamondComputationData{T}(so,si,ss,sf,sj,sa,sd,siter,snp, isempty(E), sub_structs, diamond)
        for j in inner_joins
            allowed = (j == v) ? Gset : get(incoming, j, Set{Int64}())   # self-join restricted to this group
            vec = build(j, B, allowed)
            isempty(vec) || (sub_structs[j] = vec)
        end
        DiamondsAtNode(diamond, Set{Int64}(), v)
    end

    # FACTORIZED build: the vector of independent conditioning groups for join v given E, considering only
    # parents in `allowed`. Independent groups combine by inclusion-exclusion at v (framework), so this is
    # exact and avoids joint 2^(cutset). A single correlated group reduces exactly to the non-factored
    # diamond; fanin-k splits into k O(1) groups instead of 2^k.
    function build(v, E, allowed)
        P = [p for p in get(incoming, v, Set{Int64}()) if p in allowed]
        length(P) < 2 && return DiamondsAtNode[]
        diamonds = DiamondsAtNode[]; nondiamond = Set{Int64}()
        for G in components(P, E)
            f = shared_fork(G, E)
            f === nothing ? union!(nondiamond, Set(G)) : push!(diamonds, group_diamond(v, E, G, f))
        end
        # independent (unconditioned) parents ride along as non-diamond parents on the first group
        (!isempty(nondiamond) && !isempty(diamonds)) &&
            (diamonds[1] = DiamondsAtNode(diamonds[1].diamond, nondiamond, v))
        diamonds
    end

    root_diamonds = Dict{Int64,Vector{DiamondsAtNode}}()
    for v in join_nodes
        vec = build(v, Set{Int64}(), get(incoming, v, Set{Int64}()))
        isempty(vec) || (root_diamonds[v] = vec)
    end
    (root_diamonds, unique_diamonds)
end
