# Is a contextual belief pinned to 0 or 1 (so conditioning on it is a no-op)? Used by the zero-weight
# conditioning skip in updateDiamondJoin. pbox: support ⊆[0,1] with mean bounds ml==mh==0 or ==1 ⇒ a
# point mass at 0 or 1.
_pinned01(b::Float64)  = (b == 0.0 || b == 1.0)
_pinned01(b::Interval) = (b.lower == 0.0 && b.upper == 0.0) || (b.lower == 1.0 && b.upper == 1.0)
_pinned01(b::pbox)     = (b.ml == 0.0 && b.mh == 0.0) || (b.ml == 1.0 && b.mh == 1.0)

"""
    updateDiamondJoin(...) -> T

Computes belief for a diamond join node using conditional expectation
over conditioning states. Uses bit-masking for state enumeration.

Implements: Result = sum_{s=0}^{2^n-1} P(state_s) * Belief(Join | state_s)
"""
function updateDiamondJoin(
    conditioning_nodes::Set{Int64},
    join_node::Int64,
    diamond::Diamond,
    link_probability::Dict{Tuple{Int64,Int64},T},
    node_priors::Dict{Int64,T},
    belief_dict::Dict{Int64,T},
    ancestors::Dict{Int64, Set{Int64}},
    descendants::Dict{Int64, Set{Int64}},
    iteration_sets::Vector{Set{Int64}},
    computation_lookup::Dict{UInt64, DiamondComputationData{T}},
    diamond_cache::Dict{CacheKey, DiamondCacheEntry{T}}
    ) where {T <: Union{Float64, pbox, Interval}}

    # O(1) lookup with hash key
    diamond_hash_key = Diamonds.create_diamond_hash_key(diamond)

    if !haskey(computation_lookup, diamond_hash_key)
        error("Diamond not found in computation_lookup")
    end

    computation_data = computation_lookup[diamond_hash_key]

    # Pre-computed diamond structures - skip expensive graph building
    sub_outgoing_index = computation_data.sub_outgoing_index
    sub_incoming_index = computation_data.sub_incoming_index
    fresh_sources = computation_data.sub_sources
    sub_fork_nodes = computation_data.sub_fork_nodes
    sub_join_nodes = computation_data.sub_join_nodes
    sub_ancestors = computation_data.sub_ancestors
    sub_descendants = computation_data.sub_descendants
    sub_iteration_sets = computation_data.sub_iteration_sets
    sub_diamond_structures = computation_data.sub_diamond_structures

    # Create sub_link_probability for diamond edges only
    sub_link_probability = Dict{Tuple{Int64, Int64}, T}()
    for edge in diamond.edgelist
        sub_link_probability[edge] = link_probability[edge]
    end

    # Zero-weight conditioning skip (EXACT; results identical): a conditioning node whose contextual
    # belief is already pinned to 0/1 by an OUTER conditioning has only one non-zero state, so enumerating
    # it doubles work for nothing. Treat it as an ordinary fixed source (Case 2 gives it its 0/1 belief)
    # and drop it from the enumeration. This is what lets the producer condition on every valid shared
    # fork (incl. ones already fixed upstream) without a 2^|cond| blow-up. Works for all three T:
    #   Float64  : b == 0 or 1
    #   Interval : degenerate [0,0] or [1,1]
    #   pbox     : point mass at 0 or 1  (support ⊆[0,1] with mean bounds ml==mh==0 or ==1 ⇒ point mass)
    # The dropped state carries weight complement(pinned)=0 exactly, so the skip changes no result.
    # (_pinned01 defined at module scope, above.)
    active_conditioning_nodes = Set(n for n in conditioning_nodes if !_pinned01(belief_dict[n]))

    # Build sub_node_priors with contextual beliefs
    sub_node_priors = Dict{Int64, T}()

    for node in diamond.relevant_nodes
        # Case 1: Non-source nodes - use original prior
        if node ∉ fresh_sources
            sub_node_priors[node] = node_priors[node]
            if node == join_node
                sub_node_priors[node] = one_value(T)
            end

        # Case 2: Fresh sources that are NOT (active) conditioning nodes — use contextual belief. This now
        # also covers conditioning nodes already pinned to 0/1 upstream (belief_dict is their fixed value).
        elseif node ∉ active_conditioning_nodes
            sub_node_priors[node] = belief_dict[node]

        # Case 3: Active conditioning nodes - will be set to 0 or 1 per state
        elseif node ∈ active_conditioning_nodes
            sub_node_priors[node] = one_value(T)
        end
    end

    conditioning_nodes_list = collect(unique(active_conditioning_nodes))

    # Phase 1: Compute R(s) = join belief for each conditioning state
    num_states = 2^length(conditioning_nodes_list)
    join_results = Vector{T}(undef, num_states)

    # Per-state enumeration runs SERIALLY. The Threads.@spawn path below is disabled: it has (1) a data
    # race — concurrent writes to the shared diamond_cache Dict under auto threads ("Multiple concurrent
    # writes to Dict detected!"), and (2) spawned tasks use a small stack that overflows on deep nesting.
    # Factorization keeps the conditioning set (and hence num_states = 2^|active cond|) small, so serial
    # is fast (sub-second on the corpus) and, unlike the parallel path, correct and deterministic under
    # any -t. Re-enabling would require a thread-safe cache and larger task stacks. See ROADMAP.md.
    use_parallel = false

    if use_parallel
        tasks = Vector{Task}(undef, num_states)

        for state_idx in 0:(num_states - 1)
            tasks[state_idx + 1] = Threads.@spawn begin
                # Thread-local copy to avoid race conditions
                local_sub_node_priors = copy(sub_node_priors)

                for (i, node) in enumerate(conditioning_nodes_list)
                    if (state_idx & (1 << (i-1))) != 0
                        local_sub_node_priors[node] = one_value(T)
                    else
                        local_sub_node_priors[node] = zero_value(T)
                    end
                end

                cache_key = make_cache_key(diamond.edgelist, local_sub_node_priors)

                # Cache lookup (thread-safe)
                local state_beliefs
                lock(diamond_cache_lock) do
                    if haskey(diamond_cache, cache_key)
                        cached_entry = diamond_cache[cache_key]
                        state_beliefs = cached_entry.state_beliefs
                    else
                        state_beliefs = nothing
                    end
                end

                # Recursive computation if not cached
                if state_beliefs === nothing
                    state_beliefs = update_beliefs_iterative(
                        diamond.edgelist,
                        sub_iteration_sets,
                        sub_outgoing_index,
                        sub_incoming_index,
                        fresh_sources,
                        local_sub_node_priors,
                        sub_link_probability,
                        sub_descendants,
                        sub_ancestors,
                        sub_diamond_structures,
                        sub_join_nodes,
                        sub_fork_nodes,
                        computation_lookup,
                        diamond_cache
                    )

                    lock(diamond_cache_lock) do
                        if !haskey(diamond_cache, cache_key)
                            diamond_cache[cache_key] = DiamondCacheEntry(diamond.edgelist, local_sub_node_priors, state_beliefs)
                        end
                    end
                end

                state_beliefs[join_node]
            end
        end

        # Collect results from all parallel tasks
        for idx in 1:num_states
            join_results[idx] = fetch(tasks[idx])
        end
    else
        for state_idx in 0:(num_states - 1)
            local_sub_node_priors = copy(sub_node_priors)

            for (i, node) in enumerate(conditioning_nodes_list)
                if (state_idx & (1 << (i-1))) != 0
                    local_sub_node_priors[node] = one_value(T)
                else
                    local_sub_node_priors[node] = zero_value(T)
                end
            end

            cache_key = make_cache_key(diamond.edgelist, local_sub_node_priors)

            # Lean-aware lookup. A LEAN entry holds only its own diamond's join belief, so a key
            # hit must be re-validated for THIS join_node: the same (edgelist, priors) key can in
            # principle arise from a different-join diamond (e.g. all-ones priors make the forced
            # join prior indistinguishable). get() treats that as a miss and the recomputed belief
            # is MERGED into the existing entry so both joins hit next time. Under legacy (full)
            # entries the guard never fires (a full dict covers every subgraph node incl. the
            # join), so default-OFF behavior and cache contents are identical to before.
            local join_belief
            cached_entry = get(diamond_cache, cache_key, nothing)
            jb = cached_entry === nothing ? nothing : get(cached_entry.state_beliefs, join_node, nothing)
            if jb !== nothing
                join_belief = jb
            else
                state_beliefs = update_beliefs_iterative(
                    diamond.edgelist,
                    sub_iteration_sets,
                    sub_outgoing_index,
                    sub_incoming_index,
                    fresh_sources,
                    local_sub_node_priors,
                    sub_link_probability,
                    sub_descendants,
                    sub_ancestors,
                    sub_diamond_structures,
                    sub_join_nodes,
                    sub_fork_nodes,
                    computation_lookup,
                    diamond_cache
                )
                join_belief = state_beliefs[join_node]
                if cached_entry !== nothing
                    cached_entry.state_beliefs[join_node] = join_belief
                elseif LEAN_DIAMOND_CACHE[]
                    diamond_cache[cache_key] = DiamondCacheEntry(diamond.edgelist, Dict{Int64,T}(), Dict{Int64,T}(join_node => join_belief))
                else
                    diamond_cache[cache_key] = DiamondCacheEntry(diamond.edgelist, local_sub_node_priors, state_beliefs)
                end
            end

            join_results[state_idx + 1] = join_belief
        end
    end

    # Phase 2: Combine R(s) with conditioning state probabilities.
    if T <: Interval
        m = length(conditioning_nodes_list)
        num_corners = 2^m
        min_lo = Inf
        max_hi = -Inf

        for corner_idx in 0:(num_corners - 1)
            # Fix each conditioning belief to its lower or upper bound
            corner_values = Vector{Float64}(undef, m)
            for (i, node) in enumerate(conditioning_nodes_list)
                bel = belief_dict[node]::Interval
                corner_values[i] = (corner_idx & (1 << (i-1))) != 0 ? bel.upper : bel.lower
            end

            # Compute weighted sum with scalar weights (sum to 1.0 exactly)
            lo_sum = 0.0
            hi_sum = 0.0
            for state_idx in 0:(num_states - 1)
                weight = 1.0
                for i in 1:m
                    if (state_idx & (1 << (i-1))) != 0
                        weight *= corner_values[i]
                    else
                        weight *= (1.0 - corner_values[i])
                    end
                end
                R = join_results[state_idx + 1]::Interval
                lo_sum += weight * R.lower
                hi_sum += weight * R.upper
            end

            min_lo = min(min_lo, lo_sum)
            max_hi = max(max_hi, hi_sum)
        end

        return Interval(min_lo, max_hi)
    elseif T <: pbox
        # p-box conditioning is a CONVEX COMBINATION, not a convolution. The old flat convIndep weighted
        # sum (kept below for Float64) treats it as a sum of independent RVs -> over-wide, mass>1, UNSOUND.
        # Combine ONE conditioning node at a time (nested 2-way): belief = W*A + (1-W)*B with W = node's
        # contextual belief, via the sound cvxP operator pbox_conditional_combine. m nested 2-way combos
        # reproduce the 2^m total-probability mixture. (Matches the validated reference; see
        # validation/rc_pbox_cvx.jl + memory pbox-conditioning-unsound.) m=0 -> the single join_results[1].
        m = length(conditioning_nodes_list)
        function _combine(i::Int, base_idx::Int)
            i > m && return join_results[base_idx + 1]
            bit = 1 << (i - 1)
            up   = _combine(i + 1, base_idx | bit)   # node i reachable
            down = _combine(i + 1, base_idx)          # node i not reachable
            return pbox_conditional_combine(belief_dict[conditioning_nodes_list[i]], up, down)
        end
        return _combine(1, 0)
    else
        # Float64: exact weighted sum (scalar total probability)
        final_belief = zero_value(T)
        for state_idx in 0:(num_states - 1)
            state_probability = one_value(T)
            for (i, node) in enumerate(conditioning_nodes_list)
                original_belief = belief_dict[node]
                if (state_idx & (1 << (i-1))) != 0
                    state_probability = multiply_values(state_probability, original_belief)
                else
                    state_probability = multiply_values(state_probability, complement_value(original_belief))
                end
            end
            final_belief = add_values(final_belief, multiply_values(join_results[state_idx + 1], state_probability))
        end
        return final_belief
    end
end

function calculate_diamond_groups_belief(
    diamond::DiamondsAtNode,
    belief_dict::Dict{Int64,T},
    link_probability::Dict{Tuple{Int64,Int64},T},
    node_priors::Dict{Int64,T},
    ancestors::Dict{Int64, Set{Int64}},
    descendants::Dict{Int64, Set{Int64}},
    iteration_sets::Vector{Set{Int64}},
    computation_lookup::Dict{UInt64, DiamondComputationData{T}},
    cache::Dict{CacheKey, DiamondCacheEntry{T}}
) where {T <: Union{Float64, pbox, Interval}}
    diamond_beliefs = updateDiamondJoin(
        diamond.diamond.conditioning_nodes,
        diamond.join_node,
        diamond.diamond,
        link_probability,
        node_priors,
        belief_dict,
        ancestors,
        descendants,
        iteration_sets,
        computation_lookup,
        cache
    )
    return diamond_beliefs
end
