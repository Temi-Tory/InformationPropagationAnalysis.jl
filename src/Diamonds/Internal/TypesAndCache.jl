# Struct definitions

"""
Represents a diamond structure in the network.
"""
struct Diamond
    relevant_nodes::Set{Int64}
    conditioning_nodes::Set{Int64}
    edgelist::Vector{Tuple{Int64, Int64}}
end

"""
Represents diamonds and non-diamond parents at a specific join node.
"""
struct DiamondsAtNode
    diamond::Diamond
    non_diamond_parents::Set{Int64}
    join_node::Int64
end

# Hash and equality methods for DiamondsAtNode to enable Set-based lookup tables
function Base.hash(d::DiamondsAtNode, h::UInt)
    return hash((
        d.diamond.edgelist,
        d.diamond.relevant_nodes,
        d.diamond.conditioning_nodes,
        d.non_diamond_parents,
        d.join_node
    ), h)
end

function Base.:(==)(d1::DiamondsAtNode, d2::DiamondsAtNode)
    return d1.diamond.edgelist == d2.diamond.edgelist &&
           d1.diamond.relevant_nodes == d2.diamond.relevant_nodes &&
           d1.diamond.conditioning_nodes == d2.diamond.conditioning_nodes &&
           d1.non_diamond_parents == d2.non_diamond_parents &&
           d1.join_node == d2.join_node
end

"""
Computation-ready data for a diamond - contains all pre-computed subgraph structure
"""
struct DiamondComputationData{T}
    # All pre-computed subgraph structure (replaces expensive building in updateDiamondJoin)
    sub_outgoing_index::Dict{Int64, Set{Int64}}
    sub_incoming_index::Dict{Int64, Set{Int64}}
    sub_sources::Set{Int64}
    sub_fork_nodes::Set{Int64}
    sub_join_nodes::Set{Int64}
    sub_ancestors::Dict{Int64, Set{Int64}}
    sub_descendants::Dict{Int64, Set{Int64}}
    sub_iteration_sets::Vector{Set{Int64}}
    sub_node_priors::Dict{Int64, T}
    is_rootDiamond::Bool
    # Ready-to-use inner diamonds for recursive calls. Vector per join: one DiamondsAtNode per
    # INDEPENDENT conditioning group (parents with disjoint un-conditioned ancestry). Independent groups
    # combine by inclusion-exclusion at the join (see update_beliefs_iterative), which factorizes the OR
    # and avoids joint 2^(cutset) conditioning. Degenerate case = a 1-element vector (no factorization).
    sub_diamond_structures::Dict{Int64, Vector{DiamondsAtNode}}
    diamond::Diamond

end

# Optimization statistics tracking
struct OptimizationStats
    lookups_attempted::Int
    lookups_successful::Int
    joins_looked_up::Int
    joins_computed_fresh::Int
    computation_reduction_percent::Float64
end

# Work item to replace recursive function calls
struct DiamondWorkItem
    diamond::Diamond
    join_node::Int64
    non_diamond_parents::Set{Int64}
    accumulated_excluded_nodes::Set{Int64}
    is_root_diamond::Bool
    diamond_hash::UInt64
end

# Local optimization context (instantiated per function call)
struct DiamondOptimizationContext
    # Cache repeated ancestor/descendant intersections
    ancestor_intersections::Dict{Tuple{Int64, UInt64}, Set{Int64}}
    descendant_intersections::Dict{Tuple{Int64, UInt64}, Set{Int64}}

    # Cache set operations results
    set_intersection_cache::Dict{Tuple{UInt64, UInt64}, Set{Int64}}
    set_difference_cache::Dict{Tuple{UInt64, UInt64}, Set{Int64}}

    # Cache edge filtering results
    edge_filter_cache::Dict{Tuple{UInt64, UInt64}, Vector{Tuple{Int64, Int64}}}

    # Pre-computed hash values for frequently used sets
    set_hash_cache::Dict{Set{Int64}, UInt64}
end

function DiamondOptimizationContext()
    return DiamondOptimizationContext(
        Dict{Tuple{Int64, UInt64}, Set{Int64}}(),
        Dict{Tuple{Int64, UInt64}, Set{Int64}}(),
        Dict{Tuple{UInt64, UInt64}, Set{Int64}}(),
        Dict{Tuple{UInt64, UInt64}, Set{Int64}}(),
        Dict{Tuple{UInt64, UInt64}, Vector{Tuple{Int64, Int64}}}(),
        Dict{Set{Int64}, UInt64}()
    )
end

# Cached set operations
function get_set_hash(s::Set{Int64}, ctx::DiamondOptimizationContext)::UInt64
    if haskey(ctx.set_hash_cache, s)
        return ctx.set_hash_cache[s]
    end
    h = hash(s)
    ctx.set_hash_cache[s] = h
    return h
end

function cached_intersect(set1::Set{Int64}, set2::Set{Int64}, ctx::DiamondOptimizationContext)::Set{Int64}
    h1 = get_set_hash(set1, ctx)
    h2 = get_set_hash(set2, ctx)
    cache_key = (min(h1, h2), max(h1, h2))  # Order-independent key

    if haskey(ctx.set_intersection_cache, cache_key)
        return ctx.set_intersection_cache[cache_key]
    end

    result = intersect(set1, set2)
    ctx.set_intersection_cache[cache_key] = result
    return result
end

function cached_setdiff(set1::Set{Int64}, set2::Set{Int64}, ctx::DiamondOptimizationContext)::Set{Int64}
    h1 = get_set_hash(set1, ctx)
    h2 = get_set_hash(set2, ctx)
    cache_key = (h1, h2)  # Order-dependent for setdiff

    if haskey(ctx.set_difference_cache, cache_key)
        return ctx.set_difference_cache[cache_key]
    end

    result = setdiff(set1, set2)
    ctx.set_difference_cache[cache_key] = result
    return result
end

function get_cached_ancestor_intersection(
    node::Int64,
    target_set::Set{Int64},
    ancestors::Dict{Int64, Set{Int64}},
    ctx::DiamondOptimizationContext
)::Set{Int64}
    target_hash = get_set_hash(target_set, ctx)
    cache_key = (node, target_hash)

    if haskey(ctx.ancestor_intersections, cache_key)
        return ctx.ancestor_intersections[cache_key]
    end

    node_ancestors = get(ancestors, node, Set{Int64}())
    result = cached_intersect(node_ancestors, target_set, ctx)
    ctx.ancestor_intersections[cache_key] = result
    return result
end

function cached_filter_edges(
    edgelist::Vector{Tuple{Int64, Int64}},
    relevant_nodes::Set{Int64},
    ctx::DiamondOptimizationContext
)::Vector{Tuple{Int64, Int64}}
    edgelist_hash = hash(edgelist)
    relevant_hash = get_set_hash(relevant_nodes, ctx)
    cache_key = (edgelist_hash, relevant_hash)

    if haskey(ctx.edge_filter_cache, cache_key)
        return ctx.edge_filter_cache[cache_key]
    end

    result = Vector{Tuple{Int64, Int64}}()
    for edge in edgelist
        source, target = edge
        if source ∈ relevant_nodes && target ∈ relevant_nodes
            push!(result, edge)
        end
    end

    ctx.edge_filter_cache[cache_key] = result
    return result
end
