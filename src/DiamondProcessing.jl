module DiamondProcessing

using ..InputProcessing
import ProbabilityBoundsAnalysis

# Create aliases to avoid ambiguity
const PBA = ProbabilityBoundsAnalysis
# Type aliases for convenience
const PBAInterval = ProbabilityBoundsAnalysis.Interval
const pbox = ProbabilityBoundsAnalysis.pbox
const Interval = InputProcessing.Interval

# Export all public functions and types
export DiamondsAtNode, Diamond, DiamondComputationData
export identify_and_group_diamonds
export create_diamond_hash_key
export build_unique_diamond_storage
export build_unique_diamond_storage_depth_first_parallel

# 
# STRUCT DEFINITIONS
# 

"""
Represents a diamond structure in the network.
"""
struct Diamond
    relevant_nodes::Set{Int64}
    conditioning_nodes::Set{Int64}
    edgelist::Vector{Tuple{Int64,Int64}}
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
    # All pre-computed subgraph structure (replaces expensive building in network analysis problems)
    sub_outgoing_index::Dict{Int64,Set{Int64}}
    sub_incoming_index::Dict{Int64,Set{Int64}}
    sub_sources::Set{Int64}
    sub_fork_nodes::Set{Int64}
    sub_join_nodes::Set{Int64}
    sub_ancestors::Dict{Int64,Set{Int64}}
    sub_descendants::Dict{Int64,Set{Int64}}
    sub_iteration_sets::Vector{Set{Int64}}
    sub_node_priors::Dict{Int64,T}

    sub_diamond_structures::Dict{Int64,DiamondsAtNode}
    diamond::Diamond
end


struct DiamondWorkItem
    diamond::Diamond
    join_node::Int64
    non_diamond_parents::Set{Int64}
    accumulated_excluded_nodes::Set{Int64}
    is_root_diamond::Bool
    diamond_hash::UInt64
end


struct DiamondOptimizationContext
    ancestor_intersections::Dict{Tuple{Int64,UInt64},Set{Int64}}
    descendant_intersections::Dict{Tuple{Int64,UInt64},Set{Int64}}

    set_intersection_cache::Dict{Tuple{UInt64,UInt64},Set{Int64}}
    set_difference_cache::Dict{Tuple{UInt64,UInt64},Set{Int64}}

    edge_filter_cache::Dict{Tuple{UInt64,UInt64},Vector{Tuple{Int64,Int64}}}

    set_hash_cache::Dict{Set{Int64},UInt64}
end

function DiamondOptimizationContext()
    return DiamondOptimizationContext(
        Dict{Tuple{Int64,UInt64},Set{Int64}}(),
        Dict{Tuple{Int64,UInt64},Set{Int64}}(),
        Dict{Tuple{UInt64,UInt64},Set{Int64}}(),
        Dict{Tuple{UInt64,UInt64},Set{Int64}}(),
        Dict{Tuple{UInt64,UInt64},Vector{Tuple{Int64,Int64}}}(),
        Dict{Set{Int64},UInt64}()
    )
end

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
    cache_key = (min(h1, h2), max(h1, h2)) 

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
    cache_key = (h1, h2)  

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
    ancestors::Dict{Int64,Set{Int64}},
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
    edgelist::Vector{Tuple{Int64,Int64}},
    relevant_nodes::Set{Int64},
    ctx::DiamondOptimizationContext
)::Vector{Tuple{Int64,Int64}}
    edgelist_hash = hash(edgelist)
    relevant_hash = get_set_hash(relevant_nodes, ctx)
    cache_key = (edgelist_hash, relevant_hash)

    if haskey(ctx.edge_filter_cache, cache_key)
        return ctx.edge_filter_cache[cache_key]
    end

    result = Vector{Tuple{Int64,Int64}}()
    for edge in edgelist
        source, target = edge
        if source ∈ relevant_nodes && target ∈ relevant_nodes
            push!(result, edge)
        end
    end

    ctx.edge_filter_cache[cache_key] = result
    return result
end

zero_value(::Type{Float64}) = 0.0
one_value(::Type{Float64}) = 1.0
non_fixed_value(::Type{Float64}) = 0.9
zero_value(::Type{Interval}) = Interval(0.0, 0.0)
one_value(::Type{Interval}) = Interval(1.0, 1.0)
non_fixed_value(::Type{Interval}) = Interval(0.9, 0.9)
zero_value(::Type{pbox}) = PBA.makepbox(PBA.interval(0.0, 0.0))
one_value(::Type{pbox}) = PBA.makepbox(PBA.interval(1.0, 1.1))
non_fixed_value(::Type{pbox}) = PBA.makepbox(PBA.interval(0.9, 0.9))


function create_diamond_hash_key(diamond::Diamond)::UInt64
    return hash((diamond.edgelist, diamond.conditioning_nodes))
end




function get_iteration_level(join_node::Int64, iteration_sets::Vector{Set{Int64}})::Int
    for (level, nodes) in enumerate(iteration_sets)
        if join_node in nodes
            return level
        end
    end
    return length(iteration_sets) + 1  
end

function find_highest_iteration_nodes(nodes::Set{Int64}, iteration_sets::Vector{Set{Int64}})::Set{Int64}
    highest_iter = -1
    highest_nodes = Set{Int64}()

    for (iter, set) in enumerate(iteration_sets)
        intersect_nodes = intersect(nodes, set)
        if !isempty(intersect_nodes)
            highest_iter = max(highest_iter, iter)
        end
    end

    if highest_iter > 0
        highest_nodes = intersect(nodes, iteration_sets[highest_iter])
    end

    return highest_nodes
end


function filter_irrelevant_sources(
    source_nodes::Set{Int64},
    node_priors::Union{Dict{Int64,Float64},Dict{Int64,pbox},Dict{Int64,Interval}},
    exluded_nodes::Set{Int64},
    ctx::DiamondOptimizationContext
)::Set{Int64}
    irrelevant_sources = copy(exluded_nodes)
    if !isempty(node_priors)
        first_val = first(values(node_priors))  

        if isa(first_val, pbox)
            for node in source_nodes
                prior = node_priors[node]
                if (prior.ml == 0.0 && prior.mh == 0.0) || (prior.ml == 1.0 && prior.mh == 1.0)
                    push!(irrelevant_sources, node)
                end
            end
        elseif isa(first_val, Interval)
            for node in source_nodes
                prior = node_priors[node]
                if (prior.lower == 0.0 && prior.upper == 0.0) || (prior.lower == 1.0 && prior.upper == 1.0)
                    push!(irrelevant_sources, node)
                end
            end
        else  # Float64
            for node in source_nodes
                prior_val = node_priors[node]
                if prior_val == 0.0 || prior_val == 1.0
                    push!(irrelevant_sources, node)
                end
            end
        end
    end

    return irrelevant_sources
end

function collect_shared_fork_ancestors(
    join_node::Int64,
    incoming_index::Dict{Int64,Set{Int64}},
    ancestors::Dict{Int64,Set{Int64}},
    fork_nodes::Set{Int64},
    irrelevant_sources::Set{Int64},
    ctx::DiamondOptimizationContext
)::Tuple{Set{Int64},Set{Int64},Set{Int64}}
    parents = get(incoming_index, join_node, Set{Int64}())

    if length(parents) < 2
        return Set{Int64}(), Set{Int64}(), parents
    end

    parent_fork_ancestors = Dict{Int64,Set{Int64}}()
    for parent in parents
        fork_ancestors = get_cached_ancestor_intersection(parent, fork_nodes, ancestors, ctx)
        fork_ancestors = cached_setdiff(fork_ancestors, irrelevant_sources, ctx)
        parent_fork_ancestors[parent] = fork_ancestors
    end

    ancestor_to_parents = Dict{Int64,Set{Int64}}()
    for (parent, fork_ancs) in parent_fork_ancestors
        for ancestor in fork_ancs
            if !haskey(ancestor_to_parents, ancestor)
                ancestor_to_parents[ancestor] = Set{Int64}()
            end
            push!(ancestor_to_parents[ancestor], parent)
        end
    end

    shared_fork_ancestors = Set{Int64}()
    diamond_parents = Set{Int64}()
    for (ancestor, influenced_parents) in ancestor_to_parents
        if length(influenced_parents) >= 2
            push!(shared_fork_ancestors, ancestor)
            union!(diamond_parents, influenced_parents)
        end
    end

    parent_to_parent_forks = Set{Int64}()
    for parent_a in parents
        for parent_b in parents
            if parent_a != parent_b
                parent_b_ancestors = get(ancestors, parent_b, Set{Int64}())
                if parent_a in parent_b_ancestors && parent_a ∉ irrelevant_sources
                    push!(parent_to_parent_forks, parent_a)
                    push!(shared_fork_ancestors, parent_a)
                    push!(diamond_parents, parent_a)
                    push!(diamond_parents, parent_b)
                end
            end
        end
    end

    return shared_fork_ancestors, diamond_parents, parents
end


function extract_induced_edgelist(
    shared_fork_ancestors::Set{Int64},
    join_node::Int64,
    ancestors::Dict{Int64,Set{Int64}},
    descendants::Dict{Int64,Set{Int64}},
    edgelist::Vector{Tuple{Int64,Int64}},
    ctx::DiamondOptimizationContext
)::Tuple{Vector{Tuple{Int64,Int64}},Set{Int64}}
    relevant_nodes_for_induced = copy(shared_fork_ancestors)
    push!(relevant_nodes_for_induced, join_node)

    join_ancestors = get(ancestors, join_node, Set{Int64}())

    for shared_ancestor in shared_fork_ancestors
        shared_descendants = get(descendants, shared_ancestor, Set{Int64}())
        intermediates = cached_intersect(shared_descendants, join_ancestors, ctx)
        union!(relevant_nodes_for_induced, intermediates)
    end

    induced_edgelist = cached_filter_edges(edgelist, relevant_nodes_for_induced, ctx)

    return induced_edgelist, relevant_nodes_for_induced
end

function identify_diamond_sources_and_conditioning(
    induced_edgelist::Vector{Tuple{Int64,Int64}},
    relevant_nodes_for_induced::Set{Int64},
    exluded_nodes::Set{Int64}
)::Tuple{Set{Int64},Set{Int64},Set{Int64}}
    targets_in_induced = Set{Int64}()
    for (_, target) in induced_edgelist
        push!(targets_in_induced, target)
    end
    diamond_sourcenodes = setdiff(setdiff(relevant_nodes_for_induced, targets_in_induced), exluded_nodes)

    relevant_nodes = Set{Int64}()
    for (source, target) in induced_edgelist
        push!(relevant_nodes, source)
        push!(relevant_nodes, target)
    end

    conditioning_nodes = setdiff(diamond_sourcenodes, exluded_nodes)

    return diamond_sourcenodes, relevant_nodes, conditioning_nodes
end

function identify_intermediate_nodes(
    relevant_nodes::Set{Int64},
    conditioning_nodes::Set{Int64},
    join_node::Int64
)::Set{Int64}
    return setdiff(relevant_nodes, union(conditioning_nodes, Set([join_node])))
end

function ensure_intermediate_incoming_edges(
    intermediate_nodes::Set{Int64},
    incoming_index::Dict{Int64,Set{Int64}},
    induced_edgelist::Vector{Tuple{Int64,Int64}},
    relevant_nodes_for_induced::Set{Int64}
)::Tuple{Vector{Tuple{Int64,Int64}},Set{Int64},Set{Int64}}
    final_edgelist = copy(induced_edgelist)
    final_relevant_nodes_for_induced = copy(relevant_nodes_for_induced)
    nodes_added_in_step8 = Set{Int64}()

    for intermediate_node in intermediate_nodes
        incoming_edges = get(incoming_index, intermediate_node, Set{Int64}())

        for source_node in incoming_edges
            edge = (source_node, intermediate_node)
            if edge ∉ final_edgelist
                push!(final_edgelist, edge)
                if source_node ∉ relevant_nodes_for_induced
                    push!(nodes_added_in_step8, source_node)
                end
                push!(final_relevant_nodes_for_induced, source_node)
            end
        end
    end

    return final_edgelist, final_relevant_nodes_for_induced, nodes_added_in_step8
end

function perform_subsource_analysis(
    final_edgelist::Vector{Tuple{Int64,Int64}},
    final_relevant_nodes_for_induced::Set{Int64},
    ancestors::Dict{Int64,Set{Int64}},
    descendants::Dict{Int64,Set{Int64}},
    irrelevant_sources::Set{Int64},
    join_node::Int64,
    exluded_nodes::Set{Int64},
    edgelist::Vector{Tuple{Int64,Int64}},
    ctx::DiamondOptimizationContext
)::Tuple{Vector{Tuple{Int64,Int64}},Set{Int64},Set{Int64}}
  
    subsource_analysis_depth = 0
    max_subsource_depth = 1000
    sources_changed = true

    targets_in_final = Set{Int64}()
    for (_, target) in final_edgelist
        push!(targets_in_final, target)
    end
    final_diamond_sourcenodes = setdiff(setdiff(final_relevant_nodes_for_induced, targets_in_final), exluded_nodes)

    while sources_changed && subsource_analysis_depth < max_subsource_depth
        subsource_analysis_depth += 1
        sources_changed = false

        current_sources = setdiff(final_diamond_sourcenodes, Set([join_node]))

        if length(current_sources) >= 2
            source_ancestors = Dict{Int64,Set{Int64}}()
            for source in current_sources
                source_ancs = get(ancestors, source, Set{Int64}())
                valid_ancestors = cached_setdiff(source_ancs, irrelevant_sources, ctx)
                source_ancestors[source] = valid_ancestors
            end

     
            shared_source_ancestors = Set{Int64}()
            sources_sharing_ancestors = Set{Int64}()

            sources_array = collect(current_sources)
            for i in eachindex(sources_array)
                source_i = sources_array[i]
                haskey(source_ancestors, source_i) || continue

                for j in (i+1):lastindex(sources_array)
                    source_j = sources_array[j]
                    haskey(source_ancestors, source_j) || continue

                    shared = cached_intersect(source_ancestors[source_i], source_ancestors[source_j], ctx)
                    if !isempty(shared)
                        union!(shared_source_ancestors, shared)
                        push!(sources_sharing_ancestors, source_i)
                        push!(sources_sharing_ancestors, source_j)
                    end
                end
            end

            current_sources_array = collect(current_sources)
            for i in eachindex(current_sources_array)
                source_a = current_sources_array[i]
                for j in (i+1):lastindex(current_sources_array)
                    source_b = current_sources_array[j]
                    source_b_ancestors = get(ancestors, source_b, Set{Int64}())
                    source_a_ancestors = get(ancestors, source_a, Set{Int64}())

                    if source_a in source_b_ancestors && source_a ∉ irrelevant_sources
                        push!(shared_source_ancestors, source_a)
                        push!(sources_sharing_ancestors, source_a)
                        push!(sources_sharing_ancestors, source_b)
                    elseif source_b in source_a_ancestors && source_b ∉ irrelevant_sources
                        push!(shared_source_ancestors, source_b)
                        push!(sources_sharing_ancestors, source_a)
                        push!(sources_sharing_ancestors, source_b)
                    end
                end
            end

            if !isempty(shared_source_ancestors)
                earliest_shared = shared_source_ancestors

                if !isempty(earliest_shared)
                    for ancestor in earliest_shared
                        push!(final_relevant_nodes_for_induced, ancestor)

                        if haskey(descendants, ancestor)
                            ancestor_descendants = get(descendants, ancestor, Set{Int64}())
                            join_ancestors = get(ancestors, join_node, Set{Int64}())
                            path_intermediates = cached_intersect(ancestor_descendants, join_ancestors, ctx)
                            union!(final_relevant_nodes_for_induced, path_intermediates)
                        end
                    end

                    new_edges = cached_filter_edges(edgelist, final_relevant_nodes_for_induced, ctx)
                    for edge in new_edges
                        if edge ∉ final_edgelist
                            push!(final_edgelist, edge)
                        end
                    end

                    setdiff!(final_diamond_sourcenodes, sources_sharing_ancestors)
                    union!(final_diamond_sourcenodes, earliest_shared)

                    sources_changed = true
                end
            end
        end

        targets_in_final = Set{Int64}()
        for (_, target) in final_edgelist
            push!(targets_in_final, target)
        end
        sources_in_final = Set{Int64}()
        for (source, _) in final_edgelist
            push!(sources_in_final, source)
        end
        final_diamond_sourcenodes = setdiff(setdiff(sources_in_final, targets_in_final), exluded_nodes)
    end

    return final_edgelist, final_relevant_nodes_for_induced, final_diamond_sourcenodes
end

# Step 8c: Recursive diamond completeness
function perform_recursive_diamond_completeness(
    final_edgelist::Vector{Tuple{Int64,Int64}},
    final_relevant_nodes_for_induced::Set{Int64},
    final_diamond_sourcenodes::Set{Int64},
    final_shared_fork_ancestors::Set{Int64},
    ancestors::Dict{Int64,Set{Int64}},
    descendants::Dict{Int64,Set{Int64}},
    fork_nodes::Set{Int64},
    irrelevant_sources::Set{Int64},
    incoming_index::Dict{Int64,Set{Int64}},
    join_node::Int64,
    exluded_nodes::Set{Int64},
    edgelist::Vector{Tuple{Int64,Int64}},
    ctx::DiamondOptimizationContext
)::Tuple{Vector{Tuple{Int64,Int64}},Set{Int64},Set{Int64},Set{Int64}}
    recursion_depth = 0
    max_recursion_depth = 1000

    if isempty(final_shared_fork_ancestors)
        final_highest_nodes = final_diamond_sourcenodes
    else
        final_highest_nodes = cached_intersect(final_shared_fork_ancestors, final_diamond_sourcenodes, ctx)
    end

    while recursion_depth < max_recursion_depth
        recursion_depth += 1

        diamond_source_fork_ancestors = Dict{Int64,Set{Int64}}()
        for node in final_diamond_sourcenodes
            fork_ancestors = get_cached_ancestor_intersection(node, fork_nodes, ancestors, ctx)
            fork_ancestors = cached_setdiff(fork_ancestors, irrelevant_sources, ctx)
            diamond_source_fork_ancestors[node] = fork_ancestors
        end

        source_ancestor_to_nodes = Dict{Int64,Set{Int64}}()
        for (node, fork_ancs) in diamond_source_fork_ancestors
            for ancestor in fork_ancs
                if !haskey(source_ancestor_to_nodes, ancestor)
                    source_ancestor_to_nodes[ancestor] = Set{Int64}()
                end
                push!(source_ancestor_to_nodes[ancestor], node)
            end
        end

        new_shared_fork_ancestors = Set{Int64}()
        for (ancestor, influenced_nodes) in source_ancestor_to_nodes
            if length(influenced_nodes) >= 2
                push!(new_shared_fork_ancestors, ancestor)
            end
        end
        
        source_to_source_forks = Set{Int64}()
        sources_array = collect(final_diamond_sourcenodes)
        for i in eachindex(sources_array)
            source_a = sources_array[i]
            for j in (i+1):lastindex(sources_array)
                source_b = sources_array[j]
                # Check if source_a is ancestor of source_b or vice versa
                source_b_ancestors = get(ancestors, source_b, Set{Int64}())
                source_a_ancestors = get(ancestors, source_a, Set{Int64}())

                if source_a in source_b_ancestors && source_a ∉ irrelevant_sources
                    push!(source_to_source_forks, source_a)
                elseif source_b in source_a_ancestors && source_b ∉ irrelevant_sources
                    push!(source_to_source_forks, source_b)
                end
            end
        end

        union!(new_shared_fork_ancestors, source_to_source_forks)

        if length(final_diamond_sourcenodes) <= 1
            final_shared_fork_ancestors = Set{Int64}()
            new_shared_fork_ancestors = Set{Int64}()  
        else
            new_shared_fork_ancestors = cached_setdiff(new_shared_fork_ancestors, final_shared_fork_ancestors, ctx)
        end

        current_relevant_nodes = Set{Int64}()
        for (source, target) in final_edgelist
            push!(current_relevant_nodes, source)
            push!(current_relevant_nodes, target)
        end

        current_targets = Set{Int64}()
        for (_, target) in final_edgelist
            push!(current_targets, target)
        end
        current_diamond_sources = setdiff(setdiff(final_relevant_nodes_for_induced, current_targets), exluded_nodes)

        current_intermediate_nodes = setdiff(current_relevant_nodes, union(current_diamond_sources, Set([join_node])))

        edges_added_this_iteration = false
        for intermediate_node in current_intermediate_nodes
            incoming_edges = get(incoming_index, intermediate_node, Set{Int64}())

            for source_node in incoming_edges
                edge = (source_node, intermediate_node)
                if edge ∉ final_edgelist
                    push!(final_edgelist, edge)
                    push!(final_relevant_nodes_for_induced, source_node)
                    edges_added_this_iteration = true
                end
            end
        end

        if isempty(new_shared_fork_ancestors) && !edges_added_this_iteration
            break
        end

        if !isempty(new_shared_fork_ancestors)
            union!(final_shared_fork_ancestors, new_shared_fork_ancestors)
        end

        for shared_ancestor in new_shared_fork_ancestors
            push!(final_relevant_nodes_for_induced, shared_ancestor)
            shared_descendants = get(descendants, shared_ancestor, Set{Int64}())

            for source_node in final_diamond_sourcenodes
                source_node_ancestors = get(ancestors, source_node, Set{Int64}())
                if shared_ancestor in source_node_ancestors
                    path_intermediates = cached_intersect(shared_descendants, source_node_ancestors, ctx)
                    union!(final_relevant_nodes_for_induced, path_intermediates)
                end
            end
        end

        # Extract new induced edges (only if we have new shared fork ancestors)
        if !isempty(new_shared_fork_ancestors)
            new_edges = cached_filter_edges(edgelist, final_relevant_nodes_for_induced, ctx)
            for edge in new_edges
                if edge ∉ final_edgelist
                    push!(final_edgelist, edge)
                end
            end
        end
        
        targets_in_final = Set{Int64}()
        for (_, target) in final_edgelist
            push!(targets_in_final, target)
        end
        final_diamond_sourcenodes = setdiff(setdiff(final_relevant_nodes_for_induced, targets_in_final), exluded_nodes)

        if isempty(final_shared_fork_ancestors)
            final_highest_nodes = final_diamond_sourcenodes
        else
            final_highest_nodes = cached_intersect(final_shared_fork_ancestors, final_diamond_sourcenodes, ctx)
        end

        final_relevant_nodes = Set{Int64}()
        for (source, target) in final_edgelist
            push!(final_relevant_nodes, source)
            push!(final_relevant_nodes, target)
        end

        new_intermediate_nodes = setdiff(final_relevant_nodes, union(final_highest_nodes, Set([join_node])))

        for intermediate_node in new_intermediate_nodes
            incoming_edges = get(incoming_index, intermediate_node, Set{Int64}())

            for source_node in incoming_edges
                edge = (source_node, intermediate_node)
                if edge ∉ final_edgelist
                    push!(final_edgelist, edge)
                    push!(final_relevant_nodes_for_induced, source_node)
                end
            end
        end
    end

    # Check recursion depth limit
    if recursion_depth >= max_recursion_depth
        error("Recursion depth limit ($max_recursion_depth) reached for join node $join_node")
    end

    return final_edgelist, final_relevant_nodes_for_induced, final_shared_fork_ancestors, final_highest_nodes
end



function build_final_diamond_structure(
    final_edgelist::Vector{Tuple{Int64,Int64}},
    final_relevant_nodes_for_induced::Set{Int64},
    final_shared_fork_ancestors::Set{Int64},
    final_highest_nodes::Set{Int64},
    parents::Set{Int64},
    diamond_parents::Set{Int64},
    join_node::Int64,
    exluded_nodes::Set{Int64},
    ctx::DiamondOptimizationContext
)::Tuple{Diamond,Set{Int64}}
 
    final_relevant_nodes = Set{Int64}()
    for (source, target) in final_edgelist
        push!(final_relevant_nodes, source)
        push!(final_relevant_nodes, target)
    end

    targets_in_final = Set{Int64}()
    for (_, target) in final_edgelist
        push!(targets_in_final, target)
    end
    sources_in_final = Set{Int64}()
    for (source, _) in final_edgelist
        push!(sources_in_final, source)
    end
    final_diamond_sourcenodes = setdiff(setdiff(sources_in_final, targets_in_final), exluded_nodes)
   
    if isempty(final_shared_fork_ancestors)
        final_highest_nodes = final_diamond_sourcenodes
    else
        final_highest_nodes = cached_intersect(final_shared_fork_ancestors, final_diamond_sourcenodes, ctx)
    end

    non_diamond_parents = setdiff(parents, diamond_parents)

    diamond = Diamond(final_relevant_nodes, final_highest_nodes, final_edgelist)

    return diamond, non_diamond_parents
end

function perform_hybrid_diamond_lookup(
    sub_join_nodes::Set{Int64},
    current_join::Int64,
    current_excluded_nodes::Set{Int64},
    diamond_cache::Dict{Tuple{Set{Int64},Set{Int64},UInt64},Dict{Int64,DiamondsAtNode}},
    sub_incoming_index::Dict{Int64,Set{Int64}},
    sub_ancestors::Dict{Int64,Set{Int64}},
    sub_descendants::Dict{Int64,Set{Int64}},
    sub_sources::Set{Int64},
    sub_fork_nodes::Set{Int64},
    edgelist::Vector{Tuple{Int64,Int64}},
    sub_node_priors,
    sub_iteration_sets::Vector{Set{Int64}},
    diamond_lookup_table::Dict{Int64,Set{DiamondsAtNode}},
    ctx::DiamondOptimizationContext
)

    successful_lookups = Dict{Int64,DiamondsAtNode}()
    failed_joins = Set{Int64}()

    push!(failed_joins, current_join)


    # Try lookups for non-current joins
    for join in sub_join_nodes
        if join == current_join
            continue  
        end

        if haskey(diamond_lookup_table, join)
            best_candidate = nothing
            best_score = 0

            for candidate in diamond_lookup_table[join]
                score = 0

                # Edge containment check
                candidate_edges = Set(candidate.diamond.edgelist)
                available_edges = Set(edgelist)
                if !issubset(candidate_edges, available_edges)
                    continue
                end
                score += 1

                # Conditioning conflict check
                conflicts = intersect(candidate.diamond.conditioning_nodes, current_excluded_nodes)
                if !isempty(conflicts)
                    score += 3  # Major conflict
                end

                if score < 3 && score > best_score  # Only accept if no major conflicts
                    best_candidate = candidate
                    best_score = score
                end
            end

            if best_candidate !== nothing
                successful_lookups[join] = best_candidate
            else
                push!(failed_joins, join)
            end
        end
    end

    fresh_diamonds_dict = Dict{Int64,DiamondsAtNode}()
    if !isempty(failed_joins)
        fresh_diamonds_dict = identify_and_group_diamonds(
            failed_joins,
            sub_incoming_index,
            sub_ancestors,
            sub_descendants,
            sub_sources,
            sub_fork_nodes,
            edgelist,
            sub_node_priors,
            sub_iteration_sets,
            current_excluded_nodes,
            diamond_cache,
            ctx
        )
    end

    final_diamonds_dict = merge(successful_lookups, fresh_diamonds_dict)

    return final_diamonds_dict
end


function identify_and_group_diamonds(
    join_nodes::Set{Int64},
    incoming_index::Dict{Int64,Set{Int64}},
    ancestors::Dict{Int64,Set{Int64}},
    descendants::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    fork_nodes::Set{Int64},
    edgelist::Vector{Tuple{Int64,Int64}},
    node_priors::Union{Dict{Int64,Float64},Dict{Int64,pbox},Dict{Int64,Interval}},
    iteration_sets::Vector{Set{Int64}},
    exluded_nodes::Set{Int64}=Set{Int64}(),
    DIAMOND_IDENTIFICATION_CACHE::Dict{Tuple{Set{Int64},Set{Int64},UInt64},Dict{Int64,DiamondsAtNode}}=Dict{Tuple{Set{Int64},Set{Int64},UInt64},Dict{Int64,DiamondsAtNode}}(),
    ctx::Union{DiamondOptimizationContext,Nothing}=nothing
)::Dict{Int64,DiamondsAtNode}

    #  CREATE CONTEXT ONLY IF NOT PROVIDED
    if ctx === nothing
        ctx = DiamondOptimizationContext()
    end

    result = Dict{Int64,DiamondsAtNode}()

    irrelevant_sources = filter_irrelevant_sources(source_nodes, node_priors, exluded_nodes, ctx)
    effective_fork_nodes = cached_setdiff(fork_nodes, irrelevant_sources, ctx)
    edgelist_hash = hash(sort(edgelist))
    cache_key = (join_nodes, effective_fork_nodes, edgelist_hash)

    if haskey(DIAMOND_IDENTIFICATION_CACHE, cache_key)
        return DIAMOND_IDENTIFICATION_CACHE[cache_key]
    end

    for join_node in join_nodes
        shared_fork_ancestors, diamond_parents, parents = collect_shared_fork_ancestors(
            join_node, incoming_index, ancestors, fork_nodes, irrelevant_sources, ctx
        )

        length(parents) < 2 && continue
        isempty(shared_fork_ancestors) && continue

        induced_edgelist, relevant_nodes_for_induced = extract_induced_edgelist(
            shared_fork_ancestors, join_node, ancestors, descendants, edgelist, ctx
        )

        diamond_sourcenodes, relevant_nodes, conditioning_nodes = identify_diamond_sources_and_conditioning(
            induced_edgelist, relevant_nodes_for_induced, exluded_nodes
        )

        if isempty(conditioning_nodes)
            continue
        end

        intermediate_nodes = identify_intermediate_nodes(relevant_nodes, conditioning_nodes, join_node)


        final_edgelist, final_relevant_nodes_for_induced, nodes_added_in_step8 = ensure_intermediate_incoming_edges(
            intermediate_nodes, incoming_index, induced_edgelist, relevant_nodes_for_induced
        )


        final_edgelist, final_relevant_nodes_for_induced, final_diamond_sourcenodes = perform_subsource_analysis(
            final_edgelist, final_relevant_nodes_for_induced, ancestors, descendants,
            irrelevant_sources, join_node, exluded_nodes, edgelist, ctx
        )


        final_edgelist, final_relevant_nodes_for_induced, final_shared_fork_ancestors, final_highest_nodes = perform_recursive_diamond_completeness(
            final_edgelist, final_relevant_nodes_for_induced, final_diamond_sourcenodes, shared_fork_ancestors,
            ancestors, descendants, fork_nodes, irrelevant_sources, incoming_index, join_node, exluded_nodes, edgelist, ctx
        )


        diamond, non_diamond_parents = build_final_diamond_structure(
            final_edgelist, final_relevant_nodes_for_induced, final_shared_fork_ancestors, final_highest_nodes,
            parents, diamond_parents, join_node, exluded_nodes, ctx
        )

        result[join_node] = DiamondsAtNode(diamond, non_diamond_parents, join_node)
    end

    DIAMOND_IDENTIFICATION_CACHE[cache_key] = result
    return result
end

function compute_diamond_subgraph_structure(
    diamond::Diamond,
    join_node::Int64,
    node_priors::Dict{Int64,T},
    ancestors::Dict{Int64,Set{Int64}},
    descendants::Dict{Int64,Set{Int64}},
    iteration_sets::Vector{Set{Int64}}
) where {T<:Union{Float64,pbox,Interval}}

    sub_outgoing_index = Dict{Int64,Set{Int64}}()
    sub_incoming_index = Dict{Int64,Set{Int64}}()

    for (i, j) in diamond.edgelist
        push!(get!(sub_outgoing_index, i, Set{Int64}()), j)
        push!(get!(sub_incoming_index, j, Set{Int64}()), i)
    end

    sub_sources = Set{Int64}()
    for node in keys(sub_outgoing_index)
        if !haskey(sub_incoming_index, node) || isempty(sub_incoming_index[node])
            push!(sub_sources, node)
        end
    end

    sub_fork_nodes = Set{Int64}()
    for (node, targets) in sub_outgoing_index
        if length(targets) > 1
            push!(sub_fork_nodes, node)
        end
    end

    sub_join_nodes = Set{Int64}()
    for (node, sources) in sub_incoming_index
        if length(sources) > 1
            push!(sub_join_nodes, node)
        end
    end

    sub_ancestors = Dict{Int64,Set{Int64}}()
    sub_descendants = Dict{Int64,Set{Int64}}()
    for node in diamond.relevant_nodes
        sub_ancestors[node] = Set{Int64}(intersect(ancestors[node], diamond.relevant_nodes))
        sub_descendants[node] = Set{Int64}(intersect(descendants[node], diamond.relevant_nodes))
    end

    sub_iteration_sets = Vector{Set{Int64}}()
    for iter_set in iteration_sets
        filtered_set = Set{Int64}(intersect(iter_set, diamond.relevant_nodes))
        if !isempty(filtered_set)
            push!(sub_iteration_sets, filtered_set)
        end
    end

    sub_node_priors = Dict{Int64,T}()
    for node in diamond.relevant_nodes
        if node ∉ sub_sources
            sub_node_priors[node] = node_priors[node]
            if node == join_node
                sub_node_priors[node] = one_value(T)
            end
        elseif node ∉ diamond.conditioning_nodes
            sub_node_priors[node] = non_fixed_value(T)  # Will be replaced with belief_dict in actual usage network analysis problem
        elseif node ∈ diamond.conditioning_nodes
            sub_node_priors[node] = one_value(T)
        end
    end

    return (sub_outgoing_index, sub_incoming_index, sub_sources, sub_fork_nodes,
        sub_join_nodes, sub_ancestors, sub_descendants, sub_iteration_sets, sub_node_priors)
end




function build_unique_diamond_storage(
    root_diamonds::Dict{Int64,DiamondsAtNode},
    node_priors::Dict{Int64,T},
    ancestors::Dict{Int64,Set{Int64}},
    descendants::Dict{Int64,Set{Int64}},
    iteration_sets::Vector{Set{Int64}}
) where {T<:Union{Float64,pbox,Interval}}

    DIAMOND_IDENTIFICATION_CACHE = Dict{Tuple{Set{Int64},Set{Int64},UInt64},Dict{Int64,DiamondsAtNode}}()
    ctx = DiamondOptimizationContext()
    unique_diamonds = Dict{UInt64,DiamondComputationData{T}}()

    diamond_lookup_table = Dict{Int64,Set{DiamondsAtNode}}()

    processed_diamond_hashes = Set{UInt64}()
    valid_Joins = Set(keys(root_diamonds))
    work_stack = Vector{DiamondWorkItem}()
    total_items_processed = 0
    items_skipped_early = 0

    # Group root diamonds by iteration level
    root_diamonds_by_iteration = Dict{Int64,Vector{Tuple{Int64,DiamondsAtNode}}}()
    for (join_node, diamond_at_node) in root_diamonds
        iteration_level = get_iteration_level(join_node, iteration_sets)
        if !haskey(root_diamonds_by_iteration, iteration_level)
            root_diamonds_by_iteration[iteration_level] = Vector{Tuple{Int64,DiamondsAtNode}}()
        end
        push!(root_diamonds_by_iteration[iteration_level], (join_node, diamond_at_node))
    end

    mainrootstack = sort(collect(keys(root_diamonds_by_iteration)), rev=false)
    # Initialize work stack with root diamonds - PRE-COMPUTE HASHES
    for iteration_level in mainrootstack
        for (join_node, diamond_at_node) in root_diamonds_by_iteration[iteration_level]
            diamond_hash = create_diamond_hash_key(diamond_at_node.diamond)
            push!(work_stack, DiamondWorkItem(
                diamond_at_node.diamond,
                join_node,
                diamond_at_node.non_diamond_parents,
                Set{Int64}(),
                true,
                diamond_hash  #  PRE-COMPUTED HASH
            ))
        end
    end

    while !isempty(work_stack)
        current_item = pop!(work_stack)

   
        if current_item.diamond_hash in processed_diamond_hashes
            items_skipped_early += 1
            total_items_processed += 1

            continue
        end

        push!(processed_diamond_hashes, current_item.diamond_hash)

        current_diamond = current_item.diamond
        join_node = current_item.join_node
        accumulated_excluded_nodes = current_item.accumulated_excluded_nodes
        is_root_diamond = current_item.is_root_diamond

        current_excluded_nodes = union(accumulated_excluded_nodes, current_diamond.conditioning_nodes)

        
        sub_outgoing_index, sub_incoming_index, sub_sources, sub_fork_nodes,
        sub_join_nodes, sub_ancestors, sub_descendants, sub_iteration_sets, sub_node_priors =
            compute_diamond_subgraph_structure(current_diamond, join_node, node_priors, ancestors, descendants, iteration_sets)

        if is_root_diamond
            # ROOT DIAMONDS: Always use full computation for maximal diamond discovery
            sub_diamonds_dict = identify_and_group_diamonds(
                intersect(sub_join_nodes, valid_Joins), 
                sub_incoming_index,
                sub_ancestors,
                sub_descendants,
                sub_sources,
                sub_fork_nodes,
                current_diamond.edgelist,
                sub_node_priors,
                sub_iteration_sets,
                current_excluded_nodes,
                DIAMOND_IDENTIFICATION_CACHE,
                ctx
            )
        else
            # SUB DIAMONDS: Use hybrid optimization with lookup table
            sub_diamonds_dict = perform_hybrid_diamond_lookup(
                intersect(sub_join_nodes, valid_Joins),
                join_node,
                current_excluded_nodes,
                DIAMOND_IDENTIFICATION_CACHE,
                sub_incoming_index,
                sub_ancestors,
                sub_descendants,
                sub_sources,
                sub_fork_nodes,
                current_diamond.edgelist,
                sub_node_priors,
                sub_iteration_sets,
                diamond_lookup_table,
                ctx
            )
        end

        #  FILTER OUT ALREADY PROCESSED SUB-DIAMONDS BEFORE ADDING TO STACK
        sub_diamonds_to_add = []
        filtered_sub_diamonds = Dict{Int64,DiamondsAtNode}()

        for (sub_join_node, sub_diamond_at_node) in sub_diamonds_dict
            sub_diamond = sub_diamond_at_node.diamond
            sub_diamond_hash = create_diamond_hash_key(sub_diamond)

            #  ONLY ADD TO STACK IF NOT ALREADY PROCESSED
            if sub_diamond_hash ∉ processed_diamond_hashes
                sub_work_item = DiamondWorkItem(
                    sub_diamond,
                    sub_join_node,
                    sub_diamond_at_node.non_diamond_parents,
                    current_excluded_nodes,
                    false,
                    sub_diamond_hash  #  PRE-COMPUTED HASH
                )
                push!(sub_diamonds_to_add, sub_work_item)
            end

            # Always include in filtered_sub_diamonds for the current diamond's structure
            filtered_sub_diamonds[sub_join_node] = sub_diamond_at_node
        end

        # Add NEW sub-diamonds to stack
        for sub_item in sub_diamonds_to_add
            push!(work_stack, sub_item)
        end

        # Store the computation data
        computation_data = DiamondComputationData{T}(
            sub_outgoing_index,
            sub_incoming_index,
            sub_sources,
            sub_fork_nodes,
            sub_join_nodes,
            sub_ancestors,
            sub_descendants,
            sub_iteration_sets,
            sub_node_priors,
            filtered_sub_diamonds,
            current_diamond
        )

        unique_diamonds[current_item.diamond_hash] = computation_data

        for (sub_join_node, sub_diamond_at_node) in sub_diamonds_dict
            if !haskey(diamond_lookup_table, sub_join_node)
                diamond_lookup_table[sub_join_node] = Set{DiamondsAtNode}()
            end
            push!(diamond_lookup_table[sub_join_node], sub_diamond_at_node)
        end

        total_items_processed += 1

        if total_items_processed % 500 == 0  

            cache_size = length(ctx.set_intersection_cache) + length(ctx.set_difference_cache) +
                         length(ctx.edge_filter_cache) + length(ctx.ancestor_intersections)

            if cache_size > 10000
                empty!(ctx.set_intersection_cache)
                empty!(ctx.set_difference_cache)
                empty!(ctx.edge_filter_cache)
                empty!(ctx.ancestor_intersections)
                empty!(ctx.descendant_intersections)
                empty!(ctx.set_hash_cache)
            end
         end
    end

    return unique_diamonds
end

"""
Adaptive memory management for large-scale diamond processing.
Clears caches and forces garbage collection based on dataset size and memory pressure.
"""
function manage_memory_adaptive(
    ctx::DiamondOptimizationContext,
    diamonds_processed::Int,
    thread_id::Int,
    force_gc::Bool=false
)
    cache_size = length(ctx.set_intersection_cache) + length(ctx.set_difference_cache) +
                 length(ctx.edge_filter_cache) + length(ctx.ancestor_intersections)

    # Adaptive thresholds based on processing scale
    cache_threshold = diamonds_processed > 1000 ? 500 : 1000
    gc_threshold = diamonds_processed > 2000 ? 100 : 200

    if cache_size > cache_threshold || force_gc
        empty!(ctx.set_intersection_cache)
        empty!(ctx.set_difference_cache)
        empty!(ctx.edge_filter_cache)
        empty!(ctx.ancestor_intersections)
        empty!(ctx.descendant_intersections)
        empty!(ctx.set_hash_cache)

        if diamonds_processed % gc_threshold == 0 || force_gc
            GC.gc()
        end
    end
end


function build_unique_diamond_storage_depth_first_parallel(
    root_diamonds::Dict{Int64,DiamondsAtNode},
    node_priors::Dict{Int64,T},
    ancestors::Dict{Int64,Set{Int64}},
    descendants::Dict{Int64,Set{Int64}},
    iteration_sets::Vector{Set{Int64}}
) where {T<:Union{Float64,pbox,Interval}}

    unique_diamonds = Dict{UInt64,DiamondComputationData{T}}()
    processed_hashes = Set{UInt64}()
    results_lock = Threads.SpinLock()
    processed_hashes_lock = Threads.SpinLock()

    # Shared lookup table for thread-local copying and merging
    shared_diamond_lookup_table = Dict{Int64,Set{DiamondsAtNode}}()
    lookup_table_lock = Threads.SpinLock()

    # Group root diamonds by iteration level 
    root_diamonds_by_iteration = Dict{Int64,Vector{Tuple{Int64,DiamondsAtNode}}}()
    for (join_node, diamond_at_node) in root_diamonds
        iteration_level = get_iteration_level(join_node, iteration_sets)
        if !haskey(root_diamonds_by_iteration, iteration_level)
            root_diamonds_by_iteration[iteration_level] = Vector{Tuple{Int64,DiamondsAtNode}}()
        end
        push!(root_diamonds_by_iteration[iteration_level], (join_node, diamond_at_node))
    end

    # Process iteration levels sequentially (preserves global ordering)
    mainrootstack = sort(collect(keys(root_diamonds_by_iteration)), rev=false)

    for iteration_level in mainrootstack
        level_diamonds = root_diamonds_by_iteration[iteration_level]

        # Force garbage collection between iteration levels for large datasets
        if length(level_diamonds) > 50
            GC.gc()
        end

        # Parallelize root diamonds within each level
        # Each thread gets its own LIFO subtree to process completely
        Threads.@threads for i in eachindex(level_diamonds)
            join_node, diamond_at_node = level_diamonds[i]
            thread_ctx = DiamondOptimizationContext()

            # COPY current shared lookup table as starting point
            thread_local_lookup_table = lock(lookup_table_lock) do
                deepcopy(shared_diamond_lookup_table)
            end

            # Process this entire diamond subtree in LIFO order (sequential within thread)
            thread_results = process_diamond_subtree_sequential_lifo_with_lookup(
                diamond_at_node.diamond,
                join_node,
                diamond_at_node.non_diamond_parents,
                Set{Int64}(),  # accumulated_excluded_nodes
                true,  # is_root_diamond
                thread_ctx,
                node_priors,
                ancestors,
                descendants,
                iteration_sets,
                processed_hashes,
                processed_hashes_lock,
                thread_local_lookup_table  # Pass thread-local lookup table
            )

            # Merge thread results and lookup table discoveries (thread-safe)
            lock(results_lock) do
                merge!(unique_diamonds, thread_results)
            end

            # MERGE lookup table discoveries back to shared table
            lock(lookup_table_lock) do
                for (join_node_key, diamonds_set) in thread_local_lookup_table
                    if !haskey(shared_diamond_lookup_table, join_node_key)
                        shared_diamond_lookup_table[join_node_key] = Set{DiamondsAtNode}()
                    end
                    union!(shared_diamond_lookup_table[join_node_key], diamonds_set)
                end
            end

        end

    end

    return unique_diamonds
end

"""
Process a complete diamond subtree sequentially in LIFO order within a single thread.
This preserves exact dependency resolution order while enabling parallelization across subtrees.
"""
function process_diamond_subtree_sequential_lifo_with_lookup(
    root_diamond::Diamond,
    root_join_node::Int64,
    root_non_diamond_parents::Set{Int64},
    root_accumulated_excluded_nodes::Set{Int64},
    is_root::Bool,
    ctx::DiamondOptimizationContext,
    node_priors,
    ancestors,
    descendants,
    iteration_sets,
    global_processed_hashes::Set{UInt64},
    global_processed_hashes_lock::Threads.SpinLock,
    thread_local_lookup_table::Dict{Int64,Set{DiamondsAtNode}}
)::Dict{UInt64,DiamondComputationData}

    # Each thread maintains its own LIFO stack for its subtree
    local_work_stack = Vector{DiamondWorkItem}()
    local_unique_diamonds = Dict{UInt64,DiamondComputationData}()
    local_processed_hashes = Set{UInt64}()

    # Initialize with root diamond
    root_hash = create_diamond_hash_key(root_diamond)
    push!(local_work_stack, DiamondWorkItem(
        root_diamond,
        root_join_node,
        root_non_diamond_parents,
        root_accumulated_excluded_nodes,
        is_root,
        root_hash
    ))

    # Process LIFO stack sequentially within this thread
    while !isempty(local_work_stack)
        current_item = pop!(local_work_stack) 

        # Check for duplicates (thread-local first, then global)
        if current_item.diamond_hash in local_processed_hashes
            continue
        end

        # Check global processed hashes (thread-safe)
        already_processed = false
        lock(global_processed_hashes_lock) do
            if current_item.diamond_hash in global_processed_hashes
                already_processed = true
            else
                push!(global_processed_hashes, current_item.diamond_hash)
            end
        end

        if already_processed
            continue
        end

        # Mark as locally processed
        push!(local_processed_hashes, current_item.diamond_hash)

        # Extract work item data
        current_diamond = current_item.diamond
        join_node = current_item.join_node
        accumulated_excluded_nodes = current_item.accumulated_excluded_nodes
        is_root_diamond = current_item.is_root_diamond

        current_excluded_nodes = union(accumulated_excluded_nodes, current_diamond.conditioning_nodes)

        sub_outgoing_index, sub_incoming_index, sub_sources, sub_fork_nodes,
        sub_join_nodes, sub_ancestors, sub_descendants, sub_iteration_sets, sub_node_priors =
            compute_diamond_subgraph_structure(
                current_diamond, join_node, node_priors,
                ancestors, descendants, iteration_sets
            )

        # Identify sub-diamonds using hybrid lookup for non-root diamonds
        if is_root_diamond
            # ROOT DIAMONDS: Always use full computation for maximal diamond discovery
            sub_diamonds_dict = identify_and_group_diamonds(
                sub_join_nodes,
                sub_incoming_index,
                sub_ancestors,
                sub_descendants,
                sub_sources,
                sub_fork_nodes,
                current_diamond.edgelist,
                sub_node_priors,
                sub_iteration_sets,
                current_excluded_nodes,
                Dict{Tuple{Set{Int64},Set{Int64},UInt64},Dict{Int64,DiamondsAtNode}}(),  # Fresh cache per thread
                ctx
            )
        else
            # SUB DIAMONDS: Use hybrid optimization with thread-local lookup table
            sub_diamonds_dict = perform_hybrid_diamond_lookup(
                sub_join_nodes,
                join_node,
                current_excluded_nodes,
                Dict{Tuple{Set{Int64},Set{Int64},UInt64},Dict{Int64,DiamondsAtNode}}(),  # Fresh cache per thread
                sub_incoming_index,
                sub_ancestors,
                sub_descendants,
                sub_sources,
                sub_fork_nodes,
                current_diamond.edgelist,
                sub_node_priors,
                sub_iteration_sets,
                thread_local_lookup_table,  # Use thread-local lookup table
                ctx
            )
        end

        # Create computation data
        computation_data = DiamondComputationData{eltype(values(node_priors))}(
            sub_outgoing_index,
            sub_incoming_index,
            sub_sources,
            sub_fork_nodes,
            sub_join_nodes,
            sub_ancestors,
            sub_descendants,
            sub_iteration_sets,
            sub_node_priors,
            sub_diamonds_dict,
            current_diamond
        )

        # Store result locally
        local_unique_diamonds[current_item.diamond_hash] = computation_data

        # POPULATE thread-local lookup table for future hybrid optimization
        for (sub_join_node, sub_diamond_at_node) in sub_diamonds_dict
            if !haskey(thread_local_lookup_table, sub_join_node)
                thread_local_lookup_table[sub_join_node] = Set{DiamondsAtNode}()
            end
            push!(thread_local_lookup_table[sub_join_node], sub_diamond_at_node)
        end

        # Add sub-diamonds to LOCAL LIFO stack (preserves dependency order!)
        sub_work_items = Vector{DiamondWorkItem}()
        for (sub_join_node, sub_diamond_at_node) in sub_diamonds_dict
            sub_diamond_hash = create_diamond_hash_key(sub_diamond_at_node.diamond)

            push!(sub_work_items, DiamondWorkItem(
                sub_diamond_at_node.diamond,
                sub_join_node,
                sub_diamond_at_node.non_diamond_parents,
                current_excluded_nodes,
                false,
                sub_diamond_hash
            ))
        end

        # Add to stack in reverse order to maintain LIFO processing order
        for item in reverse(sub_work_items)
            push!(local_work_stack, item)
        end

        # Adaptive memory management for large-scale processing
        thread_id = Threads.threadid()
        diamonds_processed = length(local_unique_diamonds)

        if diamonds_processed % 50 == 0  # More frequent reporting for large datasets
            manage_memory_adaptive(ctx, diamonds_processed, thread_id)
        end
    end

    thread_id = Threads.threadid()

    return local_unique_diamonds
end

end
