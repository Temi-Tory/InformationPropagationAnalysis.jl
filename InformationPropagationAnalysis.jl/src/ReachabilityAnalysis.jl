module ReachabilityAnalysis
using Combinatorics
using ..DiamondProcessing
using ..InputProcessing
import ..InputProcessing: Interval
import ProbabilityBoundsAnalysis

const PBA = ProbabilityBoundsAnalysis
const pbox = ProbabilityBoundsAnalysis.pbox


export update_beliefs_iterative, validate_network_data,
    calculate_regular_belief, inclusion_exclusion,
    updateDiamondJoin, calculate_diamond_groups_belief,
    DiamondCacheEntry, CacheKey, make_cache_key


struct DiamondCacheEntry{T}
    edgelist::Vector{Tuple{Int64,Int64}}
    current_priors::Dict{Int64,T}
    state_beliefs::Dict{Int64,T}
end

struct CacheKey
    diamond_hash::UInt64       
    priors_hash::UInt64   
end

Base.hash(k::CacheKey, h::UInt) = hash((k.diamond_hash, k.priors_hash), h)
Base.:(==)(a::CacheKey, b::CacheKey) = a.diamond_hash == b.diamond_hash && a.priors_hash == b.priors_hash

function make_cache_key(edgelist, current_priors)
    diamond_hash = hash(sort(edgelist))

    # Create a hashable representation of priors based on type
    priors_for_hash = []
    for (node, value) in current_priors
        if isa(value, Float64)
            push!(priors_for_hash, (node, value))
        elseif isa(value, pbox)
            # For pbox, use numeric bounds for hashing
            # Access the actual bounds from the pbox structure
            min_val = minimum(value.u)  # minimum of left bounds
            max_val = maximum(value.d)  # maximum of right bounds
            push!(priors_for_hash, (node, (min_val, max_val)))
        elseif isa(value, Interval)
            # For Interval, use bounds for hashing
            push!(priors_for_hash, (node, (value.lower, value.upper)))
        else
            # Fallback: convert to string
            push!(priors_for_hash, (node, string(value)))
        end
    end

    priors_hash = hash(sort(priors_for_hash))
    return CacheKey(diamond_hash, priors_hash)
end

# Helper functions for type-specific operations
# Zero and one values for different types
zero_value(::Type{Float64}) = 0.0
one_value(::Type{Float64}) = 1.0
non_fixed_value(::Type{Float64}) = 0.9
zero_value(::Type{Interval}) = Interval(0.0, 0.0)
one_value(::Type{Interval}) = Interval(1.0, 1.0)
non_fixed_value(::Type{Interval}) = Interval(0.9, 0.9)
zero_value(::Type{pbox}) = PBA.makepbox(PBA.interval(0.0, 0.0))
one_value(::Type{pbox}) = PBA.makepbox(PBA.interval(1.0, 1.0))
non_fixed_value(::Type{pbox}) = PBA.makepbox(PBA.interval(0.9, 0.9))


# Type-specific probability validation
is_valid_probability(value::Float64) = 0.0 <= value <= 1.0

function is_valid_probability(value::Interval)
    return value.lower >= 0.0 && value.upper <= 1.0
end

function is_valid_probability(value::pbox)
    min_value = PBA.minimum(value)
    max_value = PBA.maximum(value)

    # Handle the case where min/max might return intervals
    min_bound = isa(min_value, PBA.Interval) ? min_value.lo : min_value
    max_bound = isa(max_value, PBA.Interval) ? max_value.hi : max_value

    return min_bound >= 0.0 && max_bound <= 1.0
end

# Type-specific arithmetic operations
# Addition
add_values(a::Float64, b::Float64) = a + b
add_values(a::Interval, b::Interval) = Interval(a.lower + b.lower, a.upper + b.upper)
add_values(a::pbox, b::pbox) = PBA.convIndep(a, b, op=+)

# Multiplication
multiply_values(a::Float64, b::Float64) = a * b
function multiply_values(a::Interval, b::Interval)
    products = [a.lower * b.lower, a.lower * b.upper, a.upper * b.lower, a.upper * b.upper]
    return Interval(minimum(products), maximum(products))
end
multiply_values(a::pbox, b::pbox) = PBA.convIndep(a, b, op=*)

# Complement (1 - value)
complement_value(a::Float64) = 1.0 - a
complement_value(a::Interval) = Interval(1.0 - a.upper, 1.0 - a.lower)
complement_value(a::pbox) = PBA.convIndep(one_value(pbox), a, op=-)

# Subtraction
subtract_values(a::Float64, b::Float64) = a - b
subtract_values(a::Interval, b::Interval) = Interval(a.lower - b.upper, a.upper - b.lower)
subtract_values(a::pbox, b::pbox) = PBA.convIndep(a, b, op=-)

# Sum of vector
sum_values(values::Vector{Float64}) = sum(values)
function sum_values(values::Vector{Interval})
    if isempty(values)
        return zero_value(Interval)
    end
    result = values[1]
    for i in 2:length(values)
        result = add_values(result, values[i])
    end
    return result
end
function sum_values(values::Vector{pbox})
    if isempty(values)
        return zero_value(pbox)
    end
    result = values[1]
    for i in 2:length(values)
        result = add_values(result, values[i])
    end
    return result
end

# Product of vector
prod_values(values::Vector{Float64}) = prod(values)
function prod_values(values::Vector{Interval})
    if isempty(values)
        return one_value(Interval)
    end
    result = values[1]
    for i in 2:length(values)
        result = multiply_values(result, values[i])
    end
    return result
end
function prod_values(values::Vector{pbox})
    if isempty(values)
        return one_value(pbox)
    end
    result = values[1]
    for i in 2:length(values)
        result = multiply_values(result, values[i])
    end
    return result
end

function validate_network_data(
    iteration_sets::Vector{Set{Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    node_priors::Dict{Int64,T},
    link_probability::Dict{Tuple{Int64,Int64},T},
) where {T<:Union{Float64,pbox,Interval}}
    # Collect all nodes from iteration sets
    all_nodes = reduce(union, iteration_sets, init=Set{Int64}())

    # 1. Validate all nodes have priors
    nodes_without_priors = setdiff(all_nodes, keys(node_priors))
    if !isempty(nodes_without_priors)
        throw(ErrorException("The following nodes are missing priors: $nodes_without_priors"))
    end

    # 2. Validate all non-source nodes have incoming edges
    non_source_nodes = setdiff(all_nodes, source_nodes)
    for node in non_source_nodes
        if !haskey(incoming_index, node) || isempty(incoming_index[node])
            throw(ErrorException("Non-source node $node has no incoming edges"))
        end
    end

    # 3. Validate source nodes have no incoming edges
    for source in source_nodes
        if haskey(incoming_index, source) && !isempty(incoming_index[source])
            throw(ErrorException("Source node $source has incoming edges: $(incoming_index[source])"))
        end
    end

    # 4. Validate all edges have probability values
    edges = Set{Tuple{Int64,Int64}}()
    for (node, targets) in outgoing_index
        for target in targets
            push!(edges, (node, target))
        end
    end
    edges_without_probability = setdiff(edges, keys(link_probability))
    if !isempty(edges_without_probability)
        throw(ErrorException("The following edges are missing probability values: $edges_without_probability"))
    end

    # 5. Validate consistency between incoming and outgoing indices
    for (node, targets) in outgoing_index
        for target in targets
            if !haskey(incoming_index, target) || !(node in incoming_index[target])
                throw(ErrorException("Inconsistency found: edge ($node, $target) exists in outgoing_index but not in incoming_index"))
            end
        end
    end
    for (node, sources) in incoming_index
        for source in sources
            if !haskey(outgoing_index, source) || !(node in outgoing_index[source])
                throw(ErrorException("Inconsistency found: edge ($source, $node) exists in incoming_index but not in outgoing_index"))
            end
        end
    end

    # 6. Validate all prior probabilities are between 0 and 1
    invalid_priors = [(node, prior) for (node, prior) in node_priors if !is_valid_probability(prior)]
    if !isempty(invalid_priors)
        throw(ErrorException("The following nodes have invalid prior probabilities (must be between 0 and 1): $invalid_priors"))
    end

    # 7. Validate all probability values are between 0 and 1
    invalid_probabilities = [(edge, rel) for (edge, rel) in link_probability if !is_valid_probability(rel)]
    if !isempty(invalid_probabilities)
        throw(ErrorException("The following edges have invalid probability values (must be between 0 and 1): $invalid_probabilities"))
    end

    # 8. Validate iteration sets contain all nodes exactly once
    nodes_seen = Set{Int64}()
    for set in iteration_sets
        intersection = intersect(nodes_seen, set)
        if !isempty(intersection)
            throw(ErrorException("Nodes $intersection appear in multiple iteration sets"))
        end
        union!(nodes_seen, set)
    end
    if nodes_seen != all_nodes
        missing_nodes = setdiff(all_nodes, nodes_seen)
        extra_nodes = setdiff(nodes_seen, all_nodes)
        error_msg = ""
        if !isempty(missing_nodes)
            error_msg *= "Nodes missing from iteration sets: $missing_nodes. "
        end
        if !isempty(extra_nodes)
            error_msg *= "Extra nodes in iteration sets: $extra_nodes."
        end
        throw(ErrorException(error_msg))
    end
end

function update_beliefs_iterative(
    edgelist::Vector{Tuple{Int64,Int64}},
    iteration_sets::Vector{Set{Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    node_priors::Dict{Int64,T},
    link_probability::Dict{Tuple{Int64,Int64},T},
    descendants::Dict{Int64,Set{Int64}},
    ancestors::Dict{Int64,Set{Int64}},
    diamond_structures::Dict{Int64,DiamondsAtNode},
    join_nodes::Set{Int64},
    fork_nodes::Set{Int64},
    computation_lookup::Dict{UInt64,DiamondComputationData{T}},
    cache::Dict{CacheKey,DiamondCacheEntry{T}}=Dict{CacheKey,DiamondCacheEntry{T}}()  # Default empty cache
) where {T<:Union{Float64,pbox,Interval}}
    validate_network_data(iteration_sets, outgoing_index, incoming_index, source_nodes, node_priors, link_probability)
    belief_dict = Dict{Int64,T}()

    for node_set in iteration_sets
        for node in node_set
            if node in source_nodes
                belief_dict[node] = node_priors[node]
                continue
            end

            # Collect all sources of belief for this node
            all_beliefs = T[]

            # Process diamond structures if they exist
            if haskey(diamond_structures, node)
                structure = diamond_structures[node]

                # Calculate beliefs from diamond groups (now returns array of beliefs)
                diamond_beliefs = calculate_diamond_groups_belief(
                    structure,
                    belief_dict,
                    link_probability,
                    node_priors,
                    ancestors,
                    descendants,
                    iteration_sets,
                    computation_lookup,
                    cache
                )

                push!(all_beliefs, diamond_beliefs)

                # Handle non-diamond parents within the structure
                if !isempty(structure.non_diamond_parents)
                    non_diamond_beliefs = calculate_regular_belief(
                        structure.non_diamond_parents,
                        node,
                        belief_dict,
                        link_probability
                    )

                    # For simple tree paths, just take the sum
                    if !(node in join_nodes) || length(intersect(ancestors[node], source_nodes)) <= 1
                        push!(all_beliefs, sum_values(non_diamond_beliefs))
                    else
                        # For join nodes with multiple paths, use inclusion-exclusion
                        append!(all_beliefs, non_diamond_beliefs)
                    end
                end
            else
                # No diamond structures - handle regular parents
                parents = incoming_index[node]
                probability_from_parents = calculate_regular_belief(
                    parents,
                    node,
                    belief_dict,
                    link_probability
                )

                # Check if this is a join node with multiple paths from sources
                if node in join_nodes || length(intersect(ancestors[node], source_nodes)) > 1
                    # Use inclusion-exclusion for multiple paths
                    append!(all_beliefs, probability_from_parents)
                else
                    # For simple tree paths, just take the sum
                    push!(all_beliefs, sum_values(probability_from_parents))
                end
            end

            # Final combination of all belief sources
            if length(all_beliefs) == 1
                _preprior = all_beliefs[1]
                belief_dict[node] = multiply_values(node_priors[node], _preprior)
            else
                _preprior = inclusion_exclusion(all_beliefs)
                belief_dict[node] = multiply_values(node_priors[node], _preprior)
            end
        end
    end

    return belief_dict
end

function calculate_regular_belief(
    parents::Set{Int64},
    node::Int64,
    belief_dict::Dict{Int64,T},
    link_probability::Dict{Tuple{Int64,Int64},T},
) where {T<:Union{Float64,pbox,Interval}}
    combined_probability_from_parents = T[]
    for parent in parents
        if !haskey(belief_dict, parent)
            throw(ErrorException("Parent node $parent of node $node has no belief value. This indicates a processing order error."))
        end
        parent_belief = belief_dict[parent]

        if !haskey(link_probability, (parent, node))
            throw(ErrorException("No probability defined for edge ($parent, $node)"))
        end
        link_rel = link_probability[(parent, node)]

        push!(combined_probability_from_parents, multiply_values(parent_belief, link_rel))
    end

    return combined_probability_from_parents
end

function inclusion_exclusion(belief_values::Vector{T}) where {T<:Union{Float64,pbox,Interval}}
    combined_belief = zero_value(T)
    num_beliefs = length(belief_values)


    for i in 1:num_beliefs
        for combination in combinations(belief_values, i)
            intersection_probability = prod_values(collect(combination))

            if isodd(i)
                combined_belief = add_values(combined_belief, intersection_probability)
            else
                combined_belief = subtract_values(combined_belief, intersection_probability)
            end
        end
    end
    return combined_belief
end


function updateDiamondJoin(
    conditioning_nodes::Set{Int64},
    join_node::Int64,
    diamond::Diamond,
    link_probability::Dict{Tuple{Int64,Int64},T},
    node_priors::Dict{Int64,T},
    belief_dict::Dict{Int64,T},
    ancestors::Dict{Int64,Set{Int64}},
    descendants::Dict{Int64,Set{Int64}},
    iteration_sets::Vector{Set{Int64}},
    computation_lookup::Dict{UInt64,DiamondComputationData{T}},
    diamond_cache::Dict{CacheKey,DiamondCacheEntry{T}}
) where {T<:Union{Float64,pbox,Interval}}

    diamond_hash_key = create_diamond_hash_key(diamond)


    computation_data = computation_lookup[diamond_hash_key]

    sub_outgoing_index = computation_data.sub_outgoing_index
    sub_incoming_index = computation_data.sub_incoming_index
    fresh_sources = computation_data.sub_sources
    sub_fork_nodes = computation_data.sub_fork_nodes
    sub_join_nodes = computation_data.sub_join_nodes
    sub_ancestors = computation_data.sub_ancestors
    sub_descendants = computation_data.sub_descendants
    sub_iteration_sets = computation_data.sub_iteration_sets
    sub_diamond_structures = computation_data.sub_diamond_structures



    sub_link_probability = Dict{Tuple{Int64,Int64},T}()
    for edge in diamond.edgelist
        sub_link_probability[edge] = link_probability[edge]
    end

    sub_node_priors = Dict{Int64,T}()
    for node in diamond.relevant_nodes
        if node ∉ fresh_sources
            sub_node_priors[node] = node_priors[node]
            if node == join_node
                sub_node_priors[node] = one_value(T)
            end
        elseif node ∉ conditioning_nodes
            sub_node_priors[node] = belief_dict[node]
        elseif node ∈ conditioning_nodes
            sub_node_priors[node] = one_value(T)  
        end
    end

    conditioning_nodes_list = collect(unique(conditioning_nodes))


    final_belief = zero_value(T)

    for state_idx in 0:(2^length(conditioning_nodes_list)-1)
        state_probability = one_value(T)
        conditioning_state = Dict{Int64,T}()

        for (i, node) in enumerate(conditioning_nodes_list)
            original_belief = belief_dict[node]

            if (state_idx & (1 << (i - 1))) != 0
                conditioning_state[node] = one_value(T)
                state_probability = multiply_values(state_probability, original_belief)
            else
                conditioning_state[node] = zero_value(T)
                state_probability = multiply_values(state_probability, complement_value(original_belief))
            end
        end

        current_priors = copy(sub_node_priors)

        for (node, value) in conditioning_state
            current_priors[node] = value
        end

        cache_key = make_cache_key(diamond.edgelist, current_priors)


        if haskey(diamond_cache, cache_key)
            cached_entry = diamond_cache[cache_key]
            state_beliefs = cached_entry.state_beliefs
        else

            state_beliefs = update_beliefs_iterative(
                diamond.edgelist,
                sub_iteration_sets,
                sub_outgoing_index,
                sub_incoming_index,
                fresh_sources,
                current_priors,
                sub_link_probability,
                sub_descendants,
                sub_ancestors,
                sub_diamond_structures,
                sub_join_nodes,
                sub_fork_nodes,
                computation_lookup,
                diamond_cache
            )
            diamond_cache[cache_key] = DiamondCacheEntry(diamond.edgelist, current_priors, state_beliefs)
        end

        join_belief = state_beliefs[join_node]
        final_belief = add_values(final_belief, multiply_values(join_belief, state_probability))
    end


    return final_belief
end

function calculate_diamond_groups_belief(
    diamond::DiamondsAtNode,
    belief_dict::Dict{Int64,T},
    link_probability::Dict{Tuple{Int64,Int64},T},
    node_priors::Dict{Int64,T},
    ancestors::Dict{Int64,Set{Int64}},
    descendants::Dict{Int64,Set{Int64}},
    iteration_sets::Vector{Set{Int64}},
    computation_lookup::Dict{UInt64,DiamondComputationData{T}},
    cache::Dict{CacheKey,DiamondCacheEntry{T}}
) where {T<:Union{Float64,pbox,Interval}}
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


function convert_to_pbox_data(
    node_priors::Dict{Int64,Float64},
    link_probability::Dict{Tuple{Int64,Int64},Float64};
    uncertainty_type::Symbol=:none,  
    uncertainty_value::Float64=0.0
)
    pbox_node_priors = Dict{Int64,pbox}()
    for (node, value) in node_priors
        if uncertainty_type == :interval && uncertainty_value > 0.0
            min_val = max(0.0, value - uncertainty_value)
            max_val = min(1.0, value + uncertainty_value)
            pbox_node_priors[node] = PBA.makepbox(PBA.interval(min_val, max_val))
        elseif uncertainty_type == :normal && uncertainty_value > 0.0
            pbox_node_priors[node] = PBA.normal(value, uncertainty_value)
            if PBA.minimum(pbox_node_priors[node]) < 0 || PBA.maximum(pbox_node_priors[node]) > 1
                left_bound = max(0.0, PBA.minimum(pbox_node_priors[node]))
                right_bound = min(1.0, PBA.maximum(pbox_node_priors[node]))
                pbox_node_priors[node] = PBA.makepbox(PBA.interval(left_bound, right_bound))
            end
        else
            
            pbox_node_priors[node] = PBA.makepbox(PBA.interval(value, value))
        end
    end

    pbox_link_probability = Dict{Tuple{Int64,Int64},pbox}()
    for (edge, value) in link_probability
        if uncertainty_type == :interval && uncertainty_value > 0.0
            min_val = max(0.0, value - uncertainty_value)
            max_val = min(1.0, value + uncertainty_value)
            pbox_link_probability[edge] = PBA.makepbox(PBA.interval(min_val, max_val))
        elseif uncertainty_type == :normal && uncertainty_value > 0.0
            pbox_link_probability[edge] = PBA.normal(value, uncertainty_value)
            if PBA.minimum(pbox_link_probability[edge]) < 0 || PBA.maximum(pbox_link_probability[edge]) > 1
                left_bound = max(0.0, PBA.minimum(pbox_link_probability[edge]))
                right_bound = min(1.0, PBA.maximum(pbox_link_probability[edge]))
                pbox_link_probability[edge] = PBA.makepbox(PBA.interval(left_bound, right_bound))
            end
        else
            pbox_link_probability[edge] = PBA.makepbox(PBA.interval(value, value))
        end
    end

    return pbox_node_priors, pbox_link_probability
end
end