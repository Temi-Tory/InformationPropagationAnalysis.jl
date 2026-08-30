"""
    update_beliefs_iterative(...) -> Dict{Int64, T}

Main belief propagation function. Computes exact beliefs for all nodes.
Supports Float64, pbox, and Interval types.

Uses bit-masking for inclusion-exclusion (no Combinatorics dependency).
"""
function update_beliefs_iterative(
    edgelist::Vector{Tuple{Int64,Int64}},
    iteration_sets::Vector{Set{Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    node_priors::Dict{Int64,T},
    link_probability::Dict{Tuple{Int64,Int64},T},
    descendants::Dict{Int64, Set{Int64}},
    ancestors::Dict{Int64, Set{Int64}},
    diamond_structures::Dict{Int64, Vector{DiamondsAtNode}},
    join_nodes::Set{Int64},
    fork_nodes::Set{Int64},
    computation_lookup::Dict{UInt64, DiamondComputationData{T}},
    cache::Dict{CacheKey, DiamondCacheEntry{T}} = Dict{CacheKey, DiamondCacheEntry{T}}()
) where {T <: Union{Float64, pbox, Interval}}
    validate_network_data(iteration_sets, outgoing_index, incoming_index, source_nodes, node_priors, link_probability)

    belief_dict = Dict{Int64, T}()

    for node_set in iteration_sets
        for node in node_set
            if node in source_nodes
                belief_dict[node] = node_priors[node]
                continue
            end

            # Collect all sources of belief for this node
            all_beliefs = T[]

            # Process diamond structures if they exist. Each join may carry SEVERAL independent
            # conditioning groups; each contributes one belief that is inclusion-exclusion-combined below
            # (independent groups => 1 - prod(1 - b_g), the factorized OR).
            if haskey(diamond_structures, node)
                for structure in diamond_structures[node]

                    # Calculate beliefs from this diamond group
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

                        # SUM only a genuine single tree path; otherwise the independent non-diamond
                        # parents must combine by inclusion-exclusion, NOT sum (matches the regular-parent
                        # branch below). Using `&&` here: sum iff (not a join) AND (<=1 source ancestor).
                        # With `||` this wrongly summed independent parents at a join with one source
                        # ancestor — visible only once factorization emits multiple non-diamond parents.
                        if !(node in join_nodes) && length(intersect(ancestors[node], source_nodes)) <= 1
                            push!(all_beliefs, sum_values(non_diamond_beliefs))
                        else
                            # Join node / multiple source paths: independent parents -> inclusion-exclusion
                            append!(all_beliefs, non_diamond_beliefs)
                        end
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

"""
    calculate_regular_belief(parents, node, belief_dict, link_probability) -> Vector{T}

Computes belief contributions from regular (non-diamond) parent nodes.
Returns vector of contributions (one per parent).
"""
function calculate_regular_belief(
    parents::Set{Int64},
    node::Int64,
    belief_dict::Dict{Int64, T},
    link_probability::Dict{Tuple{Int64, Int64}, T},
) where {T <: Union{Float64, pbox, Interval}}
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

"""
    inclusion_exclusion(belief_values) -> T

Computes P(A1 U A2 U ... U An) for independent events using the identity:
    P(A1 U ... U An) = 1 - product_i (1 - P(Ai))

This is mathematically equivalent to the alternating-sum inclusion-exclusion
when events are independent (P(Ai ∩ Aj) = P(Ai) × P(Aj)), but each variable
appears only once, avoiding the interval dependency problem that causes
over-wide bounds (and potentially negative/above-1 values) with naive
interval arithmetic on the alternating-sum form.

Complexity: O(n) time, O(1) space (beyond input).
"""
function inclusion_exclusion(belief_values::Vector{T}) where {T <: Union{Float64, pbox, Interval}}
    survival = one_value(T)
    for bv in belief_values
        survival = multiply_values(survival, complement_value(bv))
    end
    return complement_value(survival)
end
