# Helper functions for type-specific operations
# Delegate defaults to shared uncertainty helpers
zero_value(::Type{T}) where {T} = Input.zero_value(T)
one_value(::Type{T}) where {T} = Input.one_value(T)
non_fixed_value(::Type{T}) where {T} = Input.non_fixed_value(T)

# one_value(pbox) MUST be a clean point mass at 1 (interval(1,1)). A previous [1.0, 1.1] "certainty upper
# bound" makes one_value a p-box with mean up to 1.1 > 1, so a conditioned join (whose prior is set to
# one_value) can produce belief mass ABOVE 1 -> unsound. Match Input's interval(1,1).
one_value(::Type{pbox}) = PBA.makepbox(PBA.interval(1.0, 1.0))

"""
Create a unique hash key for a diamond based on edgelist and conditioning_nodes
Much faster than using the full Sets as keys, especially for large diamonds
IMPORTANT: Sorts data to ensure consistent hashing regardless of insertion order
"""
function create_diamond_hash_key(diamond::Diamond)::UInt64
    # Sort edgelist and conditioning_nodes for consistent hashing
    sorted_edgelist = sort(diamond.edgelist)
    sorted_conditioning = sort(collect(diamond.conditioning_nodes))
    return hash((sorted_edgelist, sorted_conditioning))
end

"""
Get the topological level (iteration set index) for a join node
"""
function get_iteration_level(join_node::Int64, iteration_sets::Vector{Set{Int64}})::Int
    for (level, nodes) in enumerate(iteration_sets)
        if join_node in nodes
            return level
        end
    end
    return length(iteration_sets) + 1  # If not found, put at end
end

"""
Find highest iteration set containing any of the given nodes
Returns all nodes that appear in the highest iteration
"""
function find_highest_iteration_nodes(nodes::Set{Int64}, iteration_sets::Vector{Set{Int64}})::Set{Int64}
    highest_iter = -1
    highest_nodes = Set{Int64}()

    # First find the highest iteration
    for (iter, set) in enumerate(iteration_sets)
        intersect_nodes = intersect(nodes, set)
        if !isempty(intersect_nodes)
            highest_iter = max(highest_iter, iter)
        end
    end

    # Then collect all nodes from that iteration
    if highest_iter > 0
        highest_nodes = intersect(nodes, iteration_sets[highest_iter])
    end

    return highest_nodes
end
