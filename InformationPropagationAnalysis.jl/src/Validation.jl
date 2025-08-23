module Validation
using Combinatorics

export MC_result, path_enumeration_result, has_path



function MC_result(
    edgelist::Vector{Tuple{Int64,Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    node_priors::Dict{Int64,Float64},
    edge_probabilities::Dict{Tuple{Int64,Int64},Float64},
    N::Int=100000
)
    # Get all nodes
    all_nodes = reduce(union, values(incoming_index), init=keys(incoming_index))
    active_count = Dict{Int64,Float64}()
    for node in all_nodes
        active_count[node] = 0.0
    end

    for _ in 1:N
        # Sample node states
        node_active = Dict(
            node => rand() < node_priors[node]
            for node in all_nodes
        )

        # Sample edge states
        active_edges = Dict{Tuple{Int64,Int64},Bool}()
        for edge in edgelist
            src, dst = edge
            if node_active[src] && node_active[dst]
                active_edges[edge] = rand() < edge_probabilities[edge]
            else
                active_edges[edge] = false
            end
        end

        # Create subgraph with only active edges
        sub_outgoing = Dict{Int64,Set{Int64}}()
        for (src, dst) in edgelist
            if active_edges[(src, dst)]
                if !haskey(sub_outgoing, src)
                    sub_outgoing[src] = Set{Int64}()
                end
                push!(sub_outgoing[src], dst)
            end
        end

        # Check reachability for each node
        for node in all_nodes
            if node in source_nodes
                if node_active[node]
                    active_count[node] += 1
                end
            else
                # Check if node is reachable from any source
                reachable = false
                for source in source_nodes
                    if has_path(sub_outgoing, source, node)
                        reachable = true
                        break
                    end
                end
                if reachable
                    active_count[node] += 1
                end
            end
        end
    end

    # Convert counts to probabilities
    for node in keys(active_count)
        active_count[node] /= N
    end

    return active_count
end

# Helper function to check if there's a path between two nodes
function has_path(graph::Dict{Int64,Set{Int64}}, start::Int64, target::Int64)
    visited = Set{Int64}()
    queue = [start]

    while !isempty(queue)
        node = popfirst!(queue)
        if node == target
            return true
        end

        if haskey(graph, node)
            for neighbor in graph[node]
                if neighbor ∉ visited
                    push!(visited, neighbor)
                    push!(queue, neighbor)
                end
            end
        end
    end

    return false
end

function path_enumeration_result(
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    node_priors::Dict{Int64,Float64},
    edge_probabilities::Dict{Tuple{Int64,Int64},Float64},
    is_sink_only::Bool=false
)
    # Get all nodes
    all_nodes = reduce(union, values(incoming_index), init=keys(incoming_index))
    active_probability = Dict{Int64,Float64}()

    # For source nodes, the activation probability is simply their prior
    for node in source_nodes
        active_probability[node] = node_priors[node]
    end

    nodes_to_calc = Set{Int64}()
    if is_sink_only
        for node in all_nodes
            if !haskey(outgoing_index, node)
                push!(nodes_to_calc, node)
            end
        end
    else
        nodes_to_calc = setdiff(all_nodes, source_nodes)
    end
    # For non-source nodes, calculate using path enumeration
    for node in nodes_to_calc
        # Step 1: Find all paths from all source nodes to this node
        all_paths = Vector{Vector{Int64}}()
        for source in source_nodes
            paths_from_source = find_all_paths(outgoing_index, source, node)
            append!(all_paths, paths_from_source)
        end

        if isempty(all_paths)
            active_probability[node] = 0.0
            continue
        end

        # Step 2: Convert paths to edge sets for inclusion-exclusion
        path_edge_sets = Vector{Vector{Tuple{Int64,Int64}}}()
        for path in all_paths
            edges = Vector{Tuple{Int64,Int64}}()
            for i in 1:(length(path)-1)
                push!(edges, (path[i], path[i+1]))
            end
            push!(path_edge_sets, edges)
        end

        # Step 3: Calculate node activation probability using path enumeration
        # We need to account for both:
        # - The probability that the target node itself is active
        # - The probability that there's at least one active path to the node

        # Calculate probability that at least one path is active
        path_probability = calculate_with_inclusion_exclusion(
            path_edge_sets,
            edge_probabilities,
            node_priors,
            source_nodes
        )

        # Final probability is the product of:
        # - The probability that the target node itself is active (based on its prior)
        # - The probability that there's at least one active path to the node
        active_probability[node] = node_priors[node] * path_probability
    end

    return active_probability
end

# Find all paths from start to target using DFS
function find_all_paths(
    graph::Dict{Int64,Set{Int64}},
    start::Int64,
    target::Int64
)
    paths = Vector{Vector{Int64}}()
    visited = Set{Int64}()
    current_path = Int64[]

    function dfs(current)
        push!(visited, current)
        push!(current_path, current)

        if current == target
            # Found a path to target
            push!(paths, copy(current_path))
        else
            # Continue exploration
            if haskey(graph, current)
                for neighbor in graph[current]
                    if neighbor ∉ visited
                        dfs(neighbor)
                    end
                end
            end
        end

        # Backtrack
        pop!(current_path)
        delete!(visited, current)
    end

    dfs(start)
    return paths
end


function calculate_with_inclusion_exclusion(
    path_edge_sets::Vector{Vector{Tuple{Int64,Int64}}},
    edge_probabilities::Dict{Tuple{Int64,Int64},Float64},
    node_priors::Dict{Int64,Float64},
    source_nodes::Set{Int64}
)
    n = length(path_edge_sets)
    if n == 0
        return 0.0
    end

    # First, collect all nodes in each path (excluding the target node which is the last node in each path)
    path_node_sets = Vector{Set{Int64}}()
    for (i, path_edges) in enumerate(path_edge_sets)
        nodes = Set{Int64}()
        # Extract all nodes from edges
        for (src, dst) in path_edges
            push!(nodes, src)
            # Don't include the target node here - we'll handle it separately
            if dst != path_edges[end][2]
                push!(nodes, dst)
            end
        end
        push!(path_node_sets, nodes)
    end

    total_probability = 0.0

    # For each non-empty subset of paths (2^n - 1 subsets)
    for mask in 1:(2^n-1)
        subset_edge_sets = Vector{Vector{Tuple{Int64,Int64}}}()
        subset_node_sets = Vector{Set{Int64}}()

        # Collect paths in this subset
        for i in 0:(n-1)
            if (mask & (1 << i)) != 0
                push!(subset_edge_sets, path_edge_sets[i+1])
                push!(subset_node_sets, path_node_sets[i+1])
            end
        end

        # Calculate UNION of edges in this subset
        union_edges = union_of_edge_sets(subset_edge_sets)

        # Calculate UNION of nodes in this subset
        union_nodes = reduce(union, subset_node_sets)

        # Calculate probability for this union
        term_probability = 1.0

        # Include probability of ALL nodes in the paths (excluding target node)
        for node in union_nodes
            term_probability *= node_priors[node]
        end

        # Include probability of ALL edges in the union
        for edge in union_edges
            term_probability *= edge_probabilities[edge]
        end

        # Apply inclusion-exclusion sign rule
        subset_size = count_ones(mask)
        sign = iseven(subset_size) ? -1 : 1

        total_probability += sign * term_probability
    end

    return total_probability
end

# Find union of edge sets
function union_of_edge_sets(edge_sets::Vector{Vector{Tuple{Int64,Int64}}})
    if isempty(edge_sets)
        return Tuple{Int64,Int64}[]
    end

    if length(edge_sets) == 1
        return edge_sets[1]
    end

    # Start with an empty set
    union_edges = Set{Tuple{Int64,Int64}}()

    # Union all edge sets
    for edge_set in edge_sets
        for edge in edge_set
            push!(union_edges, edge)
        end
    end

    return collect(union_edges)
end

end
