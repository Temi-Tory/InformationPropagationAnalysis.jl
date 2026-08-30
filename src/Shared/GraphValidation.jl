module GraphValidation

export collect_edges,
       validate_non_source_nodes_have_incoming,
       validate_source_nodes_have_no_incoming,
       missing_edge_values,
    missing_values_for_edgelist,
       validate_index_consistency

function collect_edges(outgoing_index::Dict{Int64, Set{Int64}})::Set{Tuple{Int64, Int64}}
    edges = Set{Tuple{Int64, Int64}}()
    for (node, targets) in outgoing_index
        for target in targets
            push!(edges, (node, target))
        end
    end
    return edges
end

function validate_non_source_nodes_have_incoming(
    all_nodes::Set{Int64},
    source_nodes::Set{Int64},
    incoming_index::Dict{Int64, Set{Int64}}
)
    non_source_nodes = setdiff(all_nodes, source_nodes)
    for node in non_source_nodes
        if !haskey(incoming_index, node) || isempty(incoming_index[node])
            throw(ErrorException("Non-source node $node has no incoming edges"))
        end
    end
end

function validate_source_nodes_have_no_incoming(
    source_nodes::Set{Int64},
    incoming_index::Dict{Int64, Set{Int64}}
)
    for source in source_nodes
        if haskey(incoming_index, source) && !isempty(incoming_index[source])
            throw(ErrorException("Source node $source has incoming edges: $(incoming_index[source])"))
        end
    end
end

function missing_edge_values(
    outgoing_index::Dict{Int64, Set{Int64}},
    edge_values::Dict{Tuple{Int64, Int64}, T}
) where {T}
    edges = collect_edges(outgoing_index)
    return setdiff(edges, keys(edge_values))
end

function missing_values_for_edgelist(
    edgelist::Vector{Tuple{Int64, Int64}},
    edge_values::Dict{Tuple{Int64, Int64}, T}
) where {T}
    return [edge for edge in edgelist if !haskey(edge_values, edge)]
end

function validate_index_consistency(
    outgoing_index::Dict{Int64, Set{Int64}},
    incoming_index::Dict{Int64, Set{Int64}}
)
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
end

end # module GraphValidation
