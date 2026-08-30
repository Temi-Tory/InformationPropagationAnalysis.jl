module GraphTraversal

export reachable_from_sources

# Breadth-first reachability from a set of source nodes in a directed graph.
function reachable_from_sources(
    graph::Dict{Int64, Set{Int64}},
    sources::Set{Int64}
)::Set{Int64}
    reachable = Set{Int64}()
    queue = Int64[]

    for source in sources
        if source ∉ reachable
            push!(reachable, source)
            push!(queue, source)
        end
    end

    while !isempty(queue)
        node = popfirst!(queue)

        if haskey(graph, node)
            for neighbor in graph[node]
                if neighbor ∉ reachable
                    push!(reachable, neighbor)
                    push!(queue, neighbor)
                end
            end
        end
    end

    return reachable
end

end # module GraphTraversal
