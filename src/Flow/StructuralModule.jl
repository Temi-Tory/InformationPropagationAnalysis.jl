module StructuralModule

if isdefined(parentmodule(@__MODULE__), :FlowModule)
    const FlowModule = parentmodule(@__MODULE__).FlowModule
else
    include("FlowModule.jl")
end
using .FlowModule

include("_CapacityShared.jl")
if isdefined(parentmodule(@__MODULE__), :CapacityTypes)
    const CapacityTypes = parentmodule(@__MODULE__).CapacityTypes
else
    include("CapacityTypes.jl")
end
using .CapacityTypes

export StructuralResult,
       identify_spof_edges,
       identify_spof_nodes,
       enumerate_paths,
       path_flow_contributions,
       bottleneck_ranking,
       node_topological_positions,
       edge_redundancy_scores,
       analyze_structure

struct StructuralResult
    spof_edges::Vector{Tuple{Int64,Int64}}
    spof_nodes::Vector{Int64}
    paths::Vector{Vector{Int64}}
    path_flow_contributions::Vector{PathFlowContribution}
    bottleneck_ranking::Vector{BottleneckRecord}
    node_positions::Dict{Int64,Symbol}
    edge_redundancy::Dict{Tuple{Int64,Int64},Int64}
end

function _forward_reachable(
    starts::Vector{Int64},
    outgoing_index::Dict{Int64,Set{Int64}};
    banned_node::Union{Nothing,Int64}=nothing
)::Set{Int64}
    reachable = Set{Int64}()
    queue = Int64[]

    for s in starts
        if banned_node !== nothing && s == banned_node
            continue
        end
        if !(s in reachable)
            push!(reachable, s)
            push!(queue, s)
        end
    end

    head = 1
    while head <= length(queue)
        u = queue[head]
        head += 1

        for v in get(outgoing_index, u, Set{Int64}())
            if banned_node !== nothing && v == banned_node
                continue
            end
            if !(v in reachable)
                push!(reachable, v)
                push!(queue, v)
            end
        end
    end

    return reachable
end

function _backward_reachable_nodes(
    targets::Vector{Int64},
    incoming_index::Dict{Int64,Set{Int64}};
    banned_node::Union{Nothing,Int64}=nothing
)::Set{Int64}
    reachable = Set{Int64}()
    queue = Int64[]

    for t in targets
        if banned_node !== nothing && t == banned_node
            continue
        end
        if !(t in reachable)
            push!(reachable, t)
            push!(queue, t)
        end
    end

    head = 1
    while head <= length(queue)
        v = queue[head]
        head += 1

        for u in get(incoming_index, v, Set{Int64}())
            if banned_node !== nothing && u == banned_node
                continue
            end
            if !(u in reachable)
                push!(reachable, u)
                push!(queue, u)
            end
        end
    end

    return reachable
end

function _any_sink_reachable(
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    outgoing_index::Dict{Int64,Set{Int64}};
    banned_node::Union{Nothing,Int64}=nothing
)::Bool
    sink_set = Set(sink_nodes)
    visited = Set{Int64}()
    queue = Int64[]

    for s in source_nodes
        if banned_node !== nothing && s == banned_node
            continue
        end
        if s in sink_set
            return true
        end
        if !(s in visited)
            push!(visited, s)
            push!(queue, s)
        end
    end

    head = 1
    while head <= length(queue)
        u = queue[head]
        head += 1

        for v in get(outgoing_index, u, Set{Int64}())
            if banned_node !== nothing && v == banned_node
                continue
            end
            if v in sink_set
                return true
            end
            if !(v in visited)
                push!(visited, v)
                push!(queue, v)
            end
        end
    end

    return false
end

"""
    identify_spof_edges(edgelist, flow_result; tol=1e-10)

An edge is a SPOF iff it appears in every minimum cut. This uses the exact
min-cut lattice characterization and requires no additional solver calls.
"""
function identify_spof_edges(
    edgelist::Vector{Tuple{Int64,Int64}},
    flow_result::FlowSolveResult;
    tol::Float64=1e-10
)::Vector{Tuple{Int64,Int64}}
    _require_bounded_baseline(flow_result)

    original_nodes = _graph_nodes_set(edgelist)
    can_reach_sink = _backward_reachable_residual(
        flow_result.super_sink,
        flow_result.augmented_outgoing,
        flow_result.augmented_incoming,
        flow_result.augmented_capacities,
        flow_result.augmented_flow,
        tol;
        finite_caps_only=true
    )

    S_star = flow_result.mincut_S
    T_double_star = setdiff(original_nodes, can_reach_sink)

    spof_edges = Tuple{Int64,Int64}[]
    for (u, v) in edgelist
        residual = get(flow_result.residual_capacity, (u, v), Inf)
        saturated = residual <= tol
        if saturated && (u in S_star) && (v in T_double_star)
            push!(spof_edges, (u, v))
        end
    end

    sort!(spof_edges)
    return spof_edges
end

"""
    identify_spof_nodes(edgelist, outgoing_index, incoming_index, source_nodes, sink_nodes)

A node is a SPOF iff every source-to-sink path passes through it. Detection is via
reachability on the node-removed graph. No solver calls are required.
"""
function identify_spof_nodes(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64}
)::Vector{Int64}
    graph_nodes = _graph_nodes_set(edgelist)
    source_set = Set(source_nodes)
    sink_set = Set(sink_nodes)

    from_sources = _forward_reachable(source_nodes, outgoing_index)
    to_sinks = _backward_reachable_nodes(sink_nodes, incoming_index)
    on_some_st_path = intersect(from_sources, to_sinks)

    candidates = sort!(collect(setdiff(on_some_st_path, source_set, sink_set)))
    spof_nodes = Int64[]

    for node in candidates
        reachable_sink_exists = _any_sink_reachable(
            source_nodes,
            sink_nodes,
            outgoing_index;
            banned_node=node
        )
        if !reachable_sink_exists
            push!(spof_nodes, node)
        end
    end

    sort!(spof_nodes)
    return spof_nodes
end

"""
    enumerate_paths(outgoing_index, source_nodes, sink_nodes; path_limit=10_000)

Exact enumeration of all simple source-to-sink paths via DFS. Termination is guaranteed
because the input is a DAG. Results are sorted lexicographically.
"""
function enumerate_paths(
    outgoing_index::Dict{Int64,Set{Int64}},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    path_limit::Int=10_000
)::Vector{Vector{Int64}}
    path_limit <= 0 && throw(ArgumentError("path_limit must be positive."))

    sink_set = Set(sink_nodes)
    seen = Set{Tuple}()
    paths = Vector{Vector{Int64}}()

    function dfs(node::Int64, stack::Vector{Int64})
        if node in sink_set
            key = Tuple(stack)
            if !(key in seen)
                if length(paths) >= path_limit
                    throw(ArgumentError("Path enumeration exceeded path_limit=$path_limit."))
                end
                push!(seen, key)
                push!(paths, copy(stack))
            end
            return
        end

        next_nodes = sort!(collect(get(outgoing_index, node, Set{Int64}())))
        for nxt in next_nodes
            push!(stack, nxt)
            dfs(nxt, stack)
            pop!(stack)
        end
    end

    for s in sort!(unique(source_nodes))
        dfs(s, Int64[s])
    end

    sort!(paths; by=p -> Tuple(p))
    return paths
end

"""
    path_flow_contributions(paths, flow_result; tol=1e-10)

For each path, compute the bottleneck flow contribution as the minimum solved flow on
its edges. This reads directly from the solved flow dictionary and requires no solver calls.
"""
function path_flow_contributions(
    paths::Vector{Vector{Int64}},
    flow_result::FlowSolveResult;
    tol::Float64=1e-10
)::Vector{PathFlowContribution}
    _require_bounded_baseline(flow_result)
    contributions = PathFlowContribution[]

    for path in paths
        length(path) >= 2 || throw(ArgumentError("Each path must contain at least two nodes."))

        min_flow = Inf
        bottleneck_edges = Tuple{Int64,Int64}[]

        for i in 1:(length(path) - 1)
            edge = (path[i], path[i + 1])
            edge_flow = get(flow_result.flow, edge, 0.0)

            if edge_flow < min_flow - tol
                min_flow = edge_flow
                empty!(bottleneck_edges)
                push!(bottleneck_edges, edge)
            elseif abs(edge_flow - min_flow) <= tol
                push!(bottleneck_edges, edge)
            end
        end

        sort!(bottleneck_edges)
        bottleneck_edge = first(bottleneck_edges)

        push!(contributions, PathFlowContribution(copy(path), min_flow, bottleneck_edge))
    end

    sort!(contributions; by=x -> (-x.flow_contribution, Tuple(x.path)))
    return contributions
end

"""
    bottleneck_ranking(edgelist, capacities, flow_result)

Return all min-cut crossing edges ranked by ascending capacity (tightest first), with
lexicographic tie-breaking on edge tuples.
"""
function bottleneck_ranking(
    edgelist::Vector{Tuple{Int64,Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    flow_result::FlowSolveResult
)::Vector{BottleneckRecord}
    _require_bounded_baseline(flow_result)

    ranking = BottleneckRecord[]
    for edge in edgelist
        u, v = edge
        if (u in flow_result.mincut_S) && (v in flow_result.mincut_T)
            cap = capacities[edge]
            flow_value = get(flow_result.flow, edge, 0.0)
            residual = get(flow_result.residual_capacity, edge, cap - flow_value)
            push!(ranking, BottleneckRecord(edge, cap, flow_value, residual, 0))
        end
    end

    sort!(ranking; by=x -> (x.capacity, x.edge))

    ranked = BottleneckRecord[]
    for (i, item) in enumerate(ranking)
        push!(ranked, BottleneckRecord(
            item.edge,
            item.capacity,
            item.flow,
            item.residual_capacity,
            i
        ))
    end

    return ranked
end

"""
    node_topological_positions(edgelist, flow_result)

Classify each graph node relative to the solved min-cut partition as `:upstream`,
`:downstream`, or `:on_cut` (most specific for min-cut boundary nodes).
"""
function node_topological_positions(
    edgelist::Vector{Tuple{Int64,Int64}},
    flow_result::FlowSolveResult
)::Dict{Int64,Symbol}
    _require_bounded_baseline(flow_result)

    nodes = sort!(collect(_graph_nodes_set(edgelist)))
    on_cut_nodes = Set{Int64}()

    for (u, v) in edgelist
        if (u in flow_result.mincut_S) && (v in flow_result.mincut_T)
            push!(on_cut_nodes, u)
            push!(on_cut_nodes, v)
        end
    end

    positions = Dict{Int64,Symbol}()
    for node in nodes
        if node in on_cut_nodes
            positions[node] = :on_cut
        elseif node in flow_result.mincut_S
            positions[node] = :upstream
        elseif node in flow_result.mincut_T
            positions[node] = :downstream
        else
            # Invariant note: FlowModule computes mincut_T as the complement of reachable
            # (excluding super nodes) over all augmented nodes that include all original nodes,
            # so every original node should be in mincut_S ∪ mincut_T.
            throw(ArgumentError("Node $node is not classified by min-cut partition."))
        end
    end

    return positions
end

"""
    edge_redundancy_scores(edgelist, outgoing_index, incoming_index, source_nodes, sink_nodes, flow_result;
                           algorithm=:dinic, tol=1e-10, candidate_edges=nothing)

By Menger's theorem, the redundancy score equals the maximum number of edge-disjoint
bypass paths. Unit capacities are used so the integer max-flow value directly counts
disjoint paths. This is exact.

Returns a `Dict` from candidate edge to score. Dict iteration order is not guaranteed.
"""
function edge_redundancy_scores(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    candidate_edges=nothing
)::Dict{Tuple{Int64,Int64},Int64}
    _require_bounded_baseline(flow_result)

    edge_set = Set(edgelist)
    candidates = if candidate_edges === nothing
        copy(edgelist)
    else
        converted = Vector{Tuple{Int64,Int64}}(candidate_edges)
        for edge in converted
            edge in edge_set || throw(ArgumentError("candidate edge $edge is not in edgelist."))
        end
        converted
    end

    unit_caps = Dict{Tuple{Int64,Int64},Float64}(edge => 1.0 for edge in edgelist)
    scores = Dict{Tuple{Int64,Int64},Int64}()

    for edge in sort!(copy(candidates))
        modified = _copy_capacities(unit_caps)
        modified[edge] = 0.0

        rerun = _solve_with_algorithm(
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            modified,
            source_nodes,
            sink_nodes;
            tol=tol,
            validate=true
        )

        # With unit capacities on a finite graph, this should be unreachable in practice,
        # but we keep the guard for defensive correctness against malformed upstream state.
        rerun.is_unbounded && throw(ArgumentError("Redundancy score is undefined because edge-removal solve became unbounded for edge $edge."))

        score_float = rerun.max_flow
        rounded = round(score_float)
        abs(score_float - rounded) <= tol || throw(ArgumentError("Redundancy score for edge $edge is non-integer ($score_float), violating Menger/integrality assumptions."))
        scores[edge] = Int64(rounded)
    end

    return scores
end

"""
    analyze_structure(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes, flow_result;
                      algorithm=:dinic, tol=1e-10, path_limit=10_000,
                      redundancy_candidates=nothing, combination_limit=10_000)

Run the complete structural analysis pipeline and return a typed aggregate result.
"""
function analyze_structure(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    path_limit::Int=10_000,
    redundancy_candidates=nothing,
    combination_limit::Int=10_000
)::StructuralResult
    _require_bounded_baseline(flow_result)
    # API-compatibility note: combination_limit is reserved for potential future
    # combinatorial structural analyses and is intentionally unused in this module.
    combination_limit <= 0 && throw(ArgumentError("combination_limit must be positive."))

    spof_edges = identify_spof_edges(edgelist, flow_result; tol=tol)
    spof_nodes = identify_spof_nodes(edgelist, outgoing_index, incoming_index, source_nodes, sink_nodes)
    paths = enumerate_paths(outgoing_index, source_nodes, sink_nodes; path_limit=path_limit)
    path_contrib = path_flow_contributions(paths, flow_result; tol=tol)
    bottlenecks = bottleneck_ranking(edgelist, capacities, flow_result)
    node_positions = node_topological_positions(edgelist, flow_result)
    edge_redundancy = edge_redundancy_scores(
        edgelist,
        outgoing_index,
        incoming_index,
        source_nodes,
        sink_nodes,
        flow_result;
        algorithm=algorithm,
        tol=tol,
        candidate_edges=redundancy_candidates
    )

    return StructuralResult(
        spof_edges,
        spof_nodes,
        paths,
        path_contrib,
        bottlenecks,
        node_positions,
        edge_redundancy
    )
end

end
