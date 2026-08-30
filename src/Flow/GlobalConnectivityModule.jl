module GlobalConnectivityModule

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

include("NodeCapacitatedFlowModule.jl")
using .NodeCapacitatedFlowModule

export EdgeConnectivityResult,
       NodeConnectivityResult,
       GlobalMinCutResult,
       GlobalConnectivityResult,
       edge_connectivity,
       node_connectivity,
       global_min_cut,
       analyze_global_connectivity

struct EdgeConnectivityResult
    lambda::Int64
    achieving_source::Int64
    achieving_sink::Int64
    min_cut_edges::Vector{Tuple{Int64,Int64}}
    solver_calls::Int64
end

struct NodeConnectivityResult
    kappa::Int64
    achieving_source::Int64
    achieving_sink::Int64
    min_cut_nodes::Vector{Int64}
    solver_calls::Int64
end

struct GlobalMinCutResult
    min_cut_capacity::Float64
    achieving_source::Int64
    achieving_sink::Int64
    min_cut_edges::Vector{Tuple{Int64,Int64}}
    cut_S::Set{Int64}
    cut_T::Set{Int64}
    solver_calls::Int64
end

struct GlobalConnectivityResult
    edge_connectivity::EdgeConnectivityResult
    node_connectivity::NodeConnectivityResult
    global_min_cut::GlobalMinCutResult
end

function _require_valid_inputs(
    edgelist::Vector{Tuple{Int64,Int64}},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64}
)::Set{Int64}
    all_nodes = _graph_nodes_set(edgelist)
    length(all_nodes) >= 2 || throw(ArgumentError("Global connectivity requires at least 2 nodes."))
    !isempty(source_nodes) || throw(ArgumentError("source_nodes must be non-empty."))
    !isempty(sink_nodes) || throw(ArgumentError("sink_nodes must be non-empty."))

    for s in source_nodes
        s in all_nodes || throw(ArgumentError("source node $s is not present in the graph."))
    end
    for t in sink_nodes
        t in all_nodes || throw(ArgumentError("sink node $t is not present in the graph."))
    end

    return all_nodes
end

function _validate_capacities(
    edgelist::Vector{Tuple{Int64,Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64}
)::Nothing
    for e in edgelist
        haskey(capacities, e) || throw(ArgumentError("Missing capacity for edge $e."))
        c = capacities[e]
        (isnan(c) || c < 0.0) && throw(ArgumentError("Invalid capacity for edge $e: $c"))
    end
    nothing
end

function _super_sink_id(sorted_nodes::Vector{Int64})::Int64
    min_node = sorted_nodes[1]
    min_node == typemin(Int64) && throw(ArgumentError(
        "Cannot construct graph-level super sink via minimum(all_nodes)-1: minimum node ID is typemin(Int64)."
    ))
    return min_node - 1
end

function _copy_indices(
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}}
)::Tuple{Dict{Int64,Set{Int64}},Dict{Int64,Set{Int64}}}
    out_copy = Dict{Int64,Set{Int64}}()
    in_copy = Dict{Int64,Set{Int64}}()

    for (u, nbrs) in outgoing_index
        out_copy[u] = Set{Int64}(nbrs)
    end
    for (v, preds) in incoming_index
        in_copy[v] = Set{Int64}(preds)
    end

    return out_copy, in_copy
end

function _add_edge!(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing::Dict{Int64,Set{Int64}},
    incoming::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    u::Int64,
    v::Int64,
    cap::Float64
)::Nothing
    push!(edgelist, (u, v))
    if !haskey(outgoing, u); outgoing[u] = Set{Int64}(); end
    if !haskey(incoming, v); incoming[v] = Set{Int64}(); end
    if !haskey(outgoing, v); outgoing[v] = Set{Int64}(); end
    if !haskey(incoming, u); incoming[u] = Set{Int64}(); end
    push!(outgoing[u], v)
    push!(incoming[v], u)
    capacities[(u, v)] = cap
    nothing
end

function _build_indices_from_edges(
    edgelist::Vector{Tuple{Int64,Int64}}
)::Tuple{Dict{Int64,Set{Int64}},Dict{Int64,Set{Int64}}}
    nodes = union(Set(first.(edgelist)), Set(last.(edgelist)))
    outgoing = Dict{Int64,Set{Int64}}(n => Set{Int64}() for n in nodes)
    incoming = Dict{Int64,Set{Int64}}(n => Set{Int64}() for n in nodes)
    for (u, v) in edgelist
        push!(outgoing[u], v)
        push!(incoming[v], u)
    end
    return outgoing, incoming
end

function _remap_for_node_connectivity(
    edgelist::Vector{Tuple{Int64,Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source::Int64,
    super_sink::Int64
)::NamedTuple
    nodes = sort!(collect(union(Set(first.(edgelist)), Set(last.(edgelist)))))

    # Reserved mapped IDs are derived from this call's own node range (min(nodes)-1,
    # min(nodes)-2) rather than fixed constants (formerly -1/-3), so they can never
    # collide with a real node ID no matter how the caller's graph is numbered.
    # Fixed sentinels broke on a 0-indexed graph: `edge_connectivity`/`node_connectivity`
    # already inject their own synthetic super-sink at min(nodes)-1 one layer up
    # (`_super_sink_id`), which lands on exactly -1 when the caller's minimum node ID
    # is 0 -- colliding with the old hardcoded `source_mapped = -1` here.
    min_node = nodes[1]
    min_node == typemin(Int64) && throw(ArgumentError(
        "Cannot construct internal remap IDs via minimum(nodes)-1/-2: minimum node ID is typemin(Int64)."
    ))
    source_mapped = try
        Base.checked_sub(min_node, Int64(1))
    catch err
        err isa OverflowError || rethrow(err)
        throw(ArgumentError("Node ID range too large in magnitude for internal remapping (min node $min_node)."))
    end
    sink_mapped = try
        Base.checked_sub(min_node, Int64(2))
    catch err
        err isa OverflowError || rethrow(err)
        throw(ArgumentError("Node ID range too large in magnitude for internal remapping (min node $min_node)."))
    end

    orig_to_mapped = Dict{Int64,Int64}()
    mapped_to_orig = Dict{Int64,Int64}()

    orig_to_mapped[source] = source_mapped
    mapped_to_orig[source_mapped] = source
    orig_to_mapped[super_sink] = sink_mapped
    mapped_to_orig[sink_mapped] = super_sink

    next_id = Int64(1)
    for v in nodes
        if v == source || v == super_sink
            continue
        end
        orig_to_mapped[v] = next_id
        mapped_to_orig[next_id] = v
        next_id += 1
    end

    remap_edge(e::Tuple{Int64,Int64}) = (orig_to_mapped[e[1]], orig_to_mapped[e[2]])

    remapped_edgelist = Tuple{Int64,Int64}[remap_edge(e) for e in edgelist]
    remapped_caps = Dict{Tuple{Int64,Int64},Float64}()
    for e in edgelist
        remapped_caps[remap_edge(e)] = capacities[e]
    end
    remapped_out, remapped_in = _build_indices_from_edges(remapped_edgelist)

    return (
        edgelist=remapped_edgelist,
        outgoing=remapped_out,
        incoming=remapped_in,
        capacities=remapped_caps,
        source=source_mapped,
        sink=sink_mapped,
        mapped_to_orig=mapped_to_orig
    )
end

function _build_super_sink_graph(
    base_edgelist::Vector{Tuple{Int64,Int64}},
    base_outgoing::Dict{Int64,Set{Int64}},
    base_incoming::Dict{Int64,Set{Int64}},
    base_capacities::Dict{Tuple{Int64,Int64},Float64},
    sorted_nodes::Vector{Int64},
    source::Int64,
    super_sink::Int64,
    sink_edge_capacity::Float64
)::Tuple{Vector{Tuple{Int64,Int64}},Dict{Int64,Set{Int64}},Dict{Int64,Set{Int64}},Dict{Tuple{Int64,Int64},Float64}}
    edgelist = copy(base_edgelist)
    outgoing, incoming = _copy_indices(base_outgoing, base_incoming)
    capacities = Dict{Tuple{Int64,Int64},Float64}(base_capacities)

    if !haskey(outgoing, super_sink); outgoing[super_sink] = Set{Int64}(); end
    if !haskey(incoming, super_sink); incoming[super_sink] = Set{Int64}(); end

    for t in sorted_nodes
        if t == source
            continue
        end
        _add_edge!(edgelist, outgoing, incoming, capacities, t, super_sink, sink_edge_capacity)
    end

    return edgelist, outgoing, incoming, capacities
end

function _validate_integral_value(value::Float64, name::String, tol::Float64)::Int64
    rounded = round(Int64, value)
    abs(value - rounded) <= tol || throw(ArgumentError(
        "$name = $value is not integer within tol=$tol, violating integrality expectations for unit-capacity flow."
    ))
    return rounded
end

function _original_cut_edges(
    edgelist::Vector{Tuple{Int64,Int64}},
    cut_S::Set{Int64},
    cut_T::Set{Int64}
)::Vector{Tuple{Int64,Int64}}
    edges = Tuple{Int64,Int64}[]
    for (u, v) in edgelist
        if (u in cut_S) && (v in cut_T)
            push!(edges, (u, v))
        end
    end
    sort!(edges)
    return edges
end

function _pick_achieving_sink_from_super_sink_result(
    source::Int64,
    sorted_nodes::Vector{Int64},
    super_sink::Int64,
    result::FlowSolveResult,
    tol::Float64
)::Int64
    candidates = Int64[]
    for t in sorted_nodes
        if t == source
            continue
        end
        edge = (t, super_sink)
        flow_ts = get(result.flow, edge, 0.0)
        saturated = (1.0 - flow_ts) <= tol
        if saturated && (t in result.mincut_T)
            push!(candidates, t)
        end
    end

    if !isempty(candidates)
        sort!(candidates)
        return candidates[1]
    end

    t_side = sort!([t for t in sorted_nodes if t != source && (t in result.mincut_T)])
    !isempty(t_side) && return t_side[1]

    fallback = sort!([t for t in sorted_nodes if t != source])
    return fallback[1]
end

"""
    edge_connectivity(edgelist, outgoing_index, incoming_index, source_nodes, sink_nodes; algorithm=:dinic, tol=1e-10)

Compute exact directed edge connectivity using unit capacities and a graph-level
super-sink aggregation trick with one max-flow solve per source node.

Returns the minimum connectivity value `lambda`, an achieving `(source, sink)`
pair, the corresponding cut edges (restricted to original edges), and solver call count.
"""
function edge_connectivity(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10
)::EdgeConnectivityResult
    all_nodes = _require_valid_inputs(edgelist, source_nodes, sink_nodes)
    sorted_nodes = sort!(collect(all_nodes))
    super_sink = _super_sink_id(sorted_nodes)

    unit_caps = Dict{Tuple{Int64,Int64},Float64}(e => 1.0 for e in edgelist)

    solver_calls = Int64[0]
    best_lambda = typemax(Int64)
    best_source = sorted_nodes[1]
    best_sink = sorted_nodes[2]
    best_result = nothing

    for s in sorted_nodes
        aug_edgelist, aug_out, aug_in, aug_caps = _build_super_sink_graph(
            edgelist,
            outgoing_index,
            incoming_index,
            unit_caps,
            sorted_nodes,
            s,
            super_sink,
            1.0
        )

        result = _solve_with_algorithm(
            algorithm,
            aug_edgelist,
            aug_out,
            aug_in,
            aug_caps,
            Int64[s],
            Int64[super_sink];
            tol=tol,
            validate=true
        )
        solver_calls[1] += 1

        lambda_s = _validate_integral_value(result.max_flow, "edge connectivity flow", tol)
        if lambda_s < best_lambda
            best_lambda = lambda_s
            best_source = s
            best_sink = _pick_achieving_sink_from_super_sink_result(s, sorted_nodes, super_sink, result, tol)
            best_result = result
        end
    end

    cut_edges = _original_cut_edges(edgelist, best_result.mincut_S, best_result.mincut_T)

    return EdgeConnectivityResult(
        best_lambda,
        best_source,
        best_sink,
        cut_edges,
        solver_calls[1]
    )
end

"""
    node_connectivity(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes; algorithm=:dinic, tol=1e-10)

Compute exact directed node connectivity using node splitting with unit node capacities,
plus a graph-level super-sink aggregation trick with one solve per source node.

Returns `kappa`, an achieving `(source, sink)` pair, the minimum cut node set
(original node IDs only), and solver call count.
"""
function node_connectivity(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10
)::NodeConnectivityResult
    all_nodes = _require_valid_inputs(edgelist, source_nodes, sink_nodes)
    _validate_capacities(edgelist, capacities)
    sorted_nodes = sort!(collect(all_nodes))
    super_sink = _super_sink_id(sorted_nodes)

    base_caps = Dict{Tuple{Int64,Int64},Float64}(capacities)

    solver_calls = Int64[0]
    best_kappa = typemax(Int64)
    best_source = sorted_nodes[1]
    best_sink = sorted_nodes[2]
    best_cut_nodes = Int64[]

    for s in sorted_nodes
        aug_edgelist, aug_out, aug_in, aug_caps = _build_super_sink_graph(
            edgelist,
            outgoing_index,
            incoming_index,
            base_caps,
            sorted_nodes,
            s,
            super_sink,
            Inf
        )

        remapped = _remap_for_node_connectivity(aug_edgelist, aug_caps, s, super_sink)

        node_caps = Dict{Int64,Float64}()
        for v in keys(remapped.mapped_to_orig)
            if v != remapped.source && v != remapped.sink
                node_caps[v] = 1.0
            end
        end

        result = solve_node_capacitated_flow(
            remapped.edgelist,
            remapped.outgoing,
            remapped.incoming,
            remapped.capacities,
            Int64[remapped.source],
            Int64[remapped.sink],
            node_caps;
            algorithm=algorithm,
            tol=tol,
            validate=true
        )
        solver_calls[1] += 1

        kappa_s = _validate_integral_value(result.max_flow, "node connectivity flow", tol)
        if kappa_s < best_kappa
            best_kappa = kappa_s
            best_source = s

            t_side = sort!([
                remapped.mapped_to_orig[t] for t in result.mincut_T
                if haskey(remapped.mapped_to_orig, t) && remapped.mapped_to_orig[t] in all_nodes && remapped.mapped_to_orig[t] != s && remapped.mapped_to_orig[t] != super_sink
            ])
            best_sink = isempty(t_side) ? sort!([t for t in sorted_nodes if t != s])[1] : t_side[1]

            mapped_cut_nodes = [
                remapped.mapped_to_orig[v] for v in result.saturated_nodes
                if haskey(remapped.mapped_to_orig, v)
            ]
            best_cut_nodes = sort!(unique(mapped_cut_nodes))
        end
    end

    min_cut_nodes = sort!([v for v in best_cut_nodes if (v in all_nodes) && v != best_source && v != super_sink])

    return NodeConnectivityResult(
        best_kappa,
        best_source,
        best_sink,
        min_cut_nodes,
        solver_calls[1]
    )
end

"""
    global_min_cut(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes; algorithm=:dinic, tol=1e-10)

Compute the directed global minimum cut capacity with original capacities.

Runs two O(V) directional passes:
1) fixed source, varying sinks;
2) fixed sink, varying sources.

Returns min-cut capacity, achieving pair, cut edges (original edges only), cut sides,
and solver call count.
"""
function global_min_cut(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10
)::GlobalMinCutResult
    all_nodes = _require_valid_inputs(edgelist, source_nodes, sink_nodes)
    _validate_capacities(edgelist, capacities)
    sorted_nodes = sort!(collect(all_nodes))

    solver_calls = Int64[0]
    best_capacity = Inf
    best_source = sorted_nodes[1]
    best_sink = sorted_nodes[2]
    best_result = nothing

    fixed_source = sorted_nodes[1]
    for t in sorted_nodes
        if t == fixed_source
            continue
        end
        result = _solve_with_algorithm(
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            Int64[fixed_source],
            Int64[t];
            tol=tol,
            validate=true
        )
        solver_calls[1] += 1

        if result.max_flow + tol < best_capacity
            best_capacity = result.max_flow
            best_source = fixed_source
            best_sink = t
            best_result = result
        end
    end

    fixed_sink = sorted_nodes[1]
    for s in sorted_nodes
        if s == fixed_sink
            continue
        end
        result = _solve_with_algorithm(
            algorithm,
            edgelist,
            outgoing_index,
            incoming_index,
            capacities,
            Int64[s],
            Int64[fixed_sink];
            tol=tol,
            validate=true
        )
        solver_calls[1] += 1

        if result.max_flow + tol < best_capacity
            best_capacity = result.max_flow
            best_source = s
            best_sink = fixed_sink
            best_result = result
        end
    end

    cut_edges = _original_cut_edges(edgelist, best_result.mincut_S, best_result.mincut_T)

    return GlobalMinCutResult(
        best_capacity,
        best_source,
        best_sink,
        cut_edges,
        Set(best_result.mincut_S),
        Set(best_result.mincut_T),
        solver_calls[1]
    )
end

"""
    analyze_global_connectivity(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes; algorithm=:dinic, tol=1e-10)

Aggregate entry point that computes:
- edge connectivity,
- node connectivity,
- global minimum cut.

Returns a `GlobalConnectivityResult` containing all three sub-results.
"""
function analyze_global_connectivity(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10
)::GlobalConnectivityResult
    edge_result = edge_connectivity(
        edgelist,
        outgoing_index,
        incoming_index,
        source_nodes,
        sink_nodes;
        algorithm=algorithm,
        tol=tol
    )

    node_result = node_connectivity(
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes;
        algorithm=algorithm,
        tol=tol
    )

    global_cut_result = global_min_cut(
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        source_nodes,
        sink_nodes;
        algorithm=algorithm,
        tol=tol
    )

    return GlobalConnectivityResult(edge_result, node_result, global_cut_result)
end

end # module GlobalConnectivityModule
