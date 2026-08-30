module FlowDecompositionModule

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

export FlowPathComponent,
       FlowDecomposition,
       decompose_flow,
       validate_decomposition

"""
    FlowPathComponent

One path-flow component in a valid source-to-sink flow decomposition.
`bottleneck_edge` is the lexicographically smallest path edge attaining the
component bottleneck value.
"""
struct FlowPathComponent
    path::Vector{Int64}
    flow_value::Float64
    bottleneck_edge::Tuple{Int64,Int64}
end

"""
    FlowDecomposition

A validated decomposition of a solved flow into explicit path components.
The decomposition is not unique. This struct stores one valid decomposition out
of potentially many. All valid decompositions satisfy the same exactness checks.
"""
struct FlowDecomposition
    components::Vector{FlowPathComponent}
    total_flow::Float64
    is_unique::Bool
end

function _build_outgoing_index(edgelist::Vector{Tuple{Int64,Int64}})::Dict{Int64,Vector{Int64}}
    outgoing = Dict{Int64,Vector{Int64}}()
    for (u, v) in edgelist
        if !haskey(outgoing, u)
            outgoing[u] = Int64[]
        end
        push!(outgoing[u], v)
    end
    for neighbors in values(outgoing)
        sort!(neighbors)
    end
    return outgoing
end

function _source_outflow(
    source::Int64,
    outgoing_index::Dict{Int64,Vector{Int64}},
    working_flow::Dict{Tuple{Int64,Int64},Float64},
    tol::Float64
)::Float64
    total = 0.0
    for v in get(outgoing_index, source, Int64[])
        value = get(working_flow, (source, v), 0.0)
        if value > tol
            total += value
        end
    end
    return total
end

function _all_sources_exhausted(
    sorted_sources::Vector{Int64},
    outgoing_index::Dict{Int64,Vector{Int64}},
    working_flow::Dict{Tuple{Int64,Int64},Float64},
    tol::Float64
)::Bool
    for s in sorted_sources
        if _source_outflow(s, outgoing_index, working_flow, tol) > tol
            return false
        end
    end
    return true
end

function _select_source(
    sorted_sources::Vector{Int64},
    outgoing_index::Dict{Int64,Vector{Int64}},
    working_flow::Dict{Tuple{Int64,Int64},Float64},
    tol::Float64
)::Union{Nothing,Int64}
    for s in sorted_sources
        if _source_outflow(s, outgoing_index, working_flow, tol) > tol
            return s
        end
    end
    return nothing
end

function _extract_positive_flow_path(
    start_source::Int64,
    sink_set::Set{Int64},
    outgoing_index::Dict{Int64,Vector{Int64}},
    working_flow::Dict{Tuple{Int64,Int64},Float64},
    tol::Float64
)::Vector{Int64}
    path = Int64[start_source]
    current = start_source

    while !(current in sink_set)
        next_node = nothing
        for nxt in get(outgoing_index, current, Int64[])
            if get(working_flow, (current, nxt), 0.0) > tol
                next_node = nxt
                break
            end
        end
        next_node === nothing && throw(ArgumentError("Positive-flow path extraction failed at node $current; flow conservation/path-to-sink invariant violated."))
        push!(path, next_node)
        current = next_node
    end

    return path
end

function _path_edges(path::Vector{Int64})::Vector{Tuple{Int64,Int64}}
    length(path) >= 2 || return Tuple{Int64,Int64}[]
    return Tuple{Int64,Int64}[(path[i], path[i + 1]) for i in 1:(length(path) - 1)]
end

function _component_from_path(
    path::Vector{Int64},
    working_flow::Dict{Tuple{Int64,Int64},Float64},
    tol::Float64
)::FlowPathComponent
    edges = _path_edges(path)
    isempty(edges) && throw(ArgumentError("A decomposition path must contain at least one edge."))

    bottleneck = Inf
    bottleneck_edges = Tuple{Int64,Int64}[]
    for edge in edges
        value = get(working_flow, edge, 0.0)
        if value <= tol
            throw(ArgumentError("Decomposition path contains exhausted edge $edge with working flow $value."))
        end
        if value < bottleneck - tol
            bottleneck = value
            empty!(bottleneck_edges)
            push!(bottleneck_edges, edge)
        elseif abs(value - bottleneck) <= tol
            push!(bottleneck_edges, edge)
        end
    end

    sort!(bottleneck_edges)
    return FlowPathComponent(copy(path), bottleneck, first(bottleneck_edges))
end

function _subtract_component!(
    component::FlowPathComponent,
    working_flow::Dict{Tuple{Int64,Int64},Float64},
    tol::Float64
)::Nothing
    for edge in _path_edges(component.path)
        updated = get(working_flow, edge, 0.0) - component.flow_value
        if updated <= tol
            working_flow[edge] = 0.0
        else
            working_flow[edge] = updated
        end
    end
    nothing
end

function _all_integerish(values::Vector{Float64}, tol::Float64)::Bool
    for value in values
        if abs(value - round(value)) > tol
            return false
        end
    end
    return true
end

"""
    validate_decomposition(decomposition, flow_result, edgelist; tol=1e-10)

Validate a flow decomposition exactly against the solved flow. Two checks are enforced:
(1) total path flow equals the solved max flow, and (2) per-edge path accounting equals
`flow_result.flow` on every original edge to within `tol`. Throws on violation.
"""
function validate_decomposition(
    decomposition::FlowDecomposition,
    flow_result::FlowSolveResult,
    edgelist::Vector{Tuple{Int64,Int64}};
    tol::Float64=1e-10
)::Nothing
    _require_bounded_baseline(flow_result)

    total_component_flow = sum((component.flow_value for component in decomposition.components); init=0.0)
    abs(total_component_flow - flow_result.max_flow) <= tol ||
        throw(ArgumentError("Flow decomposition total mismatch: decomposition total=$total_component_flow, solved max flow=$(flow_result.max_flow)."))
    # Also verify the stored total_flow field independently —
    # catches any mismatch between components and the stored field.
    abs(decomposition.total_flow - flow_result.max_flow) <= tol ||
        throw(ArgumentError("FlowDecomposition.total_flow=$(decomposition.total_flow) does not match solved max flow=$(flow_result.max_flow)."))

    accounted = Dict{Tuple{Int64,Int64},Float64}(edge => 0.0 for edge in edgelist)
    edge_set = Set(edgelist)

    for component in decomposition.components
        component.flow_value > tol || throw(ArgumentError("Decomposition contains nonpositive component flow $(component.flow_value)."))
        path_edges = _path_edges(component.path)
        isempty(path_edges) && throw(ArgumentError("Decomposition component path must contain at least one edge."))
        component.bottleneck_edge in path_edges || throw(ArgumentError("Component bottleneck edge $(component.bottleneck_edge) is not on its path."))

        for edge in path_edges
            edge in edge_set || throw(ArgumentError("Decomposition path uses edge $edge not present in edgelist."))
            accounted[edge] = get(accounted, edge, 0.0) + component.flow_value
        end
    end

    for edge in edgelist
        expected = get(flow_result.flow, edge, 0.0)
        actual = get(accounted, edge, 0.0)
        abs(actual - expected) <= tol ||
            throw(ArgumentError("Flow decomposition edge accounting mismatch on edge $edge: accounted=$actual, solved=$expected."))
    end

    nothing
end

"""
    decompose_flow(edgelist, source_nodes, sink_nodes, flow_result; tol=1e-10)

Decompose the solved flow into explicit source-to-sink path components carrying positive
flow values. This is distinct from structural path enumeration: `StructuralModule`
functions enumerate all topological source-to-sink paths, while this function returns only
paths that actually carry flow, and the returned flow values add up exactly to the solved
max flow.

The decomposition is not unique. This function returns one deterministic canonical
valid decomposition obtained by repeatedly extracting lexicographically earliest positive-
flow source-to-sink paths from a mutable copy of `flow_result.flow`. Validation always
runs before the result is returned.
"""
function decompose_flow(
    edgelist::Vector{Tuple{Int64,Int64}},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    tol::Float64=1e-10
)::FlowDecomposition
    _require_bounded_baseline(flow_result)

    sink_set = Set(sink_nodes)
    outgoing_index = _build_outgoing_index(edgelist)
    # Precomputed once to avoid repeated sort+unique in loop
    sorted_sources = sort!(unique(copy(source_nodes)))
    working_flow = Dict{Tuple{Int64,Int64},Float64}(edge => max(0.0, get(flow_result.flow, edge, 0.0)) for edge in edgelist)
    components = FlowPathComponent[]

    while !_all_sources_exhausted(sorted_sources, outgoing_index, working_flow, tol)
        source = _select_source(sorted_sources, outgoing_index, working_flow, tol)
        source === nothing && break

        path = _extract_positive_flow_path(source, sink_set, outgoing_index, working_flow, tol)
        component = _component_from_path(path, working_flow, tol)
        push!(components, component)
        _subtract_component!(component, working_flow, tol)
    end

    sort!(components; by=c -> (-c.flow_value, Tuple(c.path)))
    decomposition = FlowDecomposition(components, sum((c.flow_value for c in components); init=0.0), false)

    positive_edge_flows = [get(flow_result.flow, edge, 0.0) for edge in edgelist if get(flow_result.flow, edge, 0.0) > tol]
    if !isempty(positive_edge_flows) && _all_integerish(positive_edge_flows, tol)
        for component in decomposition.components
            abs(component.flow_value - round(component.flow_value)) <= tol ||
                throw(ArgumentError("Integer-flow decomposition check failed: component flow $(component.flow_value) is not integral to tolerance $tol."))
        end
    end

    validate_decomposition(decomposition, flow_result, edgelist; tol=tol)
    return decomposition
end

end
