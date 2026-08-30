module NodeCapacitatedFlowModule

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

export NodeSplitGraph,
       NodeCapacitatedFlowResult,
       NodeCapacitatedAnalysisResult,
       build_node_split_graph,
       solve_node_capacitated_flow,
       node_capacitated_spof_nodes,
       analyze_node_capacitated_flow

# ─────────────────────────────────────────────────────────────────────────────
# Result structs
# ─────────────────────────────────────────────────────────────────────────────

"""
    NodeSplitGraph

Transformed graph produced by the node-splitting operation. Contains the split
graph ready to pass directly to any FlowModule solver, plus the bidirectional
mapping data needed to translate results back to original node IDs.

Node ID convention: for every split node v, v_in = 2*v and v_out = 2*v+1.
Nodes absent from node_capacities are NOT split — they appear with their
original ID unchanged in the split graph.
"""
struct NodeSplitGraph
    # Transformed graph — ready to pass directly to FlowModule
    split_edgelist    :: Vector{Tuple{Int64,Int64}}
    split_outgoing    :: Dict{Int64,Set{Int64}}
    split_incoming    :: Dict{Int64,Set{Int64}}
    split_capacities  :: Dict{Tuple{Int64,Int64},Float64}
    split_sources     :: Vector{Int64}
    split_sinks       :: Vector{Int64}

    # Mapping data — original ↔ split
    node_to_in        :: Dict{Int64,Int64}
    # original node v → split v_in ID (only for split nodes)
    node_to_out       :: Dict{Int64,Int64}
    # original node v → split v_out ID (only for split nodes)
    split_to_original :: Dict{Int64,Int64}
    # split node ID → original node ID (both v_in and v_out map to v)

    original_edgelist :: Vector{Tuple{Int64,Int64}}
    # kept for result mapping
    split_edge_to_original :: Dict{Tuple{Int64,Int64},Tuple{Int64,Int64}}
    # maps (u_out, v_in) in split graph back to (u, v) in original
    node_internal_edges :: Dict{Int64,Tuple{Int64,Int64}}
    # original node v → internal split edge (v_in, v_out); only for split nodes
end

"""
    NodeCapacitatedFlowResult

Result of solving a node-capacitated max-flow problem. All fields use original
node IDs and original edge tuples — no split node IDs appear in the public fields.
"""
struct NodeCapacitatedFlowResult
    # All fields use original node IDs and original edge tuples
    max_flow          :: Float64
    flow              :: Dict{Tuple{Int64,Int64},Float64}
    node_flow         :: Dict{Int64,Float64}
    # flow through each original node (= flow on internal edge for split nodes)
    sources           :: Vector{Int64}
    sinks             :: Vector{Int64}
    sink_flow         :: Dict{Int64,Float64}
    saturated_edges   :: Vector{Tuple{Int64,Int64}}
    # original edges whose capacity is fully used
    saturated_nodes   :: Vector{Int64}
    # nodes whose internal capacity edge is saturated
    mincut_S          :: Set{Int64}
    # original node IDs on source side of min-cut
    mincut_T          :: Set{Int64}
    # original node IDs on sink side of min-cut
    mincut_capacity   :: Float64
    is_unbounded      :: Bool
    node_split_graph  :: NodeSplitGraph
    # stored for downstream use (SPOF analysis, etc.)
end

"""
    NodeCapacitatedAnalysisResult

Aggregate result struct returned by `analyze_node_capacitated_flow`.
"""
struct NodeCapacitatedAnalysisResult
    flow_result  :: NodeCapacitatedFlowResult
    spof_nodes   :: Vector{Int64}
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

_graph_nodes_ncf(edgelist::Vector{Tuple{Int64,Int64}})::Set{Int64} = _graph_nodes_set(edgelist)

"""
Validate that node_capacities values are finite, nonneg, non-NaN, and that
all keys appear in the graph. Throws ArgumentError on violation.
"""
function _validate_node_capacities(
    node_capacities::Dict{Int64,Float64},
    all_nodes::Set{Int64}
)::Nothing
    for (v, c) in node_capacities
        if !(v in all_nodes)
            throw(ArgumentError(
                "node_capacities key $v does not appear in the graph. " *
                "All node_capacities keys must be nodes present in edgelist."
            ))
        end
        if isnan(c)
            throw(ArgumentError(
                "node_capacities[$v] is NaN. Node capacities must be finite and nonnegative."
            ))
        end
        if !isfinite(c)
            throw(ArgumentError(
                "node_capacities[$v] = $c is infinite. Node capacities must be finite. " *
                "Nodes without a capacity entry are automatically treated as infinite."
            ))
        end
        if c < 0.0
            throw(ArgumentError(
                "node_capacities[$v] = $c is negative. Node capacities must be nonnegative."
            ))
        end
    end
    nothing
end

"""
Assign overflow-safe, collision-free split IDs `(v_in, v_out)` to every node
in `split_nodes`, given the full node set `all_nodes` the split graph has to
coexist with.

Offset-based, not the old `v_in=2*v`/`v_out=2*v+1` convention: every synthetic
ID is placed strictly above `M = maximum(all_nodes)`, so it can never collide
with an UNSPLIT node's own (unchanged) ID regardless of how small or dense
the original ID range is — the old multiplicative scheme could and did
collide whenever only a SUBSET of a small, sequential ID range (1..N) was
split, which is the normal case for an optional per-node capacity list (see
the Network Model chapter: capacity is attached only where a node genuinely
has one; every other node is unconstrained, not omitted by mistake). Each
split node gets its own consecutive pair in a canonical (sorted) order, so
split-vs-split collisions are impossible too. The synthetic IDs are never
part of this function's public contract — callers translate through
`NodeSplitGraph`'s own `node_to_in`/`node_to_out`/`split_to_original` maps,
never by re-deriving the formula — so this is a behavior-preserving change
for every case that already worked, and a fix for the case that didn't.

Throws ArgumentError if the offset arithmetic would overflow Int64 (the same
failure mode the old scheme guarded against, now only reachable at genuinely
astronomical node-ID magnitudes, not at ordinary small-ID density).
"""
function _split_id_map(all_nodes::Set{Int64}, split_nodes::Set{Int64})::Dict{Int64,Tuple{Int64,Int64}}
    isempty(split_nodes) && return Dict{Int64,Tuple{Int64,Int64}}()
    m = maximum(all_nodes)
    mapping = Dict{Int64,Tuple{Int64,Int64}}()
    for (i, v) in enumerate(sort!(collect(split_nodes)))
        v_in = try
            Base.checked_add(m, Base.checked_sub(Base.checked_mul(i, 2), 1))
        catch err
            err isa OverflowError || rethrow(err)
            throw(ArgumentError(
                "Node ID range is too large in magnitude for the offset split scheme " *
                "(max original ID $m, $(length(split_nodes)) split nodes); the synthetic " *
                "ID for node $v overflows Int64."
            ))
        end
        v_out = try
            Base.checked_add(v_in, 1)
        catch err
            err isa OverflowError || rethrow(err)
            throw(ArgumentError(
                "Node ID range is too large in magnitude for the offset split scheme " *
                "(max original ID $m, $(length(split_nodes)) split nodes); the synthetic " *
                "ID for node $v overflows Int64."
            ))
        end
        mapping[v] = (v_in, v_out)
    end
    mapping
end

"""
Defense-in-depth sanity check that the proposed split-node ID space (from
`_split_id_map`) does not collide with any unsplit node's own ID, or with
itself. Mathematically unreachable given `_split_id_map`'s offset
construction (every synthetic ID exceeds every original ID by construction),
kept as an explicit assertion rather than removed silently — a correctness
check earning its keep by staying cheap, not by being load-bearing.

Throws ArgumentError with a clear description if a collision is detected.
"""
function _detect_collisions(
    all_nodes::Set{Int64},
    split_nodes::Set{Int64},
    split_id_map::Dict{Int64,Tuple{Int64,Int64}}=_split_id_map(all_nodes, split_nodes)
)::Nothing
    proposed = Dict{Int64,String}()  # ID → description of who claims it

    for v in all_nodes
        if v in split_nodes
            v_in, v_out = split_id_map[v]
            for (id, desc) in ((v_in, "v_in(node $v)"), (v_out, "v_out(node $v)"))
                if haskey(proposed, id)
                    throw(ArgumentError(
                        "Node ID collision in split graph: split ID $id is claimed by " *
                        "both $desc and $(proposed[id])."
                    ))
                end
                proposed[id] = desc
            end
        else
            # Unsplit node: uses its original ID directly
            if haskey(proposed, v)
                throw(ArgumentError(
                    "Node ID collision in split graph: original (unsplit) node $v has ID $v, " *
                    "which collides with split ID claimed by $(proposed[v])."
                ))
            end
            proposed[v] = "unsplit(node $v)"
        end
    end
    nothing
end

"""
Forward-reachability BFS from `starts`, optionally banning one node. Returns
the set of nodes reachable (starts included, unless banned).
"""
function _forward_reachable_ncf(
    starts::Vector{Int64},
    outgoing_index::Dict{Int64,Set{Int64}};
    banned_node::Union{Nothing,Int64}=nothing
)::Set{Int64}
    reachable = Set{Int64}()
    queue = Int64[]
    for s in starts
        (banned_node !== nothing && s == banned_node) && continue
        if !(s in reachable)
            push!(reachable, s)
            push!(queue, s)
        end
    end
    head = 1
    while head <= length(queue)
        u = queue[head]; head += 1
        for v in get(outgoing_index, u, Set{Int64}())
            (banned_node !== nothing && v == banned_node) && continue
            if !(v in reachable)
                push!(reachable, v)
                push!(queue, v)
            end
        end
    end
    return reachable
end

"""
Returns true iff any sink is reachable from any source when `banned_node` is removed
from the graph.
"""
function _any_sink_reachable_ncf(
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    outgoing_index::Dict{Int64,Set{Int64}};
    banned_node::Union{Nothing,Int64}=nothing
)::Bool
    reachable = _forward_reachable_ncf(source_nodes, outgoing_index; banned_node=banned_node)
    sink_set = Set(sink_nodes)
    return !isempty(intersect(reachable, sink_set))
end

"""
Backward-reachability BFS from `targets`, optionally banning one node.
"""
function _backward_reachable_ncf(
    targets::Vector{Int64},
    incoming_index::Dict{Int64,Set{Int64}};
    banned_node::Union{Nothing,Int64}=nothing
)::Set{Int64}
    reachable = Set{Int64}()
    queue = Int64[]
    for t in targets
        (banned_node !== nothing && t == banned_node) && continue
        if !(t in reachable)
            push!(reachable, t)
            push!(queue, t)
        end
    end
    head = 1
    while head <= length(queue)
        v = queue[head]; head += 1
        for u in get(incoming_index, v, Set{Int64}())
            (banned_node !== nothing && u == banned_node) && continue
            if !(u in reachable)
                push!(reachable, u)
                push!(queue, u)
            end
        end
    end
    return reachable
end

"""
Add a directed edge (u → v) with the given capacity to the mutable split graph
being constructed.
"""
function _add_split_edge!(
    split_edgelist::Vector{Tuple{Int64,Int64}},
    split_outgoing::Dict{Int64,Set{Int64}},
    split_incoming::Dict{Int64,Set{Int64}},
    split_capacities::Dict{Tuple{Int64,Int64},Float64},
    u::Int64, v::Int64, cap::Float64
)::Nothing
    push!(split_edgelist, (u, v))
    if !haskey(split_outgoing, u); split_outgoing[u] = Set{Int64}(); end
    if !haskey(split_incoming, v); split_incoming[v] = Set{Int64}(); end
    push!(split_outgoing[u], v)
    push!(split_incoming[v], u)
    split_capacities[(u, v)] = cap
    nothing
end

"""
Ensure both endpoints of every edge have entries in outgoing/incoming dicts
(possibly empty), so that downstream code can always call get(index, v, ...).
This mirrors the contract expected by _solve_with_algorithm.
"""
function _ensure_index_entries!(
    split_outgoing::Dict{Int64,Set{Int64}},
    split_incoming::Dict{Int64,Set{Int64}},
    split_edgelist::Vector{Tuple{Int64,Int64}}
)::Nothing
    for (u, v) in split_edgelist
        if !haskey(split_outgoing, u); split_outgoing[u] = Set{Int64}(); end
        if !haskey(split_incoming, v); split_incoming[v] = Set{Int64}(); end
        if !haskey(split_outgoing, v); split_outgoing[v] = Set{Int64}(); end
        if !haskey(split_incoming, u); split_incoming[u] = Set{Int64}(); end
    end
    nothing
end

"""
Map a split node's flow result back to original node IDs for node_flow.
For split nodes: node_flow[v] = flow on internal edge (v_in, v_out).
For unsplit nodes: node_flow[v] = total inflow from all split-graph predecessors.
"""
function _map_node_flow(
    all_nodes::Set{Int64},
    split_nodes::Set{Int64},
    node_to_in::Dict{Int64,Int64},
    node_to_out::Dict{Int64,Int64},
    split_incoming::Dict{Int64,Set{Int64}},
    split_result_flow::Dict{Tuple{Int64,Int64},Float64},
    node_internal_edges::Dict{Int64,Tuple{Int64,Int64}}
)::Dict{Int64,Float64}
    nf = Dict{Int64,Float64}()
    for v in all_nodes
        if v in split_nodes
            internal_edge = node_internal_edges[v]
            nf[v] = get(split_result_flow, internal_edge, 0.0)
        else
            # Unsplit: inflow = sum over all predecessors in split graph
            inflow = 0.0
            for u in get(split_incoming, v, Set{Int64}())
                inflow += get(split_result_flow, (u, v), 0.0)
            end
            nf[v] = inflow
        end
    end
    return nf
end

"""
Map split min-cut sets back to original node IDs.
Convention: a node v is in original mincut_S iff its output representative
(v_out for split nodes, v itself for unsplit) is in split_result.mincut_S.
"""
function _map_mincut(
    all_nodes::Set{Int64},
    split_nodes::Set{Int64},
    node_to_out::Dict{Int64,Int64},
    split_mincut_S::Set{Int64}
)::Tuple{Set{Int64},Set{Int64}}
    orig_S = Set{Int64}()
    orig_T = Set{Int64}()
    for v in all_nodes
        representative = v in split_nodes ? node_to_out[v] : v
        if representative in split_mincut_S
            push!(orig_S, v)
        else
            push!(orig_T, v)
        end
    end
    return orig_S, orig_T
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    build_node_split_graph(
        edgelist, outgoing_index, incoming_index,
        capacities, source_nodes, sink_nodes,
        node_capacities;
        tol=1e-10
    ) -> NodeSplitGraph

Build the node-split graph representation for a node-capacitated flow problem.

Each node v with a finite capacity c(v) is replaced by two nodes:
  v_in  = 2 * v  (receives all incoming edges)
  v_out = 2 * v + 1  (sends all outgoing edges)
  plus an internal edge (v_in, v_out) with capacity c(v)

Nodes absent from `node_capacities` are treated as having infinite capacity
and are NOT split — they appear with their original ID in the transformed graph.

Does not run any solver. Returns a NodeSplitGraph containing both the
transformed graph inputs (ready to pass to any FlowModule solver) and the
mapping data needed to translate results back to original node IDs.

# Arguments
- `node_capacities`: Dict mapping original node ID to its finite capacity.
  All values must be finite, nonnegative, and non-NaN.
  Throws `ArgumentError` for invalid values or unknown node IDs.
"""
function build_node_split_graph(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    node_capacities::Dict{Int64,Float64};
    tol::Float64=1e-10
)::NodeSplitGraph
    all_nodes = _graph_nodes_ncf(edgelist)
    split_nodes = Set{Int64}(keys(node_capacities))

    # Validate inputs
    _validate_node_capacities(node_capacities, all_nodes)
    split_id_map = _split_id_map(all_nodes, split_nodes)
    _detect_collisions(all_nodes, split_nodes, split_id_map)

    # Build ID mappings
    node_to_in  = Dict{Int64,Int64}()
    node_to_out = Dict{Int64,Int64}()
    split_to_original = Dict{Int64,Int64}()
    node_internal_edges = Dict{Int64,Tuple{Int64,Int64}}()

    for v in split_nodes
        v_in, v_out = split_id_map[v]
        node_to_in[v]  = v_in
        node_to_out[v] = v_out
        split_to_original[v_in]  = v
        split_to_original[v_out] = v
        node_internal_edges[v] = (v_in, v_out)
    end

    # Helper: get the "output representative" of a node in the split graph.
    # For edge (u, v): u uses its v_out representative, v uses its v_in representative.
    _out_rep(v::Int64) = v in split_nodes ? node_to_out[v] : v
    _in_rep(v::Int64)  = v in split_nodes ? node_to_in[v]  : v

    # Build transformed graph
    split_edgelist   = Vector{Tuple{Int64,Int64}}()
    split_outgoing   = Dict{Int64,Set{Int64}}()
    split_incoming   = Dict{Int64,Set{Int64}}()
    split_capacities = Dict{Tuple{Int64,Int64},Float64}()
    split_edge_to_original = Dict{Tuple{Int64,Int64},Tuple{Int64,Int64}}()

    # Step 1: Add internal capacity edges for all split nodes
    for v in sort!(collect(split_nodes))  # deterministic order
        v_in  = node_to_in[v]
        v_out = node_to_out[v]
        cap_v = node_capacities[v]
        _add_split_edge!(split_edgelist, split_outgoing, split_incoming, split_capacities,
                         v_in, v_out, cap_v)
        # Internal edges do not correspond to original edges (no entry in split_edge_to_original)
    end

    # Step 2: Add remapped original edges
    for (u, v) in edgelist
        u_mapped = _out_rep(u)
        v_mapped = _in_rep(v)
        cap_uv   = capacities[(u, v)]
        _add_split_edge!(split_edgelist, split_outgoing, split_incoming, split_capacities,
                         u_mapped, v_mapped, cap_uv)
        split_edge_to_original[(u_mapped, v_mapped)] = (u, v)
    end

    # Ensure all nodes have index entries (including isolated endpoints)
    _ensure_index_entries!(split_outgoing, split_incoming, split_edgelist)

    # Step 3: Remap source and sink vectors
    # Sources: flow enters at v_in (start of split), so use v_in = _in_rep
    split_sources = Int64[_in_rep(s) for s in source_nodes]
    # Sinks: flow exits at v_out (end of split), so use v_out = _out_rep
    split_sinks   = Int64[_out_rep(t) for t in sink_nodes]

    return NodeSplitGraph(
        split_edgelist,
        split_outgoing,
        split_incoming,
        split_capacities,
        split_sources,
        split_sinks,
        node_to_in,
        node_to_out,
        split_to_original,
        edgelist,
        split_edge_to_original,
        node_internal_edges
    )
end

"""
Map a solved FlowSolveResult on the split graph back to original-node-ID terms,
producing a NodeCapacitatedFlowResult.
"""
function _map_split_result(
    split_result::FlowSolveResult,
    nsg::NodeSplitGraph,
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64};
    tol::Float64=1e-10
)::NodeCapacitatedFlowResult
    all_nodes   = _graph_nodes_ncf(nsg.original_edgelist)
    split_nodes = Set{Int64}(keys(nsg.node_to_in))

    # --- flow dict: original edges ---
    # Each original edge (u, v) maps to (u_out, v_in) in the split graph.
    # u uses its v_out representative; v uses its v_in representative.
    flow_orig = Dict{Tuple{Int64,Int64},Float64}()
    for (u, v) in nsg.original_edgelist
        u_split = u in split_nodes ? nsg.node_to_out[u] : u
        v_split = v in split_nodes ? nsg.node_to_in[v]  : v
        flow_orig[(u, v)] = get(split_result.flow, (u_split, v_split), 0.0)
    end

    # --- node_flow ---
    node_flow_orig = _map_node_flow(
        all_nodes, split_nodes,
        nsg.node_to_in, nsg.node_to_out,
        nsg.split_incoming, split_result.flow,
        nsg.node_internal_edges
    )

    # --- sink_flow: original sink → flow ---
    sink_flow_orig = Dict{Int64,Float64}()
    for (i, t_orig) in enumerate(sink_nodes)
        t_split = nsg.split_sinks[i]
        sink_flow_orig[t_orig] = get(split_result.sink_flow, t_split, 0.0)
    end

    # --- saturated_edges: original edge tuples ---
    saturated_orig = Tuple{Int64,Int64}[]
    for se in split_result.saturated_edges
        if haskey(nsg.split_edge_to_original, se)
            push!(saturated_orig, nsg.split_edge_to_original[se])
        end
        # Internal node-capacity edges (v_in, v_out) intentionally have no original
        # edge counterpart, so they are excluded from saturated_edges and reported via
        # saturated_nodes below.
    end
    sort!(saturated_orig)
    unique!(saturated_orig)

    # --- saturated_nodes: nodes whose internal capacity edge is saturated ---
    saturated_nodes_set = Set{Int64}()
    for v in split_nodes
        int_edge = nsg.node_internal_edges[v]
        cap      = nsg.split_capacities[int_edge]
        flow_v   = get(split_result.flow, int_edge, 0.0)
        if cap - flow_v <= tol
            push!(saturated_nodes_set, v)
        end
    end
    saturated_nodes_orig = sort!(collect(saturated_nodes_set))

    # --- mincut_S and mincut_T mapped back to original node IDs ---
    orig_S, orig_T = _map_mincut(all_nodes, split_nodes, nsg.node_to_out, split_result.mincut_S)

    return NodeCapacitatedFlowResult(
        split_result.max_flow,
        flow_orig,
        node_flow_orig,
        copy(source_nodes),
        copy(sink_nodes),
        sink_flow_orig,
        saturated_orig,
        saturated_nodes_orig,
        orig_S,
        orig_T,
        split_result.mincut_capacity,
        split_result.is_unbounded,
        nsg
    )
end

"""
    solve_node_capacitated_flow(
        edgelist, outgoing_index, incoming_index,
        capacities, source_nodes, sink_nodes,
        node_capacities;
        algorithm=:dinic,
        tol=1e-10,
        validate=true
    ) -> NodeCapacitatedFlowResult

Solve a max-flow problem where nodes have capacity constraints in addition to
edge capacities.

The node-splitting transformation (exact bijection, no approximation) replaces
each capacity-constrained node v with two nodes v_in and v_out connected by an
internal edge of capacity c(v). Max flow is preserved exactly.

Nodes absent from `node_capacities` are treated as having infinite capacity.
All result fields use original node IDs and original edge tuples.
"""
function solve_node_capacitated_flow(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    node_capacities::Dict{Int64,Float64};
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    validate::Bool=true
)::NodeCapacitatedFlowResult
    nsg = build_node_split_graph(
        edgelist, outgoing_index, incoming_index,
        capacities, source_nodes, sink_nodes,
        node_capacities; tol=tol
    )

    split_result = _solve_with_algorithm(
        algorithm,
        nsg.split_edgelist,
        nsg.split_outgoing,
        nsg.split_incoming,
        nsg.split_capacities,
        nsg.split_sources,
        nsg.split_sinks;
        tol=tol,
        validate=validate
    )

    return _map_split_result(split_result, nsg, source_nodes, sink_nodes; tol=tol)
end

"""
    node_capacitated_spof_nodes(
        edgelist, outgoing_index, incoming_index,
        capacities, source_nodes, sink_nodes,
        node_capacities;
        algorithm=:dinic,
        tol=1e-10,
        baseline_result=nothing
    ) -> Vector{Int64}

Identify single points of COMPLETE failure (SPOFs) in the node-capacitated network.

A node v is a SPOF iff completely blocking its capacity (setting node_capacities[v] = 0)
reduces the maximum flow to zero. This is the "single point of complete failure"
semantics: v is a SPOF iff every source-to-sink path passes through v. Nodes that
reduce flow but do not stop all flow (e.g., one parallel path among several) are NOT
classified as SPOFs.

For nodes with infinite capacity (absent from `node_capacities`), structural
reachability is used: v is a structural SPOF iff removing v disconnects all
source-to-sink paths.

Returns a sorted `Vector{Int64}` of SPOF node IDs.

Keyword `baseline_result` is optional and intended for internal reuse: when
provided, it avoids re-solving the baseline node-capacitated flow.
"""
function node_capacitated_spof_nodes(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    node_capacities::Dict{Int64,Float64};
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10,
    baseline_result::Union{Nothing,NodeCapacitatedFlowResult}=nothing
)::Vector{Int64}
    source_set = Set(source_nodes)
    sink_set   = Set(sink_nodes)
    all_nodes  = _graph_nodes_ncf(edgelist)

    # Solve baseline once unless a precomputed baseline result is supplied.
    baseline = baseline_result === nothing ?
        solve_node_capacitated_flow(
            edgelist, outgoing_index, incoming_index,
            capacities, source_nodes, sink_nodes,
            node_capacities; algorithm=algorithm, tol=tol, validate=true
        ) : baseline_result
    baseline_flow = baseline.max_flow
    # If baseline is already zero, no blocking can reduce it further
    baseline_flow < tol && return Int64[]

    spof_set = Set{Int64}()

    # --- Capacity-based SPOFs: nodes with finite node capacity > 0, not source/sink ---
    # A capacity SPOF is a node whose complete blocking (cap=0) reduces max flow to zero.
    # This is the "single point of COMPLETE failure" semantics, consistent with the
    # structural SPOF definition used in StructuralModule (every s-t path passes through v).
    cap_candidates = sort!([v for (v, c) in node_capacities
                            if !(v in source_set) && !(v in sink_set) && c > tol])

    for v in cap_candidates
        perturbed_caps = copy(node_capacities)
        perturbed_caps[v] = 0.0
        result_v = solve_node_capacitated_flow(
            edgelist, outgoing_index, incoming_index,
            capacities, source_nodes, sink_nodes,
            perturbed_caps; algorithm=algorithm, tol=tol, validate=true
        )
        # SPOF iff blocking v sends max flow to (approximately) zero
        if result_v.max_flow < tol
            push!(spof_set, v)
        end
    end

    # --- Structural SPOFs: nodes with infinite capacity (not in node_capacities),
    #     not sources or sinks, that structurally disconnect all s-t paths when removed ---
    finite_cap_set = Set{Int64}(keys(node_capacities))
    inf_nodes      = sort!([v for v in all_nodes
                            if !(v in finite_cap_set) && !(v in source_set) && !(v in sink_set)])

    # Only nodes on some s-t path can be structural SPOFs
    from_sources = _forward_reachable_ncf(source_nodes, outgoing_index)
    to_sinks     = _backward_reachable_ncf(sink_nodes, incoming_index)
    on_some_path = intersect(from_sources, to_sinks)

    for v in inf_nodes
        (v in on_some_path) || continue
        if !_any_sink_reachable_ncf(source_nodes, sink_nodes, outgoing_index; banned_node=v)
            push!(spof_set, v)
        end
    end

    return sort!(collect(spof_set))
end

"""
    analyze_node_capacitated_flow(
        edgelist, outgoing_index, incoming_index,
        capacities, source_nodes, sink_nodes,
        node_capacities;
        algorithm=:dinic,
        tol=1e-10
    ) -> NodeCapacitatedAnalysisResult

Aggregate entry point for node-capacitated flow analysis. Runs max-flow
and SPOF detection in a single call.

Returns a `NodeCapacitatedAnalysisResult` with:
  - `flow_result`: full `NodeCapacitatedFlowResult`
  - `spof_nodes`: sorted vector of SPOF node IDs
"""
function analyze_node_capacitated_flow(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    node_capacities::Dict{Int64,Float64};
    algorithm::Symbol=:dinic,
    tol::Float64=1e-10
)::NodeCapacitatedAnalysisResult
    flow_result = solve_node_capacitated_flow(
        edgelist, outgoing_index, incoming_index,
        capacities, source_nodes, sink_nodes,
        node_capacities; algorithm=algorithm, tol=tol, validate=true
    )

    spof_nodes = node_capacitated_spof_nodes(
        edgelist, outgoing_index, incoming_index,
        capacities, source_nodes, sink_nodes,
        node_capacities; algorithm=algorithm, tol=tol, baseline_result=flow_result
    )

    return NodeCapacitatedAnalysisResult(flow_result, spof_nodes)
end

end  # module NodeCapacitatedFlowModule
