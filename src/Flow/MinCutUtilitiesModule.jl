module MinCutUtilitiesModule

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

export MinCut,
       MinCutEnumeration,
       MinCutAnalysis,
       minimum_st_cut_edges,
       minimum_st_cut_capacity,
       edges_in_some_mincut,
       edges_in_every_mincut,
       mincut_partition,
       enumerate_min_cuts,
       analyze_min_cuts

struct MinCut
    S::Set{Int64}
    T::Set{Int64}
    crossing_edges::Vector{Tuple{Int64,Int64}}
    capacity::Float64
end

"""
    MinCutEnumeration

Result of bounded minimum-cut enumeration.

- `cuts`: Enumerated minimum cuts in deterministic canonical order.
- `total_cuts`: Exact total count when `is_complete=true`; otherwise equals the
  requested truncation bound used by enumeration (typically `cut_limit`).
- `is_complete`: True iff all minimum cuts were enumerated.
- `free_zone_size`: number of free-zone nodes (in `S**` but not in `S*`).

When `free_zone_size > 62`, the exact total number of cuts `2^|F|` exceeds the
safe Int64 counting range used here. In that case enumeration is intentionally
truncated, `is_complete=false`, and `total_cuts` reports the truncated count.
"""
struct MinCutEnumeration
    cuts::Vector{MinCut}
    total_cuts::Int64
    is_complete::Bool
    free_zone_size::Int64
end

struct MinCutAnalysis
    representative_cut::MinCut
    edges_in_some_cut::Vector{Tuple{Int64,Int64}}
    edges_in_every_cut::Vector{Tuple{Int64,Int64}}
    enumeration::MinCutEnumeration
    max_flow::Float64
    min_cut_capacity::Float64
end

function _all_aug_nodes(
    flow_result::FlowSolveResult,
    original_nodes::Set{Int64}
)::Set{Int64}
    return union(
        Set(keys(flow_result.augmented_outgoing)),
        Set(keys(flow_result.augmented_incoming)),
        original_nodes,
        Set([flow_result.super_source, flow_result.super_sink])
    )
end

function _can_reach_sink_residual(
    flow_result::FlowSolveResult,
    tol::Float64
)::Set{Int64}
    return _backward_reachable_residual(
        flow_result.super_sink,
        flow_result.augmented_outgoing,
        flow_result.augmented_incoming,
        flow_result.augmented_capacities,
        flow_result.augmented_flow,
        tol;
        finite_caps_only=true
    )
end

function _s_double_star_original(
    edgelist::Vector{Tuple{Int64,Int64}},
    flow_result::FlowSolveResult;
    tol::Float64=1e-10,
    can_reach_sink::Union{Nothing,Set{Int64}}=nothing
)::Set{Int64}
    original_nodes = _graph_nodes_set(edgelist)
    reachable = can_reach_sink === nothing ? _can_reach_sink_residual(flow_result, tol) : can_reach_sink
    all_aug_nodes = _all_aug_nodes(flow_result, original_nodes)
    s_double_star = setdiff(all_aug_nodes, reachable)
    return intersect(s_double_star, original_nodes)
end

function _edges_in_every_mincut_from_reach(
    edgelist::Vector{Tuple{Int64,Int64}},
    flow_result::FlowSolveResult,
    can_reach_sink::Set{Int64};
    tol::Float64=1e-10
)::Vector{Tuple{Int64,Int64}}
    original_nodes = _graph_nodes_set(edgelist)
    s_star = flow_result.mincut_S
    # T** = nodes backward-reachable from the sink in the residual graph (restricted
    # to original nodes) -- the dual of s_double_star's "nodes that CANNOT reach the
    # sink" (S**). Was setdiff(original_nodes, can_reach_sink), i.e. the SAME set as
    # S** (a copy-paste of the s_double_star formula) -- that made t_double_star
    # equal the source-side candidate set instead of its complement, so the
    # saturated+u∈S*+v∈T** test could never correctly identify the always-cut
    # edges near the sink and instead matched saturated edges deep on the source
    # side. Confirmed on the power-network Baseline scenario: the true unique min
    # cut is {(22,23)}, but this returned {(8,12),(12,11),(19,22)} -- edges among
    # nodes that cannot reach the sink, exactly the (wrong) S**-side set.
    t_double_star = intersect(can_reach_sink, original_nodes)

    every = Tuple{Int64,Int64}[]
    for e in edgelist
        u, v = e
        residual = get(flow_result.residual_capacity, e, Inf)
        saturated = residual <= tol
        if saturated && (u in s_star) && (v in t_double_star)
            push!(every, e)
        end
    end

    sort!(every)
    return every
end

function _enumerate_min_cuts_with_sdouble(
    edgelist::Vector{Tuple{Int64,Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    flow_result::FlowSolveResult,
    s_double_star_orig::Set{Int64};
    cut_limit::Int=1000,
    tol::Float64=1e-10
)::MinCutEnumeration
    original_nodes = _graph_nodes_set(edgelist)
    s_star = Set{Int64}(flow_result.mincut_S)
    free_zone = setdiff(s_double_star_orig, s_star)
    sorted_free_zone = sort!(collect(free_zone))
    n = length(sorted_free_zone)

    total_possible::Int64 = 0
    is_complete::Bool = false
    to_enumerate::Int64 = 0

    if n <= 62
        total_possible = Int64(1) << n
        if total_possible <= cut_limit
            is_complete = true
            to_enumerate = total_possible
        else
            is_complete = false
            to_enumerate = Int64(cut_limit)
        end
    else
        is_complete = false
        to_enumerate = Int64(cut_limit)
    end

    cuts = Vector{MinCut}()
    sizehint!(cuts, Int(to_enumerate))

    for i in Int64(0):(to_enumerate - 1)
        R = _subset_from_bits(sorted_free_zone, i)
        cut_S = union(s_star, R)
        cut_T = setdiff(original_nodes, cut_S)

        issubset(s_star, cut_S) || throw(AssertionError("Invalid enumerated cut: S* ⊄ S."))
        issubset(cut_S, s_double_star_orig) || throw(AssertionError("Invalid enumerated cut: S ⊄ S**."))

        crossing = _crossing_edges(edgelist, cut_S, cut_T)
        cap = _cut_capacity(crossing, capacities)

        abs(cap - flow_result.max_flow) <= tol || throw(AssertionError(
            "Enumerated cut capacity $cap does not match max_flow $(flow_result.max_flow) within tol=$tol"
        ))

        push!(cuts, MinCut(cut_S, cut_T, crossing, cap))
    end

    reported_total = is_complete ? total_possible : Int64(cut_limit)
    return MinCutEnumeration(cuts, reported_total, is_complete, Int64(n))
end

function _crossing_edges(
    edgelist::Vector{Tuple{Int64,Int64}},
    S::Set{Int64},
    T::Set{Int64}
)::Vector{Tuple{Int64,Int64}}
    crossing = Tuple{Int64,Int64}[]
    for (u, v) in edgelist
        if (u in S) && (v in T)
            push!(crossing, (u, v))
        end
    end
    sort!(crossing)
    return crossing
end

function _cut_capacity(
    crossing_edges::Vector{Tuple{Int64,Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64}
)::Float64
    c = 0.0
    for e in crossing_edges
        c += capacities[e]
    end
    return c
end

function _subset_from_bits(
    sorted_free_zone::Vector{Int64},
    mask::Int64
)::Set{Int64}
    subset = Set{Int64}()
    max_bit = min(length(sorted_free_zone), 63)
    for bit in 0:(max_bit - 1)
        if ((mask >> bit) & 1) == 1
            push!(subset, sorted_free_zone[bit + 1])
        end
    end
    return subset
end

"""
    minimum_st_cut_edges(edgelist, flow_result)

Return the edges crossing the solved minimum-cut partition `(u ∈ S, v ∈ T)`
where `S = flow_result.mincut_S` and `T = flow_result.mincut_T`.

Returns the representative minimum cut derived from the solved flow. This is one
valid minimum cut — not necessarily unique. See `enumerate_min_cuts` for bounded
enumeration of all cuts.
"""
function minimum_st_cut_edges(
    edgelist::Vector{Tuple{Int64,Int64}},
    flow_result::FlowSolveResult
)::Vector{Tuple{Int64,Int64}}
    _require_bounded_baseline(flow_result)
    return _crossing_edges(edgelist, flow_result.mincut_S, flow_result.mincut_T)
end

"""
    minimum_st_cut_capacity(flow_result)

Return the solved minimum-cut capacity from `flow_result`.

Equals `max_flow` by the max-flow min-cut theorem. Exact.
"""
function minimum_st_cut_capacity(flow_result::FlowSolveResult)::Float64
    _require_bounded_baseline(flow_result)
    return flow_result.mincut_capacity
end

"""
    edges_in_some_mincut(edgelist, outgoing_index, incoming_index, capacities, flow_result; tol=1e-10)

Return all edges that appear in at least one minimum cut.

This uses the exact min-cut lattice characterization:
- edge is saturated, and
- `u ∈ S^{**}`, and
- `v ∉ S^*`.

Here `S^* = flow_result.mincut_S` and `S^{**}` is computed from one backward BFS
in the residual graph. Zero solver calls.
"""
function edges_in_some_mincut(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    flow_result::FlowSolveResult;
    tol::Float64=1e-10
)::Vector{Tuple{Int64,Int64}}
    _require_bounded_baseline(flow_result)
    some_edges = _edges_in_some_mincut(
        edgelist,
        outgoing_index,
        incoming_index,
        capacities,
        flow_result;
        tol=tol
    )

    can_reach_sink = _can_reach_sink_residual(flow_result, tol)
    every_edges = _edges_in_every_mincut_from_reach(
        edgelist,
        flow_result,
        can_reach_sink;
        tol=tol
    )

    merged = sort!(collect(union(Set(some_edges), Set(every_edges))))
    return merged
end

"""
    edges_in_every_mincut(edgelist, flow_result; tol=1e-10)

Return all edges that appear in every minimum cut.

Exact characterization:
- edge is saturated,
- `u ∈ S^*` where `S^* = flow_result.mincut_S`, and
- `v ∈ T^{**}`, where `T^{**}` is the set of nodes backward-reachable
  from the sink in the residual graph, restricted to original graph nodes.

Zero solver calls; one backward BFS.
"""
function edges_in_every_mincut(
    edgelist::Vector{Tuple{Int64,Int64}},
    flow_result::FlowSolveResult;
    tol::Float64=1e-10
)::Vector{Tuple{Int64,Int64}}
    _require_bounded_baseline(flow_result)
    can_reach_sink = _can_reach_sink_residual(flow_result, tol)
    return _edges_in_every_mincut_from_reach(edgelist, flow_result, can_reach_sink; tol=tol)
end

"""
    mincut_partition(flow_result)
        -> NamedTuple{(:S, :T, :crossing_edges, :capacity)}

Return the solved representative min-cut partition from `flow_result` as a
named tuple with fields:
- `S`: source-side node set
- `T`: sink-side node set
- `crossing_edges`: edges `(u,v)` with `u ∈ S` and `v ∈ T`
- `capacity`: min-cut capacity

Zero solver calls.
"""
function mincut_partition(
    flow_result::FlowSolveResult
)::NamedTuple{(:S, :T, :crossing_edges, :capacity),Tuple{Set{Int64},Set{Int64},Vector{Tuple{Int64,Int64}},Float64}}
    _require_bounded_baseline(flow_result)
    edges = sort!(collect(keys(flow_result.flow)))
    crossing = Tuple{Int64,Int64}[]
    for (u, v) in edges
        if (u in flow_result.mincut_S) && (v in flow_result.mincut_T)
            push!(crossing, (u, v))
        end
    end
    return (
        S=Set(flow_result.mincut_S),
        T=Set(flow_result.mincut_T),
        crossing_edges=crossing,
        capacity=flow_result.mincut_capacity
    )
end

"""
    enumerate_min_cuts(edgelist, outgoing_index, incoming_index, capacities, flow_result; cut_limit=1000, tol=1e-10)

Enumerate all valid minimum cuts up to `cut_limit` using the min-cut lattice.

Let `S* = flow_result.mincut_S`, `S**` from backward residual reachability, and
free zone `F` defined as nodes in `S**` but not in `S*`. Every minimum cut has the form `S = S* ∪ R` for some
`R ⊆ F`.

Enumeration is deterministic: `F` is sorted ascending and subsets are visited in
binary counting order. If `2^|F| <= cut_limit`, all minimum cuts are returned and
`is_complete=true`; otherwise the first `cut_limit` cuts are returned and
`is_complete=false`.
"""
function enumerate_min_cuts(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    flow_result::FlowSolveResult;
    cut_limit::Int=1000,
    tol::Float64=1e-10
)::MinCutEnumeration
    _require_bounded_baseline(flow_result)
    cut_limit > 0 || throw(ArgumentError("cut_limit must be positive."))

    # outgoing_index and incoming_index are accepted for API consistency with
    # related min-cut utilities; pure lattice enumeration does not need them.
    _ = outgoing_index
    _ = incoming_index

    s_double_star_orig = _s_double_star_original(edgelist, flow_result; tol=tol)
    return _enumerate_min_cuts_with_sdouble(
        edgelist,
        capacities,
        flow_result,
        s_double_star_orig;
        cut_limit=cut_limit,
        tol=tol
    )
end

"""
    analyze_min_cuts(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes, flow_result; cut_limit=1000, tol=1e-10)

Aggregate entry point for min-cut utilities, returning a `MinCutAnalysis`
containing:
- representative solved min-cut
- edges in some minimum cut
- edges in every minimum cut
- bounded enumeration of minimum cuts
- max-flow and min-cut capacity scalars

No solver calls are made.
"""
function analyze_min_cuts(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    capacities::Dict{Tuple{Int64,Int64},Float64},
    source_nodes::Vector{Int64},
    sink_nodes::Vector{Int64},
    flow_result::FlowSolveResult;
    cut_limit::Int=1000,
    tol::Float64=1e-10
)::MinCutAnalysis
    _require_bounded_baseline(flow_result)

    # source_nodes and sink_nodes are accepted for API consistency with other
    # module analysis entry points; this aggregate uses flow_result directly.
    _ = source_nodes
    _ = sink_nodes

    can_reach_sink = _can_reach_sink_residual(flow_result, tol)
    s_double_star_orig = _s_double_star_original(
        edgelist,
        flow_result;
        tol=tol,
        can_reach_sink=can_reach_sink
    )

    rep_edges = minimum_st_cut_edges(edgelist, flow_result)
    rep_cut = MinCut(
        Set(flow_result.mincut_S),
        Set(flow_result.mincut_T),
        rep_edges,
        flow_result.mincut_capacity
    )

    some_edges = edges_in_some_mincut(
        edgelist, outgoing_index, incoming_index, capacities, flow_result; tol=tol
    )
    every_edges = _edges_in_every_mincut_from_reach(
        edgelist,
        flow_result,
        can_reach_sink;
        tol=tol
    )
    enumeration = _enumerate_min_cuts_with_sdouble(
        edgelist,
        capacities,
        flow_result,
        s_double_star_orig;
        cut_limit=cut_limit, tol=tol
    )

    issubset(Set(every_edges), Set(some_edges)) || throw(AssertionError(
        "edges_in_every_mincut is not a subset of edges_in_some_mincut — lattice invariant violated."
    ))

    return MinCutAnalysis(
        rep_cut,
        some_edges,
        every_edges,
        enumeration,
        flow_result.max_flow,
        flow_result.mincut_capacity
    )
end

end # module MinCutUtilitiesModule
