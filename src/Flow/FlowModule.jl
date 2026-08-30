module FlowModule

export FlowSolveResult,
	   solve_max_flow_edmonds_karp,
	   solve_max_flow_dinic,
	   solve_max_flow_push_relabel,
	   sink_flows,
	   node_inflow,
	   node_outflow,
	   validate_capacity_constraints,
	   validate_flow_conservation,
	   validate_maxflow_mincut,
	   validate_exactness

struct FlowSolveResult
	max_flow::Float64
	flow::Dict{Tuple{Int64,Int64},Float64}
	augmented_flow::Dict{Tuple{Int64,Int64},Float64}
	augmented_outgoing::Dict{Int64,Set{Int64}}
	augmented_incoming::Dict{Int64,Set{Int64}}
	augmented_capacities::Dict{Tuple{Int64,Int64},Float64}
	residual_capacity::Dict{Tuple{Int64,Int64},Float64}
	node_flow::Dict{Int64,Float64}
	sources::Vector{Int64}
	sinks::Vector{Int64}
	super_source::Int64
	super_sink::Int64
	mincut_S::Set{Int64}
	mincut_T::Set{Int64}
	mincut_capacity::Float64
	saturated_edges::Vector{Tuple{Int64,Int64}}
	sink_flow::Dict{Int64,Float64}
	is_unbounded::Bool
end

function _has_infinite_augmenting_path(
	source::Int64,
	target::Int64,
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64},
	flow::Dict{Tuple{Int64,Int64},Float64},
	tol::Float64
)::Bool
	visited = Set{Int64}([source])
	queue = Int64[source]
	head = 1

	while head <= length(queue)
		u = queue[head]
		head += 1

		for v in get(outgoing_index, u, Set{Int64}())
			residual = capacities[(u, v)] - get(flow, (u, v), 0.0)
			if isinf(residual) && residual > tol && !(v in visited)
				v == target && return true
				push!(visited, v)
				push!(queue, v)
			end
		end

		for p in get(incoming_index, u, Set{Int64}())
			residual = get(flow, (p, u), 0.0)
			if isinf(residual) && residual > tol && !(p in visited)
				p == target && return true
				push!(visited, p)
				push!(queue, p)
			end
		end
	end

	return false
end

function _compute_original_residual_capacity(
	edgelist::Vector{Tuple{Int64,Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64},
	flow::Dict{Tuple{Int64,Int64},Float64}
)::Dict{Tuple{Int64,Int64},Float64}
	return Dict{Tuple{Int64,Int64},Float64}(e => capacities[e] - get(flow, e, 0.0) for e in edgelist)
end

function _compute_node_flow(
	graph_nodes::Set{Int64},
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}},
	flow::Dict{Tuple{Int64,Int64},Float64},
	source_nodes::Vector{Int64},
	sink_nodes::Vector{Int64},
	tol::Float64
)::Dict{Int64,Float64}
	sources = Set(source_nodes)
	sinks = Set(sink_nodes)
	node_flow = Dict{Int64,Float64}()

	for node in graph_nodes
		inflow = node_inflow(node, incoming_index, flow)
		outflow = node_outflow(node, outgoing_index, flow)
		if node in sources
			node_flow[node] = outflow
		elseif node in sinks
			node_flow[node] = inflow
		elseif abs(inflow - outflow) <= tol
			node_flow[node] = inflow
		else
			throw(ArgumentError("Flow conservation violated while computing node_flow at node $node: inflow=$inflow, outflow=$outflow."))
		end
	end

	return node_flow
end

function _compute_cut_capacity(
	augmented_capacities::Dict{Tuple{Int64,Int64},Float64},
	mincut_S::Set{Int64},
	mincut_T::Set{Int64},
	super_source::Int64,
	super_sink::Int64
)::Float64
	cut_capacity = 0.0
	for ((u, v), cap) in augmented_capacities
		if (u in mincut_S) && (v in mincut_T) && u != super_source && v != super_sink
			cut_capacity += cap
		end
	end
	return cut_capacity
end

function _bfs_augmenting_path(
	source::Int64,
	target::Int64,
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64},
	flow::Dict{Tuple{Int64,Int64},Float64},
	tol::Float64
)
	parent = Dict{Int64,Tuple{Int64,Symbol,Tuple{Int64,Int64}}}()
	visited = Set{Int64}([source])
	queue = Int64[source]
	head = 1

	while head <= length(queue)
		u = queue[head]
		head += 1

		for v in get(outgoing_index, u, Set{Int64}())
			edge = (u, v)
			residual = capacities[edge] - get(flow, edge, 0.0)
			if residual > tol && !(v in visited)
				parent[v] = (u, :forward, edge)
				v == target && return true, parent
				push!(visited, v)
				push!(queue, v)
			end
		end

		for p in get(incoming_index, u, Set{Int64}())
			edge = (p, u)
			residual = get(flow, edge, 0.0)
			if residual > tol && !(p in visited)
				# Backward residual step: we traverse from u -> p in the residual graph
				# using original edge (p,u). Parent map stores BFS predecessor in residual traversal,
				# so parent[p] = (u, :backward, (p,u)) is intentional.
				parent[p] = (u, :backward, edge)
				p == target && return true, parent
				push!(visited, p)
				push!(queue, p)
			end
		end
	end

	return false, parent
end

function _reachable_residual(
	source::Int64,
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64},
	flow::Dict{Tuple{Int64,Int64},Float64},
	tol::Float64
)::Set{Int64}
	reachable = Set{Int64}([source])
	queue = Int64[source]
	head = 1

	while head <= length(queue)
		u = queue[head]
		head += 1

		for v in get(outgoing_index, u, Set{Int64}())
			residual = capacities[(u, v)] - get(flow, (u, v), 0.0)
			if residual > tol && !(v in reachable)
				push!(reachable, v)
				push!(queue, v)
			end
		end

		for p in get(incoming_index, u, Set{Int64}())
			residual = get(flow, (p, u), 0.0)
			if residual > tol && !(p in reachable)
				push!(reachable, p)
				push!(queue, p)
			end
		end
	end

	return reachable
end

function _build_augmented_network(
	edgelist::Vector{Tuple{Int64,Int64}},
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64},
	source_nodes::Vector{Int64},
	sink_nodes::Vector{Int64}
)
	all_nodes = union(Set(first.(edgelist)), Set(last.(edgelist)))
	min_node = isempty(all_nodes) ? Int64(0) : minimum(all_nodes)
	min_node == typemin(Int64) && throw(ArgumentError(
		"Cannot construct super-source/super-sink IDs via minimum(node_ids)-1 and -2: minimum node ID is typemin(Int64). Remap node IDs to a safer range."
	))
	min_node <= typemin(Int64) + 1 && throw(ArgumentError(
		"Cannot construct both super-source and super-sink IDs without Int64 underflow from minimum node ID $min_node. Remap node IDs to a safer range."
	))
	super_source = min_node - 1
	super_sink = min_node - 2

	aug_out = Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in outgoing_index)
	aug_in = Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in incoming_index)
	aug_caps = Dict{Tuple{Int64,Int64},Float64}(capacities)
	aug_edges = copy(edgelist)

	for s in source_nodes
		edge = (super_source, s)
		aug_caps[edge] = Inf
		push!(aug_edges, edge)
		if !haskey(aug_out, super_source)
			aug_out[super_source] = Set{Int64}()
		end
		push!(aug_out[super_source], s)
		if !haskey(aug_in, s)
			aug_in[s] = Set{Int64}()
		end
		push!(aug_in[s], super_source)
	end

	for t in sink_nodes
		edge = (t, super_sink)
		aug_caps[edge] = Inf
		push!(aug_edges, edge)
		if !haskey(aug_out, t)
			aug_out[t] = Set{Int64}()
		end
		push!(aug_out[t], super_sink)
		if !haskey(aug_in, super_sink)
			aug_in[super_sink] = Set{Int64}()
		end
		push!(aug_in[super_sink], t)
	end

	return aug_edges, aug_out, aug_in, aug_caps, super_source, super_sink
end

"""
	solve_max_flow_edmonds_karp(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes; tol=1e-10, validate=true)

Exact max-flow solution with Edmonds-Karp on precomputed graph metadata from InputProcessing.
This API requires source/sink sets explicitly and does not perform auto-detection.

Mathematical constraints enforced:
- 0 ≤ f(u,v) ≤ c(u,v)
- Flow conservation at each non-source, non-sink node
- Max-flow = Min-cut capacity (unless unbounded flow exists)
"""
function solve_max_flow_edmonds_karp(
	edgelist::Vector{Tuple{Int64,Int64}},
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64},
	source_nodes::Vector{Int64},
	sink_nodes::Vector{Int64};
	tol::Float64=1e-10,
	validate::Bool=true
)::FlowSolveResult
	isempty(edgelist) && throw(ArgumentError("edgelist is empty."))
	isempty(source_nodes) && throw(ArgumentError("source_nodes cannot be empty."))
	isempty(sink_nodes) && throw(ArgumentError("sink_nodes cannot be empty."))
	!isempty(intersect(Set(source_nodes), Set(sink_nodes))) &&
		throw(ArgumentError("A node cannot be both source and sink."))

	graph_nodes = union(Set(first.(edgelist)), Set(last.(edgelist)))
	for s in source_nodes
		s in graph_nodes || throw(ArgumentError("source node $s is not present in graph nodes."))
	end
	for t in sink_nodes
		t in graph_nodes || throw(ArgumentError("sink node $t is not present in graph nodes."))
	end

	for e in edgelist
		haskey(capacities, e) || throw(ArgumentError("Missing capacity for edge $e"))
		c = capacities[e]
		(isnan(c) || c < 0.0) && throw(ArgumentError("Invalid capacity for edge $e: $c"))
	end

	aug_edges, aug_out, aug_in, aug_caps, super_source, super_sink =
		_build_augmented_network(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes)

	flow = Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in aug_edges)
	max_flow = 0.0
	is_unbounded = false
	_has_infinite_augmenting_path(super_source, super_sink, aug_out, aug_in, aug_caps, flow, tol) &&
		return FlowSolveResult(
			Inf,
			Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in edgelist),
			Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in aug_edges),
			Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_out),
			Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_in),
			Dict{Tuple{Int64,Int64},Float64}(aug_caps),
			Dict{Tuple{Int64,Int64},Float64}(e => capacities[e] for e in edgelist),
			# For an immediate unbounded detection, no finite flow has been pushed on original edges.
			# Using an empty original-flow dictionary yields zero inflow/outflow at all original nodes,
			# which is conservation-consistent for this sentinel return object.
			_compute_node_flow(graph_nodes, outgoing_index, incoming_index, Dict{Tuple{Int64,Int64},Float64}(), source_nodes, sink_nodes, tol),
			source_nodes,
			sink_nodes,
			super_source,
			super_sink,
			setdiff(graph_nodes, Set([super_source, super_sink])),
			Set{Int64}(),
			Inf,
			Tuple{Int64,Int64}[],
			Dict{Int64,Float64}(t => Inf for t in sink_nodes),
			true
		)

	while true
		found, parent = _bfs_augmenting_path(super_source, super_sink, aug_out, aug_in, aug_caps, flow, tol)
		!found && break

		path_edges = Tuple{Int64,Int64,Symbol,Tuple{Int64,Int64}}[]
		v = super_sink
		bottleneck = Inf

		while v != super_source
			haskey(parent, v) || throw(ArgumentError("Internal error: incomplete parent map during path reconstruction."))
			u, mode, edge = parent[v]
			residual = mode === :forward ? (aug_caps[edge] - get(flow, edge, 0.0)) : get(flow, edge, 0.0)
			bottleneck = min(bottleneck, residual)
			push!(path_edges, (u, v, mode, edge))
			v = u
		end

		if isinf(bottleneck)
			is_unbounded = true
			max_flow = Inf
			break
		end

		for (_, _, mode, edge) in path_edges
			if mode === :forward
				flow[edge] = get(flow, edge, 0.0) + bottleneck
			else
				flow[edge] = get(flow, edge, 0.0) - bottleneck
				if flow[edge] < tol
					flow[edge] = 0.0
				end
			end
		end

		max_flow += bottleneck
	end

	reachable = _reachable_residual(super_source, aug_out, aug_in, aug_caps, flow, tol)
	all_aug_nodes = union(
		Set(keys(aug_out)),
		Set(keys(aug_in)),
		graph_nodes,
		Set(source_nodes),
		Set(sink_nodes),
		Set([super_source, super_sink])
	)
	mincut_S = setdiff(reachable, Set([super_source, super_sink]))
	mincut_T = setdiff(all_aug_nodes, reachable, Set([super_source, super_sink]))
	cut_capacity = _compute_cut_capacity(aug_caps, union(mincut_S, Set([super_source])), union(mincut_T, Set([super_sink])), super_source, super_sink)

	original_flow = Dict{Tuple{Int64,Int64},Float64}(e => get(flow, e, 0.0) for e in edgelist)
	augmented_flow = Dict{Tuple{Int64,Int64},Float64}(flow)
	residual_capacity = _compute_original_residual_capacity(edgelist, capacities, original_flow)
	node_flow = _compute_node_flow(graph_nodes, outgoing_index, incoming_index, original_flow, source_nodes, sink_nodes, tol)
	saturated = [e for e in edgelist if isfinite(capacities[e]) && abs(get(original_flow, e, 0.0) - capacities[e]) <= tol]
	sink_flow = Dict{Int64,Float64}(t => get(flow, (t, super_sink), 0.0) for t in sink_nodes)

	result = FlowSolveResult(
		max_flow,
		original_flow,
		augmented_flow,
		Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_out),
		Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_in),
		Dict{Tuple{Int64,Int64},Float64}(aug_caps),
		residual_capacity,
		node_flow,
		source_nodes,
		sink_nodes,
		super_source,
		super_sink,
		mincut_S,
		mincut_T,
		cut_capacity,
		saturated,
		sink_flow,
		is_unbounded
	)

	if validate && !is_unbounded
		validate_exactness(result, edgelist, outgoing_index, incoming_index, capacities; tol=tol)
		validate_maxflow_mincut(result; tol=tol)
	end

	return result
end

function _build_level_graph_dinic(
	source::Int64,
	target::Int64,
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64},
	flow::Dict{Tuple{Int64,Int64},Float64},
	tol::Float64
)
	level = Dict{Int64,Int}()
	level[source] = 0
	queue = Int64[source]
	head = 1

	while head <= length(queue)
		u = queue[head]
		head += 1
		lu = level[u]

		for v in get(outgoing_index, u, Set{Int64}())
			residual = capacities[(u, v)] - get(flow, (u, v), 0.0)
			if residual > tol && !haskey(level, v)
				level[v] = lu + 1
				push!(queue, v)
			end
		end

		for p in get(incoming_index, u, Set{Int64}())
			residual = get(flow, (p, u), 0.0)
			if residual > tol && !haskey(level, p)
				level[p] = lu + 1
				push!(queue, p)
			end
		end
	end

	return level, haskey(level, target)
end

function _build_dinic_adjacency(
	level::Dict{Int64,Int},
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64},
	flow::Dict{Tuple{Int64,Int64},Float64},
	tol::Float64
)
	adj = Dict{Int64,Vector{Tuple{Int64,Symbol,Tuple{Int64,Int64}}}}()
	for (u, lu) in level
		moves = Tuple{Int64,Symbol,Tuple{Int64,Int64}}[]

		for v in get(outgoing_index, u, Set{Int64}())
			if get(level, v, -1) == lu + 1
				residual = capacities[(u, v)] - get(flow, (u, v), 0.0)
				residual > tol && push!(moves, (v, :forward, (u, v)))
			end
		end

		for p in get(incoming_index, u, Set{Int64}())
			if get(level, p, -1) == lu + 1
				residual = get(flow, (p, u), 0.0)
				residual > tol && push!(moves, (p, :backward, (p, u)))
			end
		end

		adj[u] = moves
	end
	return adj
end

function _dinic_dfs_blocking(
	u::Int64,
	target::Int64,
	pushed::Float64,
	adj::Dict{Int64,Vector{Tuple{Int64,Symbol,Tuple{Int64,Int64}}}},
	ptr::Dict{Int64,Int},
	capacities::Dict{Tuple{Int64,Int64},Float64},
	flow::Dict{Tuple{Int64,Int64},Float64},
	tol::Float64
)::Float64
	u == target && return pushed
	moves = get(adj, u, Tuple{Int64,Symbol,Tuple{Int64,Int64}}[])
	start_i = get(ptr, u, 1)

	for i in start_i:length(moves)
		ptr[u] = i
		v, mode, edge = moves[i]
		residual = mode === :forward ? (capacities[edge] - get(flow, edge, 0.0)) : get(flow, edge, 0.0)
		residual <= tol && continue

		candidate = min(pushed, residual)
		candidate <= tol && continue
		tr = _dinic_dfs_blocking(v, target, candidate, adj, ptr, capacities, flow, tol)

		if tr > tol
			if mode === :forward
				flow[edge] = get(flow, edge, 0.0) + tr
			else
				flow[edge] = get(flow, edge, 0.0) - tr
				if flow[edge] < tol
					flow[edge] = 0.0
				end
			end
			return tr
		end
	end

	ptr[u] = length(moves) + 1
	return 0.0
end

function _build_push_relabel_neighbors(
	aug_out::Dict{Int64,Set{Int64}},
	aug_in::Dict{Int64,Set{Int64}}
)::Dict{Int64,Vector{Int64}}
	nodes = union(Set(keys(aug_out)), Set(keys(aug_in)))
	neighbors = Dict{Int64,Vector{Int64}}()
	for u in nodes
		merged = union(get(aug_out, u, Set{Int64}()), get(aug_in, u, Set{Int64}()))
		neighbors[u] = sort!(collect(merged))
	end
	return neighbors
end

function _push_relabel_residual(
	u::Int64,
	v::Int64,
	aug_caps::Dict{Tuple{Int64,Int64},Float64},
	flow::Dict{Tuple{Int64,Int64},Float64},
	tol::Float64
)::Tuple{Float64,Symbol,Tuple{Int64,Int64}}
	forward = (u, v)
	if haskey(aug_caps, forward)
		residual = aug_caps[forward] - get(flow, forward, 0.0)
		if residual > tol
			return residual, :forward, forward
		end
	end

	backward = (v, u)
	if haskey(aug_caps, backward)
		residual = get(flow, backward, 0.0)
		if residual > tol
			return residual, :backward, backward
		end
	end

	return 0.0, :none, (u, v)
end

function _push_relabel_try_push!(
	u::Int64,
	v::Int64,
	flow::Dict{Tuple{Int64,Int64},Float64},
	aug_caps::Dict{Tuple{Int64,Int64},Float64},
	height::Dict{Int64,Int64},
	excess::Dict{Int64,Float64},
	active::Set{Int64},
	source::Int64,
	target::Int64,
	tol::Float64
)::Bool
	height[u] == height[v] + 1 || return false

	forward_edge = (u, v)
	if haskey(aug_caps, forward_edge)
		residual_forward = aug_caps[forward_edge] - get(flow, forward_edge, 0.0)
		if residual_forward > tol
			delta = min(get(excess, u, 0.0), residual_forward)
			if delta > tol
				flow[forward_edge] = get(flow, forward_edge, 0.0) + delta
				excess[u] = get(excess, u, 0.0) - delta
				excess[v] = get(excess, v, 0.0) + delta
				if v != source && v != target && excess[v] > tol
					push!(active, v)
				end
				return true
			end
		end
	end

	backward_edge = (v, u)
	if haskey(aug_caps, backward_edge)
		residual_backward = get(flow, backward_edge, 0.0)
		if residual_backward > tol
			delta = min(get(excess, u, 0.0), residual_backward)
			if delta > tol
				flow[backward_edge] = get(flow, backward_edge, 0.0) - delta
				if flow[backward_edge] < tol
					flow[backward_edge] = 0.0
				end
				excess[u] = get(excess, u, 0.0) - delta
				excess[v] = get(excess, v, 0.0) + delta
				if v != source && v != target && excess[v] > tol
					push!(active, v)
				end
				return true
			end
		end
	end

	return false
end

function _push_relabel_relabel!(
	u::Int64,
	neighbors::Dict{Int64,Vector{Int64}},
	flow::Dict{Tuple{Int64,Int64},Float64},
	aug_caps::Dict{Tuple{Int64,Int64},Float64},
	height::Dict{Int64,Int64},
	tol::Float64
)::Nothing
	min_h = typemax(Int64)
	for v in get(neighbors, u, Int64[])
		residual, _, _ = _push_relabel_residual(u, v, aug_caps, flow, tol)
		residual <= tol && continue
		min_h = min(min_h, get(height, v, 0))
	end

	min_h == typemax(Int64) && throw(ArgumentError("Push-relabel relabel failed at node $u: no residual neighbor exists while positive excess remains."))
	height[u] = min_h + 1
	nothing
end

"""
	solve_max_flow_push_relabel(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes; tol=1e-10, validate=true)

Exact max-flow solution with push-relabel preflow algorithm on precomputed graph metadata
from InputProcessing. This API requires source/sink sets explicitly.

Mathematical constraints enforced:
- 0 ≤ f(u,v) ≤ c(u,v)
- Flow conservation at each non-source, non-sink node
- Max-flow = Min-cut capacity (unless unbounded flow exists)
"""
function solve_max_flow_push_relabel(
	edgelist::Vector{Tuple{Int64,Int64}},
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64},
	source_nodes::Vector{Int64},
	sink_nodes::Vector{Int64};
	tol::Float64=1e-10,
	validate::Bool=true
)::FlowSolveResult
	isempty(edgelist) && throw(ArgumentError("edgelist is empty."))
	isempty(source_nodes) && throw(ArgumentError("source_nodes cannot be empty."))
	isempty(sink_nodes) && throw(ArgumentError("sink_nodes cannot be empty."))
	!isempty(intersect(Set(source_nodes), Set(sink_nodes))) &&
		throw(ArgumentError("A node cannot be both source and sink."))

	graph_nodes = union(Set(first.(edgelist)), Set(last.(edgelist)))
	for s in source_nodes
		s in graph_nodes || throw(ArgumentError("source node $s is not present in graph nodes."))
	end
	for t in sink_nodes
		t in graph_nodes || throw(ArgumentError("sink node $t is not present in graph nodes."))
	end

	for e in edgelist
		haskey(capacities, e) || throw(ArgumentError("Missing capacity for edge $e"))
		c = capacities[e]
		(isnan(c) || c < 0.0) && throw(ArgumentError("Invalid capacity for edge $e: $c"))
	end

	aug_edges, aug_out, aug_in, aug_caps, super_source, super_sink =
		_build_augmented_network(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes)

	flow = Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in aug_edges)
	_has_infinite_augmenting_path(super_source, super_sink, aug_out, aug_in, aug_caps, flow, tol) &&
		return FlowSolveResult(
			Inf,
			Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in edgelist),
			Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in aug_edges),
			Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_out),
			Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_in),
			Dict{Tuple{Int64,Int64},Float64}(aug_caps),
			Dict{Tuple{Int64,Int64},Float64}(e => capacities[e] for e in edgelist),
			_compute_node_flow(graph_nodes, outgoing_index, incoming_index, Dict{Tuple{Int64,Int64},Float64}(), source_nodes, sink_nodes, tol),
			source_nodes,
			sink_nodes,
			super_source,
			super_sink,
			setdiff(graph_nodes, Set([super_source, super_sink])),
			Set{Int64}(),
			Inf,
			Tuple{Int64,Int64}[],
			Dict{Int64,Float64}(t => Inf for t in sink_nodes),
			true
		)

	nodes = union(Set(keys(aug_out)), Set(keys(aug_in)), Set([super_source, super_sink]))
	height = Dict{Int64,Int64}(u => 0 for u in nodes)
	excess = Dict{Int64,Float64}(u => 0.0 for u in nodes)
	neighbors = _build_push_relabel_neighbors(aug_out, aug_in)
	total_finite_capacity = sum(c for c in values(capacities) if isfinite(c))
	@assert total_finite_capacity >= 0.0 "total_finite_capacity must be nonnegative"
	source_dispatch_bound = Dict{Int64,Float64}()
	for s in source_nodes
		bound = 0.0
		has_infinite_out = false
		for v in get(outgoing_index, s, Set{Int64}())
			cap = capacities[(s, v)]
			if isfinite(cap)
				bound += cap
			else
				has_infinite_out = true
			end
		end
		# Dispatch bound for push-relabel initialization.
		# For sources with only finite outgoing edges: exact bound is the
		# sum of outgoing capacities — the source cannot push more regardless.
		# For sources with any infinite outgoing edge: we use total_finite_capacity,
		# which is a provably sufficient bound by the following argument:
		#   max_flow = min_cut (max-flow min-cut theorem)
		#   Infinite-capacity edges cannot appear in any finite min-cut
		#   (they would make the cut capacity infinite).
		#   Therefore min_cut ≤ sum of all finite capacities = total_finite_capacity.
		#   Dispatching total_finite_capacity guarantees the algorithm is never
		#   starved of initial excess, regardless of network topology.
		#   Excess above max_flow is returned to super_source during execution.
		source_dispatch_bound[s] = has_infinite_out ? total_finite_capacity : bound
	end

	height[super_source] = Int64(length(nodes))
	for v in get(aug_out, super_source, Set{Int64}())
		edge = (super_source, v)
		residual = aug_caps[edge] - get(flow, edge, 0.0)
		dispatch = min(residual, get(source_dispatch_bound, v, 0.0))
		if dispatch > tol
			flow[edge] = get(flow, edge, 0.0) + dispatch
			excess[super_source] = get(excess, super_source, 0.0) - dispatch
			excess[v] = get(excess, v, 0.0) + dispatch
		end
	end

	# Unbounded flow note: the _has_infinite_augmenting_path pre-check
	# above guarantees no infinite augmenting path exists in the residual
	# graph before execution begins. Push-relabel saturates edges from finite
	# initial dispatch bounds and cannot generate infinite excess at any node
	# from finite capacity inputs. No mid-execution unbounded check is therefore
	# needed here, unlike augmenting-path algorithms which encounter infinite
	# bottlenecks directly.
	active = Set{Int64}()
	for u in nodes
		if u != super_source && u != super_sink && get(excess, u, 0.0) > tol
			push!(active, u)
		end
	end

	while !isempty(active)
		u = minimum(active)
		delete!(active, u)

		while get(excess, u, 0.0) > tol
			pushed = false
			for v in get(neighbors, u, Int64[])
				if _push_relabel_try_push!(u, v, flow, aug_caps, height, excess, active, super_source, super_sink, tol)
					pushed = true
					get(excess, u, 0.0) <= tol && break
				end
			end

			if !pushed
				_push_relabel_relabel!(u, neighbors, flow, aug_caps, height, tol)
			end
		end
	end

	is_unbounded = false
	reachable = _reachable_residual(super_source, aug_out, aug_in, aug_caps, flow, tol)
	all_aug_nodes = union(
		Set(keys(aug_out)),
		Set(keys(aug_in)),
		graph_nodes,
		Set(source_nodes),
		Set(sink_nodes),
		Set([super_source, super_sink])
	)
	mincut_S = setdiff(reachable, Set([super_source, super_sink]))
	mincut_T = setdiff(all_aug_nodes, reachable, Set([super_source, super_sink]))
	cut_capacity = _compute_cut_capacity(aug_caps, union(mincut_S, Set([super_source])), union(mincut_T, Set([super_sink])), super_source, super_sink)

	original_flow = Dict{Tuple{Int64,Int64},Float64}(e => get(flow, e, 0.0) for e in edgelist)
	augmented_flow = Dict{Tuple{Int64,Int64},Float64}(flow)
	residual_capacity = _compute_original_residual_capacity(edgelist, capacities, original_flow)
	node_flow = _compute_node_flow(graph_nodes, outgoing_index, incoming_index, original_flow, source_nodes, sink_nodes, tol)
	saturated = [e for e in edgelist if isfinite(capacities[e]) && abs(get(original_flow, e, 0.0) - capacities[e]) <= tol]
	sink_flow = Dict{Int64,Float64}(t => get(flow, (t, super_sink), 0.0) for t in sink_nodes)
	max_flow = sum(values(sink_flow))

	result = FlowSolveResult(
		max_flow,
		original_flow,
		augmented_flow,
		Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_out),
		Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_in),
		Dict{Tuple{Int64,Int64},Float64}(aug_caps),
		residual_capacity,
		node_flow,
		source_nodes,
		sink_nodes,
		super_source,
		super_sink,
		mincut_S,
		mincut_T,
		cut_capacity,
		saturated,
		sink_flow,
		is_unbounded
	)

	if validate && !is_unbounded
		validate_exactness(result, edgelist, outgoing_index, incoming_index, capacities; tol=tol)
		validate_maxflow_mincut(result; tol=tol)
	end

	return result
end

"""
	solve_max_flow_dinic(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes; tol=1e-10, validate=true)

Exact max-flow solution with Dinic's algorithm on precomputed graph metadata from
InputProcessing. This API requires source/sink sets explicitly.

Mathematical constraints enforced:
- 0 ≤ f(u,v) ≤ c(u,v)
- Flow conservation at each non-source, non-sink node
- Max-flow = Min-cut capacity (unless unbounded flow exists)
"""
function solve_max_flow_dinic(
	edgelist::Vector{Tuple{Int64,Int64}},
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64},
	source_nodes::Vector{Int64},
	sink_nodes::Vector{Int64};
	tol::Float64=1e-10,
	validate::Bool=true
)::FlowSolveResult
	isempty(edgelist) && throw(ArgumentError("edgelist is empty."))
	isempty(source_nodes) && throw(ArgumentError("source_nodes cannot be empty."))
	isempty(sink_nodes) && throw(ArgumentError("sink_nodes cannot be empty."))
	!isempty(intersect(Set(source_nodes), Set(sink_nodes))) &&
		throw(ArgumentError("A node cannot be both source and sink."))

	graph_nodes = union(Set(first.(edgelist)), Set(last.(edgelist)))
	for s in source_nodes
		s in graph_nodes || throw(ArgumentError("source node $s is not present in graph nodes."))
	end
	for t in sink_nodes
		t in graph_nodes || throw(ArgumentError("sink node $t is not present in graph nodes."))
	end

	for e in edgelist
		haskey(capacities, e) || throw(ArgumentError("Missing capacity for edge $e"))
		c = capacities[e]
		(isnan(c) || c < 0.0) && throw(ArgumentError("Invalid capacity for edge $e: $c"))
	end

	aug_edges, aug_out, aug_in, aug_caps, super_source, super_sink =
		_build_augmented_network(edgelist, outgoing_index, incoming_index, capacities, source_nodes, sink_nodes)

	flow = Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in aug_edges)
	max_flow = 0.0
	is_unbounded = false
	_has_infinite_augmenting_path(super_source, super_sink, aug_out, aug_in, aug_caps, flow, tol) &&
		return FlowSolveResult(
			Inf,
			Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in edgelist),
			Dict{Tuple{Int64,Int64},Float64}(e => 0.0 for e in aug_edges),
			Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_out),
			Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_in),
			Dict{Tuple{Int64,Int64},Float64}(aug_caps),
			Dict{Tuple{Int64,Int64},Float64}(e => capacities[e] for e in edgelist),
			_compute_node_flow(graph_nodes, outgoing_index, incoming_index, Dict{Tuple{Int64,Int64},Float64}(), source_nodes, sink_nodes, tol),
			source_nodes,
			sink_nodes,
			super_source,
			super_sink,
			setdiff(graph_nodes, Set([super_source, super_sink])),
			Set{Int64}(),
			Inf,
			Tuple{Int64,Int64}[],
			Dict{Int64,Float64}(t => Inf for t in sink_nodes),
			true
		)

	while true
		level, reachable_sink = _build_level_graph_dinic(super_source, super_sink, aug_out, aug_in, aug_caps, flow, tol)
		reachable_sink || break

		adj = _build_dinic_adjacency(level, aug_out, aug_in, aug_caps, flow, tol)
		ptr = Dict{Int64,Int}(u => 1 for u in keys(level))

		while true
			pushed = _dinic_dfs_blocking(super_source, super_sink, Inf, adj, ptr, aug_caps, flow, tol)
			pushed <= tol && break

			if isinf(pushed)
				is_unbounded = true
				max_flow = Inf
				break
			end

			max_flow += pushed
		end

		is_unbounded && break
	end

	reachable = _reachable_residual(super_source, aug_out, aug_in, aug_caps, flow, tol)
	all_aug_nodes = union(
		Set(keys(aug_out)),
		Set(keys(aug_in)),
		graph_nodes,
		Set(source_nodes),
		Set(sink_nodes),
		Set([super_source, super_sink])
	)
	mincut_S = setdiff(reachable, Set([super_source, super_sink]))
	mincut_T = setdiff(all_aug_nodes, reachable, Set([super_source, super_sink]))
	cut_capacity = _compute_cut_capacity(aug_caps, union(mincut_S, Set([super_source])), union(mincut_T, Set([super_sink])), super_source, super_sink)

	original_flow = Dict{Tuple{Int64,Int64},Float64}(e => get(flow, e, 0.0) for e in edgelist)
	augmented_flow = Dict{Tuple{Int64,Int64},Float64}(flow)
	residual_capacity = _compute_original_residual_capacity(edgelist, capacities, original_flow)
	node_flow = _compute_node_flow(graph_nodes, outgoing_index, incoming_index, original_flow, source_nodes, sink_nodes, tol)
	saturated = [e for e in edgelist if isfinite(capacities[e]) && abs(get(original_flow, e, 0.0) - capacities[e]) <= tol]
	sink_flow = Dict{Int64,Float64}(t => get(flow, (t, super_sink), 0.0) for t in sink_nodes)

	result = FlowSolveResult(
		max_flow,
		original_flow,
		augmented_flow,
		Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_out),
		Dict{Int64,Set{Int64}}(k => copy(v) for (k, v) in aug_in),
		Dict{Tuple{Int64,Int64},Float64}(aug_caps),
		residual_capacity,
		node_flow,
		source_nodes,
		sink_nodes,
		super_source,
		super_sink,
		mincut_S,
		mincut_T,
		cut_capacity,
		saturated,
		sink_flow,
		is_unbounded
	)

	if validate && !is_unbounded
		validate_exactness(result, edgelist, outgoing_index, incoming_index, capacities; tol=tol)
		validate_maxflow_mincut(result; tol=tol)
	end

	return result
end

"""
	sink_flows(result)

Return sink-wise throughput values.
"""
sink_flows(result::FlowSolveResult) = result.sink_flow

"""
	node_inflow(node, incoming_index, flow)

Compute total inflow at a node from directed edge flows.
"""
function node_inflow(
	node::Int64,
	incoming_index::Dict{Int64,Set{Int64}},
	flow::Dict{Tuple{Int64,Int64},Float64}
)::Float64
	total = 0.0
	for u in get(incoming_index, node, Set{Int64}())
		total += get(flow, (u, node), 0.0)
	end
	return total
end

"""
	node_outflow(node, outgoing_index, flow)

Compute total outflow at a node from directed edge flows.
"""
function node_outflow(
	node::Int64,
	outgoing_index::Dict{Int64,Set{Int64}},
	flow::Dict{Tuple{Int64,Int64},Float64}
)::Float64
	total = 0.0
	for v in get(outgoing_index, node, Set{Int64}())
		total += get(flow, (node, v), 0.0)
	end
	return total
end

"""
	validate_capacity_constraints(result, edgelist, capacities; tol=1e-10)

Verify 0 ≤ f(u,v) ≤ c(u,v) for every original edge.
"""
function validate_capacity_constraints(
	result::FlowSolveResult,
	edgelist::Vector{Tuple{Int64,Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64};
	tol::Float64=1e-10
)::Nothing
	for e in edgelist
		f = get(result.flow, e, 0.0)
		c = capacities[e]
		if f < -tol
			throw(ArgumentError("Capacity constraint violated: flow on edge $e is negative ($f)."))
		end
		if isfinite(c) && f - c > tol
			throw(ArgumentError("Capacity constraint violated: flow($e)=$f exceeds capacity=$c."))
		end
	end
	nothing
end

"""
	validate_flow_conservation(result, outgoing_index, incoming_index; tol=1e-10)

Verify inflow equals outflow for each non-source and non-sink original node.
"""
function validate_flow_conservation(
	result::FlowSolveResult,
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}};
	tol::Float64=1e-10
)::Nothing
	all_nodes = union(Set(keys(outgoing_index)), Set(keys(incoming_index)))
	terminals = union(Set(result.sources), Set(result.sinks))

	for node in all_nodes
		node in terminals && continue
		inflow = node_inflow(node, incoming_index, result.flow)
		outflow = node_outflow(node, outgoing_index, result.flow)
		if abs(inflow - outflow) > tol
			throw(ArgumentError("Flow conservation violated at node $node: inflow=$inflow, outflow=$outflow."))
		end
	end
	nothing
end

"""
	validate_maxflow_mincut(result; tol=1e-10)

Verify max-flow min-cut theorem numerically using stored min-cut capacity in the result.
"""
function validate_maxflow_mincut(
	result::FlowSolveResult;
	tol::Float64=1e-10
)::Nothing
	if !(isinf(result.max_flow) && isinf(result.mincut_capacity)) &&
		abs(result.max_flow - result.mincut_capacity) > tol
		throw(ArgumentError("Max-flow/min-cut mismatch: max_flow=$(result.max_flow), cut_capacity=$(result.mincut_capacity)."))
	end
	nothing
end

"""
	validate_exactness(result, edgelist, outgoing_index, incoming_index, capacities; tol=1e-10)

Run exactness checks (capacity constraints + flow conservation).
"""
function validate_exactness(
	result::FlowSolveResult,
	edgelist::Vector{Tuple{Int64,Int64}},
	outgoing_index::Dict{Int64,Set{Int64}},
	incoming_index::Dict{Int64,Set{Int64}},
	capacities::Dict{Tuple{Int64,Int64},Float64};
	tol::Float64=1e-10
)::Nothing
	validate_capacity_constraints(result, edgelist, capacities; tol=tol)
	validate_flow_conservation(result, outgoing_index, incoming_index; tol=tol)
	nothing
end

end
