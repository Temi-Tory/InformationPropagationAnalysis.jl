# DiamondProcessing Module

## Overview

The `DiamondProcessing` module implements diamond structure identification and preprocessing for directed acyclic graphs (DAGs). It contains three main algorithms: the 8-step diamond identification process, sequential diamond storage building, and parallel diamond processing with thread-safe optimization.

## Core Data Structures

### Diamond
Represents a diamond structure in the network:
- `relevant_nodes`: All nodes within the diamond subgraph
- `conditioning_nodes`: Source nodes that control diamond behavior
- `edgelist`: All edges within the diamond structure

### DiamondsAtNode
Associates diamond structures with join nodes:
- `diamond`: The diamond structure
- `non_diamond_parents`: Parents not part of diamond structures
- `join_node`: The convergence point

### DiamondComputationData
Precomputed subgraph structure for efficient reuse:
- `sub_outgoing_index`, `sub_incoming_index`: Edge indices within diamond
- `sub_sources`, `sub_fork_nodes`, `sub_join_nodes`: Node classifications
- `sub_ancestors`, `sub_descendants`: Precomputed reachability
- `sub_iteration_sets`: Topological processing order
- `sub_node_priors`: Probability assignments
- `sub_diamond_structures`: Nested diamond structures within this diamond

### DiamondOptimizationContext
Caching system for expensive operations:
- `set_intersection_cache`, `set_difference_cache`: Set operation results
- `edge_filter_cache`: Filtered edge lists
- `ancestor_intersections`, `descendant_intersections`: Reachability intersections

## Algorithm 1: `identify_and_group_diamonds`

The core 8-step diamond identification algorithm that processes each join node:

### Step 1: Collect Shared Fork Ancestors
```julia
shared_fork_ancestors, diamond_parents, parents = collect_shared_fork_ancestors(
    join_node, incoming_index, ancestors, fork_nodes, irrelevant_sources, ctx
)
```
Identifies fork nodes that influence multiple parents of the join node.

### Step 2: Extract Induced Edgelist
```julia
induced_edgelist, relevant_nodes_for_induced = extract_induced_edgelist(
    shared_fork_ancestors, join_node, ancestors, descendants, edgelist, ctx
)
```
Builds subgraph containing all paths between shared forks and join node.

### Step 3: Identify Diamond Sources and Conditioning Nodes
```julia
diamond_sourcenodes, relevant_nodes, conditioning_nodes = identify_diamond_sources_and_conditioning(
    induced_edgelist, relevant_nodes_for_induced, exluded_nodes
)
```
Determines which nodes are sources (no incoming edges) and will serve as conditioning nodes.

### Step 4: Identify Intermediate Nodes
```julia
intermediate_nodes = identify_intermediate_nodes(relevant_nodes, conditioning_nodes, join_node)
```
Finds nodes that are neither conditioning nodes nor the join node.

### Step 5: Ensure Intermediate Incoming Edges
```julia
final_edgelist, final_relevant_nodes_for_induced, nodes_added_in_step8 = ensure_intermediate_incoming_edges(
    intermediate_nodes, incoming_index, induced_edgelist, relevant_nodes_for_induced
)
```
Adds missing incoming edges for intermediate nodes to complete the subgraph.

### Step 6: Perform Subsource Analysis
```julia
final_edgelist, final_relevant_nodes_for_induced, final_diamond_sourcenodes = perform_subsource_analysis(
    final_edgelist, final_relevant_nodes_for_induced, ancestors, descendants,
    irrelevant_sources, join_node, exluded_nodes, edgelist, ctx
)
```
Iteratively finds shared ancestors among current sources and replaces them with earlier common ancestors.

### Step 7: Recursive Diamond Completeness
```julia
final_edgelist, final_relevant_nodes_for_induced, final_shared_fork_ancestors, final_highest_nodes = perform_recursive_diamond_completeness(
    final_edgelist, final_relevant_nodes_for_induced, final_diamond_sourcenodes, shared_fork_ancestors,
    ancestors, descendants, fork_nodes, irrelevant_sources, incoming_index, join_node, exluded_nodes, edgelist, ctx
)
```
Recursively ensures diamond completeness by finding additional shared fork ancestors among diamond sources.

### Step 8: Build Final Diamond Structure
```julia
diamond, non_diamond_parents = build_final_diamond_structure(
    final_edgelist, final_relevant_nodes_for_induced, final_shared_fork_ancestors, final_highest_nodes,
    parents, diamond_parents, join_node, exluded_nodes, ctx
)
```
Constructs the final Diamond object with all identified components.

## Algorithm 2: `build_unique_diamond_storage`

Sequential algorithm that builds precomputed structures for all diamonds using a LIFO work stack:

### Process Flow:
1. **Initialize Work Stack**: Add root diamonds grouped by iteration level
2. **Process Each Diamond**: Pop from stack, compute subgraph structure
3. **Recursive Discovery**: 
   - For root diamonds: Use full `identify_and_group_diamonds`
   - For sub-diamonds: Use `perform_hybrid_diamond_lookup` with lookup table optimization
4. **Store Computation Data**: Create `DiamondComputationData` with precomputed indices
5. **Add Sub-diamonds**: Push newly discovered diamonds onto work stack
6. **Memory Management**: Clear caches periodically to manage memory usage

### Key Features:
- **LIFO Processing**: Ensures proper dependency resolution order
- **Hybrid Optimization**: Uses lookup table for sub-diamonds to avoid recomputation  
- **Duplicate Detection**: Skips already processed diamonds using hash-based tracking
- **Cache Management**: Adaptive cache clearing based on memory pressure

## Algorithm 3: `build_unique_diamond_storage_depth_first_parallel`

Parallel version that processes diamond subtrees concurrently:

### Parallelization Strategy:
1. **Sequential Iteration Levels**: Process iteration levels in order to preserve dependencies
2. **Parallel Subtrees**: Within each level, process diamond subtrees in parallel threads
3. **Thread-Local Processing**: Each thread maintains its own LIFO stack and lookup table
4. **Thread-Safe Merging**: Results merged with locks to shared storage

### Thread Management:
```julia
Threads.@threads for i in eachindex(level_diamonds)
    thread_ctx = DiamondOptimizationContext()
    thread_local_lookup_table = deepcopy(shared_diamond_lookup_table)
    
    thread_results = process_diamond_subtree_sequential_lifo_with_lookup(...)
    
    # Thread-safe merging
    lock(results_lock) do
        merge!(unique_diamonds, thread_results)
    end
end
```

### Memory Management:
- **Adaptive Thresholds**: Cache clearing frequency based on processing scale
- **Garbage Collection**: Forced GC between iteration levels for large datasets
- **Thread-Local Caches**: Each thread maintains separate optimization context

## Supporting Functions

### Caching Operations
- `cached_intersect`, `cached_setdiff`: Cached set operations
- `cached_filter_edges`: Edge filtering with memoization
- `get_cached_ancestor_intersection`: Ancestor intersection with caching

### Memory Management
- `manage_memory_adaptive`: Clears caches and forces GC based on processing scale
- `DiamondOptimizationContext`: Maintains all caching structures

### Hybrid Optimization
- `perform_hybrid_diamond_lookup`: Tries lookup table first, falls back to full computation
- Scores candidates based on edge containment and conditioning conflicts

## Performance Characteristics

### Sequential Algorithm (`build_unique_diamond_storage`)
- **Time**: O(d × 2^c) where d = diamonds, c = avg conditioning nodes
- **Space**: O(d × s) where s = avg diamond size
- **Use Case**: Smaller networks, deterministic processing order

### Parallel Algorithm (`build_unique_diamond_storage_depth_first_parallel`)
- **Time**: O((d × 2^c) / t) where t = threads
- **Space**: O(d × s × t) due to thread-local storage
- **Use Case**: Large networks, multi-core systems

### Memory Optimizations
- **Cache Clearing**: Prevents memory growth on large datasets
- **Hash-Based Deduplication**: O(1) duplicate detection
- **LIFO Processing**: Minimizes memory footprint vs BFS

## Integration

The module works with:
- **InputProcessing**: Uses graph structures and probability data
- **ReachabilityAnalysis**: Provides precomputed `DiamondComputationData` for efficient belief propagation
- **ProbabilityBoundsAnalysis**: Supports pbox probability types

The diamond identification and preprocessing separate network structure analysis from reachability computation, allowing for efficient reuse of structural information across multiple analysis runs.