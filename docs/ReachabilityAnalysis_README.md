# ReachabilityAnalysis Module

## Overview

The `ReachabilityAnalysis` module computes exact belief propagation for network reachability analysis in directed acyclic graphs (DAGs). It handles complex dependency structures using diamond identification, precomputed subgraph structures, and caching.

## Core Components

### Probability Type Support
Supports three probability representations:
- **Float64**: Standard floating-point probabilities
- **Interval**: Interval arithmetic for uncertainty bounds  
- **pbox**: Probability boxes for distributional uncertainty

### Precomputed Diamond Structures
Uses `DiamondComputationData` from the DiamondProcessing module:
- Prebuilt subgraph indices (`sub_outgoing_index`, `sub_incoming_index`)
- Identified diamond structures within diamonds (`sub_diamond_structures`)
- Topological orderings (`sub_iteration_sets`)

### Caching System
- **Diamond Cache**: Stores diamond computation results by structure + conditioning states
- **Computation Lookup**: Maps diamond hashes to precomputed subgraph data

## Main Algorithm: `update_beliefs_iterative`

Processes nodes in topological order:

### Source Nodes
```julia
belief_dict[node] = node_priors[node]
```

### Diamond Nodes
Uses precomputed structures:
```julia
if haskey(diamond_structures, node)
    diamond_beliefs = calculate_diamond_groups_belief(...)
end
```

### Regular Join Nodes
Applies inclusion-exclusion for multiple independent paths:
```julia
_preprior = inclusion_exclusion(all_beliefs)
belief_dict[node] = multiply_values(node_priors[node], _preprior)
```

## Diamond Processing: `updateDiamondJoin`

### Algorithm Steps:

1. **Get Precomputed Structure**:
   ```julia
   computation_data = computation_lookup[diamond_hash_key]
   ```

2. **Build Subgraph Context**:
   - Extract subgraph indices and diamond structures from precomputed data
   - Create localized probability dictionaries

3. **Conditional Enumeration**:
   ```julia
   for state_idx in 0:(2^length(conditioning_nodes)-1)
       # Calculate state probability P(conditioning_state)
       state_probability = ∏ᵢ belief[cᵢ]^bit_i × (1-belief[cᵢ])^(1-bit_i)
       
       # Check cache
       if haskey(diamond_cache, cache_key)
           state_beliefs = cached_entry.state_beliefs
       else
           # Recursive belief propagation on diamond subgraph
           state_beliefs = update_beliefs_iterative(diamond_subgraph, ...)
           diamond_cache[cache_key] = DiamondCacheEntry(...)
       end
       
       # Weight and accumulate
       final_belief += state_beliefs[join_node] × state_probability
   end
   ```

### Mathematical Foundation
For conditioning nodes C₁, C₂, ..., Cₙ:

```
Belief(J) = Σₛ P(state_s) × Belief(J | state_s)
```

Where:
- `state_s` ∈ {0,1}ⁿ represents conditioning node states
- `P(state_s) = ∏ᵢ Belief(Cᵢ)^sᵢ × (1-Belief(Cᵢ))^(1-sᵢ)`
- `Belief(J | state_s)` computed via recursive belief propagation

## Inclusion-Exclusion for Independent Paths

For multiple independent belief sources B₁, B₂, ..., Bₖ:

```
P(≥1 signal) = Σᵢ Bᵢ - Σᵢ<ⱼ BᵢBⱼ + Σᵢ<ⱼ<ₖ BᵢBⱼBₖ - ... + (-1)^(k+1) B₁B₂...Bₖ
```

Implementation:
```julia
function inclusion_exclusion(belief_values::Vector{T})
    combined_belief = zero_value(T)
    for i in 1:length(belief_values)
        for combination in combinations(belief_values, i)
            intersection_probability = prod_values(collect(combination))
            if isodd(i)
                combined_belief = add_values(combined_belief, intersection_probability)
            else
                combined_belief = subtract_values(combined_belief, intersection_probability)
            end
        end
    end
    return combined_belief
end
```

## Performance Features

### Precomputed Subgraph Structures
- Diamond subgraphs built once during preprocessing
- Avoids expensive graph operations during belief propagation
- Allows parallel diamond processing

### Multi-Level Caching
- **Structure Level**: `computation_lookup` maps diamond hashes to subgraph data
- **State Level**: `diamond_cache` stores computation results for specific conditioning states

### Type-Agnostic Operations
Unified interface for all probability types:
```julia
add_values(a, b)        # Addition
multiply_values(a, b)   # Multiplication  
complement_value(a)     # 1 - a
inclusion_exclusion(values)  # Complex probability combinations
```

## Complexity Analysis

| Structure Type | Time Complexity | Space Complexity |
|---------------|----------------|------------------|
| Tree paths | O(n) | O(n) |
| Independent joins | O(2^k × n) | O(n) |
| Single diamond | O(2^c × d) | O(2^c) |
| Nested diamonds | O(2^Σcᵢ × Πdᵢ) | O(2^max(cᵢ)) |

Where:
- n = total nodes
- k = maximum parents per join node  
- c = conditioning nodes per diamond
- d = diamond subgraph size

## Network Validation

The module validates:
- DAG structure (no cycles)
- Complete probability specifications
- Valid probability bounds [0,1] for all types
- Consistent incoming/outgoing edge indices
- Proper source node identification

## Integration with DiamondProcessing

Uses precomputed structures from DiamondProcessing:
- `DiamondComputationData`: Complete subgraph metadata
- `Diamond`: Identified diamond structures with conditioning nodes
- `DiamondsAtNode`: Diamond groupings at specific join nodes

This provides exact probabilistic inference with computational efficiency through caching and separating network structure preprocessing from analysis problem.