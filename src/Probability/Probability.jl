"""
    Probability

Exact belief propagation algorithm for DAG networks with diamond structures.
Supports Float64, pbox, and Interval types via parametric polymorphism.

# Mathematical Foundation (see Probability_Maths.md for details)

This module computes exact beliefs using:
    P(N) = Prior(N) × P(N receives ≥1 signal from sources)

For nodes with multiple incoming paths, uses inclusion-exclusion principle:
    P(A ∪ B ∪ C) = Σᵢ P(Aᵢ) - Σᵢ<ⱼ P(Aᵢ ∩ Aⱼ) + Σᵢ<ⱼ<ₖ P(Aᵢ ∩ Aⱼ ∩ Aₖ) - ...

For diamonds (convergent path structures with conditioning nodes), uses conditional expectation:
    Result = Σ_{states} P(state) × P(Join | state)
where P(state) = ∏ᵢ [Belief(cᵢ)]^{bit_i} × [1-Belief(cᵢ)]^{1-bit_i}

For nested diamonds, implements nested conditional expectation:
    E[E[...E[Belief(Outer_Join) | Inner_Layer] ... | Outer_Layer]]

# Key Implementation Features

1. **Thread-Safe Parallelism**: Uses Threads.@spawn for parallel state enumeration
   - Each diamond state is mathematically independent
   - Thread-local copies prevent race conditions

2. **Contextual Belief Mechanism**: Critical for nested diamonds
   - Non-conditioning source nodes in diamond subgraphs receive contextual beliefs
     from outer computation
   - This makes each diamond computation dependent on outer context
   - Results in low cache hit rate but maintains exactness

3. **Recursive Processing**: Recursive call to update_beliefs_iterative
   - Each diamond spawns complete belief propagation on its subgraph
   - Nested diamonds create call stack 50+ levels deep for complex networks

4. **Optimizations**:
   - Bit-masking for state enumeration (no Combinatorics library needed)
   - Stream hashing for cache keys (no intermediate arrays)
   - Pre-computed diamond structures (computation_lookup) - O(1) retrieval
   - Diamond result caching (though low hit rate due to contextual beliefs)

# Type Support

All functions are parametric in T <: Union{Float64, pbox, Interval}.
Arithmetic uses Input helpers (multiply_values, add_values, etc.)
so that pbox/Interval arithmetic is handled transparently.
"""
module Probability

    # No Combinatorics needed — using bit-masking instead
    using ..Diamonds
    using ..Input
    using ..GraphValidation

    # Import all uncertainty operations from Input
    import ..Input: Interval, pbox, PBA,
           zero_value, one_value, non_fixed_value,
           is_valid_probability, add_values, multiply_values,
           complement_value, subtract_values, sum_values, prod_values,
           pbox_conditional_combine, PBOX_COND_BLEND

    # Public API. `update_beliefs_iterative` is the entry point; `validate_network_data`
    # checks a network before propagation.
    export update_beliefs_iterative, validate_network_data
    # Secondary surface — reusable computational primitives.
    public calculate_regular_belief, inclusion_exclusion
    # Internal (reach as `Probability.x` if ever needed): updateDiamondJoin,
    # calculate_diamond_groups_belief, DiamondCacheEntry, CacheKey, make_cache_key

    include(joinpath(@__DIR__, "Internal", "TypesAndCache.jl"))
    include(joinpath(@__DIR__, "Internal", "Validation.jl"))
    include(joinpath(@__DIR__, "Internal", "CorePropagation.jl"))
    include(joinpath(@__DIR__, "Internal", "DiamondPropagation.jl"))
    include(joinpath(@__DIR__, "Internal", "Conversions.jl"))
end
