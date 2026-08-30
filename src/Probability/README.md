# ProbabilityPropagationModule

`ProbabilityPropagationModule` performs exact belief propagation for DAG networks with support for nested diamond structures and uncertainty-aware numeric types.

## Main Entry Points

- `update_beliefs_iterative`
- `validate_network_data`
- `calculate_regular_belief`
- `inclusion_exclusion`
- `updateDiamondJoin`
- `calculate_diamond_groups_belief`

## Caching Types

- `DiamondCacheEntry`
- `CacheKey`
- `make_cache_key`

## Dependencies

- `InputProcessingModule`
- `DiamondDecompositionModule`

## Internal Layout

- `ProbabilityPropagationModule.jl`
	- Public module wrapper, imports, exports, and include wiring.
- `Internal/TypesAndCache.jl`
	- `DiamondCacheEntry`, `CacheKey`, `make_cache_key`, and cache lock.
- `Internal/Validation.jl`
	- Graph/probability contract validation before propagation starts.
- `Internal/CorePropagation.jl`
	- Main iterative pass (`update_beliefs_iterative`) plus regular-path helpers.
- `Internal/DiamondPropagation.jl`
	- Diamond-state conditional expectation and nested recursion path.
- `Internal/Conversions.jl`
	- Float64-to-pbox conversion helper utilities.

## Input Contracts (High Level)

- Graph topology: `edgelist`, `incoming_index`, `outgoing_index`
- Structural metadata: `iteration_sets`, `ancestors`, `descendants`, join/fork sets
- Probabilistic data: node priors and edge probabilities for `Float64`, `Interval`, or `pbox`

## Packaging Guidance

- Keep this module pure from server concerns.
- Preserve the parametric type signatures for uncertainty support.
- Keep internal includes focused by responsibility (types/cache, validation, core pass, diamond logic, conversion helpers).
