# DiamondDecompositionModule

`DiamondDecompositionModule` detects and groups convergent path structures (diamonds) and prepares subgraph metadata used by propagation analyses.

## Main Types

- `Diamond`
- `DiamondsAtNode`
- `DiamondComputationData`

## Primary Entry Points

- `identify_and_group_diamonds`
- `build_unique_diamond_storage`
- `build_unique_diamond_storage_depth_first_parallel`
- `create_diamond_hash_key`

## Dependencies

- `InputProcessingModule`
- `ProbabilityBoundsAnalysis` (for uncertainty-type compatibility)

## Internal Layout

- `DiamondDecompositionModule.jl`
	- Public module wrapper and include wiring.
- `Internal/TypesAndCache.jl`
	- Core types (`Diamond`, `DiamondsAtNode`, `DiamondComputationData`) and cache-context helpers.
- `Internal/UtilityFunctions.jl`
	- Hash-key, iteration helper functions, and module-local typed defaults.
- `Internal/Pipeline.jl`
	- Detection and precomputation pipeline functions extracted from the former monolith body.

## Notes for Packaging

- This module is stateful only through local computation contexts; no global mutable pipeline state is required.
- Keep dependency direction one-way: this module can depend on `Shared/`, but higher modules should not reimplement diamond logic.
