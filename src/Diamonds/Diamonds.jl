#= The diamond detection automatically identifies:

Common failure points (convergent subsystems)
Natural organizational boundaries
Places where delays propagate through multiple paths =#
module Diamonds

    using ..Input 
    import ProbabilityBoundsAnalysis
    
    # Create aliases to avoid ambiguity
    const PBA = ProbabilityBoundsAnalysis
    # Type aliases for convenience
    const PBAInterval = ProbabilityBoundsAnalysis.Interval
    const pbox = ProbabilityBoundsAnalysis.pbox
    const Interval = Input.Interval

    # Public API. `new_identify` is the sole diamond producer (factorized,
    # correct-by-construction); `Diamond` / `DiamondsAtNode` / `DiamondComputationData`
    # are its result types.
    export Diamond, DiamondsAtNode, DiamondComputationData, new_identify
    # Secondary surface — importable and documented, not pulled in by `using`.
    public create_diamond_hash_key
    # RETIRED (buggy hybrid-reuse + completeness loop): identify_and_group_diamonds,
    # build_unique_diamond_storage[_depth_first_parallel] from Pipeline*.jl. Replaced by new_identify,
    # which emits root_diamonds + unique_diamonds together. See ROADMAP.md / PIPELINE_REWRITE_STATUS.md.

    include(joinpath(@__DIR__, "Internal", "TypesAndCache.jl"))
    include(joinpath(@__DIR__, "Internal", "UtilityFunctions.jl"))
    include(joinpath(@__DIR__, "Internal", "NewIdentify.jl"))

end
