#
# CriticalPathV2 — mode-based rebuild of the critical path toolkit.
#
# Design: validation/CPM_REBUILD_DESIGN.md. Every shipped mode carries forward AND
# backward semantics; the sum family's backward object is the adjoint pass. Value
# types are generic; interval results are computed by scheme (corner runs), not by
# operator overloading — see Internal/IntervalScheme.jl once the interval phase lands.
#

module CriticalPath

# Public API: the modes, the two entry points (`analyze` / `accumulation_analysis`),
# the interval variants, and the result types.
export AnalysisMode, LONGEST_PATH, SHORTEST_PATH, MAX_SCALING,
       PathResult, AccumulationResult,
       analyze, accumulation_analysis,
       ValueInterval, IntervalPathResult, interval_analyze, interval_analyze_exact,
       interval_analyze_split, SplitDeclined,
       width, is_degenerate

# Secondary surface — the low-level fold primitives `analyze` is built on.
public forward_fold, reverse_fold

include(joinpath(@__DIR__, "Internal", "Modes.jl"))
include(joinpath(@__DIR__, "Internal", "Results.jl"))
include(joinpath(@__DIR__, "Internal", "Kernels.jl"))
include(joinpath(@__DIR__, "Internal", "IntervalScheme.jl"))
include(joinpath(@__DIR__, "Internal", "DominationSplit.jl"))

end # module
