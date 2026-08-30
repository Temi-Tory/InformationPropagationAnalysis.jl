module InformationPropagationAnalysis

# ─────────────────────────────────────────────────────────────────────────────
# Shared infrastructure (internal — not re-exported)
# ─────────────────────────────────────────────────────────────────────────────
include(joinpath(@__DIR__, "Input", "Input.jl"))
include(joinpath(@__DIR__, "Shared", "GraphValidation.jl"))
include(joinpath(@__DIR__, "Shared", "GraphTraversal.jl"))

# ─────────────────────────────────────────────────────────────────────────────
# Toolkits
# ─────────────────────────────────────────────────────────────────────────────
include(joinpath(@__DIR__, "Diamonds", "Diamonds.jl"))
include(joinpath(@__DIR__, "CriticalPath", "CriticalPath.jl"))
include(joinpath(@__DIR__, "Probability", "Probability.jl"))
include(joinpath(@__DIR__, "Flow", "Flow.jl"))

# ─────────────────────────────────────────────────────────────────────────────
# Public API
#
# The five toolkits are the surface — `IPA.Input`, `IPA.Diamonds`,
# `IPA.Probability`, `IPA.CriticalPath`, `IPA.Flow` — each `using`-able on its own
# (e.g. `using InformationPropagationAnalysis.Flow`). At the top level we
# re-export the toolkit names, the two universal value types, and one blessed
# entry point per toolkit; everything else is reached through its toolkit.
# ─────────────────────────────────────────────────────────────────────────────
using .Input: Interval, pbox, PBA
using .Diamonds: new_identify
using .Probability: update_beliefs_iterative
using .Flow: analyze_all

"""
    critical_path(edgelist, source_values; mode = CriticalPath.LONGEST_PATH, kwargs...)

Top-level alias for `CriticalPath.analyze`.
"""
const critical_path = CriticalPath.analyze

export Input, Diamonds, Probability, CriticalPath, Flow
export Interval, pbox
export new_identify, update_beliefs_iterative, critical_path, analyze_all

end # module InformationPropagationAnalysis
