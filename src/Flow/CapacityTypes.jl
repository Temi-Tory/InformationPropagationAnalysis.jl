# CapacityTypes.jl
# Concrete record types for Flow results.
# Included directly by modules that need these types.
# Canonical Julia module wrapper.

module CapacityTypes

export CriticalEdgeRecord,
       SingleEdgeFailureRecord,
       KEdgeFailureRecord,
       DegradationScenarioRecord,
       PathFlowContribution,
       BottleneckRecord

"""
    CriticalEdgeRecord

Concrete critical-edge ranking record produced by `SensitivityModule`
(`critical_edge_ranking` and `analyze_sensitivity`).
"""
struct CriticalEdgeRecord
    edge           :: Tuple{Int64,Int64}
    baseline_flow  :: Float64
    perturbed_flow :: Float64
    drop           :: Float64
end

"""
    SingleEdgeFailureRecord

Concrete single-edge failure impact record produced by `FailureImpactModule`
(`analyze_single_edge_failures` and `analyze_failure_impact`).
"""
struct SingleEdgeFailureRecord
    edge           :: Tuple{Int64,Int64}
    baseline_flow  :: Float64
    perturbed_flow :: Float64
    drop           :: Float64
    is_critical    :: Bool
    is_unbounded   :: Bool
end

"""
    KEdgeFailureRecord

Concrete k-edge failure impact record produced by `FailureImpactModule`
(`analyze_k_edge_failures` and `analyze_failure_impact`).
"""
struct KEdgeFailureRecord
    edges          :: Tuple
    baseline_flow  :: Float64
    perturbed_flow :: Float64
    drop           :: Float64
    is_unbounded   :: Bool
end

"""
    DegradationScenarioRecord

Concrete capacity-degradation scenario record produced by `FailureImpactModule`
(`analyze_capacity_degradation` and `analyze_failure_impact`).
"""
struct DegradationScenarioRecord
    scenario_id         :: Int64
    scenario_capacities :: Dict{Tuple{Int64,Int64},Float64}
    max_flow            :: Float64
    sink_flow           :: Dict{Int64,Float64}
    saturated_edges     :: Vector{Tuple{Int64,Int64}}
    drop_from_baseline  :: Float64
    is_unbounded        :: Bool
end

"""
    PathFlowContribution

Concrete per-path flow contribution record produced by `StructuralModule`
(`path_flow_contributions` and `analyze_structure`).
"""
struct PathFlowContribution
    path              :: Vector{Int64}
    flow_contribution :: Float64
    bottleneck_edge   :: Tuple{Int64,Int64}
end

"""
    BottleneckRecord

Concrete bottleneck ranking record produced by `StructuralModule`
(`bottleneck_ranking` and `analyze_structure`).
"""
struct BottleneckRecord
    edge              :: Tuple{Int64,Int64}
    capacity          :: Float64
    flow              :: Float64
    residual_capacity :: Float64
    rank              :: Int64
end

end # module CapacityTypes
