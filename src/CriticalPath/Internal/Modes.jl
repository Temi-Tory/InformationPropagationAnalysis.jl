#
# ANALYSIS MODES
#
# A mode bundles the operator pair WITH its backward semantics. Only bundles whose
# backward story closes are defined here; the sum/+ family lives in the linear kernel
# (accumulation_analysis) because its backward object is an adjoint, not a residuation.
#

struct AnalysisMode
    name::Symbol
    combine::Function       # reduce incoming propagated values (maximum / minimum)
    propagate::Function     # compound a value along one edge: (value, edge_value) -> value
    apply_node::Function    # fold the node's own value in: (combined, node_value) -> value
    through::Function       # join forward and reverse completions: (F, R) -> value through node
    margin::Function        # (project_value, through) -> nonnegative margin; critical iff ~ 0
    project::Function       # reduce sink forward values to the project value
    neutral::Function       # identity element constructor: Type -> value (zero / one)
    additive_schedule::Bool # ES/LS/LF are meaningful (additive modes only)
    margin_name::Symbol
end

const LONGEST_PATH = AnalysisMode(
    :longest_path, maximum, +, +, +,
    (P, t) -> P - t,
    maximum, zero, true, :slack)

const SHORTEST_PATH = AnalysisMode(
    :shortest_path, minimum, +, +, +,
    (P, t) -> t - P,
    minimum, zero, true, :margin)

const MAX_SCALING = AnalysisMode(
    :max_scaling, maximum, *, *, *,
    (P, t) -> P / t - one(P / t),
    maximum, one, false, :ratio_slack)
