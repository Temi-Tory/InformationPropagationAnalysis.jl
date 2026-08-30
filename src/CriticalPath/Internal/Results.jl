#
# RESULT TYPES
#

"""
Result of a path-mode analysis (LongestPath / ShortestPath / MaxScaling).
`method` records how the numbers were computed (:exact_scalar for Float64;
interval methods carry their own tags so exact and conservative can never be confused).
Schedule fields (early_start / late_finish / late_start) are populated for additive
modes only and empty otherwise.
"""
struct PathResult{T}
    mode::Symbol
    method::Symbol
    forward::Dict{Int64,T}             # F: best value reaching each node (inclusive of the node)
    reverse_completion::Dict{Int64,T}  # R: best completion from each node to a sink (exclusive)
    project_value::T
    through::Dict{Int64,T}             # best value of a complete path constrained through the node
    margin::Dict{Int64,T}
    margin_name::Symbol
    critical::Vector{Int64}
    early_start::Dict{Int64,T}
    late_finish::Dict{Int64,T}
    late_start::Dict{Int64,T}
end

"""
Result of the linear (sum-family) analysis. The backward object is the adjoint:
multiplicity = number of directed paths from the node to the target, sensitivity
= d(total)/d(node value), contribution = node value x multiplicity. Allowance is
present only when a budget was supplied.
"""
struct AccumulationResult{T}
    mode::Symbol
    method::Symbol
    forward::Dict{Int64,T}
    target::Int64
    total::T
    multiplicity::Dict{Int64,Int64}
    sensitivity::Dict{Int64,T}
    contribution::Dict{Int64,T}
    allowance::Dict{Int64,T}
    ranking::Vector{Int64}             # nodes by contribution, largest first
end
