#
# KERNELS
#
# Idempotent kernel: one forward fold and one reverse fold over the layered iteration
# sets, both driven entirely by the mode's operator bundle. The reverse fold is the
# forward fold on the reversed graph with the same operators, which is what makes the
# backward pass exist for every shipped mode by construction.
#
# Linear kernel: forward accumulation plus the adjoint (path-multiplicity) pass.
#

function forward_fold(
    iteration_sets::Vector{Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    node_values::Dict{Int64,T},
    edge_values::Dict{Tuple{Int64,Int64},T},
    mode::AnalysisMode,
    initial::T
)::Dict{Int64,T} where T
    F = Dict{Int64,T}()
    neu = mode.neutral(T)
    for layer in iteration_sets
        for node in layer
            dv = get(node_values, node, neu)
            if node in source_nodes
                F[node] = mode.apply_node(initial, dv)
            else
                parents = incoming_index[node]
                vals = T[mode.propagate(F[p], get(edge_values, (p, node), neu)) for p in parents]
                F[node] = mode.apply_node(mode.combine(vals), dv)
            end
        end
    end
    return F
end

function reverse_fold(
    iteration_sets::Vector{Set{Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    node_values::Dict{Int64,T},
    edge_values::Dict{Tuple{Int64,Int64},T},
    mode::AnalysisMode
)::Dict{Int64,T} where T
    R = Dict{Int64,T}()
    neu = mode.neutral(T)
    for i in length(iteration_sets):-1:1
        for node in iteration_sets[i]
            succs = get(outgoing_index, node, Set{Int64}())
            if isempty(succs)
                R[node] = neu
            else
                vals = T[mode.propagate(mode.apply_node(R[s], get(node_values, s, neu)),
                                        get(edge_values, (node, s), neu)) for s in succs]
                R[node] = mode.combine(vals)
            end
        end
    end
    return R
end

"""
Full analysis for one path mode: forward and reverse folds, project value, per-node
through-values and margins, critical set, and (for additive modes) the classical
schedule quantities.
"""
function analyze(
    iteration_sets::Vector{Set{Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    node_values::Dict{Int64,T},
    edge_values::Dict{Tuple{Int64,Int64},T};
    mode::AnalysisMode = LONGEST_PATH,
    initial = nothing,
    atol::Float64 = 1e-9
)::PathResult{T} where T
    init = initial === nothing ? mode.neutral(T) : convert(T, initial)
    F = forward_fold(iteration_sets, incoming_index, source_nodes, node_values, edge_values, mode, init)
    R = reverse_fold(iteration_sets, outgoing_index, node_values, edge_values, mode)

    sinks = Int64[n for n in keys(F) if isempty(get(outgoing_index, n, Set{Int64}()))]
    P = mode.project(T[F[s] for s in sinks])

    through = Dict{Int64,T}(n => mode.through(F[n], R[n]) for n in keys(F))
    margin = Dict{Int64,T}(n => mode.margin(P, through[n]) for n in keys(F))
    critical = sort!(Int64[n for (n, m) in margin if isapprox(m, zero(T); atol = atol)])

    es = Dict{Int64,T}(); lf = Dict{Int64,T}(); ls = Dict{Int64,T}()
    if mode.additive_schedule
        neu = mode.neutral(T)
        for n in keys(F)
            dv = get(node_values, n, neu)
            es[n] = F[n] - dv
            lf[n] = P - R[n]
            ls[n] = lf[n] - dv
        end
    end

    return PathResult{T}(mode.name, :exact_scalar, F, R, P, through, margin,
                         mode.margin_name, critical, es, lf, ls)
end

"""
Linear (sum-family) analysis. Forward: every node accumulates the sum of its parents'
accumulated values plus edge values plus its own value, so shared ancestry is counted
once per path (path-multiplicity semantics). Backward: the adjoint pass, computing the
number of directed paths from each node to the target, which IS d(total)/d(node value).
"""
function accumulation_analysis(
    iteration_sets::Vector{Set{Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    node_values::Dict{Int64,T},
    edge_values::Dict{Tuple{Int64,Int64},T};
    target::Union{Int64,Nothing} = nothing,
    budget::Union{T,Nothing} = nothing,
    initial = nothing
)::AccumulationResult{T} where T
    init = initial === nothing ? zero(T) : convert(T, initial)
    F = Dict{Int64,T}()
    for layer in iteration_sets
        for node in layer
            dv = get(node_values, node, zero(T))
            if node in source_nodes
                F[node] = init + dv
            else
                acc = zero(T)
                for p in incoming_index[node]
                    acc += F[p] + get(edge_values, (p, node), zero(T))
                end
                F[node] = acc + dv
            end
        end
    end

    t = target
    if t === nothing
        sinks = Int64[n for n in keys(F) if isempty(get(outgoing_index, n, Set{Int64}()))]
        best = sinks[1]
        for s in sinks
            if F[s] > F[best] || (F[s] == F[best] && s < best)
                best = s
            end
        end
        t = best
    end

    m = Dict{Int64,Int64}(t => 1)
    for i in length(iteration_sets):-1:1
        for node in iteration_sets[i]
            node == t && continue
            cnt = 0
            for s in get(outgoing_index, node, Set{Int64}())
                cnt += get(m, s, 0)
            end
            m[node] = cnt
        end
    end

    sensitivity = Dict{Int64,T}(n => convert(T, c) for (n, c) in m)
    contribution = Dict{Int64,T}(n => get(node_values, n, zero(T)) * convert(T, c) for (n, c) in m)

    allowance = Dict{Int64,T}()
    if budget !== nothing
        head = budget - F[t]
        for (n, c) in m
            c > 0 && (allowance[n] = head / convert(T, c))
        end
    end

    ranking = sort!(Int64[n for (n, c) in m if c > 0]; by = n -> (-contribution[n], n))

    return AccumulationResult{T}(:accumulation, :exact_scalar, F, t, F[t], m,
                                 sensitivity, contribution, allowance, ranking)
end
