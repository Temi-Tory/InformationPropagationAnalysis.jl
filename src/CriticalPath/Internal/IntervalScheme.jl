#
# INTERVAL SCHEME
#
# Interval results are computed by SCHEME, never by running the scalar algorithm on an
# overloaded interval type (the failure mode of the old module: partial-order minimum
# selection and dependency widening).
#
#   Forward quantities (F, R, through, P) are monotone in every input, so two crisp
#   corner runs give EXACT interval bounds (method :exact_corners).
#
#   Margins are differences of dependent monotone quantities, so they get either a
#   sound conservative enclosure (tier 1, method :conservative_enclosure) or exact
#   bounds by corner enumeration over the interval-valued inputs (tier 2, method
#   :exact_corners_exhaustive), together with the necessary / possible criticality
#   classification of the interval-PERT literature.
#

struct ValueInterval
    lo::Float64
    hi::Float64
    function ValueInterval(lo::Float64, hi::Float64)
        lo <= hi || throw(ArgumentError("interval bounds out of order: [$lo, $hi]"))
        new(lo, hi)
    end
end

width(v::ValueInterval) = v.hi - v.lo
is_degenerate(v::ValueInterval) = v.lo == v.hi

struct IntervalPathResult
    mode::Symbol
    method::Symbol                       # margin/criticality method actually used
    forward::Dict{Int64,ValueInterval}   # exact (corner pair)
    through::Dict{Int64,ValueInterval}   # exact (corner pair)
    project_value::ValueInterval         # exact (corner pair)
    margin::Dict{Int64,ValueInterval}    # tagged by `method`
    margin_name::Symbol
    necessarily_critical::Vector{Int64}
    possibly_critical::Vector{Int64}
    corner_count::Int                    # corners enumerated (0 for tier 1)
end

lo_values(d::Dict{Int64,ValueInterval}) = Dict{Int64,Float64}(k => v.lo for (k, v) in d)
hi_values(d::Dict{Int64,ValueInterval}) = Dict{Int64,Float64}(k => v.hi for (k, v) in d)
lo_values(d::Dict{Tuple{Int64,Int64},ValueInterval}) = Dict{Tuple{Int64,Int64},Float64}(k => v.lo for (k, v) in d)
hi_values(d::Dict{Tuple{Int64,Int64},ValueInterval}) = Dict{Tuple{Int64,Int64},Float64}(k => v.hi for (k, v) in d)

"""
Tier 1: exact forward bounds by corner pair, conservative margin enclosure, sound
criticality classification (necessarily-critical certificates are certain; the
possibly-critical list is a superset that tier 2 can tighten).
"""
function interval_analyze(
    iteration_sets::Vector{Set{Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    node_values::Dict{Int64,ValueInterval},
    edge_values::Dict{Tuple{Int64,Int64},ValueInterval};
    mode::AnalysisMode = LONGEST_PATH,
    atol::Float64 = 1e-9
)::IntervalPathResult
    rlo = analyze(iteration_sets, outgoing_index, incoming_index, source_nodes,
                  lo_values(node_values), lo_values(edge_values); mode = mode, atol = atol)
    rhi = analyze(iteration_sets, outgoing_index, incoming_index, source_nodes,
                  hi_values(node_values), hi_values(edge_values); mode = mode, atol = atol)

    F = Dict{Int64,ValueInterval}(n => ValueInterval(rlo.forward[n], rhi.forward[n]) for n in keys(rlo.forward))
    through = Dict{Int64,ValueInterval}(n => ValueInterval(rlo.through[n], rhi.through[n]) for n in keys(rlo.through))
    P = ValueInterval(rlo.project_value, rhi.project_value)

    margin = Dict{Int64,ValueInterval}()
    for (n, t) in through
        mlo = max(0.0, mode.margin(P.lo, t.hi))
        mhi = mode.margin(P.hi, t.lo)
        margin[n] = ValueInterval(min(mlo, mhi), max(mlo, mhi))
    end
    necessarily = sort!([n for (n, m) in margin if m.hi <= atol])
    possibly = sort!([n for (n, m) in margin if m.lo <= atol])

    return IntervalPathResult(mode.name, :conservative_enclosure, F, through, P, margin,
                              mode.margin_name, necessarily, possibly, 0)
end

"""
Tier 2: exact margins and criticality by enumerating every corner assignment of the
interval-valued inputs (degenerate intervals stay fixed). Margin extremes are attained
at endpoint configurations, so the result is exact. Cost 2^k for k interval inputs;
refuses beyond `max_corners` unless forced. The forward/through/project fields keep
their exact corner-pair values.
"""
function interval_analyze_exact(
    iteration_sets::Vector{Set{Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}},
    source_nodes::Set{Int64},
    node_values::Dict{Int64,ValueInterval},
    edge_values::Dict{Tuple{Int64,Int64},ValueInterval};
    mode::AnalysisMode = LONGEST_PATH,
    atol::Float64 = 1e-9,
    max_corners::Int = 1 << 16
)::IntervalPathResult
    base = interval_analyze(iteration_sets, outgoing_index, incoming_index, source_nodes,
                            node_values, edge_values; mode = mode, atol = atol)

    var_nodes = sort!([k for (k, v) in node_values if !is_degenerate(v)])
    var_edges = sort!([k for (k, v) in edge_values if !is_degenerate(v)])
    k = length(var_nodes) + length(var_edges)
    ncorners = 1 << k
    ncorners <= max_corners ||
        throw(ArgumentError("2^$k corners exceeds max_corners=$max_corners; raise the cap or use the diamond-guided driver"))

    d = lo_values(node_values)
    w = lo_values(edge_values)
    mlo = Dict{Int64,Float64}(); mhi = Dict{Int64,Float64}()
    ncrit = Dict{Int64,Bool}(); pcrit = Dict{Int64,Bool}()
    first_corner = true
    for mask in 0:(ncorners - 1)
        for (i, n) in enumerate(var_nodes)
            d[n] = (mask >> (i - 1)) & 1 == 1 ? node_values[n].hi : node_values[n].lo
        end
        off = length(var_nodes)
        for (i, e) in enumerate(var_edges)
            w[e] = (mask >> (off + i - 1)) & 1 == 1 ? edge_values[e].hi : edge_values[e].lo
        end
        r = analyze(iteration_sets, outgoing_index, incoming_index, source_nodes, d, w;
                    mode = mode, atol = atol)
        critset = Set(r.critical)
        for (n, m) in r.margin
            if first_corner
                mlo[n] = m; mhi[n] = m
                ncrit[n] = n in critset; pcrit[n] = n in critset
            else
                mlo[n] = min(mlo[n], m); mhi[n] = max(mhi[n], m)
                ncrit[n] &= n in critset; pcrit[n] |= n in critset
            end
        end
        first_corner = false
    end

    # Precautionary, not a fix for a live bug here the way the split's own
    # snap is (see DominationSplit.jl): mlo/mhi track min/max of the SAME
    # per-corner value stream in the SAME loop, so mlo[n] <= mhi[n] holds by
    # construction — there is no second, independently-summed sweep for
    # floating-point noise to diverge against. Snapping near-zero margins to
    # exactly 0 anyway keeps this tier consistent with the split's and with
    # the atol the criticality classification below already uses.
    for n in keys(mlo)
        abs(mlo[n]) <= atol && (mlo[n] = 0.0)
        abs(mhi[n]) <= atol && (mhi[n] = 0.0)
    end
    margin = Dict{Int64,ValueInterval}(n => ValueInterval(mlo[n], mhi[n]) for n in keys(mlo))
    necessarily = sort!([n for (n, b) in ncrit if b])
    possibly = sort!([n for (n, b) in pcrit if b])

    return IntervalPathResult(mode.name, :exact_corners_exhaustive, base.forward, base.through,
                              base.project_value, margin, mode.margin_name,
                              necessarily, possibly, ncorners)
end
