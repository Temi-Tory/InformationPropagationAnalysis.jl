"""
    DiamondCacheEntry{T}

Stores the result of a diamond computation for potential reuse.
Parametric in T to support Float64, pbox, and Interval types.
"""
struct DiamondCacheEntry{T}
    edgelist::Vector{Tuple{Int64,Int64}}
    current_priors::Dict{Int64,T}
    state_beliefs::Dict{Int64,T}
end

# Opt-in lean cache mode. Default OFF = legacy behavior byte-identical (including the server's
# cache serialization in Handlers/AnalysisCommon.jl). When ON, updateDiamondJoin stores each entry
# as (edgelist, empty priors, Dict(join_node => belief)) instead of two full per-subgraph dicts —
# ~100x smaller per entry, results EXACT: the compute path reads cached entries only at
# [join_node], and the rare same-key/different-join collision is handled by a get()-miss guard
# that recomputes and merges (see updateDiamondJoin). Motivation: interval-mode context explosion
# (link-bnlearn: 2.05M entries = 24GB OOM with full entries vs 2,982 entries in Float64).
const LEAN_DIAMOND_CACHE = Ref(false)

"""
    CacheKey

Compact hash-based key for diamond cache lookups.
Uses hashes instead of full data structures for fast equality checks.
"""
struct CacheKey
    diamond_hash::UInt64
    priors_hash::UInt64
end

Base.hash(k::CacheKey, h::UInt) = hash((k.diamond_hash, k.priors_hash), h)
Base.:(==)(a::CacheKey, b::CacheKey) = a.diamond_hash == b.diamond_hash && a.priors_hash == b.priors_hash

"""
    make_cache_key(edgelist, current_priors) -> CacheKey

Creates cache key using stream hashing for performance.
Handles Float64, pbox, and Interval types for hashing.
"""
function make_cache_key(edgelist, current_priors::Dict{Int64, T}) where {T}
    diamond_hash = hash(sort(edgelist))

    # Stream hashing - no intermediate array allocation
    priors_hash = UInt64(0)
    sorted_nodes = sort(collect(keys(current_priors)))

    for node in sorted_nodes
        value = current_priors[node]
        if isa(value, Float64)
            priors_hash = hash((node, value), priors_hash)
        elseif isa(value, pbox)
            # For pbox, use numeric bounds for hashing
            min_val = minimum(value.u)
            max_val = maximum(value.d)
            priors_hash = hash((node, (min_val, max_val)), priors_hash)
        elseif isa(value, Interval)
            priors_hash = hash((node, (value.lower, value.upper)), priors_hash)
        else
            priors_hash = hash((node, string(value)), priors_hash)
        end
    end

    return CacheKey(diamond_hash, priors_hash)
end

# Thread-safe lock for diamond cache access in parallel execution
const diamond_cache_lock = ReentrantLock()
