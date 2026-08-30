module Input
    # Dependency audit (port to InformationPropagationAnalysis.jl): the only third-party
    # symbols this module actually touches are JSON.parsefile, DelimitedFiles.readdlm, and
    # DataStructures' Queue/enqueue!/dequeue!. Random, DataFrames, Distributions,
    # SparseArrays and Combinatorics were imported but unused — dropped.
    using DelimitedFiles, JSON, DataStructures

    # Import ProbabilityBoundsAnalysis for pbox construction
    import ProbabilityBoundsAnalysis
    const PBA = ProbabilityBoundsAnalysis
    const pbox = ProbabilityBoundsAnalysis.pbox

    export Interval, pbox, PBA,
           # Value types + the generic uncertainty layer the toolkits build on
           # (kept exported: `Diamonds` extends `zero_value` / `one_value` via `using`)
           zero_value, one_value, non_fixed_value, is_valid_probability,
           add_values, multiply_values, min_values, max_values, sum_values,
           complement_value, subtract_values, prod_values, divide_values,
           pbox_conditional_combine, PBOX_COND_BLEND,
           # Graph structure
           read_graph_to_dict, read_complete_network,
           identify_fork_and_join_nodes, find_iteration_sets,
           # Priors / edge probabilities — generic (auto-detect type)
           read_node_priors_from_json, read_edge_probabilities_from_json,
           # Priors / edge probabilities — type-specific (guaranteed return types)
           read_node_priors_from_json_pbox, read_edge_probabilities_from_json_pbox,
           read_node_priors_from_json_interval, read_edge_probabilities_from_json_interval,
           read_node_priors_from_json_float64, read_edge_probabilities_from_json_float64,
           # Capacities
           read_edge_capacities_from_json, read_node_capacities_from_json, read_capacities_input

    # Basic interval type
    struct Interval
        lower::Float64
        upper::Float64
        
        function Interval(lower::Float64, upper::Float64)
            if lower > upper
                throw(ArgumentError("Lower bound must be ≤ upper bound"))
            end
            new(lower, upper)
        end
    end

     # Min operation (crucial for capacity constraints)
    min_values(a::Float64, b::Float64) = min(a, b)
    min_values(a::Interval, b::Interval) = Interval(min(a.lower,
    b.lower), min(a.upper, b.upper))
    min_values(a::pbox, b::pbox) = PBA.convIndep(a, b, op = min)

    # Max operation (for flow aggregation)
    max_values(a::Float64, b::Float64) = max(a, b)
    max_values(a::Interval, b::Interval) = Interval(max(a.lower,
    b.lower), max(a.upper, b.upper))
    max_values(a::pbox, b::pbox) = PBA.convIndep(a, b, op = max)

    # Sum over collections
    sum_values(values::Vector{T}) where T = reduce(add_values, values;      
    init=zero_value(T))

    # Helper functions for type-specific operations
    # Zero and one values for different types
    zero_value(::Type{Float64}) = 0.0
    one_value(::Type{Float64}) = 1.0    
    non_fixed_value(::Type{Float64}) = 0.9
    zero_value(::Type{Interval}) = Interval(0.0, 0.0)
    one_value(::Type{Interval}) = Interval(1.0, 1.0)    
    non_fixed_value(::Type{Interval}) = Interval(0.9, 0.9)  
    zero_value(::Type{pbox}) = PBA.makepbox(PBA.interval(0.0, 0.0))
    one_value(::Type{pbox}) = PBA.makepbox(PBA.interval(1.0, 1.0))   
    non_fixed_value(::Type{pbox}) = PBA.makepbox(PBA.interval(0.9, 0.9)) 
    

    # Type-specific probability validation
    is_valid_probability(value::Float64) = 0.0 <= value <= 1.0

    function is_valid_probability(value::Interval)
        return value.lower >= 0.0 && value.upper <= 1.0
    end

    function is_valid_probability(value::pbox)
        # PBA.minimum(value)/PBA.maximum(value) do NOT reliably return a scalar or an
        # IntervalArithmetic.Interval to bound against -- for an IMPRECISE pbox (built from
        # interval-valued parameters, e.g. a parametric_interval construction), both return
        # ANOTHER pbox object. The old isa(..., PBA.Interval) check then fails silently (it's
        # a pbox, not a PBA.Interval), min_bound/max_bound fall through as the raw pbox, and
        # `pbox >= 0.0` returns an IntervalArithmetic.Interval{Float64} (an interval-valued
        # truth value), not a Bool -- which throws downstream wherever the caller does `!` or
        # `if` on this function's result (confirmed: MethodError: no method matching
        # !(::IntervalArithmetic.Interval{Float64}), from Validation.jl's edge/node prior
        # checks, on any non-degenerate pbox link probability). A precise pbox (built from a
        # plain scalar or interval, not a parametric_interval) never hit this, which is why it
        # went unnoticed until a genuinely imprecise pbox was propagated.
        #
        # Fixed: read the true support bounds directly off the pbox's own discretized bounding
        # CDFs (`.u` upper bound, `.d` lower bound, both plain Float64 arrays regardless of
        # whether the pbox is precise or imprecise) -- minimum(value.u) is the leftmost point
        # the support can reach, maximum(value.d) the rightmost, matching the bounds the pbox's
        # own `range=[...]` summary reports. No PBA.minimum/maximum call needed at all.
        min_bound = minimum(value.u)
        max_bound = maximum(value.d)

        return min_bound >= 0.0 && max_bound <= 1.0
    end

    # Type-specific arithmetic operations
    # Addition
    add_values(a::Float64, b::Float64) = a + b
    add_values(a::Interval, b::Interval) = Interval(a.lower + b.lower, a.upper + b.upper)
    add_values(a::pbox, b::pbox) = PBA.convIndep(a, b, op = +)

    # Multiplication
    multiply_values(a::Float64, b::Float64) = a * b
    function multiply_values(a::Interval, b::Interval)
        products = [a.lower * b.lower, a.lower * b.upper, a.upper * b.lower, a.upper * b.upper]
        return Interval(minimum(products), maximum(products))
    end
    multiply_values(a::pbox, b::pbox) = PBA.convIndep(a, b, op = *)

    # Complement (1 - value)
    complement_value(a::Float64) = 1.0 - a
    complement_value(a::Interval) = Interval(1.0 - a.upper, 1.0 - a.lower)
    complement_value(a::pbox) = PBA.convIndep(one_value(pbox), a, op = -)

    # Subtraction
    subtract_values(a::Float64, b::Float64) = a - b
    subtract_values(a::Interval, b::Interval) = Interval(a.lower - b.upper, a.upper - b.lower)
    subtract_values(a::pbox, b::pbox) = PBA.convIndep(a, b, op = -)

    # Division
    divide_values(a::Float64, b::Float64) = a / b
    function divide_values(a::Interval, b::Interval)
        if b.lower <= 0.0 && b.upper >= 0.0
            throw(ArgumentError("Division by interval containing zero"))
        end
        quotients = [a.lower / b.lower, a.lower / b.upper, a.upper / b.lower, a.upper / b.upper]
        return Interval(minimum(quotients), maximum(quotients))
    end
    divide_values(a::pbox, b::pbox) = PBA.convIndep(a, b, op = /)

    # ---- Conditioning recombination (diamond total-probability step) ----------------------------------
    # belief = W*A + (1-W)*B for ONE conditioning node (weight W = its contextual belief) is a CONVEX
    # COMBINATION, NOT a convolution. Float64: exact scalar. pbox: the OLD convIndep weighted-sum treated
    # this as a convolution -> over-wide, mass>1, UNSOUND. Correct operator (cvxP): integrate over W's own
    # distribution (mixture over its discretisation levels), blend the two branches with a POSITIVE-
    # DEPENDENCE bound env(convIndep,convPerfect) (branches are monotone-increasing => positively
    # dependent), and envelope over W's imprecision (u/d). Validated sound vs Monte Carlo across the corpus
    # (validation/rc_pbox_cvx.jl, cvx_sound.jl: 20/20). Set PBOX_COND_BLEND[]=:frechet for the guaranteed-
    # sound (but conservative) convFrechet branch blend. Interval uses corner enumeration in
    # updateDiamondJoin (exact), so it is NOT routed through here.
    const PBOX_COND_BLEND = Ref(:positive)   # :positive (cvxP, tight) | :frechet (cvxF, guaranteed-sound)
    _pbox_branch_blend(x::pbox, y::pbox) = PBOX_COND_BLEND[] == :frechet ?
        PBA.convFrechet(x, y, op = +) :
        PBA.env(PBA.convIndep(x, y, op = +), PBA.convPerfect(x, y, op = +))
    # belief is a PROBABILITY (in [0,1] by construction); PBA's quantile arithmetic uses plain Float64
    # (no directed rounding) so ULP-level leakage past 0/1 is possible. imp() with the [0,1] box is a SOUND
    # projection (true belief in [0,1] always) that can only tighten, never widen, the result. Built fresh
    # per call (not module-level) so its discretisation always matches the CURRENT PBA.setSteps() level —
    # a cached const captured the step count active at module load time, mismatching later setSteps calls.
    #
    # GUARDED, not silent: a clamp that swallows an excursion of ANY size would risk masking a REAL
    # conditioning bug (the old convIndep-as-convolution bug leaked mass by up to 0.34 -- see
    # pbox-conditioning-unsound memory) behind an innocuous-looking "SOUND, [0,1]" result. Genuine Float64
    # rounding leakage from ~steps sequential ops is of order steps*eps (<<1e-9 even at steps=800); an
    # excursion above EXCURSION_TOL is therefore a correctness bug in the operator, not FP noise, and must
    # fail loudly instead of being trimmed away.
    const EXCURSION_TOL = 1e-6

    # Unit-box [0,1] clamp target, cached KEYED BY the operand's discretisation length (== steps
    # at its construction). Rebuilding per call was the safe-but-slow fix for the original
    # module-load-const DimensionMismatch bug; keying by length keeps correctness under any
    # later setSteps() while removing per-combine construction cost. Exactness-identical: the
    # cached box has the same u=zeros/d=ones content makepbox(interval(0,1)) produces at that
    # steps count. (User-approved performance repair, 2026-08-16.)
    const _UNIT_BOX_N = Ref(0)
    const _UNIT_BOX = Ref{Any}(nothing)
    function _unit_box(n::Int)
        _UNIT_BOX_N[] == n && return _UNIT_BOX[]::pbox
        ub = PBA.makepbox(PBA.interval(0.0, 1.0))
        length(ub.u) == n || (ub = PBA.pbox(fill(0.0, n), fill(1.0, n)))
        _UNIT_BOX_N[] = n; _UNIT_BOX[] = ub
        ub
    end

    function pbox_conditional_combine(W::pbox, A::pbox, B::pbox)
        n = length(W.u); ps = fill(1.0 / n, n)
        Mu = PBA.mixture([_pbox_branch_blend(W.u[i] * A, (1.0 - W.u[i]) * B) for i in 1:n], ps)
        Md = PBA.mixture([_pbox_branch_blend(W.d[i] * A, (1.0 - W.d[i]) * B) for i in 1:n], ps)
        raw = PBA.env(Mu, Md)
        excursion = max(maximum(raw.u), maximum(raw.d)) - 1.0
        excursion = max(excursion, -min(minimum(raw.u), minimum(raw.d)))
        excursion > EXCURSION_TOL && error(
            "pbox_conditional_combine: [0,1] excursion of $excursion exceeds the FP-noise tolerance " *
            "($EXCURSION_TOL) -- this indicates a REAL soundness bug in the conditioning operator " *
            "(cf. the pre-cvxP convIndep bug, which leaked up to 0.34), not floating-point rounding. " *
            "Refusing to silently clamp it away; investigate before trusting this result.")
        return PBA.imp(raw, _unit_box(length(raw.u)))
    end
    pbox_conditional_combine(W::Float64, A::Float64, B::Float64) = W * A + (1.0 - W) * B

    # Sum of vector
    sum_values(values::Vector{Float64}) = sum(values)
    function sum_values(values::Vector{Interval})
        if isempty(values)
            return zero_value(Interval)
        end
        result = values[1]
        for i in 2:length(values)
            result = add_values(result, values[i])
        end
        return result
    end
    function sum_values(values::Vector{pbox})
        if isempty(values)
            return zero_value(pbox)
        end
        result = values[1]
        for i in 2:length(values)
            result = add_values(result, values[i])
        end
        return result
    end

    # Product of vector
    prod_values(values::Vector{Float64}) = prod(values)
    function prod_values(values::Vector{Interval})
        if isempty(values)
            return one_value(Interval)
        end
        result = values[1]
        for i in 2:length(values)
            result = multiply_values(result, values[i])
        end
        return result
    end
    function prod_values(values::Vector{pbox})
        if isempty(values)
            return one_value(pbox)
        end
        result = values[1]
        for i in 2:length(values)
            result = multiply_values(result, values[i])
        end
        return result
    end
    
    # Constructor for creating interval from single value (deterministic case)
    Interval(value::Float64) = Interval(value, value)

    # Base overloads for Interval — required by GeneralizedCriticalPathModule and CapacityAnalysisModule
    Base.zero(::Type{Interval}) = Interval(0.0, 0.0)
    Base.typemax(::Type{Interval}) = Interval(typemax(Float64), typemax(Float64))
    Base.minimum(v::AbstractVector{Interval}) = Interval(minimum(x.lower for x in v), minimum(x.upper for x in v))
    Base.maximum(v::AbstractVector{Interval}) = Interval(maximum(x.lower for x in v), maximum(x.upper for x in v))
    Base.:(>)(a::Interval, b::Interval) = a.lower > b.upper
    Base.:(<)(a::Interval, b::Interval) = a.upper < b.lower
    Base.:(+)(a::Interval, b::Interval) = Interval(a.lower + b.lower, a.upper + b.upper)
    Base.:(-)(a::Interval, b::Interval) = Interval(a.lower - b.upper, a.upper - b.lower)
    Base.:(*)(a::Interval, b::Interval) = begin
        products = (a.lower * b.lower, a.lower * b.upper, a.upper * b.lower, a.upper * b.upper)
        Interval(min(products...), max(products...))
    end
    Base.:(/)(a::Interval, b::Real) = Interval(a.lower / b, a.upper / b)
    Base.:(/)(a::Interval, b::Interval) = begin
        (b.lower <= 0.0 && b.upper >= 0.0) && throw(ArgumentError("Division by interval containing zero"))
        quotients = (a.lower / b.lower, a.lower / b.upper, a.upper / b.lower, a.upper / b.upper)
        Interval(min(quotients...), max(quotients...))
    end

    """
        read_graph_to_dict(filename::String)

        Reads a graph from either CSV adjacency matrix (integers 0/1) or edge list file (.edge).
        Auto-detects format based on file extension and content.
    """  
    function read_graph_to_dict(filename::String)::Tuple{Vector{Tuple{Int64,Int64}}, Dict{Int64,Set{Int64}}, Dict{Int64,Set{Int64}}, Set{Int64}}
        isfile(filename) || throw(SystemError("File not found: $filename"))
        
        # Determine file type by extension or content
        if endswith(filename, ".edge") || endswith(filename, ".EDGES")
            return read_graph_from_edgelist(filename)
        else
            # Try to detect if it's an edge list by reading first few lines
            try
                lines = readlines(filename)
                if length(lines) >= 2
                    header = strip(lines[1])
                    # Check if first line looks like edge list header
                    if occursin("source", lowercase(header)) && occursin("destination", lowercase(header))
                        return read_graph_from_edgelist(filename)
                    end
                    
                    # Check if second line has comma-separated integers (edge list format)
                    second_line = strip(lines[2])
                    if occursin(',', second_line) && !occursin(' ', second_line)
                        parts = split(second_line, ',')
                        if length(parts) == 2
                            try
                                parse(Int, parts[1])
                                parse(Int, parts[2])
                                return read_graph_from_edgelist(filename)
                            catch
                                # Not edge list format, continue to adjacency matrix
                            end
                        end
                    end
                end
            catch
                # If reading fails, fall back to adjacency matrix
            end
            
            # Default to adjacency matrix
            return read_graph_from_adjacency_matrix(filename)
        end
    end

    """
        read_graph_from_edgelist(filename::String)

        Reads a graph from edge list file with format:
        source,destination
        1,2
        1,3
        2,4
    """
    function read_graph_from_edgelist(filename::String)::Tuple{Vector{Tuple{Int64,Int64}}, Dict{Int64,Set{Int64}}, Dict{Int64,Set{Int64}}, Set{Int64}}
        isfile(filename) || throw(SystemError("File not found: $filename"))
        
        edgelist = Vector{Tuple{Int64,Int64}}()
        outgoing_index = Dict{Int64,Set{Int64}}()
        incoming_index = Dict{Int64,Set{Int64}}()
        all_nodes = Set{Int64}()
        
        open(filename, "r") do file
            lines = readlines(file)
            
            # Skip header if present
            start_line = 1
            if length(lines) > 0
                header = strip(lines[1])
                if occursin("source", lowercase(header)) || occursin("destination", lowercase(header))
                    start_line = 2
                end
            end
            
            # Process edge data
            for i in start_line:length(lines)
                line = strip(lines[i])
                isempty(line) && continue
                startswith(line, "#") && continue
                
                # Parse edge
                parts = split(line, ',')
                if length(parts) != 2
                    throw(ArgumentError("Invalid edge format at line $i: '$line'. Expected 'source,destination'"))
                end
                
                try
                    source = parse(Int64, strip(parts[1]))
                    target = parse(Int64, strip(parts[2]))
                    
                    # Check for self-loops
                    if source == target
                        throw(ArgumentError("Self-loop detected at node $source (line $i)"))
                    end
                    
                    # Add to edge list
                    push!(edgelist, (source, target))
                    
                    # Track all nodes
                    push!(all_nodes, source, target)
                    
                    # Update outgoing index
                    if !haskey(outgoing_index, source)
                        outgoing_index[source] = Set{Int64}()
                    end
                    push!(outgoing_index[source], target)
                    
                    # Update incoming index
                    if !haskey(incoming_index, target)
                        incoming_index[target] = Set{Int64}()
                    end
                    push!(incoming_index[target], source)
                    
                catch e
                    throw(ArgumentError("Invalid integer format at line $i: '$line'. Error: $e"))
                end
            end
        end
        
        # Validate DAG property
        if has_cycle(outgoing_index)
            throw(ArgumentError("Graph contains cycles - must be a DAG"))
        end

        # Find source nodes (nodes with no incoming edges)
        source_nodes = setdiff(all_nodes, keys(incoming_index))
        
        # Initialize incoming index for source nodes
        for node in source_nodes
            incoming_index[node] = Set{Int64}()
        end
        
        return edgelist, outgoing_index, incoming_index, source_nodes
    end

    """
        read_graph_from_adjacency_matrix(filename::String)

        Reads a graph from adjacency matrix CSV file (integers 0/1 only).
        This is the original functionality.
    """
    function read_graph_from_adjacency_matrix(filename::String)::Tuple{Vector{Tuple{Int64,Int64}}, Dict{Int64,Set{Int64}}, Dict{Int64,Set{Int64}}, Set{Int64}}
        isfile(filename) || throw(SystemError("File not found: $filename"))
        
        # Read adjacency matrix (integers only)
        adj_matrix = readdlm(filename, ',', Int)
        
        # Validate square matrix
        n_rows, n_cols = size(adj_matrix)
        if n_rows != n_cols
            throw(ArgumentError("Adjacency matrix must be square, got $(n_rows)x$(n_cols)"))
        end
        
        # Validate values are 0 or 1 only
        if !all(x -> x in [0, 1], adj_matrix)
            throw(ArgumentError("Adjacency matrix must contain only 0 and 1 values"))
        end
        
        edgelist = Vector{Tuple{Int64,Int64}}()
        outgoing_index = Dict{Int64,Set{Int64}}()
        incoming_index = Dict{Int64,Set{Int64}}()
        all_nodes = Set{Int64}(1:n_rows)
        
        # Build graph structure from adjacency matrix
        for i in 1:n_rows, j in 1:n_cols
            if adj_matrix[i, j] == 1
                # Check for self-loops
                if i == j
                    throw(ArgumentError("Self-loop detected at node $i"))
                end
                
                push!(edgelist, (i, j))
                
                # Update outgoing index
                if !haskey(outgoing_index, i)
                    outgoing_index[i] = Set{Int64}()
                end
                push!(outgoing_index[i], j)
                
                # Update incoming index
                if !haskey(incoming_index, j)
                    incoming_index[j] = Set{Int64}()
                end
                push!(incoming_index[j], i)
            end
        end

        # Validate DAG property
        if has_cycle(outgoing_index)
            throw(ArgumentError("Graph contains cycles - must be a DAG"))
        end

        # Find source nodes (nodes with no incoming edges)
        source_nodes = setdiff(all_nodes, keys(incoming_index))
        
        # Initialize incoming index for source nodes
        for node in source_nodes
            incoming_index[node] = Set{Int64}()
        end
        
        return edgelist, outgoing_index, incoming_index, source_nodes
    end

    """
        has_cycle(graph::Dict{Int64,Set{Int64}})

        Helper function to detect cycles in the graph using DFS.
    """
    function has_cycle(graph::Dict{Int64,Set{Int64}})
        visited = Set{Int64}()
        temp_visited = Set{Int64}()
        
        function dfs(node::Int64)
            if node in temp_visited
                return true  # Cycle detected
            end
            if node in visited
                return false
            end
            push!(temp_visited, node)
            
            if haskey(graph, node)
                for neighbor in graph[node]
                    if dfs(neighbor)
                        return true
                    end
                end
            end
            
            delete!(temp_visited, node)
            push!(visited, node)
            return false
        end
        
        for node in keys(graph)
            if dfs(node)
                return true
            end
        end
        return false
    end

    """
        deserialize_probability_value(data::Any)

        Deserialize probability values from JSON, returning actual pbox objects.
        Supports all ProbabilityBoundsAnalysis constructors.
    """
    function deserialize_probability_value(data::Any)
        # Handle simple numeric values
        if isa(data, Real)
            return Float64(data)
        end
        
        # Handle Dict (complex types)
        if !isa(data, AbstractDict)
            throw(ArgumentError("Invalid probability data format: $(typeof(data))"))
        end
        
        if data["type"] == "interval"
            return Interval(Float64(data["lower"]), Float64(data["upper"]))
            
        elseif data["type"] == "pbox"
            construction_type = data["construction_type"]
            
            if construction_type == "scalar"
                # pbox(value) -> Create precise pbox using makepbox(interval(value, value))
                value = Float64(data["value"])
                return PBA.makepbox(PBA.interval(value, value))
                
            elseif construction_type == "interval"
                # pbox(lower, upper) -> Create interval pbox
                lower = Float64(data["lower"])
                upper = Float64(data["upper"])
                return PBA.makepbox(PBA.interval(lower, upper))
                
            elseif construction_type == "parametric"
                # normal(mean, std), uniform(a, b), etc.
                return create_parametric_pbox(data)
                
            elseif construction_type == "parametric_interval"
                # normal(interval(0,1), 1), uniform(interval(0,1), interval(2,3))
                return create_parametric_interval_pbox(data)
                
            elseif construction_type == "envelope"
                # env(d1, d2, ...)
                return create_envelope_pbox(data)
                
            elseif construction_type == "distribution_free"
                # meanVar(ml, mh, vl, vh), meanMin(ml, mh, min_val), etc.
                return create_distribution_free_pbox(data)
                
            elseif construction_type == "complex"
                # Fallback - create using moments
                ml = get(data, "ml", 0.0)
                mh = get(data, "mh", 1.0) 
                vl = get(data, "vl", 0.0)
                vh = get(data, "vh", 1.0)
                return PBA.meanVar(ml, mh, vl, vh)
                
            else
                throw(ArgumentError("Unknown pbox construction type: $construction_type"))
            end
            
        else
            throw(ArgumentError("Unknown probability type: $(data["type"])"))
        end
    end

    """
        parse_capacity_value(data::Any)::Float64

        Parse capacity values from JSON. Supports finite numeric values and infinite
        values represented as "Inf", "+Inf", "Infinity", "+Infinity", or "∞".
    """
    function parse_capacity_value(data::Any)::Float64
        value = if isa(data, Real)
            Float64(data)
        elseif isa(data, String)
            token = lowercase(strip(data))
            if token in ("inf", "+inf", "infinity", "+infinity", "∞")
                Inf
            else
                try
                    parse(Float64, token)
                catch
                    throw(ArgumentError("Invalid capacity value string: '$data'"))
                end
            end
        else
            throw(ArgumentError("Unsupported capacity type: $(typeof(data))"))
        end

        if isnan(value)
            throw(ArgumentError("Capacity cannot be NaN"))
        end
        if value < 0.0
            throw(ArgumentError("Capacity must be non-negative"))
        end

        return value
    end

    """
        read_edge_capacities_from_json(filename::String)

        Read capacity inputs from JSON using schema:

        {
            "edges": [
                {"source": 1, "destination": 2, "capacity": 10.0},
                {"source": 2, "destination": 3, "capacity": "Inf"}
            ]
        }

        Returns Dict{Tuple{Int64, Int64}, Float64}.
    """
    function read_edge_capacities_from_json(filename::String)::Dict{Tuple{Int64, Int64}, Float64}
        isfile(filename) || throw(SystemError("File not found: $filename"))

        data = JSON.parsefile(filename)
        haskey(data, "edges") || throw(ArgumentError("JSON file must contain 'edges' key"))
        isa(data["edges"], AbstractVector) || throw(ArgumentError("'edges' must be an array"))

        capacities = Dict{Tuple{Int64, Int64}, Float64}()

        for (i, edge_data) in enumerate(data["edges"])
            isa(edge_data, AbstractDict) || throw(ArgumentError("Invalid edge entry at index $i: expected object"))
            haskey(edge_data, "source") || throw(ArgumentError("Missing 'source' in edges[$i]"))
            haskey(edge_data, "destination") || throw(ArgumentError("Missing 'destination' in edges[$i]"))
            haskey(edge_data, "capacity") || throw(ArgumentError("Missing 'capacity' in edges[$i]"))

            source = Int64(edge_data["source"])
            target = Int64(edge_data["destination"])
            source == target && throw(ArgumentError("Self-loop detected at node $source (edges[$i])"))

            capacities[(source, target)] = parse_capacity_value(edge_data["capacity"])
        end

        isempty(capacities) && throw(ArgumentError("No edges found in 'edges' array"))
        return capacities
    end

    """
        read_node_capacities_from_json(filename::String)
            -> Dict{Int64, Float64}

    Read per-node capacity constraints from a JSON file.

    Schema:
      {
        "nodes": [
          {"node": 5, "capacity": 4.0},
          {"node": 6, "capacity": "Inf"}
        ]
      }

    Node IDs must be positive integers present in the graph.
    Capacity values must be finite, nonnegative, and non-NaN.
    Supports "Inf" string for unconstrained nodes.
    Returns Dict{Int64,Float64} mapping node ID to capacity.
    Throws ArgumentError for invalid values.
    Throws SystemError if file not found.
    """
    function read_node_capacities_from_json(
        filename::String
    )::Dict{Int64,Float64}
        isfile(filename) || throw(SystemError(
            "File not found: $filename"))

        data = JSON.parsefile(filename)
        haskey(data, "nodes") || throw(ArgumentError(
            "JSON must contain 'nodes' key"))
        isa(data["nodes"], AbstractVector) || throw(ArgumentError(
            "'nodes' must be an array"))

        node_capacities = Dict{Int64,Float64}()

        for (i, entry) in enumerate(data["nodes"])
            isa(entry, AbstractDict) || throw(ArgumentError(
                "Invalid node entry at index $i: expected object"))
            haskey(entry, "node") || throw(ArgumentError(
                "Missing 'node' key in nodes[$i]"))
            haskey(entry, "capacity") || throw(ArgumentError(
                "Missing 'capacity' key in nodes[$i]"))

            node_id = Int64(entry["node"])
            node_id > 0 || throw(ArgumentError(
                "Node ID must be positive, got $node_id at nodes[$i]"))

            node_capacities[node_id] =
                parse_capacity_value(entry["capacity"])
        end

        isempty(node_capacities) && throw(ArgumentError(
            "No nodes found in 'nodes' array"))

        return node_capacities
    end

    """
        read_capacities_input(filename::String)

        Alias for capacity input parsing entry point.
    """
    read_capacities_input(filename::String) = read_edge_capacities_from_json(filename)

    """
        create_parametric_pbox(data::Dict)

        Create parametric pbox using PBA constructors: normal(mean, std), uniform(a, b), etc.
    """
    function create_parametric_pbox(data::AbstractDict)
        shape = data["shape"]
        params = data["params"]
        
        try
            if shape == "normal"
                if length(params) >= 2
                    return PBA.normal(params[1], params[2])
                else
                    return PBA.normal(params[1], 1.0)  # Default std = 1
                end
                
            elseif shape == "uniform"
                if length(params) >= 2
                    return PBA.uniform(params[1], params[2])
                else
                    return PBA.uniform(0.0, params[1])  # Default min = 0
                end
                
            elseif shape == "beta"
                if length(params) >= 2
                    return PBA.beta(params[1], params[2])
                else
                    return PBA.beta(params[1], 1.0)
                end
                
            elseif shape == "exponential"
                return PBA.exponential(params[1])
                
            elseif shape == "erlang"
                if length(params) >= 2
                    return PBA.erlang(Int(params[1]), params[2])
                else
                    return PBA.erlang(Int(params[1]), 1.0)
                end
                
            elseif shape == "cauchy"
                if length(params) >= 2
                    return PBA.cauchy(params[1], params[2])
                else
                    return PBA.cauchy(params[1], 1.0)
                end
                
            elseif shape == "chi"
                return PBA.chi(Int(params[1]))
                
            elseif shape == "chisq"
                return PBA.chisq(Int(params[1]))
                
            elseif shape == "cosine"
                if length(params) >= 2
                    return PBA.cosine(params[1], params[2])
                else
                    return PBA.cosine(params[1], 1.0)
                end
                
            # Add more distributions as needed
            else
                @warn "Unknown parametric distribution: $shape, creating normal with mean $(params[1])"
                return PBA.normal(params[1], length(params) >= 2 ? params[2] : 1.0)
            end
            
        catch e
            @warn "Error creating $shape pbox: $e, falling back to interval"
            if length(params) >= 2
                return PBA.makepbox(PBA.interval(params[1], params[2]))
            else
                val = params[1]
                return PBA.makepbox(PBA.interval(val, val))
            end
        end
    end

    """
        create_parametric_interval_pbox(data::Dict)

        Create parametric pbox with interval parameters: normal(interval(0,1), 1)
    """
    function create_parametric_interval_pbox(data::AbstractDict)
        shape = data["shape"]
        params = data["params"]  # Array of parameter specifications
        
        try
            if shape == "normal"
                # params could be [{"type": "interval", "lower": 0, "upper": 1}, 1.0]
                mean_param = params[1]
                std_param = length(params) >= 2 ? params[2] : 1.0
                
                if isa(mean_param, Dict) && mean_param["type"] == "interval"
                    mean_interval = PBA.interval(mean_param["lower"], mean_param["upper"])
                    if isa(std_param, Dict) && std_param["type"] == "interval"
                        std_interval = PBA.interval(std_param["lower"], std_param["upper"])
                        return PBA.normal(mean_interval, std_interval)
                    else
                        return PBA.normal(mean_interval, std_param)
                    end
                else
                    return PBA.normal(mean_param, std_param)
                end
                
            elseif shape == "uniform"
                a_param = params[1]
                b_param = length(params) >= 2 ? params[2] : 1.0
                
                if isa(a_param, Dict) && a_param["type"] == "interval"
                    a_interval = PBA.interval(a_param["lower"], a_param["upper"])
                    if isa(b_param, Dict) && b_param["type"] == "interval"
                        b_interval = PBA.interval(b_param["lower"], b_param["upper"])
                        return PBA.uniform(a_interval, b_interval)
                    else
                        return PBA.uniform(a_interval, b_param)
                    end
                else
                    return PBA.uniform(a_param, b_param)
                end
                
            # Add more interval-parametric distributions as needed
            else
                @warn "Unknown interval-parametric distribution: $shape"
                return create_parametric_pbox(data)  # Fallback to regular parametric
            end
            
        catch e
            @warn "Error creating interval-parametric $shape pbox: $e"
            return create_parametric_pbox(data)  # Fallback
        end
    end

    """
        create_envelope_pbox(data::Dict)

        Create envelope pbox: env(d1, d2, ...)
    """
    function create_envelope_pbox(data::AbstractDict)
        components = data["components"]
        
        try
            # Recursively deserialize each component
            pbox_components = []
            for component in components
                push!(pbox_components, deserialize_probability_value(component))
            end
            
            # Create envelope
            if length(pbox_components) >= 2
                result = pbox_components[1]
                for i in 2:length(pbox_components)
                    result = PBA.env(result, pbox_components[i])
                end
                return result
            else
                return pbox_components[1]
            end
            
        catch e
            @warn "Error creating envelope pbox: $e"
            # Fallback to simple interval
            return PBA.makepbox(PBA.interval(0.0, 1.0))
        end
    end

    """
        create_distribution_free_pbox(data::Dict)

        Create distribution-free pbox: meanVar(ml, mh, vl, vh), etc.
    """
    function create_distribution_free_pbox(data::AbstractDict)
        method = data["method"]
        params = data["params"]
        
        try
            if method == "meanVar"
                if length(params) >= 4
                    return PBA.meanVar(params[1], params[2], params[3], params[4])
                else
                    @warn "meanVar requires 4 parameters, got $(length(params))"
                    return PBA.meanVar(params[1], params[2], 0.0, 1.0)
                end
                
            elseif method == "meanMin"
                if length(params) >= 3
                    return PBA.meanMin(params[1], params[2], params[3])
                else
                    return PBA.meanMin(params[1], params[2], 0.0)
                end
                
            elseif method == "meanMax"
                if length(params) >= 3
                    return PBA.meanMax(params[1], params[2], params[3])
                else
                    return PBA.meanMax(params[1], params[2], 1.0)
                end
                
            elseif method == "meanMinMax"
                if length(params) >= 4
                    return PBA.meanMinMax(params[1], params[2], params[3], params[4])
                else
                    @warn "meanMinMax requires 4 parameters"
                    return PBA.meanMinMax(params[1], params[2], 0.0, 1.0)
                end
                
            elseif method == "minMaxMeanVar"
                if length(params) >= 4
                    return PBA.minMaxMeanVar(params[1], params[2], params[3], params[4])
                else
                    @warn "minMaxMeanVar requires 4 parameters"
                    return PBA.minMaxMeanVar(0.0, 1.0, params[1], params[2])
                end
                
            else
                @warn "Unknown distribution-free method: $method"
                return PBA.meanVar(params[1], params[2], 0.0, 1.0)
            end
            
        catch e
            @warn "Error creating distribution-free pbox ($method): $e"
            # Fallback
            return PBA.makepbox(PBA.interval(params[1], params[2]))
        end
    end

    """
        deserialize_pbox_value(data::Any)

        Deserialize probability values from JSON, returning only pbox objects.
        Throws ArgumentError for non-pbox compatible data.
    """
    function deserialize_pbox_value(data::Any)::pbox
        # Handle Dict (pbox types only)
        if !isa(data, AbstractDict)
            throw(ArgumentError("pbox deserializer requires Dict format, got $(typeof(data))"))
        end
        
        if !haskey(data, "type") || data["type"] != "pbox"
            throw(ArgumentError("pbox deserializer requires type='pbox', got type='$(get(data, "type", "missing"))'"))
        end
        
        construction_type = data["construction_type"]
        
        if construction_type == "scalar"
            # pbox(value) -> Create precise pbox using makepbox(interval(value, value))
            value = Float64(data["value"])
            return PBA.makepbox(PBA.interval(value, value))
            
        elseif construction_type == "interval"
            # pbox(lower, upper) -> Create interval pbox
            lower = Float64(data["lower"])
            upper = Float64(data["upper"])
            return PBA.makepbox(PBA.interval(lower, upper))
            
        elseif construction_type == "parametric"
            # normal(mean, std), uniform(a, b), etc.
            return create_parametric_pbox(data)
            
        elseif construction_type == "parametric_interval"
            # normal(interval(0,1), 1), uniform(interval(0,1), interval(2,3))
            return create_parametric_interval_pbox(data)
            
        elseif construction_type == "envelope"
            # env(d1, d2, ...)
            return create_envelope_pbox(data)
            
        elseif construction_type == "distribution_free"
            # meanVar(ml, mh, vl, vh), meanMin(ml, mh, min_val), etc.
            return create_distribution_free_pbox(data)
            
        elseif construction_type == "complex"
            # Fallback - create using moments
            ml = get(data, "ml", 0.0)
            mh = get(data, "mh", 1.0)
            vl = get(data, "vl", 0.0)
            vh = get(data, "vh", 1.0)
            return PBA.meanVar(ml, mh, vl, vh)
            
        else
            throw(ArgumentError("Unknown pbox construction type: $construction_type"))
        end
    end

    """
        read_node_priors_from_json_pbox(filename::String)
        
        Read pbox node priors from JSON file. Returns Dict{Int64, pbox}.
    """
    function read_node_priors_from_json_pbox(filename::String)::Dict{Int64, pbox}
        isfile(filename) || throw(SystemError("File not found: $filename"))
        
        data = JSON.parsefile(filename)
        if !haskey(data, "nodes")
            throw(ArgumentError("JSON file must contain 'nodes' key"))
        end
        
        result = Dict{Int64, pbox}()
        for (node_str, node_data) in data["nodes"]
            node_id = parse(Int, node_str)
            result[node_id] = deserialize_pbox_value(node_data)
        end
        return result
    end

    """
        read_edge_probabilities_from_json_pbox(filename::String)
        
        Read pbox edge probabilities from JSON file. Returns Dict{Tuple{Int64, Int64}, pbox}.
    """
    function read_edge_probabilities_from_json_pbox(filename::String)::Dict{Tuple{Int64, Int64}, pbox}
        isfile(filename) || throw(SystemError("File not found: $filename"))
        
        data = JSON.parsefile(filename)
        if !haskey(data, "links")
            throw(ArgumentError("JSON file must contain 'links' key"))
        end
        
        result = Dict{Tuple{Int64, Int64}, pbox}()
        for (edge_str, edge_data) in data["links"]
            edge_match = match(r"\((\d+),(\d+)\)", edge_str)
            if edge_match !== nothing
                source = parse(Int, edge_match.captures[1])
                target = parse(Int, edge_match.captures[2])
                result[(source, target)] = deserialize_pbox_value(edge_data)
            end
        end
        return result
    end

    """
        read_node_priors_from_json_interval(filename::String)
        
        Read Interval node priors from JSON file. Returns Dict{Int64, Interval}.
    """
    function read_node_priors_from_json_interval(filename::String)::Dict{Int64, Interval}
        isfile(filename) || throw(SystemError("File not found: $filename"))
        
        data = JSON.parsefile(filename)
        if !haskey(data, "nodes")
            throw(ArgumentError("JSON file must contain 'nodes' key"))
        end
        
        result = Dict{Int64, Interval}()
        for (node_str, node_data) in data["nodes"]
            node_id = parse(Int, node_str)
            result[node_id] = deserialize_probability_value(node_data)::Interval
        end
        return result
    end

    """
        read_edge_probabilities_from_json_interval(filename::String)
        
        Read Interval edge probabilities from JSON file. Returns Dict{Tuple{Int64, Int64}, Interval}.
    """
    function read_edge_probabilities_from_json_interval(filename::String)::Dict{Tuple{Int64, Int64}, Interval}
        isfile(filename) || throw(SystemError("File not found: $filename"))
        
        data = JSON.parsefile(filename)
        if !haskey(data, "links")
            throw(ArgumentError("JSON file must contain 'links' key"))
        end
        
        result = Dict{Tuple{Int64, Int64}, Interval}()
        for (edge_str, edge_data) in data["links"]
            edge_match = match(r"\((\d+),(\d+)\)", edge_str)
            if edge_match !== nothing
                source = parse(Int, edge_match.captures[1])
                target = parse(Int, edge_match.captures[2])
                result[(source, target)] = deserialize_probability_value(edge_data)::Interval
            end
        end
        return result
    end

    """
        read_node_priors_from_json_float64(filename::String)
        
        Read Float64 node priors from JSON file. Returns Dict{Int64, Float64}.
    """
    function read_node_priors_from_json_float64(filename::String)::Dict{Int64, Float64}
        isfile(filename) || throw(SystemError("File not found: $filename"))
        
        data = JSON.parsefile(filename)
        if !haskey(data, "nodes")
            throw(ArgumentError("JSON file must contain 'nodes' key"))
        end
        
        result = Dict{Int64, Float64}()
        for (node_str, node_data) in data["nodes"]
            node_id = parse(Int, node_str)
            result[node_id] = Float64(node_data)
        end
        return result
    end

    """
        read_edge_probabilities_from_json_float64(filename::String)
        
        Read Float64 edge probabilities from JSON file. Returns Dict{Tuple{Int64, Int64}, Float64}.
    """
    function read_edge_probabilities_from_json_float64(filename::String)::Dict{Tuple{Int64, Int64}, Float64}
        isfile(filename) || throw(SystemError("File not found: $filename"))
        
        data = JSON.parsefile(filename)
        if !haskey(data, "links")
            throw(ArgumentError("JSON file must contain 'links' key"))
        end
        
        result = Dict{Tuple{Int64, Int64}, Float64}()
        for (edge_str, edge_data) in data["links"]
            edge_match = match(r"\((\d+),(\d+)\)", edge_str)
            if edge_match !== nothing
                source = parse(Int, edge_match.captures[1])
                target = parse(Int, edge_match.captures[2])
                result[(source, target)] = Float64(edge_data)
            end
        end
        return result
    end

    """
        read_node_priors_from_json(filename::String)

        Generic function that auto-detects type and calls the appropriate specific function.
        Returns properly typed dictionaries based on JSON data_type field.
    """
    function read_node_priors_from_json(filename::String)
        isfile(filename) || throw(SystemError("File not found: $filename"))
        
        data = JSON.parsefile(filename)
        data_type = get(data, "data_type", "Float64")
        
        if data_type == "Float64"
            return read_node_priors_from_json_float64(filename)::Dict{Int64, Float64}
        elseif data_type == "Interval"
            return read_node_priors_from_json_interval(filename)::Dict{Int64, Interval}
        elseif data_type == "pbox" || data_type == "ProbabilityBoundsAnalysis.pbox"
            return read_node_priors_from_json_pbox(filename)::Dict{Int64, pbox}
        else
            throw(ArgumentError("Unknown data_type: $data_type. Expected 'Float64', 'Interval', or 'pbox'"))
        end
    end

    """
        read_edge_probabilities_from_json(filename::String)

        Generic function that auto-detects type and calls the appropriate specific function.
        For guaranteed return types, use the type-specific functions directly.
    """
    function read_edge_probabilities_from_json(filename::String)
        isfile(filename) || throw(SystemError("File not found: $filename"))
        
        data = JSON.parsefile(filename)
        data_type = get(data, "data_type", "Float64")
        
        if data_type == "Float64"
            return read_edge_probabilities_from_json_float64(filename)
        elseif data_type == "Interval"
            return read_edge_probabilities_from_json_interval(filename)
        elseif data_type == "pbox" || data_type == "ProbabilityBoundsAnalysis.pbox"
            return read_edge_probabilities_from_json_pbox(filename)
        else
            # Generic fallback
            if !haskey(data, "links")
                throw(ArgumentError("JSON file must contain 'links' key"))
            end
            result = Dict{Tuple{Int64, Int64}, Any}()
            for (edge_str, edge_data) in data["links"]
                edge_match = match(r"\((\d+),(\d+)\)", edge_str)
                if edge_match !== nothing
                    source = parse(Int, edge_match.captures[1])
                    target = parse(Int, edge_match.captures[2])
                    result[(source, target)] = deserialize_probability_value(edge_data)
                end
            end
            return result
        end
    end

    """
        read_complete_network(adj_matrix_file::String, node_priors_file::String, edge_probs_file::String)

        Convenience function to read complete network from separate files.
    """
    function read_complete_network(adj_matrix_file::String, node_priors_file::String, edge_probs_file::String)
        # Read graph structure
        edgelist, outgoing_index, incoming_index, source_nodes = read_graph_to_dict(adj_matrix_file)
        
        # Read probabilities
        node_priors = read_node_priors_from_json(node_priors_file)
        edge_probabilities = read_edge_probabilities_from_json(edge_probs_file)
        
        # Validate that all edges in graph have corresponding probabilities
        for (source, target) in edgelist
            if !haskey(edge_probabilities, (source, target))
                throw(ArgumentError("Missing probability data for edge ($source,$target)"))
            end
        end
        
        # Validate that all nodes have prior probabilities
        all_nodes = union(Set(first.(edgelist)), Set(last.(edgelist)))
        for node in all_nodes
            if !haskey(node_priors, node)
                throw(ArgumentError("Missing prior probability for node $node"))
            end
        end
        
        return edgelist, outgoing_index, incoming_index, source_nodes, node_priors, edge_probabilities
    end

    # Include other functions (identify_fork_and_join_nodes, find_iteration_sets) as before...
    """
        identify_fork_and_join_nodes(outgoing_index, incoming_index)
    """
    function identify_fork_and_join_nodes(
        outgoing_index::Dict{Int64,Set{Int64}},
        incoming_index::Dict{Int64,Set{Int64}}
        )::Tuple{Set{Int64},Set{Int64}}
        
        fork_nodes = Set{Int64}()
        join_nodes = Set{Int64}()
    
        # Identify fork nodes
        for (node, children) in outgoing_index
            if length(children) > 1
                push!(fork_nodes, node)
            end
        end
    
        # Identify join nodes
        for (node, parents) in incoming_index
            if length(parents) > 1
                push!(join_nodes, node)
            end
        end
    
        return fork_nodes, join_nodes
    end
    
    """
        find_iteration_sets(edgelist, outgoing_index, incoming_index)
    """
    function find_iteration_sets(
        edgelist::Vector{Tuple{Int64,Int64}},
        outgoing_index::Dict{Int64,Set{Int64}},
        incoming_index::Dict{Int64,Set{Int64}}
        )::Tuple{Vector{Set{Int64}}, Dict{Int64, Set{Int64}}, Dict{Int64, Set{Int64}}}
        
        isempty(edgelist) && return (Vector{Set{Int64}}(), Dict{Int64,Set{Int64}}(), Dict{Int64,Set{Int64}}())
        
        # Find the maximum node id
        n = maximum(max(first(edge), last(edge)) for edge in edgelist)
        
        in_degree = zeros(Int, n)
        all_nodes = Set{Int64}()
        
        # Calculate initial in-degrees and collect all nodes
        for (source, target) in edgelist
            in_degree[target] += 1
            push!(all_nodes, source, target)
        end
        
        ancestors = Dict(node => Set{Int64}([node]) for node in all_nodes)
        descendants = Dict(node => Set{Int64}() for node in all_nodes)
    
        queue = Queue{Int64}()
        for node in all_nodes
            if !haskey(incoming_index, node) || isempty(incoming_index[node])
                enqueue!(queue, node)
            end
        end
        
        iteration_sets = Vector{Set{Int64}}()
        
        while !isempty(queue)
            current_set = Set{Int64}()
            
            # Process all nodes in the current level
            level_size = length(queue)
            for _ in 1:level_size
                node = dequeue!(queue)
                push!(current_set, node)
                
                # Process outgoing edges
                for target in get(outgoing_index, node, Set{Int64}())
                    # Update ancestors efficiently
                    if !issubset(ancestors[node], ancestors[target])
                        union!(ancestors[target], ancestors[node])
                    end
                    
                    # Update descendants efficiently
                    new_descendants = setdiff(descendants[target], descendants[node])
                    if !isempty(new_descendants)
                        union!(descendants[node], new_descendants, Set([target]))
                        # Propagate new descendants to all ancestors of the current node
                        for ancestor in ancestors[node]
                            if ancestor != node
                                union!(descendants[ancestor], new_descendants, Set([target]))
                            end
                        end
                    elseif !(target in descendants[node])
                        push!(descendants[node], target)
                        # Propagate new descendant to all ancestors of the current node
                        for ancestor in ancestors[node]
                            if ancestor != node
                                push!(descendants[ancestor], target)
                            end
                        end
                    end
                    
                    in_degree[target] -= 1
                    if in_degree[target] == 0
                        enqueue!(queue, target)
                    end
                end
            end
            
            push!(iteration_sets, current_set)
        end
        
        return (iteration_sets, ancestors, descendants)
    end

end
