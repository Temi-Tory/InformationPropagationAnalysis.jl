module InputProcessing

using Random, DataFrames, DelimitedFiles, Distributions, JSON,
    DataStructures, SparseArrays, Combinatorics

# Import ProbabilityBoundsAnalysis for pbox construction
import ProbabilityBoundsAnalysis
const PBA = ProbabilityBoundsAnalysis
const pbox = ProbabilityBoundsAnalysis.pbox

export Interval, read_graph_to_dict,
    # Generic functions (auto-detect type)
    read_node_priors_from_json, read_edge_probabilities_from_json,
    read_node_priors_from_json_pbox, read_edge_probabilities_from_json_pbox,
    read_node_priors_from_json_interval, read_edge_probabilities_from_json_interval,
    read_node_priors_from_json_float64, read_edge_probabilities_from_json_float64,
    identify_fork_and_join_nodes, find_iteration_sets, read_complete_network


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

Interval(value::Float64) = Interval(value, value)

function read_graph_to_dict(filename::String)::Tuple{Vector{Tuple{Int64,Int64}},Dict{Int64,Set{Int64}},Dict{Int64,Set{Int64}},Set{Int64}}
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

function read_graph_from_edgelist(filename::String)::Tuple{Vector{Tuple{Int64,Int64}},Dict{Int64,Set{Int64}},Dict{Int64,Set{Int64}},Set{Int64}}
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


function read_graph_from_adjacency_matrix(filename::String)::Tuple{Vector{Tuple{Int64,Int64}},Dict{Int64,Set{Int64}},Dict{Int64,Set{Int64}},Set{Int64}}
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


function deserialize_probability_value(data::Any)
    # Handle simple numeric values
    if isa(data, Real)
        return Float64(data)
    end

    # Handle Dict (complex types)
    if !isa(data, Dict)
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


function create_parametric_pbox(data::Dict)
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


function create_parametric_interval_pbox(data::Dict)
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


function create_envelope_pbox(data::Dict)
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
            for i in 2:eachindex(pbox_components)
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


function create_distribution_free_pbox(data::Dict)
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


function deserialize_pbox_value(data::Any)::pbox
    # Handle Dict (pbox types only)
    if !isa(data, Dict)
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


function read_node_priors_from_json_pbox(filename::String)::Dict{Int64,pbox}
    isfile(filename) || throw(SystemError("File not found: $filename"))

    data = JSON.parsefile(filename)
    if !haskey(data, "nodes")
        throw(ArgumentError("JSON file must contain 'nodes' key"))
    end

    result = Dict{Int64,pbox}()
    for (node_str, node_data) in data["nodes"]
        node_id = parse(Int, node_str)
        result[node_id] = deserialize_pbox_value(node_data)
    end
    return result
end


function read_edge_probabilities_from_json_pbox(filename::String)::Dict{Tuple{Int64,Int64},pbox}
    isfile(filename) || throw(SystemError("File not found: $filename"))

    data = JSON.parsefile(filename)
    if !haskey(data, "links")
        throw(ArgumentError("JSON file must contain 'links' key"))
    end

    result = Dict{Tuple{Int64,Int64},pbox}()
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


function read_node_priors_from_json_interval(filename::String)::Dict{Int64,Interval}
    isfile(filename) || throw(SystemError("File not found: $filename"))

    data = JSON.parsefile(filename)
    if !haskey(data, "nodes")
        throw(ArgumentError("JSON file must contain 'nodes' key"))
    end

    result = Dict{Int64,Interval}()
    for (node_str, node_data) in data["nodes"]
        node_id = parse(Int, node_str)
        result[node_id] = deserialize_probability_value(node_data)::Interval
    end
    return result
end


function read_edge_probabilities_from_json_interval(filename::String)::Dict{Tuple{Int64,Int64},Interval}
    isfile(filename) || throw(SystemError("File not found: $filename"))

    data = JSON.parsefile(filename)
    if !haskey(data, "links")
        throw(ArgumentError("JSON file must contain 'links' key"))
    end

    result = Dict{Tuple{Int64,Int64},Interval}()
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


function read_node_priors_from_json_float64(filename::String)::Dict{Int64,Float64}
    isfile(filename) || throw(SystemError("File not found: $filename"))

    data = JSON.parsefile(filename)
    if !haskey(data, "nodes")
        throw(ArgumentError("JSON file must contain 'nodes' key"))
    end

    result = Dict{Int64,Float64}()
    for (node_str, node_data) in data["nodes"]
        node_id = parse(Int, node_str)
        result[node_id] = Float64(node_data)
    end
    return result
end


function read_edge_probabilities_from_json_float64(filename::String)::Dict{Tuple{Int64,Int64},Float64}
    isfile(filename) || throw(SystemError("File not found: $filename"))

    data = JSON.parsefile(filename)
    if !haskey(data, "links")
        throw(ArgumentError("JSON file must contain 'links' key"))
    end

    result = Dict{Tuple{Int64,Int64},Float64}()
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


function read_node_priors_from_json(filename::String)
    isfile(filename) || throw(SystemError("File not found: $filename"))

    data = JSON.parsefile(filename)
    data_type = get(data, "data_type", "Float64")

    if data_type == "Float64"
        return read_node_priors_from_json_float64(filename)::Dict{Int64,Float64}
    elseif data_type == "Interval"
        return read_node_priors_from_json_interval(filename)::Dict{Int64,Interval}
    elseif data_type == "pbox" || data_type == "ProbabilityBoundsAnalysis.pbox"
        return read_node_priors_from_json_pbox(filename)::Dict{Int64,pbox}
    else
        throw(ArgumentError("Unknown data_type: $data_type. Expected 'Float64', 'Interval', or 'pbox'"))
    end
end


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
        result = Dict{Tuple{Int64,Int64},Any}()
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

function find_iteration_sets(
    edgelist::Vector{Tuple{Int64,Int64}},
    outgoing_index::Dict{Int64,Set{Int64}},
    incoming_index::Dict{Int64,Set{Int64}}
)::Tuple{Vector{Set{Int64}},Dict{Int64,Set{Int64}},Dict{Int64,Set{Int64}}}

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
