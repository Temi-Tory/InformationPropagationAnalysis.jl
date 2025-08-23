
# Check if this is the first run of the script for this julia repl session
# This is useful to avoid re-initializing the environment multiple times
if !@isdefined(script_initialized)
    println("First run - initializing...")

    import Fontconfig
    using DataFrames, DelimitedFiles, Distributions,
        DataStructures, SparseArrays, BenchmarkTools,
        Combinatorics, Dates

    # Ensure we're running from the project root directory
    current_dir = pwd()
    # Include the InformationPropagationAnalysis module
    include("../src/InformationPropagationAnalysis.jl")
    using .InformationPropagationAnalysis

    # Mark as initialized
    global script_initialized = true
    println("Initialization complete!")
else
    println("Subsequent run - skipping initialization")
end

network_name = "power-network"

#more examples available in data folder such as:
# network_name = "mlgw-gas-network"
# network_name = "single-mission-drone-network"
# network_name = "drone-medical-delivery-network"





data_type = "float"


# Construct file paths using new folder structure
base_path = joinpath(@__DIR__, "data", network_name)

# Option 1: Use edge file (recommended)
filepath_graph = joinpath(base_path, network_name * ".EDGES");
json_network_name = replace(network_name, "_" => "-")  # Convert underscores to hyphens for JSON files
filepath_node_json = joinpath(base_path, data_type, json_network_name * "-nodepriors.json")
filepath_edge_json = joinpath(base_path, data_type, json_network_name * "-linkprobabilities.json")



if !isfile(filepath_graph)
    error("Graph file not found: $filepath_graph")
end
if !isfile(filepath_node_json)
    error("Node priors file not found: $filepath_node_json")
end
if !isfile(filepath_edge_json)
    error("Edge probabilities file not found: $filepath_edge_json")
end

# Read the graph and node priors

# Option 1: Separate calls (gives you more control)
edgelist, outgoing_index, incoming_index, source_nodes = read_graph_to_dict(filepath_graph)

allnodes = # Get all nodes from the outgoing index
    collect(keys(incoming_index));
sink_nodes = #nodes with no keys in outgoing_index or with empty outgoing_index
    filter(node -> !haskey(outgoing_index, node) || isempty(outgoing_index[node]), allnodes);

node_priors = read_node_priors_from_json(filepath_node_json)

edge_probabilities = read_edge_probabilities_from_json(filepath_edge_json)


# Identify network structure
fork_nodes, join_nodes = identify_fork_and_join_nodes(outgoing_index, incoming_index)
iteration_sets, ancestors, descendants = find_iteration_sets(edgelist, outgoing_index, incoming_index)


println(" finding root diamonds");
# Diamond structure analysis (if you have this function)
root_diamonds = identify_and_group_diamonds(
    join_nodes,
    incoming_index,
    ancestors,
    descendants,
    source_nodes,
    fork_nodes,
    edgelist,
    node_priors,
    iteration_sets
);
l_root_diamonds = length(root_diamonds);
diamond_joins = keys(root_diamonds);
println("Found $l_root_diamonds root_diamonds");
println("Starting build unique diamond storage");
unique_diamonds = build_unique_diamond_storage_depth_first_parallel(
    root_diamonds,
    node_priors,
    ancestors,
    descendants,
    iteration_sets
);
l_unique_diamonds = length(unique_diamonds)
println("Found $l_unique_diamonds unique_diamonds");
# show(keys(root_diamonds))

println("Starting iterative belief update");
start_time = time()
output = InformationPropagationAnalysis.update_beliefs_iterative(
    edgelist,
    iteration_sets,
    outgoing_index,
    incoming_index,
    source_nodes,
    node_priors,
    edge_probabilities,
    descendants,
    ancestors,
    root_diamonds,
    join_nodes,
    fork_nodes,
    unique_diamonds
);

# Calculate computation time
computation_time = time() - start_time

println("Starting exact_computation");
exact_start_time = time()

exact_results = (path_enumeration_result(
    outgoing_index,
    incoming_index,
    source_nodes,
    node_priors,
    edge_probabilities,
   # true
));


exact_computation_time = time() - exact_start_time 