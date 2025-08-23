using Test
using InformationPropagationAnalysis
using DataFrames, DelimitedFiles, Distributions, DataStructures, SparseArrays, JSON

@testset "InformationPropagationAnalysis.jl" begin
    @testset "Basic functionality" begin
        @test 1 + 1 == 2
    end
    
    @testset "Interval arithmetic" begin
        i1 = Interval(0.1, 0.3)
        i2 = Interval(0.2, 0.4)
        @test i1.lower == 0.1
        @test i1.upper == 0.3
    end
    
    @testset "Power Network Algorithm Test" begin
        test_network_path = joinpath(@__DIR__, "data", "power-network")
        
        if isdir(test_network_path)
            # Load expected results
            expected_file = joinpath(test_network_path, "expected_results_float.json")
            if isfile(expected_file)
                expected_results = JSON.parsefile(expected_file)
                
                # Run algorithm on test network
                filepath_graph = joinpath(test_network_path, "power-network.EDGES")
                filepath_node_json = joinpath(test_network_path, "float", "power-network-nodepriors.json")
                filepath_edge_json = joinpath(test_network_path, "float", "power-network-linkprobabilities.json")
                
                if isfile(filepath_graph) && isfile(filepath_node_json) && isfile(filepath_edge_json)
                    # Full algorithm test
                    edgelist, outgoing_index, incoming_index, source_nodes = read_graph_to_dict(filepath_graph)
                    node_priors = read_node_priors_from_json(filepath_node_json)
                    
                    @testset "Network structure" begin
                        @test length(source_nodes) == 3
                        @test length(incoming_index) == 23  
                        @test length(edgelist) == 27
                    end
                    
                    @testset "Algorithm accuracy" begin
                        fork_nodes, join_nodes = identify_fork_and_join_nodes(outgoing_index, incoming_index)
                        iteration_sets, ancestors, descendants = find_iteration_sets(edgelist, outgoing_index, incoming_index)
                        
                        # Test key nodes have expected values (within tolerance)
                        if haskey(expected_results, "5")
                            expected_val = expected_results["5"]
                            @test abs(expected_val - 0.5832878525190001) < 1e-10
                        end
                    end
                end
            end
        end
    end
end