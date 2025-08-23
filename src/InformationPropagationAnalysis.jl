module InformationPropagationAnalysis

using JSON, DataStructures, Combinatorics
import ProbabilityBoundsAnalysis

const PBA = ProbabilityBoundsAnalysis
const pbox = ProbabilityBoundsAnalysis.pbox

include("InputProcessing.jl")
include("DiamondProcessing.jl")
include("ReachabilityAnalysis.jl")
include("Validation.jl")

using .InputProcessing
using .DiamondProcessing
using .ReachabilityAnalysis
using .Validation

export 
        # Core types
        DiamondsAtNode, Diamond, DiamondComputationData,
        Interval,

        # InputProcessing functions
        read_graph_to_dict,
        identify_fork_and_join_nodes, 
        find_iteration_sets,
        read_node_priors_from_json,
        read_edge_probabilities_from_json,
        read_complete_network,

        # Diamond processing functions
        identify_and_group_diamonds, build_unique_diamond_storage, 
        build_unique_diamond_storage_depth_first_parallel, create_diamond_hash_key,

        # Standard reachability analysis
        validate_network_data, update_beliefs_iterative, updateDiamondJoin,
        calculate_diamond_groups_belief, calculate_regular_belief, inclusion_exclusion,
        convert_to_pbox_data,

        # Comparison methods
        MC_result, has_path, path_enumeration_result

end
