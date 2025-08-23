# InformationPropagationAnalysis.jl

A Julia package for exact reachability analysis and belief propagation in directed acyclic graphs (DAGs) with diamond structure optimization.

## Overview

InformationPropagationAnalysis.jl implements exact probabilistic inference for network reachability problems. The package identifies diamond structures in DAGs to optimize belief propagation computation, supporting multiple uncertainty representations including standard probabilities, interval arithmetic, and probability boxes.

## Key Features

- **Exact Belief Propagation**: Computes exact reachability probabilities without approximation
- **Diamond Structure Optimization**: 8-step algorithm identifies convergent path structures for computational efficiency
- **Multi-Type Uncertainty**: Supports Float64, Interval arithmetic, and probability box (pbox) representations  
- **Parallel Processing**: Thread-safe diamond processing with adaptive memory management
- **Network Validation**: Comprehensive input validation and error checking
- **Multiple File Formats**: Automatic detection of edge lists and adjacency matrices

## Installation

```julia
using Pkg
Pkg.add("InformationPropagationAnalysis")
```

## Quick Start

```julia
using InformationPropagationAnalysis

# Load network data
edgelist, outgoing_index, incoming_index, source_nodes = read_graph_to_dict("network.edges")
node_priors = read_node_priors_from_json("nodepriors.json")
edge_probabilities = read_edge_probabilities_from_json("linkprobs.json")

# Analyze network structure
fork_nodes, join_nodes = identify_fork_and_join_nodes(outgoing_index, incoming_index)
iteration_sets, ancestors, descendants = find_iteration_sets(edgelist, outgoing_index, incoming_index)

# Identify diamond structures
root_diamonds = identify_and_group_diamonds(
    join_nodes, incoming_index, ancestors, descendants,
    source_nodes, fork_nodes, edgelist, node_priors, iteration_sets
)

# Build computation structures
unique_diamonds = build_unique_diamond_storage_depth_first_parallel(
    root_diamonds, node_priors, ancestors, descendants, iteration_sets
)

# Compute reachability probabilities
results = update_beliefs_iterative(
    edgelist, iteration_sets, outgoing_index, incoming_index,
    source_nodes, node_priors, edge_probabilities, descendants,
    ancestors, root_diamonds, join_nodes, fork_nodes, unique_diamonds
)

# View results
for (node, probability) in results
    println("Node $node: $(round(probability, digits=4))")
end
```

## Supported File Formats

### Network Structure
- **Edge Lists**: CSV format with source,destination columns
- **Adjacency Matrices**: Square CSV matrices with 0/1 values
- **Auto-detection**: Automatically determines file format

### Probability Data (JSON)
```julia
# Float64 probabilities
{"1": 0.8, "2": 0.7, "3": 0.5}

# Interval probabilities  
{"1": {"type": "interval", "lower": 0.7, "upper": 0.9}}

# Probability boxes
{"1": {"type": "pbox", "construction": "parametric", "shape": "normal", "params": [0.8, 0.1]}}
```

## Core Algorithms

### Diamond Identification
8-step algorithm that identifies convergent path structures:
1. Collect shared fork ancestors
2. Extract induced subgraph
3. Identify conditioning nodes
4. Find intermediate nodes
5. Ensure edge completeness
6. Perform subsource analysis
7. Recursive diamond completeness
8. Build final diamond structure

### Belief Propagation
Handles three network cases:
- **Tree paths**: Direct probability propagation
- **Independent joins**: Inclusion-exclusion principle
- **Diamond structures**: Conditional enumeration over diamond states

### Parallel Processing
- Thread-safe diamond processing across iteration levels
- Adaptive memory management for large networks
- Hybrid optimization using lookup tables

## Examples

The package includes several example networks:

```julia
# Run power network example
include("examples/examples.jl")

# Available networks:
# - power-network
# - mlgw-gas-network  
# - drone-medical-delivery-network
# - munin-dag
# - water
```

## Validation Methods

Built-in validation against ground truth:

```julia
# Monte Carlo validation
mc_results = MC_result(edgelist, incoming_index, source_nodes, node_priors, edge_probabilities)

# Exact path enumeration
exact_results = path_enumeration_result(outgoing_index, incoming_index, source_nodes, node_priors, edge_probabilities)

# Compare results
max_difference = maximum(abs(results[node] - exact_results[node]) for node in keys(results))
```

## Documentation

Detailed documentation for each module:
- [DiamondProcessing](docs/DiamondProcessing_README.md): Diamond identification and preprocessing
- [ReachabilityAnalysis](docs/ReachabilityAnalysis_README.md): Belief propagation algorithms
- [InputProcessing](docs/InputProcessing_README.md): Network data loading and validation

## Performance

The package handles networks with:
- Thousands of nodes and edges
- Complex diamond structures with nested dependencies
- Multiple probability types with uncertainty quantification
- Parallel processing on multi-core systems

Computational complexity scales with diamond structure complexity rather than raw network size, providing significant speedups for networks with convergent path patterns.

## Applications

- **Network Reliability**: Infrastructure failure analysis
- **System Reachability**: Communication and transportation networks
- **Uncertainty Propagation**: Analysis with imprecise probabilities
- **Risk Assessment**: Multi-path dependency modeling

## Citation

If you use InformationPropagationAnalysis.jl in your research, please cite:

```bibtex
@software{InformationPropagationAnalysis,
  title={InformationPropagationAnalysis.jl: Exact Reachability Analysis with Diamond Structure Optimization},
  author={T. Ohiani and E. Patelli},
  year={2024},
  url={https://github.com/Temi-Tory/IPM_Working/tree/ipa.jlV1deployment}
}
```

### Related Publications

The theoretical foundation and diamond structure concepts are described in:

```bibtex
@inproceedings{ohiani2023information,
  title={The Information Propagation Method for Efficient Network Reliability Analysis},
  author={T. Ohiani and E. Patelli},
  booktitle={2023 7th International Conference on System Reliability and Safety (ICSRS)},
  pages={580--584},
  year={2023},
  organization={IEEE},
  address={Bologna, Italy},
  doi={10.1109/ICSRS59833.2023.10381157},
  keywords={Directed acyclic graph, Monte Carlo methods, Heuristic algorithms, Redundancy, Diamonds, Approximation algorithms, Safety, System Reliability, Probability Propagation, Network Graphs, Message Passing, Simulation}
}
```

## Acknowledgments

T. Ohiani is supported by the UK National Decommissioning Authority Bursary Studentship "Developing a resilience framework for decommissioning plan of a nuclear facility".

## License

MIT License. See [LICENSE](LICENSE) file for details.

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality  
4. Submit a pull request

## Requirements

- Julia ≥ 1.6
- Dependencies: JSON, DataStructures, Combinatorics, DelimitedFiles, Random, DataFrames, Distributions, SparseArrays, ProbabilityBoundsAnalysis