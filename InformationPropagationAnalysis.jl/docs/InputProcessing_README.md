# InputProcessing Module

## Overview

The `InputProcessing` module handles network data import, validation, and preprocessing for graph-based analysis. It supports multiple file formats and probability types, with automatic format detection and comprehensive validation.

## Core Components

### Graph Input Formats

#### Edge List Format (.EDGES)
```
source,destination
1,3
2,3
3,4
```

#### Adjacency Matrix Format (.csv)
```
0,1,0,1
0,0,1,0
0,0,0,1
0,0,0,0
```

### Probability Type Support

#### Float64 (Standard)
```json
{
  "1": 0.8,
  "2": 0.7,
  "3": 0.5
}
```

#### Interval Arithmetic
```json
{
  "1": {"type": "interval", "lower": 0.7, "upper": 0.9},
  "2": {"type": "interval", "lower": 0.6, "upper": 0.8}
}
```

#### Probability Boxes (pbox)
```json
{
  "1": {
    "type": "pbox",
    "construction": "parametric",
    "shape": "normal",
    "params": [0.8, 0.1]
  },
  "2": {
    "type": "pbox", 
    "construction": "interval",
    "lower": 0.6,
    "upper": 0.8
  }
}
```

## Main Functions

### Graph Processing

#### `read_graph_to_dict(filename)`
Auto-detects file format and returns:
- `edgelist`: Vector of (source, target) tuples
- `outgoing_index`: Dict mapping nodes to their children
- `incoming_index`: Dict mapping nodes to their parents  
- `source_nodes`: Set of nodes with no incoming edges

Features:
- Automatic format detection (edge list vs adjacency matrix)
- Cycle detection for DAG validation
- Self-loop detection and prevention
- Comprehensive error handling

#### `identify_fork_and_join_nodes(outgoing_index, incoming_index)`
Returns:
- `fork_nodes`: Nodes with multiple outgoing edges
- `join_nodes`: Nodes with multiple incoming edges

#### `find_iteration_sets(edgelist, outgoing_index, incoming_index)`
Computes topological ordering and graph analysis:
- `iteration_sets`: Vector of node sets in processing order
- `ancestors`: Dict mapping each node to its ancestor set
- `descendants`: Dict mapping each node to its descendant set

Uses breadth-first topological sort with ancestor/descendant propagation.

### Probability Data Loading

#### Generic Functions (Auto-detect Type)
```julia
# Automatically detects probability type from JSON structure
node_priors = read_node_priors_from_json("nodepriors.json")
edge_probs = read_edge_probabilities_from_json("linkprobs.json")
```

#### Type-Specific Functions
```julia
# Explicit type specification
node_priors_float = read_node_priors_from_json_float64("nodepriors.json")
node_priors_interval = read_node_priors_from_json_interval("nodepriors.json")
node_priors_pbox = read_node_priors_from_json_pbox("nodepriors.json")
```

### Probability Box Construction

Supports multiple pbox construction methods:

#### Parametric Distributions
- Normal: `{"shape": "normal", "params": [mean, std]}`
- Uniform: `{"shape": "uniform", "params": [min, max]}`
- Beta: `{"shape": "beta", "params": [α, β]}`
- Exponential: `{"shape": "exponential", "params": [λ]}`

#### Interval-Based
```json
{
  "type": "pbox",
  "construction": "interval", 
  "lower": 0.6,
  "upper": 0.8
}
```

#### Envelope Method
```json
{
  "type": "pbox",
  "construction": "envelope",
  "distributions": [
    {"shape": "normal", "params": [0.7, 0.1]},
    {"shape": "uniform", "params": [0.6, 0.8]}
  ]
}
```

## Network Validation

The module performs comprehensive validation:

### Graph Structure
- DAG property (no cycles)
- No self-loops
- Consistent edge indices
- Square adjacency matrices
- Valid integer node IDs

### Probability Data
- All nodes have prior probabilities
- All edges have probability values
- Probability bounds [0,1] for all types
- Consistent JSON structure

### Error Handling
```julia
try
    edgelist, outgoing, incoming, sources = read_graph_to_dict("network.edges")
catch ArgumentError as e
    println("Graph validation failed: ", e)
end
```

## Usage Examples

### Basic Network Loading
```julia
# Load complete network
edgelist, outgoing, incoming, sources, node_priors, edge_probs = 
    read_complete_network("graph.edges", "nodepriors.json", "linkprobs.json")

# Analyze structure
fork_nodes, join_nodes = identify_fork_and_join_nodes(outgoing, incoming)
iteration_sets, ancestors, descendants = find_iteration_sets(edgelist, outgoing, incoming)
```

### Working with Different Probability Types
```julia
# Float64 probabilities
float_priors = read_node_priors_from_json_float64("float_priors.json")

# Interval probabilities  
interval_priors = read_node_priors_from_json_interval("interval_priors.json")

# Probability box uncertainties
pbox_priors = read_node_priors_from_json_pbox("pbox_priors.json")
```

## File Format Requirements

### Edge List Files
- CSV format with source,destination columns
- Integer node IDs only
- Optional header row
- No self-loops or duplicate edges

### Adjacency Matrix Files
- CSV format with 0/1 values only
- Square matrix (n×n)
- Row i, Column j = 1 indicates edge from i to j

### JSON Probability Files
- Valid JSON structure
- String keys (node IDs or edge tuples)
- Probability values between 0 and 1
- Type-specific structure for intervals and pboxes

## Performance Features

- **Memory Efficient**: Sparse representation using sets and dictionaries
- **Format Detection**: Automatic detection avoids user format specification
- **Incremental Validation**: Early error detection during parsing
- **Batch Processing**: Efficient topological sort with ancestor/descendant computation

The module provides a robust foundation for network analysis with comprehensive error handling and support for uncertainty quantification through different probability representations.