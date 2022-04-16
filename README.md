# Structure Network
## Introduction
Residue interaction networks (RINs) is an abstraction of protein 3D structure as a graph. 
RINs represent every residue of the protein as node (v) and every non-bonded residue interaction as edge (E).
We weight the E using a CHARMM-based interaction energy. 
Single point mutations on the protein can be represented as a mutant RIN and compared to wild-type RIN.  
Change in residue priorities upon mutation can be evaluated by a comparative RIN analysis. 

## Project Structure
Current Project uses multiple concepts from protein structures, geometry, and network analysis.
All these concepts are organized as submodule in the project.

*Module Name* |  *Description*
------------- |  -------------
`ds` | Data structure 
`geometry` | Geometry operations 
`networks` | Network analysis
`structure` | Protein structure parsing

## Installation
The `structure_network` is a Python library developed on Python __3.8__. 
We recommend a virtual environment setup with the following steps : 
```commandline
conda create -n <virtname> python=3.8
conda activate <virtname>
pip -r <install directory>/requirements/common.txt
```
Prior to library usage, please ensure the root directory (installation directory) is in the python path.

```commandline
# for linux environments
export PYTHONPATH="<install directory>":$PYTHONPATH
```

## Usage 
### Example notebook : NetworkAnalysisDemo
This notebook demonstrates 
 - how to import a target PDB and view it ![formula](https://render.githubusercontent.com/render/math?math=C_{\alpha}-C_{\alpha}) network.
 - how to compute network properties for the pdb-derived structure network 
 - how to derive strongly connected components on the undirected structure network (Markov based clustering)
 - how to create a point mutation of choice and estimate ts effect on network properties 
 - how to define sites (groups of residues) within the structure network and analyze network communication properties

#### Build structure network
The network generation uses various parameters, viz. `distance_threshold`, and `energy_threshold`.
Distance threshold is the maximum distance bound in Angstrom to draw an edge in the network. 
Energy threshold is the minimum contact energy value needed to be a valid contact. 
`threshold_type` dynamically changes the thresholding conditions. 
The `potential` parameter provides a choice for pair potential selection, supported potentials are `MJ` and `CHARMM`.

The Derived network represents ![formula](https://render.githubusercontent.com/render/math?math=C_{\alpha})  network.
The edges of the structure network is weighted by the pair potential score.
Each node of the network represents a residue of the PDB structure. The library
allows query the amino acid properties for each node.

#### Compute Structure Network Properties
`structure_network` library allow computation of the various network properties on the 
pdb-derived network. Typically the pair potential is a positive score. Unlike energetics 
higher score represents favourable residue interaction. For such reason to derive a meanigful and
consistent analysis based on shortest path, prior to analysis each edge weights inverted by 
```max_value + min_value - current_value```. Different network analysis supported by the libraries are
node and edge centrality computation, network flow analysis, connected component analysis, and markov 
clustering. 

#### Markov Clustering
Markov clustering is carried out the weighted protein structure network, where the weights are pair potentials/interaction energy. 
Unlike the shortest path analysis, higher edge weights in markov clustering suggest higher chances of random walk through the edge.
Our main objective is to perform markov clustering on protein structure network to discover strongly connected structural components. 

#### Perform Sequence Mutation
`MutationAnalyzer` creates an internal representation of point mutation on structure. It assumes the point mutation does not affect residue contact pattern, but affects the side-chain contacts. 
It allows comparison of a single point residue change with respect to the base structure. 
`register_mutation` uses base structure as input, and a valid residue id and valid amino-acid change. Multiple single point mutations can be prepared in one-shot. A comparative analysis of network properties can be carried out among the point mutations using this library. 
In situation of improper mutation specification the model silently ignores the specified mutation. `registered_mutations` API allows access to all prepared mutations.

#### Group Analysis
Library allows group definitions by specifying a set of residues in the structure. 
Residue groups can indicate source (ligand binding) and target (effector) residues, or sequence-based motifs. 
Further analysis of structure-driven communication between these residue groups can be carried out with the library. 
The library allows the assessment of mutant effect on the structure network using `MutationGroupAnalyzer`.

