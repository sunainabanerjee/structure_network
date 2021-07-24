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




