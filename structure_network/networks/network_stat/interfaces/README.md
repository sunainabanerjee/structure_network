# Interfaces

This package defines the basic interfaces required to carry out detailed distribution sampling via network perturbation. 
This module is designed to analyze the effect of network flow property as a function of restricted variation in the edge attribute.
Various components that is needed for the sampling algorithm to work are, 


1. __Selector__
   
   A selector interface perform network substructure selection.
   The selection can be a single network nodes, or network edges (pair of network nodes), 
   or any size of nodes group. The `n_pair` API allows to access the cardinality of the 
   group selection. This interface works like an iterator object.
   
2. __Perturbation__
   
   Perturbation interface carry out the network changes. The base API is designed very 
   generic. The `Perturbation` interface accepts selector, which restrincts the network
   substructure selection on which the network perturbation will be executed. Network 
   perturbation has two steps, update attributes/parameter associated with the network 
   substructure, and then update the network property globally.
   
   `Perturbation` interface can be extended to support a preferential search, viz., MCMC,
    importance sampling, etc.

3. __Lookup__
   
   This Interface is an auxiliary interface to carry out controlled parameter update. 
   `Perturbation` and `Updater` interfaces often uses `Lookup` object for attribute 
   updation.
   
4. __Support__
   
   `Support` interface is an auxiliary interface, often required to update second 
   order global context dependent update. `Lookup` interface updates are generally
   agnostic of global graph structure, and associated with local information only,
   where `Support` updates are global context update, like computing shortest path 
   or cluster coefficient of the local structure.

5. __Updater__
   
   `Updater` interface updates all network attributes, post perturbation. Updater 
    updates the network inplace, creating a copy of the original network is handled
    by `Perturbation` interface, as it internally uses `Updater` to finally  return 
   the updated network. There can be a list of `Updater` object applied to final 
   network. Logic of updater sequence is left on the designer and the usecase. Updaters 
   are applied to the network in the same order as it is registered.
