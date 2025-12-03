# causaliflower

The causaliflower package provides simple ways of working with causal graphs in R.

It makes use of existing functions for analysing and visualising directed acylic graphs (DAGs), including 'dagitty' and 'ggdag'. 


## Overview

###  What is causaliflower?

>Causality + cauliflower = 'causaliflower' (working title!)

The central focus of causaliflower is to provide open-source software for building causal graphs and analysing their structures in a cohesive workflow.

The code below gives an idea, however an introductory vignette provided on release will include more in-depth examples.


### In This Version of causaliflower


- build a directed (hopefully acyclic) graph from scratch with flexible inputs (as few or many as needed).

```R

dag <- build_graph(type = 'ordered', variables, treatments, outcomes) # required inputs

```


- draw edges between nodes in a dagitty objectto create a saturated or fully connected graph

```R
saturated_graph <- build_graph(type = c('full', 'saturated'), # choose a type
                               variables = dag) # existing dagitty object inputted
```


- assess the implied causal relationships using causal criteria (deciding which edges to keep)

```R
edges <- assess_edges(saturated_graph, edges_to_keep = dag) 

edges_to_keep <- assess_edges(saturated_graph, edges_to_keep = dag, 
                      assess_causal_criteria = TRUE) # guided causal criteria sequence

dag <- keep_edges(saturated_graph, edges_to_keep)
```


- update graph with new coordinates 

```R
dag <- add_coords(dag) # also called in build_graph()
```


- merge two dagitty objects

```R
new_dag <- build_graph(type = 'ordered', variables, treatments, outcomes, # required inputs
                       mediators, latent_variables, instrumental_variables, # more inputs
                       coords_spec = 2) # higher spec increases volatility in node placement
            
new_dag <- merge_graphs(dag, new_dag)
```


- simplify graph using markov properties (removes unnecessary edges)

```R
markov_dag <- markov_graph(dag) # simple algorithm, may produce some unexpected results
```


- obtain minimally sufficient adjustment sets (returns smallest 5 sets by default)

```R
minimal_sets(dag, effect = "direct")
```


- display dagitty graphs using ggdag (causaliflower presets)

```R
ggdagitty(new_dag)
```


- copy existing nodes to new nodes at different time-points

```R
copy_nodes(dag, existing_nodes)
```


- add nodes, or the saturate connection of existing nodes (grouped by role, e.g., add confounder/treatment/outcome nodes)

```R
add_nodes(dag, new_nodes)

saturate_nodes(dag, new_nodes)
```


- functions for extracting edges and node roles in a dagitty object    

```R
get_edges(dag)
get_roles(dag)  

node_roles(dag)
node_structure(dag)
```


- other dag utility functions

```R
variables(dag)  
colliders(dag)  
competing_exposures(dag) 
mediator_outcome_confounders(dag) 
confounders(dag)  
mediators(dag)  
instrumental_variables(dag) 
proxies(dag)
```


In the near future, causaliflower will be made public so that anyone interested in causal analysis in R can use the code in this package.

If you have any questions or suggestions, please get in touch and/or contribute to the development of causaliflower.

