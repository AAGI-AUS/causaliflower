# causaliflower

The causaliflower package provides simple ways of working with causal graphs in R.

It makes use of existing functions for analysing and visualising directed acylic graphs (DAGs), including 'dagitty' and 'ggdag'. 

>What do you get when you can't think of a name for a causality + farming related R package? 

Causality + cauliflower = 'causaliflower' (working title - please help!)


## Overview
The central focus of causaliflower is to provide open-source software for building causal graphs and analysing their structures in a cohesive workflow.

The code below gives some idea of what this looks like, however an introductory vignette will be provided on release including an in-depth example.


### In This Version of causaliflower


- build a directed (hopefully acyclic) graph from scratch with flexible inputs (as few or many as needed).

```R

dag <- build_graph(type = 'ordered', variables, treatments, outcomes, # required inputs
                   mediators, latent_variables, instrumental_variables, mediator_outcome_confounders, competing_exposures, colliders) # optional inputs

```


- draw edges between nodes in a dagitty objectto create a saturated or fully connected graph

```R
saturated_graph <- build_graph(type = c('full', 'saturated'), # choose a type
                               variables = dag) # input existing dagitty object
```


- assess the implied causal relationships using causal criteria (deciding which edges to keep)

```R
edges <- assess_edges(saturated_graph, edges_to_keep = dag, ) 

edges_to_keep <- assess_edges(saturated_graph, edges_to_keep = dag, assess_causal_criteria = TRUE)

dag <- keep_edges(saturated_graph, edges_to_keep)
```


- update graph with new coordinates 

```R
dag <- add_coords(dag) # also called in build_graph()
```


- merge two dagitty objects

```R
new_dag <- build_graph(type = 'ordered', variables, treatments, outcomes, # required inputs
                       mediators, latent_variables, instrumental_variables, coords_spec = 2) # optional inputs
            
new_dag <- merge_graphs(dag, new_dag)
```


- simplify representation using markov principles to remove edges (a simple algorithm, but may produce some unexpected results)

```R
markov_dag <- markov_graph(dag)
```


- obtain minimally sufficient adjustment sets (returns smallest 5 sets by default)

```R
minimal_sets(dag, effect = "direct")
minimal_sets(new_dag, effect = "direct")
minimal_sets(markov_dag, effect = "direct")
```


- display dagitty graphs using ggdag (causaliflower presets)

```R
ggdagitty(dag)
ggdagitty(new_dag)
ggdagitty(markov_dag)
```


- add nodes at different time-points using existing nodes in a dagitty object 

```R
add_nodes(dag, existing_nodes)
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

