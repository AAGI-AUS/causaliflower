# causaliflower

The causaliflower package provides functions for working with causal graph objects in R.

It extends 'dagitty' and 'ggdag' for analysing and visualising directed acylic graphs (DAGs) to provide repeatable workflows to simplify the process of building a DAG from prior knowledge.


## Overview

###  What is causaliflower?

>Causality + cauliflower = 'causaliflower' (working title!)

The aim of causaliflower is to enable repeatable workflows for building and analysing causal graphs in R.

Example code is included below, however an upcoming vignette will provide a more in-depth guide to some of the functions in this package.


### In This Version


- Build a dagitty object, inputting node roles (if known):

```R

dag <- build_graph(type = 'ordered', variables, treatments, outcomes) # required inputs

```


- Saturate dagitty object nodes, up to a fully connected graph:

```R
saturated_graph <- build_graph(type = c('full', 'saturated'), # choose a type
                               variables = dag) # existing dagitty object inputted
```


- Assess implied causal relationships to remove edges using a set of causal criteria:

```R
edges <- assess_edges(saturated_graph, edges_to_keep = dag) 

edges_to_keep <- assess_edges(saturated_graph, edges_to_keep = dag, 
                      assess_causal_criteria = TRUE) # guided causal criteria sequence

dag <- keep_edges(saturated_graph, edges_to_keep)
```


- Add or update node coordinates:

```R
dag <- add_coords(dag) # also called in build_graph()
```


- Join two dagitty objects, keeping the coordinates of the first:

```R
new_dag <- build_graph(type = 'ordered', variables, treatments, outcomes, # required inputs
                       mediators, latent_variables, instrumental_variables, # more inputs
                       coords_spec = 2) # higher spec increases volatility in node placement
            
new_dag <- join_graphs(dag, new_dag)
```


- Output minimal sufficient adjustment sets (returns smallest 5 sets by default):

```R
minimal_sets(dag, effect = "direct")
```


- Display dagitty objects using ggdag (wrapper function with custom causaliflower presets):

```R
plot_dagitty(new_dag)
```


- Add nodes to an existing graph:

```R
add_nodes(dag, new_nodes)

saturate_nodes(dag, new_nodes)
```



- Create copies of nodes, e.g. occurring at different time-points:

```R
copy_nodes(dag, existing_nodes)
```


- Functions to get edges and node structure information from a dagitty object:

```R
get_edges(dag)

get_roles(dag)  

get_roles(dag)

get_structure(dag)
```


- Other utility functions:

```R
colliders(dag)  

competing_exposures(dag) 

mediator_outcome_confounders(dag) 

confounders(dag)  

mediators(dag)  

instrumental_variables(dag) 

proxies(dag)
```

If you have any questions, suggestions, or would like to contribute, please let me know!
