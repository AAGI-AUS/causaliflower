# causaliflower

The causaliflower package aims to support causal analysis in R. It extends 'dagitty' and 'ggdag' functions for analysing and visualising directed acylic graphs (DAGs) to enable building and assessing causal DAGs using expert knowledge.

## Overview

###  What is causaliflower?

>Causality + cauliflower = 'causaliflower' (working title!)

The goal of this package is to provide functions for reproducible causal analytical workflows in R.

Some examples are included below, however an in-depth tutorial will be provided in an upcoming vignette.


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
