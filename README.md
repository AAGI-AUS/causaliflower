# causaliflower

The causaliflower package provides simple ways of working with causal graphs in R.

It makes use of existing functions for analysing and visualising directed acylic graphs (DAGs), including 'dagitty' and 'ggdag'. 

>What do you get when you can't think of a name for a causality + farming related R package? 

Causality + cauliflower = 'causaliflower' (working title - please help!)


## Overview
The central focus of causaliflower is to provide open-source software for building causal graphs and analysing their structures in a cohesive workflow.

The code below gives some idea of what this looks like, however an introductory vignette will be provided on release including an in-depth example.


### In This Version of causaliflower


- build a directed acyclic graph from scratch, with minimal input    

```R
buildGraph(type = 'ordered', variables, treatments, outcomes)
```

- saturate (fully connect) a graph before using causal criteria to assess each implied relationship and decide which edges to keep (or remove)

```R
buildGraph(type = c('full', 'saturated'), variables, treatments, outcomes)

assessEdges(dag)

keepEdges(dag, edges_to_keep)
```

- generate coordinates

```R
getCoords(dag) - also automatically generated with buildGraph()
```

- add nodes to a dagitty object at different time-points, based on existing node names

```R
addNodes(dag, reference_nodes)
```

- extract edges and node roles in a dagitty object    

```R
getEdges(dag)
```

- other functions for extracting node names and roles 

```R
getFeatureMap(dag)  
variables(dag)  
colliders(dag)  
competing_exposures(dag) 
mediator_outcome_confounders(dag) 
confounders(dag)  
mediators(dag)  
instrumental_variables(dag) 
proxies(dag)
```


In the near future causaliflower will be made public so that anyone interested in causal analysis can use the package.

If you have any questions or suggestions, please get in touch and/or contribute to the development of causaliflower.

