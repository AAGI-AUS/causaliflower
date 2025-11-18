# causaliflower

This package aims to provide intuitive methods for building causal graphs, assessing edges, and creating useful workflows for causal analysis in R.

The causaliflower package makes use of existing functions for analysing and visualising directed acylic graphs (DAGs). 

Wrapper functions in causaliflower extend the already impressive functionality of dagitty and ggdag packages. 

Others assist with causal inference tasks by working directly with the underlying causal structure of causal diagrams in R. 

This aims to allow researchers, applied statisticians, students, and anyone else interested in causal diagrams to integrate a causal analysis as part of their workflows, or where it is useful in a variety of causal inference tasks.


## Overview

>What do you get when you can't think of a name for a causality + farming related R package? 

Causality + cauliflower - 'causaliflower' (working title - please help!) - provides causal analysis functions, built on top of 'dagitty' and 'ggdag' in R. 

The below examples provide a brief introduction to some of the functions available in causaliflower (working title - please help!) for assessing causal relationships in causal diagrams. For a more detailed example, see the following section for using causaliflower in a causal analytical workflow.


## In This Version of causaliflower (working title - please help!)

Currently, causaliflower gives the ability to:

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

## More About causaliflower

The aim of causaliflower is to assist with causal inference tasks by providing ways of working directly with causal diagrams in R, and integration with existing workflows.

If you have any suggestions about how it can do this, or do it better, please get in touch and/or contribute to its development.

