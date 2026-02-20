# causaliflower

The causaliflower package aims to support causal analysis in R. It extends 'dagitty' and 'ggdag' functions for analysing and visualising directed acylic graphs (DAGs) to enable building and assessing causal DAGs using expert knowledge.

### Installation

Causaliflower is currently in development. The most recent version can be installed from [GitHub](https://github.com/AAGI-AUS/causaliflower) with:

```R
o <- options() # store original options

options(pkg.build_vignettes = TRUE)

if (!require("pak")) {
  install.packages("pak")
}

pak::pak("AAGI-AUS/causaliflower")
options(o) # reset options

```

###  What is causaliflower?

>Causality + cauliflower = 'causaliflower' (working title!)

The goal of this package is to provide functions for reproducible causal analytical workflows in R.

Some examples are included below, however an in-depth tutorial will be provided in an upcoming vignette.


### Example code


- Build and plot a basic graph (dagitty object)

```R
variables <- c("Z3", "Z2", "Z1")
treatments <- "X"
outcomes <- "Y"

dag <- build_graph(variables = variables,
                   treatments = treatments,
                   outcomes = outcomes)
```

- Plot dagitty objects (ggdag wrapper with custom causaliflower presets):

```R

plot_dagitty(dag)

```

- Create a fully connected graph 

```R

fc_graph <- saturate_nodes(dag)
plot_dagitty(fc_graph)

```


- Assess implied causal relationships to remove edges using a set of causal criteria:

```R

edges <- assess_edges(fc_graph, edges_to_keep = dag,
                      assess_causal_criteria = TRUE) # guided causal criteria sequence

dag <- keep_edges(fc_graph, edges$edges)
plot_dagitty(dag)

```

- Join two dagitty objects, keeping the coordinates of the first:

```R

mediators <- "M"

new_dag <- build_graph(treatments = treatments, 
                       outcomes = outcomes, 
                       mediators = mediators)
plot_dagitty(new_dag)

merged_dag <- join_graphs(dag, new_dag)
plot_dagitty(merged_dag)

```


- Output minimal sufficient adjustment sets (returns smallest 5 sets by default):

```R

minimal_sets(dag, effect = "direct")

```


- Add nodes to an existing graph:

```R

new_nodes <- c("Z4", "Z5")
descendants <- names(merged_dag)

new_dag <- add_nodes(merged_dag, new_nodes, descendants = descendants)
plot_dagitty(new_dag)

```


- Functions to get edges and node structure information from a dagitty object:

```R
get_edges(dag)

get_nodes(dag)  

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
