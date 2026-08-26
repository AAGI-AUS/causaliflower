# causaliflower

>causality + cauliflower = 'causaliflower'!

Graph-based causal analysis in R.
The `causaliflower` package extends `dagitty` and `ggdag` with functions for building and assessing directed acyclic graphs (DAGs).

The goal of `causaliflower` is to support reproducible causal analytical workflows,
with a focus on applied causal inference techniques in agriculture.

### Installation

The most recent version can be installed from [GitHub](https://github.com/AAGI-AUS/causaliflower):

```R
if (!require("pak")) {
  install.packages("pak")
}

pak::pak("AAGI-AUS/causaliflower")

```

Some examples are included below, and an introductory vignette will be provided in an upcoming release.


### Example code

- Build a basic graph (dagitty object)

```R

variables <- c("Z3", "Z2", "Z1")
treatments <- "X"
outcomes <- "Y"


dag <- build_graph(variables = variables,
                   treatments = treatments,
                   outcomes = outcomes,
                   type = "ordered")

```

- Plot dagitty objects:

```R

plot_dagitty(dag)

# Add coordinates to a dagitty object (more detail below)
dag <- add_coords(dag)

```

- Generate new coordinates:

```R

# Changing the input parameters affects how coordinates are generated
dag <- add_coords(dag,
                  x_step = 2, # default horizontal spacing between temporal ranks
                  y_step = 1) # default vertical spacing between nodes within a rank

```

- Connect graph edges in a fully connected or saturated graph:

```R

fc_graph <- connect_nodes(dag) # default connects all nodes in both directions (type = "full")
fc_graph |> plot_dagitty()

saturated_graph <- connect_nodes(dag, type = "saturated", print_edges = TRUE) # saturated graph connects earlier nodes to later nodes

```


- Assess graph edges using causal criteria:

```R
## Assess edges to keep and build a new graph
new_graph <- fc_graph |>
  assess_edges(edges_to_keep = dag,
               assess_causal_criteria = TRUE) |> # guided causal criteria sequence
  keep_edges(dag = fc_graph)


plot_dagitty(new_graph)

## Or, save answers to causal criteria in a list
edges_list <- fc_graph |>
  assess_edges(edges_to_keep = dag,
               assess_causal_criteria = TRUE,
               save_answers = TRUE) # saves answers to causal criteria, output becomes a list

edges_list$edges
edges_list$answers

```


- Join two dagitty objects, keeping the coordinates of the first:

```R

mediators <- "M"

new_dag <- build_graph(treatments = treatments,
                       outcomes = outcomes,
                       mediators = mediators)

new_dag <- add_coords(new_dag)
plot_dagitty(new_dag)

dag <- join_graphs(dag, new_dag)
plot_dagitty(dag)

```


- Output minimal sufficient adjustment sets (returns smallest 5 sets by default):

```R

minimal_sets(dag, effect = "direct")

```


- Add nodes to an existing graph:

```R

new_nodes <- c("Z4", "Z5")
descendants <- names(dag)

dag <- add_nodes(dag, new_nodes, descendants = descendants)
plot_dagitty(dag)

```

- Place new nodes by causal role and temporal position:

```R

dag <- add_nodes(
  dag,
  new_nodes = "W",
  node_role = "confounder",
  type = "saturated",
  position = first()
)
plot_dagitty(dag)

```

`position` takes four placement helpers, used alone or combined with `c()`:

```R

first()       # before every existing node
last()        # after every existing node
before("X")   # immediately before an existing anchor node
after("Z2")   # immediately after an existing anchor node

# Name the new nodes assigned to each position clause.
dag <- add_nodes(
  dag,
  new_nodes = c("W1", "W2", "U"),
  node_role = "confounder",
  type = "saturated",
  position = c(first("W1"), after("Z2", "W2"))
)

```

- Get edges and node structure information from a dagitty object:

```R
get_edges(dag)

get_ancestor_edges(dag)

get_structure(dag)

get_nodes(dag)

get_roles(dag)

get_diff_roles(dag, new_dag)

get_diff_edges(dag, new_dag)
```


- Other utility functions:

```R
confounders(dag)

mediators(dag)

instruments(dag)

colliders(dag)

competing_causes(dag)

mediator_outcome_confounders(dag)

proxies(dag)
```

If you have any questions, suggestions, or would like to contribute, please let me know!
