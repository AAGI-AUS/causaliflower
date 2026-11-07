# causaliflower

> causality + cauliflower = 'causaliflower'!

Graph-based causal analysis in R. The `causaliflower` package extends `dagitty` and `ggdag` with functions for building and assessing directed acyclic graphs (DAGs).

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

Development version 0.0.2.9000 requires `dagitty` exactly 0.3.4.
Installing from source also requires a C++17 compiler.

Examples are included below. If the introductory vignette is installed,
open it with `vignette("causaliflower-introduction")`.


### DAG-by-query

`get_dag()` constructs a DAG by querying the bundled structural causal model
derived from APSIM. The slice contains 576 quantities with readable names
across barley, wheat and maize. Supply the crop, treatment and outcome, and state
whether the treatment was randomised.

```R
library(causaliflower)

crop_dag <- get_dag(
  crop = "barley",
  treatments = "nitrogen_applied",
  outcomes = "grain_protein",
  randomised_treatment = FALSE,
  latents = NULL,
  quiet = TRUE
)

plot_dagitty(crop_dag)
```

Here `latents = NULL` assumes all quantities are observed for this example.
The ranked search reports its step budget and completion
limits, including when larger adjustment options may remain. `quiet = TRUE`
suppresses routine progress, while these warnings remain visible.

Inspect the saved query, ranked options and construction checks with
`original_query(crop_dag)`, `dag_options(crop_dag)` and `dag_checks(crop_dag)`.
`dag_history(crop_dag)` records the transformations. The selected adjustment
sets and the catalogue of alternatives are checked separately.


### Example code

- Build a basic graph (dagitty object)

```R

variables <- c("X3", "X2", "X1")
treatments <- "T"
outcomes <- "Y"


dag <- build_graph(variables = variables,
                   treatments = treatments,
                   outcomes = outcomes,
                   type = "ordered")

```

`full` uses bidirected edges following the ESC-DAGs Mapping Stage in
Ferguson et al. (2020). `chain` connects adjacent temporal positions within
each role except instrumental variables and adds no arrows between roles
apart from explicit latent targets and instrument edges. `last` adds arrows
between roles from the last position of the sending role to the first
position of the receiving role, and `first` sends them from the first
position instead. `ordered` connects every earlier position to every later
position within each role except instrumental variables, and `saturated`
also adds the directed arrows allowed by the package's order of causal
roles. `first` and `last` are `build_graph()` types. `add_nodes()` and
`connect_nodes()` use `full`, `chain`, `ordered`, and `saturated`.

Use a list element to place nodes at the same temporal position:

```R

dag <- build_graph(
  variables = variables,
  treatments = treatments,
  outcomes = outcomes,
  mediators = list("M1", c("M2", "M3"), "M4"),
  type = "last"
)

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
                  x_step = 2, # horizontal spacing between node columns
                  y_step = 1) # vertical spacing between nodes in the same column

```

- Add edges to existing nodes:

```R

full_graph <- connect_nodes(dag, type = "full")
full_graph |> plot_dagitty()

saturated_graph <- connect_nodes(dag, type = "saturated", print_edges = TRUE)

```


- Assess graph edges using causal criteria:

```R
## Assess edges to keep and build a new graph
new_graph <- full_graph |>
  assess_edges(keep_edges = dag,
               assess_causal_criteria = TRUE) |> # guided causal criteria sequence
  keep_edges(dag = full_graph)


plot_dagitty(new_graph)

## Or, save answers to causal criteria in a list
edges_list <- full_graph |>
  assess_edges(keep_edges = dag,
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
                       mediators = mediators,
                       type = "last")

new_dag <- add_coords(new_dag)
plot_dagitty(new_dag)

dag <- join_graphs(dag, new_dag)
plot_dagitty(dag)

```


- Check if the effects are identifiable:

```R

get_identification(dag)

```

The result is grouped by effect. Identification verdicts are shown separately
from the strategies. Each strategy shows its method, assumptions, variables
and identification status. Use `id_status = "all"` to include unsuccessful
strategy checks. Identification depends on the supplied graph and the
assumptions stated with each result.


- Find single door sets for the direct path coefficient in a linear model (returns up to five sets, ordered by size):

```R

minimal_sets(dag, effect = "direct")

```


- Add nodes to an existing graph:

```R

nodes <- c("X4", "X5")
children <- names(dag)

dag <- add_nodes(dag, nodes, children = children)
plot_dagitty(dag)

```

- Place new nodes by causal role and temporal position:

```R

dag <- add_nodes(
  dag,
  nodes = "W",
  role = "confounder",
  type = "saturated",
  position = first()
)
plot_dagitty(dag)

```

`position` takes four placement helpers, used alone or combined with `c()`.
Their names also work as values, as in `position = "first"` or
`position = c(after = "X2")`:

```R

first()       # before every existing node
last()        # after every existing node
before("T")   # immediately before an existing anchor node
after("X2")   # immediately after an existing anchor node

# Choose a position for each named new node.
dag <- add_nodes(
  dag,
  nodes = c("W1", "W2", "U"),
  role = "confounder",
  type = "saturated",
  position = c(first("W1"), after("X2", "W2"))
)

```

- Choose which nodes receive and send arrows generated between roles:

```R

dag <- build_graph(
  variables = variables,
  treatments = treatments,
  outcomes = outcomes,
  mediators = c(incoming("M1"), "M2", outgoing("M3")),
  type = "ordered"
)

```

- Copy a node into successive positions:

```R

dag <- copy_nodes(dag, "M2", copy_suffix = c("survey", "harvest"))

```

- Get edges and node structure information from a dagitty object:

```R
get_edges(dag)

get_node_edges(dag)

get_paths(dag)

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
