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


- Build and a basic graph (dagitty object)

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

- Generate new coordinates

```R

# Changing the input parameters affects how coordinates are generated
dag <- add_coords(dag,
                  coords_spec = 0.2, # >1 increases node placement volatility
                  threshold = 0.9) # controls the spacing/closeness allowed between nodes

```

- Connect graph edges in a fully connected or saturated graph:

```R

fc_graph <- connect_nodes(dag) # default connects all nodes in both directions (type = "full")
fc_graph |> plot_dagitty()

saturated_graph <- connect_nodes(dag, type = "saturated", print_edges = TRUE) # saturated graph connects earlier nodes to default connects all nodes in both directions

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
               
edges_list$edges_to_keep
edges_list$answers
  
```


- Join two dagitty objects, keeping the coordinates of the first:

```R

mediators <- "M"

new_dag <- build_graph(treatments = treatments, 
                       outcomes = outcomes, 
                       mediators = mediators)
                       
new_dag <- add_coords(new_dag, coords_spec = 0.9)
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


- Functions to get edges and node structure information from a dagitty object:

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

competing_exposures(dag) 

mediator_outcome_confounders(dag) 

proxies(dag)

get_nodes_from_treatment_to_outcome(dag)
```


- Obtain remote sensing data from Sentinel 2:

```R
bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, 
                      ymin = -34.902, ymax = -34.906),
                    crs = sf::st_crs(4326))
## General function for obtaining time series data, NDVI example:
ndvi_df <- getS2_data(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001,
                      asset_names = s2_index$NDVI$assets,
                      index_function = s2_index$NDVI$fun)
                      
## You can also use a convenience function:
ndvi_df <- getNDVI(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
```


- Raster composite of Sentinel 2 data:

```R
bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, 
                      ymin = -34.902, ymax = -34.906),
                    crs = sf::st_crs(4326))
                    
## General function, NDVI:
ndvi_rast <- getS2_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001,
                          asset_names = s2_index$NDVI$assets,
                          index_function = s2_index$NDVI$fun)
## Convenience function, NDVI:
ndvi_rast <- getNDVI_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)

## Plot using terra:
terra::plot(ndvi_rast)

```

If you have any questions, suggestions, or would like to contribute, please let me know!
