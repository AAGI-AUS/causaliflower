# causaliflower 0.0.1

Causal graph core release.

## Public interface

- Adds the status accessors `treatments()`, `unobserved()`, `observed()`, and
  `outcomes()` (with an `outcomes(dag) <- value` replacement form delegating to
  dagitty).
- Renames `competing_exposures()` to `competing_causes()`. The new name is
  used throughout: the role label in `get_roles()`, `get_edges()`,
  `get_nodes()`, and `get_structure()` output, plot legends, and
  `build_graph()`'s `competing_causes` argument. Calling
  `competing_exposures()` stops with a pointer to the new name.
- Removes `get_nodes_from_treatment_to_outcome()` and makes
  `nodes_between_treatment_and_outcome()` internal; `mediators()` is the
  supported way to ask which nodes carry the treatment effect.

## Behaviour

- `connect_nodes()` now honours its `nodes` selection: only edges involving a
  selected node are added.
- Repairs the `node_role` pathway in `add_nodes()`, which previously errored
  for every role; a `node_role` cannot be combined with explicit
  `ancestors`/`descendants`.
- `assess_edges(edges_to_assess = "bidirectional")` now filters to bidirected
  ("<->") and undirected ("--") edges in every mode, not only within the
  causal criteria sequence.
- `minimal_sets()` validates `num_sets` as a positive whole number, gives a
  clearer error for an invalid graph, and warns only when a same length
  adjustment set has actually been excluded.
- Role classification uses directed reachability: mediators are nodes on
  directed treatment-to-outcome paths, and confounders include ancestral
  common causes; partially directed graphs report an `undetermined` role.
- `join_graphs()` warns when the two graphs disagree on a shared node's
  exposure/outcome/latent declarations and when the merged graph contains a
  directed cycle. Under `type = "ordered"`, edge templates that involve the
  outcome now use every declared outcome rather than only the first.
- Preserves isolated nodes, role declarations, coordinates, and native edge
  types when graphs are rebuilt; reciprocal directed edges are no longer
  collapsed into a bidirected ("<->") edge. Uses one parser and serialisation
  path for graph text and retained edges, and supports node names containing
  spaces, hyphens, and percent signs. Adds `string_to_dag()` and the temporal
  position helpers `first()`, `last()`, `before()`, and `after()`.

## Coordinates

- Replaces the random coordinate engine with a deterministic layout: nodes are
  ranked left to right by longest directed path and spaced within each rank.
  The same graph always receives the same coordinates. Cyclic graphs are laid
  out on an acyclic subset of their edges, with a warning.
- `add_coords()` and `plot_dagitty()` take `x_step`, `y_step`, and
  `layout_shape` in place of `coords_spec` and `threshold`;
  `add_coords_helper()` and the layout time limits are removed. Joining graphs
  keeps the first graph's coordinates and places only genuinely new nodes.

## Package maintenance

- Removes `remote_sensing_funcs.R` (`get_rs_data()`, `get_rs_raster()`, `crop_search()`,
  `crop_search_defaults()`, `point_to_bbox()`, `metres_to_degrees()`,
  and `validate_data_fusion()`, now in a dedicated package: <https://github.com/AAGI-AUS/ocular>)
  and an `rstac` dependency.
- Adds automated correctness and regression tests plus a cross-platform
  R CMD check workflow.
- Corrects package and MIT licence metadata, including the ESC-DAG Protocol CC BY 4.0 attribution
  (causal criteria citation corrected to Glass et al., 2013).
- R 4.1 or later required.
