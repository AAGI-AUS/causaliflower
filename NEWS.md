# causaliflower 0.0.2

- Adds `collapse_nodes()` to marginalise named nodes, replacing their paths
  with directed or bidirected edges as required by the graph.
- Adds `merge_nodes()` to combine an explicitly selected group into one named
  node.
- Adds `recombine_nodes()` to find and merge groups with matching connections
  in one pass. Its `preview` argument returns the proposed groups and names
  for inspection or selective use with `merge_nodes()`.
- `print_edges` now defaults to `FALSE` in `add_nodes()`, `connect_nodes()`,
  and `copy_nodes()`, as it does in the new collapse and merge functions.
- Adds `remove_nodes()` to remove named nodes and all their incident edges.
  Other nodes are kept, including neighbours left without edges.
- `remove_edges()` now removes nodes left without edges by the current call,
  with a warning naming them. Nodes that were already isolated are kept.
  `keep_edges()` continues to retain every node.
- Graph construction now distinguishes `full`, `chain`, `first`, `last`,
  `ordered`, and `saturated`. Vectors give successive temporal positions,
  and each list element can group nodes at the same position. `chain`
  connects adjacent positions within each role except instrumental variables
  and adds no arrows between roles apart from explicit latent targets and
  instrument edges. `last` adds arrows between roles from the last position
  of the sending role to the first position of the receiving role, and
  `first` sends them from the first position instead. `ordered` connects
  every earlier position to every later position within each role except
  instrumental variables, and `saturated` also adds the directed arrows
  allowed by the package's order of causal roles. `first` and `last` are
  `build_graph()` types. `add_nodes()` and `connect_nodes()` use `full`,
  `chain`, `ordered`, and `saturated`. Instrumental variable inputs and
  generated edges are unchanged.
- With `type = "full"`, `connect_nodes()` adds bidirected edges between
  nodes other than instruments when at least one endpoint is selected. In
  `add_nodes()`, a `role` other than `"instrument"` connects the new nodes
  to each other and to existing nodes other than instruments unless
  `parents` or `children` selects the connections. Without a `role` or
  explicit targets, the new nodes connect to every existing node.
- In `add_nodes()`, `parents` and `children` can be combined with the
  `"treatment"`, `"outcome"`, or `"latent"` role to set the node status
  without generating edges from that role. Other roles cannot be combined
  with these inputs. For example, a new latent node can be a common cause of
  two specified nodes without an arrow to the outcome.
- Adds `incoming()` and `outgoing()` so role inputs can select which nodes
  receive or send generated arrows between roles.
- Role inputs also accept `incoming` and `outgoing` as names, and `position`
  accepts `"first"`, `"last"`, `c(before = "M2")`, and `c(after = "M2")`.
- `add_nodes()` now uses `nodes`, `parents`, `children`, and `role` as its main
  argument names. The older `new_nodes`, `ancestors`, `descendants`, and
  `node_role` names remain as deprecated aliases. `connect_nodes()` now takes
  `dag`, `nodes`, `type`, and `print_edges`, in that order. Named `role`,
  `node_role`, and `position` arguments warn and are ignored.
- Adds `set_treatment()`, `set_outcome()`, and `set_latent()`, which add a node
  status and change nothing else.
- `build_graph()` now uses `latents` for latent node names. `latent_variables`
  remains as a deprecated alias.
- Deprecation messages state the version in which the deprecation took effect.
- `competing_exposures()` is a deprecated alias of `competing_causes()`: it
  warns that the name is deprecated as of 0.0.1 and returns the same result
  instead of stopping.
- Adds `copy_nodes()` with sequential and parallel copies, optional
  suffixes, inherited edges, and placement relative to an anchor. The anchor
  joins one end of the copies, and `copy_nodes(type = "chain")` places the
  copies in the anchor's path instead of inheriting the original's edges.
  Its `renew_coords` argument defaults to `TRUE`, which regenerates every
  node's coordinates. `FALSE` places only the copies.
- `add_nodes()` and `join_graphs()` also accept `renew_coords`, which defaults
  to `FALSE` to keep existing coordinates and place only the new nodes. In
  all three functions, `FALSE` falls back to a full layout with a warning if
  a new node coincides with another node.
- Node placement keeps the minimum vertical gap when it levels the ends of a
  bidirected or undirected edge, so a new node is no longer placed on an
  existing one.
- `add_coords()` now applies `x_step` to graphs without edges.
- Adds `get_paths()` and `get_node_edges()`.
- `get_structure()` now returns a nested list with each node's role, parents,
  children, ancestors, descendants, undirected and bidirected connections,
  and directed paths. Its new `max_paths = 100L` limit applies separately to
  each node's incoming and outgoing paths. `max_length = Inf` leaves path
  length unlimited. Truncation is recorded and reported when printed.
- `get_ancestor_edges()` is deprecated in favour of `get_node_edges()`. Use
  `direction = "any"` to include undirected and bidirected edges.
- Graph names and bounding box values are no longer mistaken for node names
  by status accessors and graph inspection functions.
- Adds `get_identification()` to check whether causal effects are identifiable
  in `dagitty` DAGs.
- With `effect = "direct"`, `get_identification()` also checks the standard
  mediation formula for natural direct and indirect effects. Use `mediator`
  to select a mediator or joint mediator set. These checks report their
  additional assumptions separately from the controlled direct effect verdict.
- Adds `add_edges()` for adding specified edges between nodes already present
  in a graph, and `remove_edges()` for removing specified edges while
  preserving the remaining graph.
- `build_graph()` and `add_nodes()` now warn when the generated edges contain
  a directed cycle, as `connect_nodes()` does. `add_edges()` warns when the
  supplied graph already contains a cycle and rejects additions that would
  create a cycle in an acyclic graph.
- Checks for directed cycles in the graph construction and editing functions
  and in `get_identification()` no longer follow every directed path, which
  took minutes on dense graphs. `get_identification()` also narrows its
  front door and mediation searches before enumerating candidate sets.
- `add_nodes()`, `connect_nodes()`, `add_edges()`, `remove_edges()`,
  `keep_edges()`, and `copy_nodes(type = "chain")` now keep the graph name
  and bounding box, and the labels, notes, and other native attributes of the
  remaining nodes and edges.
- `add_edges()`, `keep_edges()`, and `remove_edges()` use `edges` for their
  edge input. `assess_edges()` uses `keep_edges` and `remove_edges` for edges
  retained or discarded without assessment. The older `edges_to_keep`
  argument remains as a deprecated alias in both `keep_edges()` and
  `assess_edges()`.
- When selection or filtering leaves no edges to assess, `assess_edges()`
  reports "No edges remain to assess." and returns normally. An input graph
  with no edges still raises an error.
- `add_edges()`, `keep_edges()`, and `remove_edges()` warn and ignore missing
  edge strings or incomplete rows in edge tables. An input with no complete
  edges still raises an error.
- Adds `add_curves()` for curving directed arrows in a `dagitty` object, with
  relative or absolute height and automatic or specified placement.
- `plot_dagitty()` automatically curves directed arrows to avoid nodes and
  reduce crossings. Use `curve = FALSE` to retain existing edge paths
  without recalculating them.
- `plot_dagitty()` has a new presentation. Labels sit just above their nodes,
  or below them, on the roomier side (`"auto"`) or inside them, set by
  `label_placement`. Latent nodes are hollow, every role present is listed in
  the legend, edges are lighter and bidirected edges are dashed. Node, label,
  arrow and legend sizes use physical units. Use `node_size`, `label_size`
  and `arrow_size` to change node, label and arrowhead sizes, and `edge_gap`
  to adjust the space between edges and nodes. The margins grow to fit the
  labels.
- `plot_dagitty()` places the legend with one argument, `legend_position`,
  such as `"top"` (the default), `"right-top"` or `"none"`.
- `plot_dagitty()` can save plots as SVG (the default) or PNG using `file`.
  File dimensions are fitted to the graph by default. svglite and ragg are
  used when installed.
- Automatic curves in `plot_dagitty()` have bounded heights, matching
  `add_curves(mode = "absolute")` with a height set by the gaps between
  nodes. Generated layouts default to wider, evenly spread columns
  (`x_step = 3`, `y_step = 1.5`, `layout_shape = 1`).
- `plot_dagitty()`'s `seed` argument and `label_placement = "text_box"` are
  deprecated as of 0.0.2.
- Adds an introductory vignette covering graph construction, inspection, and
  identification.
- Adds runnable examples for graph construction, inspection, `get_labels()`,
  and `string_to_dag()`.
- Revises wording in the help pages, README and introductory vignette.

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
  directed cycle. Under `type = "ordered"`, edge rules that involve the
  outcome now use every declared outcome rather than only the first.
- Preserves isolated nodes, role declarations, coordinates, and native edge
  types when graphs are rebuilt; reciprocal directed edges are no longer
  collapsed into a bidirected ("<->") edge. Uses one parser and serialisation
  path for graph text and retained edges, and supports node names containing
  spaces, hyphens, and percent signs. Adds `string_to_dag()` and the temporal
  position helpers `first()`, `last()`, `before()`, and `after()`.

## Coordinates

- Replaces the random coordinate engine with a deterministic layout. Root nodes
  are placed on the left. Every other node is placed one step to the right of
  its latest parent, and nodes sharing a horizontal position are spaced
  vertically. The same graph always receives the same coordinates. Cyclic
  graphs are laid out using an acyclic subset of their edges, with a warning.
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
