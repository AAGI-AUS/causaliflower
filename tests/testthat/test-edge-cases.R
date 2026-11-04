test_that("public graph helpers handle edgeless graphs", {
  dag <- dagitty::dagitty("dag { T [exposure] Y [outcome] }")

  added <- suppressMessages(
    add_nodes(dag, "W", role = "confounder", type = "saturated",
              print_edges = FALSE)
  )
  expect_identical(edge_keys(added), c("W -> T", "W -> Y"))

  connected <- suppressMessages(connect_nodes(dag, print_edges = FALSE))
  expect_identical(edge_keys(connected), "T <-> Y")

  joined <- suppressMessages(join_graphs(string_to_dag("A"),
                                         string_to_dag("B")))
  expect_setequal(names(joined), c("A", "B"))
  expect_length(edge_keys(joined), 0L)

  expect_warning(kept <- keep_edges("T -> Y", dag), "not present")
  expect_setequal(names(kept), c("T", "Y"))
  expect_identical(treatments(kept), "T")
  expect_identical(causaliflower:::get_outcome_names(kept), "Y")
  expect_length(edge_keys(kept), 0L)

  expect_warning(
    ancestor_edges <- get_ancestor_edges(dag),
    paste0("get_ancestor_edges(): 'get_ancestor_edges()' is deprecated as of ",
           "0.0.2. Use 'get_node_edges(direction = \"any\")' instead."),
    fixed = TRUE
  )
  expect_length(ancestor_edges, 0L)
  structure <- get_structure(dag)
  expect_setequal(names(structure), c("T", "Y"))
  expect_true(all(vapply(structure, function(node) {
    length(node$parents) == 0L && length(node$children) == 0L &&
      length(node$ancestors) == 0L && length(node$descendants) == 0L
  }, logical(1L))))
})

test_that("rebuilds preserve partial coordinates and node declarations", {
  dag <- string_to_dag(c("A -> B", "C"))
  dagitty::coordinates(dag) <- list(x = c(A = 1, C = 3),
                                    y = c(A = 4, C = 6))
  kept <- keep_edges("A -> B", dag)
  coordinates <- dagitty::coordinates(kept)

  expect_equal(coordinates$x[c("A", "C")], c(A = 1, C = 3))
  expect_equal(coordinates$y[c("A", "C")], c(A = 4, C = 6))
  expect_true(is.na(coordinates$x[["B"]]))
  expect_true(is.na(coordinates$y[["B"]]))

  pdag <- dagitty::dagitty(
    "pdag { A [exposure,adjusted] B [outcome,selected]
            C [latent] A -- B; C -> A }"
  )
  rebuilt <- keep_edges(c("A -- B", "C -> A"), pdag)
  expect_match(as.character(rebuilt), "^pdag ")
  expect_identical(dagitty::exposures(rebuilt), "A")
  expect_identical(dagitty::outcomes(rebuilt), "B")
  expect_identical(dagitty::latents(rebuilt), "C")
  expect_identical(dagitty::adjustedNodes(rebuilt), "A")
  expect_identical(causaliflower:::selected_nodes(rebuilt), "B")

  mag <- dagitty::dagitty(
    "mag { U; T [exposure]; Y [outcome]
           U -> T; U <-> Y; T -> Y }"
  )
  rebuilt_mag <- keep_edges(c("T -> Y", "U -> T", "U <-> Y"), mag)
  expect_match(as.character(rebuilt_mag), "^mag ")
  expect_identical(edge_keys(rebuilt_mag),
                   c("T -> Y", "U -> T", "U <-> Y"))
  expect_s3_class(suppressMessages(plot_dagitty(rebuilt_mag)), "ggplot")
})

test_that("reciprocal directed edges remain a directed cycle", {
  dag <- string_to_dag(c("A -> B", "B -> A"))

  expect_identical(edge_keys(dag), c("A -> B", "B -> A"))
  expect_false(dagitty::isAcyclic(dag))
  expect_identical(sort(paste(get_edges(dag)$ancestor,
                              get_edges(dag)$edge,
                              get_edges(dag)$descendant)),
                   c("A -> B", "B -> A"))
})

test_that("invalid add and connect boundaries fail clearly", {
  dag <- build_graph(c("X1", "X2"), "T", "Y", type = "ordered")

  expect_error(add_nodes(dag, character(), print_edges = FALSE),
               "at least one node name")
  expect_error(add_nodes(dag, NA_character_, print_edges = FALSE),
               "non-missing, non-empty")
  expect_error(add_nodes(dag, "", print_edges = FALSE),
               "non-missing, non-empty")
  expect_error(
    add_nodes(dag, "X1", role = "confounder", print_edges = FALSE),
    "already in the graph"
  )
  expect_warning(
    mixed <- suppressMessages(
      add_nodes(dag, c("X1", "W"), role = "confounder",
                print_edges = FALSE)
    ),
    "skipping"
  )
  expect_true("W" %in% names(mixed))
  expect_error(connect_nodes(dag, "W", print_edges = FALSE),
               "not in the graph")
})

test_that("legacy placement inputs remain compatible", {
  dag <- build_graph(c("X1", "X2"), "T", "Y", type = "ordered")

  expect_warning(
    first_graph <- suppressMessages(
      add_nodes(dag, c("W1", "W2"), role = "confounder",
                type = "first", print_edges = FALSE)
    ),
    paste0("add_nodes(): 'type = \"first\"' is deprecated as of 0.0.1. ",
           "Use 'position = first()' instead."),
    fixed = TRUE
  )
  expect_true(all(c("W1 -> X1", "W2 -> X1") %in% edge_keys(first_graph)))
  expect_false("W1 -> W2" %in% edge_keys(first_graph))

  expect_warning(
    last_graph <- suppressMessages(
      add_nodes(dag, c("W1", "W2"), role = "confounder",
                type = "last", print_edges = FALSE)
    ),
    paste0("add_nodes(): 'type = \"last\"' is deprecated as of 0.0.1. ",
           "Use 'position = last()' instead."),
    fixed = TRUE
  )
  expect_true(all(c("X2 -> W1", "X2 -> W2") %in% edge_keys(last_graph)))
  expect_false("W1 -> W2" %in% edge_keys(last_graph))

  expect_warning(
    expect_error(
      add_nodes(dag, "C", node_role = "competing_exposure",
                print_edges = FALSE),
      "Invalid node_role"
    ),
    "add_nodes(): 'node_role' is deprecated as of 0.0.2. Use 'role' instead.",
    fixed = TRUE
  )
})

test_that("assessment input uses the shared edge parser", {
  dag <- build_graph(c("X1", "X2"), "T", "Y", type = "saturated")

  capture.output(
    reverse <- suppressMessages(
      assess_edges(dag, keep_edges = "T <- X1")
    )
  )
  capture.output(
    compact <- suppressMessages(
      assess_edges(dag, keep_edges = "X1->T")
    )
  )
  expect_identical(reverse, compact)
  expect_false(any(reverse$v == "X1" & reverse$e == "->" &
                   reverse$w == "T"))

  expect_error(assess_edges(dag, keep_edges = "A => B"),
               "could not parse")
  expect_error(assess_edges(dag, edges_to_assess = c("all", "bidirected")),
               "must be one of")
  testthat::local_mocked_bindings(assessment_is_interactive = function() FALSE)
  expect_error(assess_edges(dag, assess_causal_criteria = TRUE),
               "interactive R session")
})

test_that("assessment pre-keeps and pre-removes edges without changing output names", {
  dag <- string_to_dag(c("A -> B", "B -> C", "A <-> C", "C -> D"))

  invisible(capture.output(
    result <- suppressMessages(assess_edges(
      dag,
      keep_edges = "A -> B",
      remove_edges = c("B -> C", "C <-> A")
    ))
  ))
  expect_identical(names(result), c("edges_to_assess", "edges_to_keep"))
  expect_identical(unname(result$edges_to_keep), "A -> B")
  expect_identical(
    paste(result$edges_to_assess$v, result$edges_to_assess$e,
          result$edges_to_assess$w),
    "C -> D"
  )

  expect_warning(
    invisible(capture.output(
      legacy <- suppressMessages(assess_edges(dag, edges_to_keep = "A -> B"))
    )),
    paste0("assess_edges(): 'edges_to_keep' is deprecated as of 0.0.2. ",
           "Use 'keep_edges' instead."),
    fixed = TRUE
  )
  expect_identical(unname(legacy$edges_to_keep), "A -> B")
  expect_error(
    assess_edges(dag, keep_edges = "A -> B", edges_to_keep = "B -> C"),
    "supply either 'keep_edges' or its deprecated alias 'edges_to_keep'",
    fixed = TRUE
  )
  expect_error(
    assess_edges(dag, keep_edges = "C <-> A", remove_edges = "A <-> C"),
    "cannot be supplied in both"
  )
})

test_that("assessment selectors validate membership and complete triage", {
  dag <- string_to_dag(c("A -> B", "B -> C"))

  expect_warning(
    invisible(capture.output(
      missing <- suppressMessages(assess_edges(dag, keep_edges = "X -> Y"))
    )),
    "not present in 'dag'"
  )
  expect_type(missing, "character")

  expect_message(
    complete <- assess_edges(
      dag, keep_edges = "A -> B", remove_edges = "B -> C"
    ),
    "No edges remain to assess"
  )
  expect_length(complete$edges_to_assess, 3L)
  expect_identical(nrow(complete$edges_to_assess), 0L)
  expect_identical(unname(complete$edges_to_keep), "A -> B")

  expect_warning(
    expect_error(
      assess_edges(dag, keep_edges = data.frame(v = NA, e = NA, w = NA)),
      "no complete edges"
    ),
    "missing values"
  )
})

test_that("causal assessment logs are complete and replayable", {
  criteria <- data.frame(
    name = "temporality",
    question = "question",
    description = "description",
    source = "source",
    required = "yes",
    stringsAsFactors = FALSE
  )
  dag <- string_to_dag("A -> B")

  invisible(capture.output(
    all_yes <- testthat::with_mocked_bindings(
      suppressMessages(assess_edges(
        dag,
        assess_causal_criteria = TRUE,
        save_answers = TRUE,
        causal_criteria = criteria
      )),
      assessment_is_interactive = function() TRUE,
      assessment_readline = function(...) "y",
      .package = "causaliflower"
    )
  ))
  expect_identical(names(all_yes$answers),
                   c("v", "e", "w", "temporality"))

  invisible(capture.output(
    replayed <- suppressMessages(assess_edges(
      dag,
      assess_causal_criteria = TRUE,
      save_answers = all_yes$answers,
      causal_criteria = criteria
    ))
  ))
  expect_identical(sort(paste(replayed$v, replayed$e, replayed$w)),
                   edge_keys(dag))

  multi <- string_to_dag(c("A -> B", "B -> C", "C -> D"))
  kept <- data.frame(v = c("A", "B"), e = "->", w = c("B", "C"))
  replay_log <- data.frame(v = "C", e = "->", w = "D",
                           temporality = "y")
  invisible(capture.output(
    replayed_multi <- suppressMessages(assess_edges(
      multi,
      keep_edges = kept,
      assess_causal_criteria = TRUE,
      save_answers = replay_log,
      causal_criteria = criteria
    ))
  ))
  expect_identical(nrow(replayed_multi), 3L)

  two_criteria <- rbind(criteria, transform(criteria, name = "theory"))
  rejected_log <- data.frame(v = "A", e = "->", w = "B",
                             temporality = "n", theory = NA_character_)
  invisible(capture.output(
    rejected <- suppressMessages(assess_edges(
      dag,
      assess_causal_criteria = TRUE,
      save_answers = rejected_log,
      causal_criteria = two_criteria
    ))
  ))
  expect_identical(nrow(rejected), 0L)
})

test_that("internal empty and overlapping-role paths are stable", {
  empty_edges <- data.frame(v = character(), e = character(), w = character())
  expect_identical(causaliflower:::print_edges_helper(empty_edges), list())

  dag <- string_to_dag("A -> B")
  dag_df <- data.frame(name = c("A", "B"))
  labels <- c(A = "A", B = "B")
  roles <- list(confounder = "A", observed = c("A", "B"))
  labelled <- causaliflower:::add_labels(dag, dag_df, labels, roles)
  expect_identical(labelled$role, c("confounder", "observed"))
})

test_that("instruments retains its legacy result shape with optional details", {
  dag <- build_graph("X", "T", "Y", instrumental_variables = "IV",
                     type = "saturated")

  legacy <- instruments(dag)
  expect_identical(legacy, list(T = list(Y = "IV")))

  details <- instruments(dag, details = TRUE)
  expect_s3_class(details, "data.table")
  expect_identical(names(details),
                   c("treatment", "outcome", "instrument",
                     "conditioning_set"))
  expect_identical(details$instrument, "IV")

  no_instrument <- build_graph("X", "T", "Y", type = "saturated")
  expect_null(instruments(no_instrument)$T$Y)
})

test_that("add_nodes preserves every declared node status", {
  dag <- build_graph("X", "T", "Y", type = "ordered")

  treatment <- suppressMessages(
    add_nodes(dag, "T2", role = "treatment", print_edges = FALSE)
  )
  outcome <- suppressMessages(
    add_nodes(dag, "Y2", role = "outcome", print_edges = FALSE)
  )
  latent <- suppressMessages(
    add_nodes(dag, "U", role = "latent", print_edges = FALSE)
  )
  observed <- suppressMessages(
    add_nodes(dag, "O", role = "observed", print_edges = FALSE)
  )

  expect_true("T2" %in% treatments(treatment))
  expect_true("Y2" %in% causaliflower:::get_outcome_names(outcome))
  expect_true("U" %in% unobserved(latent))
  expect_true("O" %in% names(observed))
  expect_true("O" %in% get_roles(observed)$observed)
})


test_that("partially directed graphs retain reportable roles and can plot", {
  dag <- dagitty::dagitty(
    "pdag { T -- Z; Z -- W; T -> Y; T [exposure]; Y [outcome] }"
  )
  roles <- get_roles(dag)

  expect_setequal(roles$undetermined, c("W", "Z"))
  expect_true("competing_cause" %in% names(roles))
  expect_false("competing_exposure" %in% names(roles))
  expect_identical(get_nodes(dag)$W, "undetermined")
  expect_identical(get_nodes(dag)$Z, "undetermined")
  expect_equal(nrow(get_edges(dag)), 3L)
  expect_equal(nrow(get_edges(dag, "undetermined")), 2L)

  expect_warning(
    rebuilt <- build_graph(variables = dag, type = "saturated"),
    "constructs a new DAG"
  )
  expect_setequal(names(rebuilt), names(dag))
  expect_match(as.character(rebuilt), "^dag ")
  expect_s3_class(suppressMessages(plot_dagitty(dag)), "ggplot")

  pag <- dagitty::dagitty("pag { T @-@ Y }")
  expect_error(get_edges(pag), "does not currently support")
})

test_that("coordinate helpers handle isolated graphs and no-op joins", {
  empty <- string_to_dag("")
  expect_identical(add_coords(empty), empty)
  expect_identical(add_coords(empty, x_step = 4), empty)

  isolated <- string_to_dag(c("A", "B", "C"))
  expect_silent(positioned <- add_coords(isolated))
  coordinates <- dagitty::coordinates(positioned)
  expect_equal(coordinates$x, c(A = 1, B = 3, C = 5))
  expect_equal(coordinates$y, c(A = 0, B = 0, C = 0))

  spaced <- dagitty::coordinates(add_coords(isolated, x_step = 4))
  expect_equal(spaced$x, c(A = 1, B = 5, C = 9))
  expect_equal(spaced$y, c(A = 0, B = 0, C = 0))
  close <- dagitty::coordinates(add_coords(isolated, x_step = 0.25))
  expect_equal(close$x, c(A = 1, B = 1.25, C = 1.5))

  single <- add_coords(string_to_dag("A"), x_step = 4)
  expect_equal(dagitty::coordinates(single), list(x = c(A = 1), y = c(A = 0)))

  dag <- string_to_dag("A -> B")
  dagitty::coordinates(dag) <- list(x = c(A = 1, B = 2),
                                    y = c(A = 3, B = 4))
  joined <- join_graphs(dag, dag)
  expect_equal(dagitty::coordinates(joined), dagitty::coordinates(dag))
})

test_that("repeated relative-position clauses retain clause order", {
  position <- c(after("Z2", "V"), after("Z2", "W"))
  resolved <- causaliflower:::resolve_position(
    position, c("V", "W"), c("Z1", "Z2", "Z3")
  )
  expect_identical(resolved, c("Z1", "Z2", "V", "W", "Z3"))
})

test_that("causal criteria replay handles optional-only criteria and non-tabular answers", {
  criteria <- data.frame(
    name = c("temporality", "theory"),
    question = "question",
    description = "description",
    source = "source",
    required = "no",
    stringsAsFactors = FALSE
  )
  dag <- string_to_dag(c("A -> B", "B -> C", "A -> C"))
  replay_log <- data.frame(v = c("A", "B"), e = "->", w = c("C", "C"),
                           temporality = c("y", "n"), theory = c("n", "y"))

  invisible(capture.output(
    optional <- suppressMessages(assess_edges(
      dag,
      keep_edges = "A -> B",
      assess_causal_criteria = TRUE,
      save_answers = replay_log,
      causal_criteria = criteria
    ))
  ))
  expect_identical(sort(paste(optional$v, optional$e, optional$w)), edge_keys(dag))

  required_criteria <- criteria
  required_criteria$required <- "yes"
  invisible(capture.output(
    dropped <- suppressMessages(assess_edges(
      dag,
      keep_edges = "A -> B",
      assess_causal_criteria = TRUE,
      save_answers = replay_log,
      causal_criteria = required_criteria
    ))
  ))
  expect_identical(nrow(dropped), 1L)

  expect_error(
    assess_edges(dag, keep_edges = "A -> B", assess_causal_criteria = TRUE,
                 save_answers = c("y", "y", "y")),
    "Length of answers does not match"
  )
})
