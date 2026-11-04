test_that("get_node_edges() keeps focal direction separate from symmetric edges", {
  dag <- dagitty::dagitty(
    "dag { A -> B; B -> C; B <-> D; B -- E; F }"
  )

  any_edges <- get_node_edges(dag, nodes = c("B", "D", "F"))
  expect_identical(names(any_edges), c("B", "D", "F"))
  expect_identical(nrow(any_edges$B), 4L)
  expect_identical(nrow(any_edges$D), 1L)
  expect_identical(nrow(any_edges$F), 0L)
  expect_setequal(any_edges$B$e, c("->", "--", "<->"))

  incoming <- get_node_edges(dag, "B", direction = "incoming")$B
  outgoing <- get_node_edges(dag, "B", direction = "outgoing")$B
  expect_identical(paste(incoming$v, incoming$e, incoming$w), "A -> B")
  expect_identical(paste(outgoing$v, outgoing$e, outgoing$w), "B -> C")
  expect_false(any(c("--", "<->") %in% c(incoming$e, outgoing$e)))

  expect_error(get_node_edges(dag, "missing"), "not in the graph")
  expect_error(get_node_edges(dag, direction = "sideways"), "must be one of")
})


test_that("get_paths() indexes maximal directed paths by focal and other node", {
  dag <- dagitty::dagitty(
    "dag { Y [outcome]; C -> A; D -> A; D -> B; A -> Y; B -> Y }"
  )
  paths <- get_paths(dag, nodes = "Y")

  expect_s3_class(paths, "causaliflower_paths")
  expect_identical(
    attr(paths$Y, "incoming"),
    list(
      a = c("C", "A", "Y"),
      b = c("D", "A", "Y"),
      c = c("D", "B", "Y")
    )
  )
  expect_identical(paths$Y$A$incoming, attr(paths$Y, "incoming")[c("a", "b")])
  expect_identical(paths$Y$D$starts_at,
                   attr(paths$Y, "incoming")[c("b", "c")])
  expect_length(paths$Y$A$starts_at, 0L)
  expect_length(paths$Y$A$ends_at, 0L)
  expect_identical(paths$Y$A$any, paths$Y$A$incoming)
  expect_null(paths$Y$Y)

  prefix_only <- get_paths(string_to_dag("AA -> Y"), nodes = "Y")
  expect_null(causaliflower:::`$.causaliflower_path_index`(
    prefix_only[["Y"]], "A"
  ))
  expect_identical(prefix_only$Y$AA$incoming$a, c("AA", "Y"))

  focal_prefix <- get_paths(string_to_dag("AA -> Y"))
  expect_null(causaliflower:::`$.causaliflower_paths`(focal_prefix, "A"))
})


test_that("get_paths() keeps the focal node as a path endpoint", {
  paths <- get_paths(string_to_dag("A -> F -> B"), nodes = "F")

  expect_identical(attr(paths$F, "incoming")$a, c("A", "F"))
  expect_identical(attr(paths$F, "outgoing")$b, c("F", "B"))
  expect_false(any(vapply(
    c(attr(paths$F, "incoming"), attr(paths$F, "outgoing")),
    identical, logical(1L), c("A", "F", "B")
  )))
})


test_that("get_paths() supports non-syntactic names and outgoing lookup", {
  dag <- string_to_dag(c(
    '"soil N" -> "crop growth"',
    '"crop growth" -> "grain yield"'
  ))
  paths <- get_paths(dag, nodes = "crop growth")
  focal <- paths[["crop growth"]]

  expect_identical(
    focal[["soil N"]]$incoming$a,
    c("soil N", "crop growth")
  )
  expect_identical(
    focal[["grain yield"]]$outgoing$b,
    c("crop growth", "grain yield")
  )
  expect_identical(focal[["grain yield"]]$ends_at,
                   focal[["grain yield"]]$outgoing)
})


test_that("get_paths() records finite path and length truncation", {
  star <- string_to_dag(c("A -> Y", "B -> Y", "C -> Y"))
  limited <- get_paths(star, nodes = "Y", max_paths = 2L)
  unlimited <- get_paths(star, nodes = "Y", max_paths = Inf)

  expect_identical(length(attr(limited$Y, "incoming")), 2L)
  expect_true(attr(limited$Y, "truncated")$incoming[["max_paths"]])
  expect_identical(length(attr(unlimited$Y, "incoming")), 3L)
  expect_false(attr(unlimited$Y, "truncated")$incoming[["max_paths"]])

  both_directions <- string_to_dag(c(
    "A -> F", "B -> F", "F -> C", "F -> D"
  ))
  both_limited <- get_paths(both_directions, nodes = "F", max_paths = 1L)
  expect_identical(length(attr(both_limited$F, "incoming")), 1L)
  expect_identical(length(attr(both_limited$F, "outgoing")), 1L)
  expect_true(attr(both_limited$F, "truncated")$incoming[["max_paths"]])
  expect_true(attr(both_limited$F, "truncated")$outgoing[["max_paths"]])

  chain <- string_to_dag("A -> B -> C -> Y")
  short <- get_paths(chain, nodes = "Y", max_length = 2L)
  full <- get_paths(chain, nodes = "Y", max_length = Inf)
  expect_length(attr(short$Y, "incoming"), 0L)
  expect_true(attr(short$Y, "truncated")$incoming[["max_length"]])
  expect_identical(attr(full$Y, "incoming")$a,
                   c("A", "B", "C", "Y"))

  exact <- get_paths(string_to_dag("A -> B -> Y"), nodes = "Y",
                     max_length = 2L)
  expect_identical(attr(exact$Y, "incoming")$a, c("A", "B", "Y"))
  expect_false(attr(exact$Y, "truncated")$incoming[["max_length"]])

  many <- get_paths(
    string_to_dag(paste0("N", seq_len(27L), " -> Y")),
    nodes = "Y", max_paths = Inf
  )
  expect_identical(names(attr(many$Y, "incoming")), c(letters, "aa"))

  printed <- capture.output(print.causaliflower_paths(limited))
  expect_match(paste(printed, collapse = "\n"),
               "Incoming paths truncated at max_paths = 2")

  for( bad in list(0, -1, NA_real_, 1.5, c(1, 2), "2", -Inf) ){
    expect_error(get_paths(star, max_paths = bad),
                 "positive whole number or Inf")
    expect_error(get_paths(star, max_length = bad),
                 "positive whole number or Inf")
  }
})


test_that("get_structure() reports node relationships and minimum distances", {
  dag <- dagitty::dagitty(
    "dag { Y [outcome]; C -> A; D -> A; D -> B; A -> Y; B -> Y }"
  )
  structure <- get_structure(dag)
  y <- structure$Y

  expect_s3_class(structure, "causaliflower_structure")
  expect_identical(names(y), c(
    "role", "parents", "children", "undirected", "bidirected",
    "ancestors", "descendants", "paths"
  ))
  expect_identical(y$role, "outcome")
  expect_identical(y$parents, c("A", "B"))
  expect_identical(y$children, character())
  expect_identical(y$ancestors$d1,
                   list(A = c("a", "b"), B = "c"))
  expect_identical(y$ancestors$d2,
                   list(C = "a", D = c("b", "c")))
  expect_identical(y$descendants, list())
  expect_identical(unname(lapply(y$paths, identity)), list(
    c("C", "A", "Y"), c("D", "A", "Y"), c("D", "B", "Y")
  ))

  printed <- capture.output(print.causaliflower_structure(structure["Y"]))
  printed <- paste(printed, collapse = "\n")
  expect_match(printed, "1\\. \\[A \\(a\\)\\(b\\)\\] \\[B \\(c\\)\\]")
  expect_match(printed, "2\\. \\[C \\(a\\)\\] \\[D \\(b\\)\\(c\\)\\]")
  expect_match(printed, "\\(a\\) C -> A -> Y")
})


test_that("get_structure() lists each node at its closest distance and keeps every path key", {
  dag <- dagitty::dagitty(
    "dag { Y [outcome]; D -> Y; D -> A; A -> Y }"
  )
  y <- get_structure(dag)$Y

  expect_identical(names(y$ancestors), "d1")
  expect_identical(y$ancestors$d1$A, "a")
  expect_identical(y$ancestors$d1$D, c("a", "b"))
})


test_that("get_structure() separates direct mixed-edge relationships", {
  dag <- dagitty::dagitty(
    "dag { A -> B; B -> C; B <-> D; B -- E; F }"
  )
  structure <- get_structure(dag)

  expect_identical(structure$B$parents, "A")
  expect_identical(structure$B$children, "C")
  expect_identical(structure$B$bidirected, "D")
  expect_identical(structure$B$undirected, "E")
  expect_identical(structure$F$parents, character())
  expect_identical(structure$F$ancestors, list())
  expect_length(structure$F$paths, 0L)
})


test_that("get_structure() passes through finite and infinite path limits", {
  dag <- string_to_dag(c("A -> Y", "B -> Y", "C -> Y"))
  limited <- get_structure(dag, max_paths = 2L)
  unlimited <- get_structure(dag, max_paths = Inf, max_length = Inf)

  expect_identical(length(limited$Y$paths), 2L)
  expect_true(attr(limited$Y$paths, "truncated")$incoming[["max_paths"]])
  expect_identical(length(unlimited$Y$paths), 3L)
  expect_false(attr(unlimited$Y$paths, "truncated")$incoming[["max_paths"]])

  long_chain <- string_to_dag("A -> B -> C -> Y")
  length_limited <- get_structure(long_chain, max_length = 2L)
  expect_length(length_limited$Y$paths, 0L)
  expect_true(
    attr(length_limited$Y$paths, "truncated")$incoming[["max_length"]]
  )
  expect_identical(names(length_limited$Y$ancestors), c("d1", "d2", "d3"))
  expect_identical(names(length_limited$Y$ancestors$d1), "C")

  printed <- capture.output(print.causaliflower_structure(limited["Y"]))
  expect_match(paste(printed, collapse = "\n"),
               "Incoming paths truncated at max_paths = 2")
})


test_that("get_ancestor_edges() warns and groups edges under their left-hand node", {
  dag <- string_to_dag(c("A -> B", "B -> C", "B <-> D"))
  expect_warning(
    old <- get_ancestor_edges(dag),
    paste0("get_ancestor_edges(): 'get_ancestor_edges()' is deprecated as of ",
           "0.0.2. Use 'get_node_edges(direction = \"any\")' instead."),
    fixed = TRUE
  )

  expect_identical(names(old), c("A", "B"))
  expect_true(all(vapply(old, data.table::is.data.table, logical(1L))))
  expect_identical(paste(old$A$v, old$A$e, old$A$w), "A -> B")
  expect_setequal(paste(old$B$v, old$B$e, old$B$w), c("B -> C", "B <-> D"))

  empty <- dagitty::dagitty("dag { T [exposure]; Y [outcome] }")
  expect_warning(
    empty_old <- get_ancestor_edges(empty),
    paste0("get_ancestor_edges(): 'get_ancestor_edges()' is deprecated as of ",
           "0.0.2. Use 'get_node_edges(direction = \"any\")' instead."),
    fixed = TRUE
  )
  expect_length(empty_old, 0L)
})


test_that("inspection keeps node names that contain square brackets", {
  dag <- dagitty::dagitty(
    'dag { "yield [t/ha]" [outcome]; "N rate [kg/ha]" [exposure];
           "N rate [kg/ha]" -> "yield [t/ha]" }'
  )
  node_edges <- get_node_edges(dag)

  expect_identical(causaliflower:::get_node_names(dag),
                   as.character(names(dag)))
  expect_identical(observed(dag), as.character(names(dag)))
  expect_identical(names(node_edges), c("N rate [kg/ha]", "yield [t/ha]"))
  expect_identical(
    sort(paste(node_edges[[1L]]$v, node_edges[[1L]]$e, node_edges[[1L]]$w)),
    edge_keys(dag)
  )
  expect_identical(get_structure(dag)[["yield [t/ha]"]]$parents,
                   "N rate [kg/ha]")
  expect_identical(get_paths(dag)[["yield [t/ha]"]][["N rate [kg/ha]"]]$incoming,
                   list(a = c("N rate [kg/ha]", "yield [t/ha]")))
})


test_that("inspection ignores the graph name and bounding box of a native graph", {
  dag <- dagitty::dagitty(
    'dag named { bb="-1,-1,3,3"
       T [exposure,pos="0,0"] M [pos="1,0"] Y [outcome,pos="2,0"]
       T -> M M -> Y T -> Y }'
  )
  edited <- remove_edges(dag, "T -> Y")

  for( graph in list(dag, edited) ){
    expect_identical(causaliflower:::get_node_names(graph),
                     as.character(names(graph)))
    expect_identical(observed(graph), c("M", "T", "Y"))
    expect_identical(names(get_structure(graph)), c("M", "T", "Y"))
    expect_identical(names(get_paths(graph)), c("M", "T", "Y"))
    expect_identical(names(get_node_edges(graph)), c("M", "T", "Y"))
  }

  expect_identical(edge_keys(edited), sort(c("M -> Y", "T -> M")))
  expect_identical(get_structure(edited)$Y$parents, "M")
})
