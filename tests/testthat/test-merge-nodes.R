test_that("merge_nodes() combines explicit groups and deduplicates typed edges", {
  dag <- dagitty::dagitty("dag { T -> M1 -> Y; T -> M2 -> Y; M1 <-> Y; I }")
  before <- dag
  merged <- merge_nodes(dag, c("M1", "M2"), "M")
  expected <- dagitty::dagitty("dag { T -> M -> Y; M <-> Y; I }")
  expect_setequal(names(merged), names(expected))
  expect_setequal(edge_key(dag_edges(merged)), edge_key(dag_edges(expected)))
  expect_identical(dag, before)

  dag <- dagitty::dagitty("dag { A -> M1; M2 -> B }")
  merged <- merge_nodes(dag, c("M1", "M2"), "M")
  expected <- dagitty::dagitty("dag { A -> M -> B }")
  expect_setequal(edge_key(dag_edges(merged)), edge_key(dag_edges(expected)))
})

test_that("merge_nodes() removes internal edges and retains all other nodes", {
  dag <- dagitty::dagitty("dag { A -> B; A <-> B; I }")
  merged <- merge_nodes(dag, c("A", "B"), "M")
  expect_setequal(names(merged), c("I", "M"))
  expect_equal(nrow(dag_edges(merged)), 0L)
  expect_identical(dag_graph_type(merged), "dag")

  dag <- dagitty::dagitty("dag { A; B }")
  merged <- merge_nodes(dag, c("A", "B"), "A")
  expect_identical(names(merged), "A")
  expect_equal(nrow(dag_edges(merged)), 0L)
})

test_that("merge_nodes() keeps the first node's metadata when reusing a name", {
  dag <- structure(paste0(
    'dag named { bb="0,0,9,9"\n',
    'A [pos="1,2",label="first",note="keep"]\n',
    'B [pos="3,4",label="second",note="discard"]\n',
    'Y [pos="7,8",note="unchanged"]\n',
    'I [pos="8,1"]\n',
    'A -> Y [pos="4,6",note="first edge"]\n',
    'B -> Y [pos="6,5",note="second edge"]\n',
    'I -> Y [note="unaffected edge"]\n}'
  ), class = "dagitty")
  attr(dag, "labels") <- c(A = "alpha", B = "beta", Y = "outcome")
  attr(dag, "custom") <- list(nodes = c("A", "B"), data = 3L)
  before <- dag
  merged <- merge_nodes(dag, c("A", "B"), "B")
  text <- as.character(merged)
  expect_match(text, 'dag named { bb="0,0,9,9"', fixed = TRUE)
  expect_match(text, '"B" [pos="1,2",label="first",note="keep"]', fixed = TRUE)
  expect_match(text, 'Y [pos="7,8",note="unchanged"]', fixed = TRUE)
  expect_match(text, '"B" -> Y [pos="4,6",note="first edge"]', fixed = TRUE)
  expect_match(text, 'I -> Y [note="unaffected edge"]', fixed = TRUE)
  expect_false(grepl("second", text, fixed = TRUE))
  expect_identical(unname(dagitty::coordinates(merged)$x["B"]),
                   unname(dagitty::coordinates(dag)$x["A"]))
  expect_identical(unname(dagitty::coordinates(merged)$y["B"]),
                   unname(dagitty::coordinates(dag)$y["A"]))
  expect_identical(attr(merged, "labels"), c(B = "alpha", Y = "outcome"))
  expect_identical(attr(merged, "custom"), attr(dag, "custom"))
  expect_identical(dag, before)
})

test_that("merge_nodes() makes edges inherited from other members straight", {
  dag <- dagitty::dagitty(paste0(
    'dag { A [pos="0,0"]; B [pos="0,2"]; Y [pos="2,0"]; ',
    'B -> Y [pos="1,3",note="removed owner"] }'
  ))
  merged <- merge_nodes(dag, c("A", "B"), "M")
  expect_setequal(edge_key(dag_edges(merged)), "M\r->\rY")
  expect_equal(nrow(edge_controls(merged, "test")), 0L)
  expect_false(grepl("removed owner", as.character(merged), fixed = TRUE))
  expect_identical(unname(dagitty::coordinates(merged)$x["M"]),
                   unname(dagitty::coordinates(dag)$x["A"]))
  expect_identical(unname(dagitty::coordinates(merged)$y["M"]),
                   unname(dagitty::coordinates(dag)$y["A"]))
})

test_that("merge_nodes() requires agreement for every native status", {
  readers <- list(
    exposure = dagitty::exposures, outcome = dagitty::outcomes,
    latent = dagitty::latents, adjusted = dagitty::adjustedNodes,
    selected = selected_nodes
  )
  for( status in names(readers) ){
    dag <- dagitty::dagitty(paste0(
      "dag { A [", status, "]; B [", status, "]; A -> Y; B -> Y }"
    ))
    merged <- merge_nodes(dag, c("A", "B"), "M")
    parsed <- dagitty::dagitty(as.character(merged))
    expect_setequal(readers[[status]](parsed), "M")
    mixed <- dagitty::dagitty(paste0(
      "dag { A [", status, "]; B; A -> Y; B -> Y }"
    ))
    before <- mixed
    expect_error(merge_nodes(mixed, c("A", "B"), "M"),
                 "matching native statuses", fixed = TRUE)
    expect_identical(mixed, before)
  }
})

test_that("merge_nodes() validates selections, names and the resulting graph", {
  dag <- dagitty::dagitty("dag { A -> B -> C; I }")
  before <- dag
  expect_error(merge_nodes(dag, c("A", "C"), "M"), "directed cycle")
  expect_error(merge_nodes(dag, c("A", "B"), "I"), "outside the group")
  expect_error(merge_nodes(dag, c("A", "absent"), "M"), "not present")
  expect_error(merge_nodes(dag, "A", "M"), "at least two")
  expect_error(merge_nodes(dag, c("A", "A"), "M"), "at least two")
  for( nodes in list(character(), NA_character_, "", 1L) ){
    expect_error(merge_nodes(dag, nodes, "M"), "node names")
  }
  for( name in list(character(), NA_character_, "", c("M", "N"), 1L) ){
    expect_error(merge_nodes(dag, c("A", "B"), name), "one non-missing")
  }
  for( value in list(NA, NULL, "yes", c(TRUE, FALSE)) ){
    expect_error(merge_nodes(dag, c("A", "B"), "M", print_edges = value),
                 "TRUE or FALSE")
  }
  merged <- merge_nodes(dag, c("A", "B", "A"), "M")
  expect_setequal(names(merged), c("M", "C", "I"))
  expect_identical(dag, before)
})

test_that("merge functions require acyclic native DAGs with supported edges", {
  for( function_name in c("merge_nodes", "recombine_nodes") ){
    operation <- get(function_name)
    for( text in c("pdag { A -> B }", "mag { A -> B }") ){
      args <- list(dag = dagitty::dagitty(text))
      if( function_name == "merge_nodes" ){
        args <- c(args, list(nodes = c("A", "B"), name = "M"))
      }
      expect_error(do.call(operation, args), "dag graph type")
    }
    for( text in c("dag { A -- B }", "dag { A -> B -> A }",
                   "dag { A; B; I <-> I }") ){
      args <- list(dag = dagitty::dagitty(text))
      if( function_name == "merge_nodes" ){
        args <- c(args, list(nodes = c("A", "B"), name = "M"))
      }
      expect_error(do.call(operation, args), "supported|acyclic")
    }
  }
})

test_that("recombine_nodes() previews candidates using the whole graph", {
  dag <- dagitty::dagitty("dag { T -> M1 -> Y; T -> M2 -> Y; A -> B; I }")
  expect_identical(unname(recombine_nodes(dag, preview = TRUE)),
                   list(c("M1", "M2")))
  expect_identical(unname(recombine_nodes(dag, c("M2", "M1"), TRUE)),
                   list(c("M1", "M2")))
  expect_identical(recombine_nodes(dag, "M1", TRUE), list())
  different <- dagitty::dagitty("dag { T -> M1 -> Y; X -> M2 -> Y }")
  expect_identical(recombine_nodes(different, c("M1", "M2"), TRUE), list())
  expect_error(recombine_nodes(dag, c("M1", "absent"), TRUE), "not present")
  expect_error(recombine_nodes(dag, NA_character_, TRUE), "node names")
  expect_identical(recombine_nodes(dagitty::dagitty("dag {}"), preview = TRUE),
                   list())
  expect_identical(recombine_nodes(dagitty::dagitty("dag { A; B }"),
                                    preview = TRUE), list())
})

test_that("recombine_nodes() compares open and closed bidirected neighbours", {
  dag <- dagitty::dagitty("dag { A <-> B; A <-> C; B <-> C }")
  expect_identical(unname(recombine_nodes(dag, preview = TRUE)),
                   list(c("A", "B", "C")))
  dag <- dagitty::dagitty("dag { A <-> U; B <-> U }")
  expect_identical(unname(recombine_nodes(dag, preview = TRUE)),
                   list(c("A", "B")))
  dag <- dagitty::dagitty("dag { A <-> B; B <-> C; C <-> D }")
  expect_identical(recombine_nodes(dag, preview = TRUE), list())
  dag <- dagitty::dagitty("dag { T -> A -> Y; T -> B -> Y; A <-> U }")
  expect_identical(recombine_nodes(dag, c("A", "B"), TRUE), list())
})

test_that("recombine_nodes() protects status endpoints and separates latent nodes", {
  dag <- dagitty::dagitty(paste0(
    "dag { T [exposure]; Y [outcome]; S [selected]; X [adjusted]; ",
    "U1 [latent]; U2 [latent]; ",
    "T -> K; Y -> K; S -> K; X -> K; U1 -> K; U2 -> K; A -> K; B -> K }"
  ))
  expect_identical(unname(recombine_nodes(dag, preview = TRUE)),
                   list(c("A", "B"), c("U1", "U2")))
})

test_that("merge functions preserve quoted names and produce no hidden state", {
  dag <- dagitty::dagitty('dag { T -> "M 1" -> Y; T -> "M[2]" -> Y }')
  expect_identical(unname(recombine_nodes(dag, preview = TRUE)),
                   list(c("M 1", "M[2]")))
  merged <- merge_nodes(dag, c("M 1", "M[2]"), "M [combined]")
  expected <- dagitty::dagitty('dag { T -> "M [combined]" -> Y }')
  expect_setequal(names(merged), names(expected))
  expect_setequal(edge_key(dag_edges(merged)), edge_key(dag_edges(expected)))
  expect_identical(attributes(merged), attributes(dag))
  expect_error(merge_nodes(dag, c("M 1", "M[2]"), 'M "combined"'),
                 "cannot contain quotes", fixed = TRUE)
})

test_that("merge_nodes() printing is optional and reports only retained edges", {
  dag <- dagitty::dagitty("dag { A -> B -> Y; I }")
  expect_output(quiet <- merge_nodes(dag, c("A", "B"), "M"), NA)
  output <- capture.output(loud <- merge_nodes(dag, c("A", "B"), "M", TRUE))
  expect_identical(loud, quiet)
  expect_match(paste(output, collapse = "\n"), "redirected edges", fixed = TRUE)
  expect_true(any(grepl("M", output, fixed = TRUE)))
  expect_true(any(grepl("Y", output, fixed = TRUE)))
  expect_false(any(grepl("I", output, fixed = TRUE)))
})

test_that("existing edge-printing arguments default to FALSE", {
  for( operation in list(add_nodes, connect_nodes, copy_nodes, merge_nodes) ){
    expect_identical(formals(operation)$print_edges, FALSE)
  }
  dag <- dagitty::dagitty("dag { T -> M -> Y; X }")
  expect_output(added <- add_nodes(dag, "N", children = "Y"), NA)
  expect_output(connected <- connect_nodes(dag, c("T", "X"), type = "full"), NA)
  expect_output(copied <- copy_nodes(dag, "M"), NA)
  expect_output(add_nodes(dag, "N", children = "Y", print_edges = TRUE), "Y")
  expect_output(connect_nodes(dag, c("T", "X"), type = "full", print_edges = TRUE), "T")
  expect_output(copy_nodes(dag, "M", print_edges = TRUE), "M")
})
