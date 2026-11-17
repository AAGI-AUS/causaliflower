test_that("constructor retains the source and prints only the graph", {
  q <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; U -> T; U -> M; U -> Y; T -> M -> Y }")
  attr(q, "author_note") <- "preserve"
  before <- serialize(q, NULL)
  for (size in list(Inf, 0)) {
    g <- suppressWarnings(construct_dag(q, size, quiet = TRUE))
    expect_identical(serialize(original_query(g), NULL), before)
    expect_identical(paste(capture.output(print(g)), collapse = "\n"), as.character(g))
  }
  expect_identical(serialize(q, NULL), before)
})


test_that("construction on a restricted emitted query uses it as its reference", {
  q <- query_graph("barley", treatments = "N_application", outcomes = "grain_N",
    randomised_treatment = TRUE, quiet = TRUE, params = query_params(query_depth = 1))
  g <- suppressWarnings(construct_dag(q, 0, quiet = TRUE))
  expect_true(any(grepl("supplied query_depth", dag_options(g)$limitations)))
  expect_identical(original_query(g), q)
  expect_true(all(names(temporal_dag(g)) %in% names(q)))
})
