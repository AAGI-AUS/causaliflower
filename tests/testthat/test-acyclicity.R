test_that("the acyclicity check agrees with dagitty", {
  graphs <- c(
    "dag { }",
    "dag { A B C }",
    "dag { A -> B B -> C }",
    "dag { A -> B B -> C C -> A }",
    "dag { A -> B B -> A }",
    "dag { A <-> B B <-> C C <-> A }",
    "dag { A -> B B <-> A }",
    "dag { \"node one\" -> \"node two\" \"node two\" -> \"node one\" }",
    "pdag { A -- B B -> C C -> A }",
    "pdag { A -- B B -- C C -- A }"
  )
  for (text in graphs) {
    dag <- dagitty::dagitty(text)
    expect_identical(dag_is_acyclic(dag), dagitty::isAcyclic(dag), info = text)
  }

  set.seed(20261006)
  for (index in seq_len(200L)) {
    nodes <- paste0("V", seq_len(sample(2:8, 1L)))
    pairs <- expand.grid(from = nodes, to = nodes, stringsAsFactors = FALSE)
    pairs <- pairs[pairs$from != pairs$to, , drop = FALSE]
    keep <- stats::runif(nrow(pairs)) < 0.2
    type <- ifelse(stats::runif(nrow(pairs)) < 0.8, "->", "<->")
    lines <- paste(pairs$from[keep], type[keep], pairs$to[keep])
    dag <- dagitty::dagitty(paste("dag {", paste(c(nodes, lines),
                                                 collapse = "\n"), "}"))
    expect_identical(dag_is_acyclic(dag), dagitty::isAcyclic(dag),
                     info = as.character(dag))
  }
})

test_that("package code does not call dagitty's path-enumerating cycle check", {
  # dagitty::isAcyclic() follows every directed path and can take minutes on
  # dense graphs built by build_graph().
  namespace <- asNamespace("causaliflower")
  callers <- Filter(function(name) {
    value <- get(name, envir = namespace)
    is.function(value) &&
      any(grepl("isAcyclic", deparse(body(value)), fixed = TRUE))
  }, ls(namespace, all.names = TRUE))
  expect_identical(callers, character())
})
