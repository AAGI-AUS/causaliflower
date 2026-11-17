reuse_query <- function(columns = NULL, policy = NULL, depth = Inf) {
  graph <- dagitty::dagitty(paste("dag { T [exposure]; Y [outcome];",
    "A -> B -> T; A -> M; T -> M -> Y; T -> Y; V -> Y }"))
  nodes <- names(graph)
  if (is.null(policy)) policy <- data.frame(node_id = nodes, recovered_node = "yes")
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "barley", columns = columns,
      params = list(water_model = "SoilWater", query_depth = depth)),
    requests = list(columns = list(resolved = stats::setNames(columns, columns))),
    node_map = data.frame(node = nodes, quantity = nodes, member_id = nodes, occurrence = nodes),
    policy = policy)
  graph
}


test_that("CTYPE changes the prepared search identity only when names sort differently", {
  ctype <- Sys.getlocale("LC_CTYPE")
  on.exit(invisible(Sys.setlocale("LC_CTYPE", ctype)), add = TRUE)
  if (!nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))))
    skip("A second CTYPE locale is unavailable")
  graphs <- list(reuse_query(), dagitty::dagitty(paste0(
    "dag { T [exposure]; Y [outcome]; T -> Y; Z -> T; Z -> z -> Y; ",
    "\"Å\" -> T; \"Å\" -> \"é\" -> Y }")))
  for (i in seq_along(graphs)) {
    graph <- graphs[[i]]
    invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
    before <- .cd_search_identity(.cd_prepare(graph, NULL, NULL))
    invisible(Sys.setlocale("LC_CTYPE", "C"))
    after <- .cd_search_identity(.cd_prepare(graph, NULL, NULL))
    if (i == 1L) expect_identical(after, before) else expect_false(identical(after, before))
  }
})
