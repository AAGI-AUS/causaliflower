variable_query <- function(edges, ids = NULL, exposures = "T") {
  graph <- dagitty::dagitty(paste("dag {", paste(exposures, "[exposure];", collapse = " "),
    "Y [outcome];", edges, "}"))
  nodes <- names(graph)
  if (is.null(ids)) ids <- stats::setNames(nodes, nodes)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = NULL, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = list())),
    node_map = data.frame(node = nodes, quantity = unname(ids[nodes]),
      member_id = unname(ids[nodes]), occurrence = nodes),
    policy = data.frame(node_id = unique(unname(ids)), recovered_node = "yes"))
  graph
}

same_size_query <- function() {
  variable_query("T -> Y; T -> P_t0 -> P_t1",
    c(T = "T", Y = "Y", P_t0 = "P", P_t1 = "P"))
}

test_that("literal plain node names are separate variables and exact requests stay exact", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; P_t0 -> T; P_t0 -> Y; P_t1 -> T; P_t1 -> Y; T -> Y }")
  expect_identical(.cd_variable_count(graph), 4L)
  expect_identical(.cd_match_nodes(graph, "P_t0"), "P_t0")
  expect_error(.cd_match_nodes(graph, "P"), "absent")
  built <- suppressWarnings(dagby_construct(graph, quiet = TRUE, size = "largest"))
  expect_setequal(names(built), names(graph))
  expect_identical(dagby_options(built)$level$variables, 4L)
  expect_setequal(dagby_options(built)$quantities$total$members, c("P_t0", "P_t1"))
})

test_that("source maps count canonical members and refuse fabricated identities", {
  graph <- same_size_query()
  expect_identical(.cd_variable_count(graph), 3L)
  record <- attr(graph, "query_record")
  record$requests$columns$resolved <- c(alias = "P")
  attr(graph, "query_record") <- record
  expect_setequal(.cd_match_nodes(graph, "alias"), c("P_t0", "P_t1"))
  record$node_map$quantity[record$node_map$member_id == "P"] <- "printed-combination"
  extra <- record$node_map[record$node_map$node == "P_t0", ]
  extra$member_id <- "other-source"
  record$node_map <- rbind(record$node_map, extra)
  attr(graph, "query_record") <- record
  expect_identical(.cd_variable_count(graph), 4L)
  record$node_map$member_id[[1L]] <- NA_character_
  attr(graph, "query_record") <- record
  expect_error(dagby_construct(graph, quiet = TRUE, size = "largest"), "node map")
  expect_error(.cd_variable_map(graph), "identity map")
  record$node_map <- record$node_map[record$node_map$node != record$node_map$node[[1L]], ]
  attr(graph, "query_record") <- record
  expect_error(.cd_variable_map(graph), "identity map")
})

test_that("same-variable alternatives retain exact temporal membership and stable numbers", {
  graph <- same_size_query()
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  full <- dagby_construct(graph, quiet = TRUE, size = "largest")
  choices <- dagby_options(full)$alternatives
  expect_identical(dagby_options(full)$level$variables, 3L)
  expect_identical(choices$select, 1:2)
  expect_identical(choices$selected, c(TRUE, FALSE))
  expect_true(dagby_options(full)$alternatives_complete)
  expect_true(all(choices$static_available))
  expect_setequal(names(dagby_temporal(full)), c("T", "Y", "P_t0", "P_t1"))
  other <- dagby_select(full, select = 2, quiet = TRUE)
  expect_setequal(names(dagby_temporal(other)), c("T", "Y", "P_t0"))
  expect_identical(dagby_options(other)$level$variables, 3L)
  expect_identical(dagby_options(other)$level$select, 2L)
  expect_identical(.gd_signature(dagby_select(other, quiet = TRUE)), .gd_signature(other))
  expect_identical(.gd_signature(dagby_temporal(dagby_select(other, quiet = TRUE))),
    .gd_signature(dagby_temporal(other)))
  expect_identical(.gd_signature(dagby_temporal(eval(dagby_call(other)))), .gd_signature(dagby_temporal(other)))
  expect_identical(dagby_options(dagby_construct(graph, Inf, quiet = TRUE))$candidate,
    dagby_options(full)$candidate)
  .cf_cache$clear()
  expect_identical(dagby_options(dagby_construct(graph, quiet = TRUE, size = "largest"))$alternatives, choices)
  withr::local_options(causaliflower.accelerate = FALSE)
  expect_identical(dagby_options(dagby_construct(graph, quiet = TRUE, size = "largest"))$alternatives, choices)
})

test_that("size and strict resize operate in variables with deterministic snapping", {
  graph <- variable_query("T -> A -> B -> Y")
  small <- dagby_construct(graph, 0, quiet = TRUE)
  large <- dagby_construct(graph, quiet = TRUE, size = "largest")
  expect_identical(dagby_options(small)$level$variables, 4L)
  expect_identical(dagby_prune(small, 1, quiet = TRUE), small)
  expect_identical(dagby_grow(large, 1, quiet = TRUE), large)
  graph <- same_size_query()
  large <- dagby_construct(graph, quiet = TRUE, size = "largest")
  small <- dagby_prune(large, 1, quiet = TRUE)
  expect_identical(dagby_options(small)$level$variables, 2L)
  expect_identical(dagby_options(dagby_grow(small, 1, quiet = TRUE))$level$variables, 3L)
  expect_warning(tied <- dagby_construct(graph, 2.5, quiet = TRUE), "equally close")
  expect_identical(dagby_options(tied)$level$variables, 3L)
  other <- dagby_select(large, select = 2, quiet = TRUE)
  expect_identical(dagby_grow(other, 1, quiet = TRUE), other)
  expect_identical(dagby_options(dagby_select(other, 2, quiet = TRUE))$level$variables, 2L)
  expect_no_warning(dagby_construct(graph, 1e300, quiet = TRUE))
  expect_no_warning(dagby_prune(large, 1e300, quiet = TRUE))
  expect_no_warning(dagby_grow(small, 1e300, quiet = TRUE))
})

test_that("view changes preserve candidate identity and defaults are static", {
  graph <- same_size_query()
  static <- dagby_construct(graph, quiet = TRUE, size = "largest")
  temporal <- dagby_select(static, unroll = TRUE, quiet = TRUE)
  expect_identical(dagby_options(temporal)$candidate, dagby_options(static)$candidate)
  expect_identical(dagby_options(temporal)$level$variables, 3L)
  expect_identical(dagby_options(temporal)$level$nodes, 4L)
  expect_true(dagby_options(temporal)$arguments$unroll)
  expect_identical(.gd_signature(temporal), .gd_signature(dagby_temporal(temporal)))
  expect_true(dagby_options(dagby_construct(temporal, quiet = TRUE, size = "largest"))$arguments$unroll)
  expect_true(dagby_options(dagby_grow(dagby_prune(temporal, quiet = TRUE), quiet = TRUE))$arguments$unroll)
  expect_identical(.gd_signature(dagby_prune(temporal, 0, unroll = FALSE, quiet = TRUE)), .gd_signature(static))
  expect_identical(.gd_signature(eval(dagby_call(temporal))), .gd_signature(temporal))
  expect_true(dagby_options(dagby_grow(static, 1, unroll = TRUE, quiet = TRUE))$arguments$unroll)
  endpoint <- variable_query("T_t0 -> T_t1 -> Y",
    c(T_t0 = "T", T_t1 = "T", Y = "Y"), exposures = c("T_t0", "T_t1"))
  expect_error(suppressWarnings(dagby_construct(endpoint, quiet = TRUE, size = "largest")), "unroll = TRUE")
  allowed <- suppressWarnings(dagby_construct(endpoint, unroll = TRUE, quiet = TRUE, size = "largest"))
  expect_identical(dagby_options(allowed)$level$variables, 2L)
  expect_false(dagby_options(allowed)$level$static_available)
  attr(endpoint, "cf_node_map") <- attr(endpoint, "query_record")$node_map
  attr(endpoint, "query_record") <- NULL
  expect_error(suppressWarnings(dagby_construct(endpoint, quiet = TRUE, size = "largest")), "unroll = TRUE")
})

test_that("selectors are validated and old objects migrate only by exact identity", {
  graph <- same_size_query()
  dag <- dagby_construct(graph, quiet = TRUE, size = "largest")
  for (bad in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(dagby_select(dag, select = bad, quiet = TRUE), "positive whole")
  }
  expect_error(dagby_select(dag, select = 3, quiet = TRUE), "Available alternatives.*1, 2")
  expect_error(dagby_select(graph), "dagby_construct")
  plan <- dagby_construct(graph, dry_run = TRUE, quiet = TRUE, size = "largest")
  expect_s3_class(dagby_select(plan, select = 2, quiet = TRUE), "cf_dag")
  expect_s3_class(eval(dagby_call(plan)), "cf_dag_plan")
  old <- .gd_state(dag)
  old$construction$metric <- NULL
  old$construction$candidate <- NULL
  old$construction$level$nodes <- 999L
  old <- .gd_attach(dag, old)
  expect_identical(.gd_signature(dagby_temporal(dagby_select(old, quiet = TRUE))), .gd_signature(dagby_temporal(dag)))
  expect_identical(.gd_signature(eval(dagby_call(old))), .gd_signature(dag))
  call <- dagby_call(dag)
  call$candidate$binding$version <- "unavailable"
  expect_error(eval(call), "unsupported construction version")
})

test_that("plain latent declarations are inherited unless explicitly overridden", {
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; Z [latent]; Z -> T; Z -> Y; T -> Y }")
  inherited <- suppressWarnings(dagby_construct(graph, quiet = TRUE, size = "largest"))
  observed <- suppressWarnings(dagby_construct(graph, latents = NULL, quiet = TRUE, size = "largest"))
  expect_identical(dagitty::latents(inherited), "Z")
  expect_length(dagitty::latents(observed), 0L)
  expect_identical(dagby_options(inherited)$level$variables, dagby_options(observed)$level$variables)
  query <- variable_query("Z -> T; Z -> Y; T -> Y")
  record <- attr(query, "query_record")
  dagitty::latents(query) <- "Z"
  attr(query, "query_record") <- record
  omitted <- dagby_construct(query, quiet = TRUE, size = "largest")
  cleared <- dagby_construct(query, latents = NULL, quiet = TRUE, size = "largest")
  expect_length(dagitty::latents(omitted), 0L)
  expect_identical(.gd_signature(omitted), .gd_signature(cleared))
  expect_identical(dagitty::latents(dagby_source(omitted)), "Z")
})

test_that("replay keeps an unused ranked budget independent of later options", {
  withr::local_options(causaliflower.ranked_steps = "unused")
  dag <- dagby_construct(same_size_query(), quiet = TRUE, size = "largest")
  saved <- dagby_call(dag)
  options(causaliflower.ranked_steps = 20000L)
  expect_identical(.gd_signature(dagby_temporal(eval(saved))), .gd_signature(dagby_temporal(dag)))
})

test_that("refused rule families never claim complete alternatives or session-dependent maxima", {
  graph <- variable_query("R -> T; R -> Y; T -> Y")
  record <- attr(graph, "query_record")
  row <- record$node_map[record$node_map$node == "R", ]
  row$member_id <- "R-other"
  record$node_map <- rbind(record$node_map, row)
  attr(graph, "query_record") <- record
  .cf_cache$clear()
  expect_error(dagby_construct(graph, quiet = TRUE, size = "largest"), "several source members")
  dag <- dagby_construct(graph, 0, quiet = TRUE)
  expect_false(dagby_options(dag)$alternatives_complete)
  expect_false(dagby_options(dag)$sizes_complete)
  expect_match(dagby_options(dag)$alternatives_limit, "several source members")
  expect_identical(dagby_options(dag)$level$variables, 4L)
  expect_error(dagby_construct(graph, quiet = TRUE, size = "largest"), "several source members")
  withr::local_options(causaliflower.accelerate = FALSE)
  expect_error(dagby_construct(graph, quiet = TRUE, size = "largest"), "several source members")
})
