diagnostic_choice_query <- function() {
  graph <- dagitty::dagitty("dag {
    T [exposure]; Y [outcome];
    Q1_t2 -> Q3_t2; Q2_t0 -> Q9_t1; Q2_t0 -> S2; Q2_t0 -> Y;
    Q2_t1 -> Q4_t1; Q2_t1 -> Q6_t1; Q2_t1 -> Y; Q3_t2 <-> T;
    Q4_t1 -> Q1_t2; Q5_t1 -> Q3_t2; Q5_t1 -> Q4_t1; Q5_t1 -> Q6_t1;
    Q5_t1 -> Q8_t1; Q6_t0 -> Q3_t2; Q6_t0 -> Q8_t1; Q6_t0 -> Q9_t1;
    Q6_t1 -> Q3_t2; Q7_t0 -> Q3_t2; Q7_t0 -> Q5_t1; Q7_t0 -> Q6_t0;
    Q7_t0 -> Q6_t1; Q7_t0 -> Q9_t1; Q7_t0 -> S2; Q7_t0 -> Y;
    Q8_t0 -> Q1_t2; Q8_t0 -> Q2_t1; Q8_t0 -> Q3_t2; Q8_t0 -> Q5_t1;
    Q8_t0 -> Q6_t1; Q8_t0 -> Q8_t1; Q8_t0 -> S2; Q8_t1 -> Q4_t1;
    Q8_t1 -> Q9_t1; Q8_t1 -> Y; Q8_t1 <-> T; Q9_t1 -> Q3_t2; Q9_t1 -> Y;
    S1 -> Q1_t2; S1 -> Q2_t0; S1 -> Q3_t2; S1 -> Q4_t1; S1 -> Q5_t1;
    S1 -> Q6_t0; S1 -> Q6_t1; S1 -> Q7_t0; S1 -> S2; S1 -> T;
    S2 -> Q1_t2; S2 -> Q2_t1; S2 -> Q3_t2; S2 -> Q4_t1; S2 -> Q6_t0;
    S2 -> Q6_t1; S2 -> Y; T -> Q2_t0; T -> Q2_t1; T -> Q3_t2;
    T -> Q4_t1; T -> Q5_t1; T -> Q6_t0; T -> Q7_t0; T -> Q9_t1; T -> S2;
    Y -> Q1_t2; Y -> Q3_t2
  }")
  members <- c(Q1_t2 = "Q1", Q2_t0 = "Q2", Q2_t1 = "Q2", Q3_t2 = "Q3",
    Q4_t1 = "Q4", Q5_t1 = "Q5", Q6_t0 = "Q6", Q6_t1 = "Q6", Q7_t0 = "Q7",
    Q8_t0 = "Q8", Q8_t1 = "Q8", Q9_t1 = "Q9", S1 = "S1", S2 = "S2", T = "T", Y = "Y")
  nodes <- names(graph)
  attr(graph, "query_record") <- list(version = 1L,
    arguments = list(crop = "fixture", columns = NULL, params = list(query_depth = Inf)),
    requests = list(columns = list(resolved = list())),
    node_map = data.frame(node = nodes, quantity = unname(members[nodes]),
      member_id = unname(members[nodes]), occurrence = nodes),
    policy = data.frame(node_id = unique(unname(members)), recovered_node = "yes"))
  class(graph) <- c("dag_by_query", class(graph))
  graph
}

diagnostic_choice_temporal <- function() {
  graph <- dagitty::dagitty("dag {
    T [exposure]; Y [outcome];
    Q2_t0 -> Q9_t1; Q2_t0 -> S2; Q2_t0 -> Y; Q2_t0 <-> Q5_t1;
    Q2_t0 <-> Q6_t0; Q2_t0 <-> Q7_t0; Q2_t0 <-> S2; Q2_t0 <-> T;
    Q2_t1 -> Y; Q2_t1 <-> Q5_t1; Q2_t1 <-> Q9_t1; Q2_t1 <-> S2; Q2_t1 <-> Y;
    Q5_t1 -> Q9_t1; Q5_t1 -> Y; Q5_t1 <-> Q6_t0; Q5_t1 <-> Q7_t0;
    Q5_t1 <-> Q9_t1; Q5_t1 <-> S2; Q5_t1 <-> T; Q5_t1 <-> Y;
    Q6_t0 -> Q9_t1; Q6_t0 -> Y; Q6_t0 <-> Q7_t0; Q6_t0 <-> S2; Q6_t0 <-> T;
    Q7_t0 -> Q5_t1; Q7_t0 -> Q6_t0; Q7_t0 -> Q9_t1; Q7_t0 -> S2; Q7_t0 -> Y;
    Q7_t0 <-> S2; Q7_t0 <-> T; Q9_t1 -> Y; Q9_t1 <-> S2; Q9_t1 <-> T; Q9_t1 <-> Y;
    S2 -> Q2_t1; S2 -> Q6_t0; S2 -> Y; S2 <-> T; S2 <-> Y;
    T -> Q2_t0; T -> Q2_t1; T -> Q5_t1; T -> Q6_t0; T -> Q7_t0; T -> Q9_t1; T -> S2;
    T <-> Y
  }")
  attr(graph, "input_labels") <- c(Q2 = "Q2", Q2 = "Q2", Q5 = "Q5", Q6 = "Q6",
    Q7 = "Q7", Q9 = "Q9", S2 = "S2", T = "T", Y = "Y")
  graph
}

test_that("advertised diagnostic sizes remain reachable across cache states", {
  withr::local_options(causaliflower.ranked_steps = 10L)
  graph <- diagnostic_choice_query()
  expected <- .gd_signature(diagnostic_choice_temporal())
  on.exit(.cf_cache$clear(), add = TRUE)
  for (mode in c("cold", "warm", "disabled")) {
    withr::local_options(causaliflower.accelerate = mode != "disabled")
    .cf_cache$clear()
    run <- function(expr) {
      if (mode == "cold") .cf_cache$clear()
      suppressWarnings(expr)
    }
    if (mode == "warm") suppressWarnings(dagby_construct(graph, quiet = TRUE, unroll = TRUE, size = "largest"))
    small <- run(dagby_construct(graph, 0, quiet = TRUE, unroll = TRUE))
    expect_identical(dagby_options(small)$level$variables, 7L)
    expect_identical(dagby_options(small)$sizes$variables, c(7L, 8L, 9L, 12L, 13L))
    expect_true(dagby_options(small)$sizes_complete)
    exact <- run(dagby_construct(graph, 8, quiet = TRUE, unroll = TRUE))
    above <- run(dagby_construct(graph, 9, quiet = TRUE, unroll = TRUE))
    calls <- list(exact,
      run(dagby_select(small, size = 8, quiet = TRUE)),
      run(dagby_grow(small, 1, quiet = TRUE)),
      run(dagby_prune(above, 1, quiet = TRUE)),
      run(dagby_select(exact, quiet = TRUE)),
      run(dagby_select(exact, select = 1, quiet = TRUE)),
      run(eval(dagby_call(exact))))
    for (dag in calls) {
      expect_identical(dagby_options(dag)$level$variables, 8L)
      expect_identical(.gd_signature(dagby_temporal(dag)), expected)
      expect_identical(dagby_options(dag)$level$cap, 0)
      expect_false(dagby_options(dag)$level$certified)
      expect_identical(dagby_options(dag)$arguments$steps, 10L)
      expect_identical(dagby_options(dag)$candidate, dagby_options(exact)$candidate)
    }
    expect_identical(dagby_options(run(dagby_prune(exact, 1, quiet = TRUE)))$level$variables, 7L)
    expect_identical(dagby_options(run(dagby_grow(exact, 1, quiet = TRUE)))$level$variables, 9L)
  }
  snapped <- function(size) {
    .cf_cache$clear()
    warnings <- character()
    dag <- withCallingHandlers(suppressMessages(dagby_construct(graph, size, unroll = TRUE)),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      })
    list(dag = dag, warnings = warnings)
  }
  tie <- snapped(10.5)
  expect_true(any(grepl("equally close", tie$warnings)))
  expect_identical(dagby_options(tie$dag)$level$variables, 12L)
  expect_true(dagby_options(tie$dag)$sizes_complete)
  near <- snapped(10.6)
  expect_false(any(grepl("equally close", near$warnings)))
  expect_true(any(grepl("No eligible total adjustment set", near$warnings)))
  expect_identical(dagby_options(near$dag)$level$variables, 12L)
  expect_false(dagby_options(near$dag)$sizes_complete)
  middle <- snapped(8.5)
  expect_true(any(grepl("Sizes 8 and 9 are equally close", middle$warnings)))
  expect_identical(dagby_options(middle$dag)$level$variables, 9L)
  .cf_cache$clear()
  expect_message(suppressWarnings(dagby_construct(graph, 8.7, quiet = NULL, unroll = TRUE)),
    "Requested 8.7 variables\\. Returned 9.*Next smaller size: 8 variables\\. Use dagby_prune\\(dag, 1\\)")
  static <- suppressWarnings(dagby_construct(graph, 8, quiet = TRUE))
  expect_false(dagby_options(static)$arguments$unroll)
  expect_true(dagby_options(static)$level$static_available)
  expect_identical(.gd_signature(dagby_temporal(static)), expected)
  expect_identical(.gd_signature(dagby_temporal(suppressWarnings(dagby_select(static, unroll = TRUE, quiet = TRUE)))), expected)
})

test_that("recipe diagnostic choices match the constructor across cache states", {
  withr::local_options(causaliflower.ranked_steps = 10L)
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  testthat::local_mocked_bindings(.dagby_query_source = function(...) diagnostic_choice_query(), .package = "causaliflower")
  direct <- suppressWarnings(dagby_construct(diagnostic_choice_query(), unroll = TRUE, quiet = TRUE, size = "largest"))
  for (mode in c("cold", "warm", "disabled")) {
    withr::local_options(causaliflower.accelerate = mode != "disabled")
    if (mode != "warm") .cf_cache$clear()
    guided <- suppressWarnings(dagby("fixture", "T", "Y", FALSE,
      size = "largest", unroll = TRUE) |> dagby_query(quiet = TRUE) |>
      dagby_construct(quiet = TRUE))
    expect_identical(dagby_options(guided)$level$variables, 13L)
    expect_false(dagby_options(guided)$level$smallest)
    expect_identical(.gd_signature(dagby_temporal(guided)), .gd_signature(dagby_temporal(direct)))
    expect_identical(dagby_options(guided), dagby_options(direct))
    expect_identical(dagby_checks(guided), dagby_checks(direct))
    expect_identical(.gd_signature(dagby_temporal(suppressWarnings(eval(dagby_call(guided))))),
      .gd_signature(dagby_temporal(guided)))
    explicit <- suppressWarnings(dagby("fixture", "T", "Y", FALSE,
      size = "smallest", unroll = TRUE) |> dagby_query(quiet = TRUE) |>
      dagby_construct(quiet = TRUE))
    expect_identical(dagby_options(explicit)$level$variables, 7L)
  }
})

test_that("safe discovery skips do not claim an unproved smaller neighbour", {
  withr::local_options(causaliflower.ranked_steps = "unused")
  .cf_cache$clear()
  on.exit(.cf_cache$clear(), add = TRUE)
  graph <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; I -> T; V -> Y; T -> Y }")
  messages <- character()
  dag <- withCallingHandlers(suppressWarnings(dagby_construct(graph, 3.7)), message = function(message) {
    messages <<- c(messages, conditionMessage(message))
    invokeRestart("muffleMessage")
  })
  expect_identical(dagby_options(dag)$level$variables, 4L)
  expect_false(dagby_options(dag)$sizes_complete)
  expect_null(dagby_options(dag)$sizes_steps)
  expect_true(any(grepl("Smaller, if any: dagby_prune(dag, 1).", messages, fixed = TRUE)))
  expect_false(any(grepl("Next smaller size:", messages, fixed = TRUE)))
})
