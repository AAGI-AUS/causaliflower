test_that("construction retains observed treatment parents when adjustment is unavailable", {
  g <- dagitty::dagitty('dag { X1 [exposure] X2 [exposure] Y [outcome]
    X1 -> L -> X2 -> M -> Y; L -> M }')
  result <- suppressWarnings(construct_dag(g, quiet = TRUE))
  guard <- dag_options(result)$identification
  expect_true(guard$applied)
  expect_identical(guard$parents, 'L')
  expect_true('L' %in% names(result))
  expect_true('L' %in% guard$added_nodes)
  expect_identical(guard$status, "preserved by retaining exact treatment parents")
  expect_length(dagitty::adjustmentSets(result, effect = 'total'), 0L)
  expect_true(attr(get_identification(result, quiet = TRUE), 'analyses')$total$identifiable)
  expect_true(isTRUE(.gd_state(result)$controls$simplify == FALSE))
  expect_identical(as.character(original_query(result)), as.character(g))
})

test_that("the safeguard preserves exact temporal parents and observation assumptions", {
  g <- dagitty::dagitty('dag { X1 [exposure] X2 [exposure] Y [outcome]
    X1 -> L_t0 -> X2 -> M_t1 -> M_t2 -> Y; L_t0 -> M_t2; L_t1 -> Y }')
  result <- suppressWarnings(construct_dag(g, quiet = TRUE))
  expect_setequal(dag_options(result)$identification$parents, 'L_t0')
  expect_true('L_t0' %in% names(result))
  expect_false('L' %in% names(result))
  expect_false(.gd_state(result)$summary$applied)
  expect_match(.gd_state(result)$summary$reason, "joint total identification safeguard")
  expect_length(dagitty::latents(result), 0L)
  restored <- suppressWarnings(restore_nodes(result, 'L_t1'))
  expect_true(all(c('L_t0', 'L_t1') %in% names(restored)))
  expect_length(dagitty::latents(restored), 0L)
})

test_that("ordinary adjustment keeps the established selection policy", {
  g <- dagitty::dagitty('dag { X [exposure] Y [outcome] U -> X; U -> Y; X -> M -> Y }')
  result <- suppressWarnings(construct_dag(g, quiet = TRUE))
  expect_false(dag_options(result)$identification$applied)
  expect_setequal(names(result), c('X','Y','U','M'))
  expect_true(dagitty::isAdjustmentSet(result, 'U'))
})

test_that("latent and bidirected references do not receive an unsupported guarantee", {
  g <- dagitty::dagitty('dag { X [exposure] Y [outcome] X -> Y; X <-> Y }')
  result <- suppressWarnings(construct_dag(g, quiet = TRUE))
  expect_false(dag_options(result)$identification$applied)
  expect_match(dag_options(result)$identification$status, 'not certified')
  expect_length(dagitty::adjustmentSets(result, effect = 'total'), 0L)
  # A declared unmeasured confounder no longer removes the total option {U}, so
  # no retention is needed. The selected set's need for U is reported instead.
  g <- dagitty::dagitty('dag { X [exposure] Y [outcome] U [latent] U -> X; U -> Y; X -> Y }')
  result <- suppressWarnings(construct_dag(g, latents = 'U', quiet = TRUE))
  expect_false(dag_options(result)$identification$applied)
  expect_identical(dag_options(result)$identification$status, 'not required: a total adjustment option was selected')
  expect_identical(dag_checks(result)$total$availability, 'needs unmeasured: U')
  # A retained treatment parent declared unmeasured keeps the exact occasions
  # but leaves the joint total effect guarantee uncertified.
  g <- dagitty::dagitty('dag { X1 [exposure] X2 [exposure] Y [outcome]
    X1 -> L -> X2 -> M -> Y; L -> M }')
  result <- suppressWarnings(construct_dag(g, latents = 'L', quiet = TRUE))
  guard <- dag_options(result)$identification
  expect_true(guard$applied)
  expect_identical(guard$unmeasured_parents, 'L')
  expect_match(guard$status, '^not certified: treatment parent\\(s\\) L are declared unmeasured')
  expect_identical(dagitty::latents(result), 'L')
  expect_false(.gd_state(result)$summary$applied)
})

test_that("construction rejects overlapping endpoint declarations before certification", {
  g <- dagitty::dagitty('dag { X [exposure,outcome] Y [outcome] X -> Y }')
  expect_error(construct_dag(g, keep_nodes = c('X','Y'), quiet = TRUE),
               'requires disjoint endpoint sets')
})
