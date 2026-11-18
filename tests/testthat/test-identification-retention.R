test_that("construction retains observed treatment parents when adjustment is unavailable", {
  g <- dagitty::dagitty('dag { X1 [exposure] X2 [exposure] Y [outcome]
    X1 -> L -> X2 -> M -> Y; L -> M }')
  result <- suppressWarnings(dagby_construct(g, quiet = TRUE, size = "largest"))
  guard <- dagby_options(result)$identification
  expect_true(guard$applied)
  expect_identical(guard$parents, 'L')
  expect_true('L' %in% names(result))
  expect_false('L' %in% guard$added_nodes)
  expect_identical(guard$status, "preserved by retaining exact treatment parents")
  expect_length(dagitty::adjustmentSets(result, effect = 'total'), 0L)
  expect_true(attr(get_identification(result, quiet = TRUE), 'analyses')$total$identifiable)
  expect_identical(as.character(dagby_source(result)), as.character(g))
  expect_true('L' %in% names(suppressWarnings(dagby_construct(g, 0, quiet = TRUE))))
})

test_that("the safeguard preserves exact temporal parents and observation assumptions", {
  g <- dagitty::dagitty('dag { X1 [exposure] X2 [exposure] Y [outcome]
    X1 -> L_t0 -> X2 -> M_t1 -> M_t2 -> Y; L_t0 -> M_t2; L_t1 -> Y }')
  nodes <- names(g)
  attr(g, "cf_node_map") <- data.frame(node = nodes, quantity = .strip_slice(nodes),
    member_id = .strip_slice(nodes), occurrence = nodes)
  expect_error(suppressWarnings(dagby_construct(g, quiet = TRUE, size = "largest")),
    "joint total identification safeguard", class = "cd_static_unavailable")
  result <- suppressWarnings(dagby_construct(g, quiet = TRUE, unroll = TRUE, size = "largest"))
  expect_setequal(dagby_options(result)$identification$parents, 'L_t0')
  expect_true('L_t0' %in% names(result))
  expect_false('L' %in% names(result))
  expect_false(dagby_options(result)$summary$applied)
  expect_match(dagby_options(result)$level$static_reason, "joint total identification safeguard")
  expect_length(dagitty::latents(result), 0L)
})

test_that("ordinary adjustment needs no retention", {
  g <- dagitty::dagitty('dag { X [exposure] Y [outcome] U -> X; U -> Y; X -> M -> Y }')
  result <- suppressWarnings(dagby_construct(g, quiet = TRUE, size = "largest"))
  expect_false(dagby_options(result)$identification$applied)
  expect_setequal(names(result), c('X', 'Y', 'U', 'M'))
  expect_true(dagitty::isAdjustmentSet(result, 'U'))
})

test_that("latent and bidirected references do not receive an unsupported guarantee", {
  g <- dagitty::dagitty('dag { X [exposure] Y [outcome] X -> Y; X <-> Y }')
  result <- suppressWarnings(dagby_construct(g, quiet = TRUE, size = "largest"))
  expect_false(dagby_options(result)$identification$applied)
  expect_match(dagby_options(result)$identification$status, 'not certified')
  expect_length(dagitty::adjustmentSets(result, effect = 'total'), 0L)
  g <- dagitty::dagitty('dag { X [exposure] Y [outcome] U [latent] U -> X; U -> Y; X -> Y }')
  result <- suppressWarnings(dagby_construct(g, latents = 'U', quiet = TRUE, size = "largest"))
  expect_false(dagby_options(result)$identification$applied)
  expect_identical(dagby_options(result)$identification$status, 'not required: a total adjustment option was selected')
  expect_identical(dagby_checks(result)$total$availability, 'needs an unmeasured variable')
  g <- dagitty::dagitty('dag { X1 [exposure] X2 [exposure] Y [outcome]
    X1 -> L -> X2 -> M -> Y; L -> M }')
  result <- suppressWarnings(dagby_construct(g, latents = 'L', quiet = TRUE, size = "largest"))
  guard <- dagby_options(result)$identification
  expect_true(guard$applied)
  expect_identical(guard$unmeasured_parents, 'L')
  expect_match(guard$status, '^not certified: treatment parent\\(s\\) L are declared unmeasured')
  expect_identical(dagitty::latents(result), 'L')
  expect_false(dagby_options(result)$summary$applied)
})

test_that("construction rejects overlapping endpoint declarations", {
  g <- dagitty::dagitty('dag { X [exposure,outcome] Y [outcome] X -> Y }')
  expect_error(dagby_construct(g, keep_nodes = c('X', 'Y'), quiet = TRUE, size = "largest"), 'requires disjoint endpoint sets')
  expect_error(dagby_construct(g, 0, quiet = TRUE), 'requires disjoint endpoint sets')
})
