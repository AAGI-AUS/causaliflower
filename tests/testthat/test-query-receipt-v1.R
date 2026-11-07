receipt_query <- function(...) {
  query_graph('barley', treatments = 'N_application',
              outcomes = c('grain_N', 'grain_wt'), quiet = TRUE, ...)
}

check_receipt_identity <- function(g) {
  r <- attr(g, 'query_record')
  expect_identical(r$version, 1L)
  expect_setequal(r$node_map$node, names(g))
  expect_false(anyNA(r$node_map))
  expect_true(all(r$node_map$member_id %in% r$policy$node_id))
  expect_true(all(r$node_map$member_id %in% names(r$aliases)))
  expect_true(all(vapply(split(r$node_map$quantity, r$node_map$node),
                         function(q) length(unique(q)) == 1L, logical(1L))))
  expect_identical(r$policy_scope, 'Barley recovery review plus explicit non-APSIM provenance across crops')
  expect_true(r$policy_applied)
}

test_that('query receipt pins original requests and their actual resolution', {
  columns <- c(nitrate = 'soil_NO3')
  g <- receipt_query(columns = columns, randomised_treatment = TRUE, label = NULL)
  r <- attr(g, 'query_record')
  check_receipt_identity(g)
  expect_identical(r$input$columns, columns)
  expect_identical(r$requests$columns$raw, columns)
  expect_identical(r$requests$columns$resolved, c(nitrate = 'soil_no3'))
  expect_length(r$requests$columns$unresolved$unmatched, 0L)
  expect_identical(r$requests$columns$present, 'soil_no3')
  expect_length(r$requests$columns$absent, 0L)
  expect_true(r$time_unroll_defaulted)
  expect_identical(r$arguments$time_unroll, 'all')
  expect_identical(names(r$arguments$params), names(r$input$params))
  expect_true(r$scope$full_query_defaults)
  expect_setequal(r$requests$treatments$resolved_ids, 'nitrogen_applied')
  expect_setequal(r$requests$outcomes$resolved_ids, c('grain_n', 'grain_wt'))
  expect_setequal(dagitty::exposures(g), attr(g, 'treatments'))
  expect_setequal(dagitty::outcomes(g), attr(g, 'outcomes'))
  expect_true(attr(g, 'randomised_treatment'))
  expect_error(receipt_query(columns = c(columns, absent = 'not_a_real_column')),
               'These columns match no node.*not_a_real_column')
})

test_that('external table receipts retain unmatched inputs without invoking slice guards', {
  graph <- slice_external_query(columns = c(left = 'toy.soil.quantity.left', absent = 'not_a_real_column'))
  receipt <- attr(graph, 'query_record')
  check_receipt_identity(graph)
  expect_identical(receipt$requests$columns$unresolved$unmatched, 'not_a_real_column')
  expect_identical(receipt$requests$columns$resolved, c(left = 'toy.soil.quantity.left'))
  expect_null(receipt$source_order)
  expect_false("source_order" %in% names(receipt))
})

test_that('all label modes retain every canonical node identity', {
  plain <- receipt_query(label = NULL)
  plain_map <- attr(plain, 'query_record')$node_map
  for (label in c('APSIM', 'input', 'readable')) {
    g <- receipt_query(label = label)
    check_receipt_identity(g)
    mapped <- attr(g, 'query_record')$node_map
    expect_equal(mapped[-1], plain_map[-1])
    expect_identical(attr(g, 'query_record')$arguments$label, label)
    edges <- attr(g, 'edges')
    reverse <- setNames(mapped$occurrence, mapped$node)
    canonical_keys <- paste(reverse[edges$from], reverse[edges$to])
    original_edges <- attr(plain, 'edges')
    original_reverse <- setNames(plain_map$occurrence, plain_map$node)
    expect_setequal(canonical_keys,
      paste(original_reverse[original_edges$from], original_reverse[original_edges$to]))
  }
})

test_that('manual scope is recorded without substituting a full query', {
  g <- receipt_query(label = NULL, time_unroll = 'pre_post',
                     params = query_params(query_depth = 1L, show_hidden_nodes = FALSE))
  check_receipt_identity(g)
  r <- attr(g, 'query_record')
  expect_identical(r$arguments$time_unroll, 'pre_post')
  expect_false(r$time_unroll_defaulted)
  expect_identical(r$scope$query_depth, 1L)
  expect_false(r$scope$full_query_defaults)
  expect_lt(length(names(g)), length(attr(g, 'full_vertices')))
  expect_identical(r$arguments$params$show_hidden_nodes, FALSE)
  expect_true(any(grepl('_pre_treatment|_post_treatment|_t[0-9]+$', r$node_map$occurrence)))
})

test_that('recombination and detail contraction retain original memberships', {
  g <- receipt_query(label = 'readable',
    params = query_params(node_recombination = TRUE))
  check_receipt_identity(g)
  r <- attr(g, 'query_record')
  expect_true(r$scope$node_recombination)
  expect_true(any(lengths(split(r$node_map$member_id, r$node_map$node)) > 1L))
  coarse <- slice_external_query(label = NULL, time_unroll = 'pre_post',
    params = query_params(apsim_nodes_detail = 'subnode', force_acyclic = TRUE))
  check_receipt_identity(coarse)
  expect_true(any(attr(coarse, 'query_record')$node_map$quantity !=
                    attr(coarse, 'query_record')$node_map$member_id))
})

test_that('readable_dag retains a construction receipt and endpoint identity', {
  g <- receipt_query(label = NULL)
  h <- readable_dag(g)
  check_receipt_identity(h)
  expect_identical(attr(h, 'query_record')$node_map[-1],
                   attr(g, 'query_record')$node_map[-1])
  expect_identical(attr(h, 'query_record')$input, attr(g, 'query_record')$input)
  expect_identical(attr(h, 'query_record')$display_transform[[1]]$operation, 'readable_dag')
  expect_setequal(dagitty::exposures(h), attr(h, 'treatments'))
  expect_setequal(dagitty::outcomes(h), attr(h, 'outcomes'))
})

test_that('manual queries use package tables despite similarly named workspace objects', {
  keys <- c('apsim_dag_nodes', 'node_aliases', 'management_defaults', 'crop_phases', 'phase_of')
  baseline <- receipt_query(label = NULL)
  before <- lapply(keys, function(k) if (exists(k, .GlobalEnv, inherits = FALSE)) get(k, .GlobalEnv) else NULL)
  existed <- vapply(keys, exists, logical(1L), envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    for (i in seq_along(keys)) {
      if (existed[i]) assign(keys[i], before[[i]], .GlobalEnv)
      else if (exists(keys[i], .GlobalEnv, inherits = FALSE)) rm(list = keys[i], envir = .GlobalEnv)
    }
  })
  for (k in keys) assign(k, 'must_not_be_used', .GlobalEnv)
  current <- receipt_query(label = NULL)
  expect_identical(attr(current, 'dag_text'), attr(baseline, 'dag_text'))
  expect_equal(attr(current, 'query_record')$node_map, attr(baseline, 'query_record')$node_map)
})
