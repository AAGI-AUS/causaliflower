construction_check_status <- function(source, final = source, max_sets = 10L) {
  prepared <- .cd_prepare(dagitty::dagitty(source), NULL, NULL)
  effects <- c("direct", "total")
  search <- stats::setNames(lapply(effects, function(effect)
    .cd_search(prepared, effect, 10L)), effects)
  choices <- stats::setNames(lapply(effects, function(effect)
    .cd_select(search[[effect]]$table, "auto", effect)), effects)
  .cd_checks(prepared, dagitty::dagitty(final), search, choices, max_sets)
}

test_that("a preserved selected set does not certify unseen final alternatives", {
  source <- "dag { T [exposure]; Y [outcome]; B -> T; B -> Y; T -> Y; A }"
  final <- "dag { T [exposure]; Y [outcome]; A -> T; A -> B -> Y; T -> Y }"
  checked <- construction_check_status(source, final, max_sets = 1L)
  for (effect in c("direct", "total")) {
    check <- checked[[effect]]
    expect_true(check$reference_complete)
    expect_false(check$final_complete)
    expect_identical(check$required, list("B"))
    expect_identical(check$final_sets, list("B"))
    expect_identical(check$records[[1L]]$status, "preserved")
    expect_length(check$introduced, 0L)
    expect_identical(check$catalogue_status, "not fully checked")
    expect_identical(check$status, "not fully checked")
    full <- .gd_adjustment_sets(dagitty::dagitty(final), effect, 10L)
    expect_true(full$complete)
    expect_true(any(vapply(full$sets, identical, logical(1), "A")))
    expect_false(.gd_validate_sets(dagitty::dagitty(source), list("A"), effect)[[1L]])
  }
  complete <- construction_check_status(source, final)
  expect_identical(vapply(complete, `[[`, character(1), "status"),
                   c(direct = "changed", total = "changed"))
})

test_that("a missing reference family does not hide an invalid introduced empty set", {
  source <- "dag { T [exposure]; Y [outcome]; T <-> Y; T -> Y }"
  final <- "dag { T [exposure]; Y [outcome]; T -> Y }"
  checked <- construction_check_status(source, final)
  for (check in checked) {
    expect_true(check$reference_complete)
    expect_true(check$final_complete)
    expect_length(check$required, 0L)
    expect_length(check$reference_sets, 0L)
    expect_identical(check$final_sets, list(character()))
    expect_identical(check$introduced[[1L]]$status, "not valid on reference")
    expect_identical(check$catalogue_status, "changed")
    expect_identical(check$status, "changed")
  }
})

test_that("unchanged empty and absent families retain their distinct statuses", {
  empty <- construction_check_status("dag { T [exposure]; Y [outcome]; T -> Y }")
  absent <- construction_check_status("dag { T [exposure]; Y [outcome]; T <-> Y; T -> Y }")
  for (effect in c("direct", "total")) {
    expect_identical(empty[[effect]]$required, list(character()))
    expect_identical(empty[[effect]]$status, "preserved")
    expect_length(absent[[effect]]$required, 0L)
    expect_length(absent[[effect]]$introduced, 0L)
    expect_identical(absent[[effect]]$status, "no eligible reference set")
  }
})
