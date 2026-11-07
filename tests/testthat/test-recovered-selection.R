selection_fixture <- function(text, max_sets = 100L) {
  g <- dagitty::dagitty(text)
  attr(g, "fixture_attribute") <- list(meaning = "preserve exactly", time = 3L)
  found <- .gd_adjustment_sets(g, "direct", max_sets)
  candidates <- .gd_candidates(found$sets, g, list())
  views <- lapply(seq_len(nrow(candidates)), function(i) {
    focused <- .gd_focus(g, candidates$concepts[[i]])
    list(graph = focused$graph, temporal = focused$graph, focus = focused$ledger,
         introduced_edges = focused$introduced_edges, summary = NULL)
  })
  names(views) <- candidates$id
  recovery <- .gd_recover(g, character(), apply_policy = FALSE)
  list(query = g, reference = g, analysis = g, selection = g, candidates = candidates,
    sets = found$sets, views = views, dags = lapply(views, `[[`, "graph"), recovery = recovery,
    availability = list(mode = "scm", statement = "fixture", ledger = data.frame()),
    search = list(effect = "direct", complete = found$complete, limit = max_sets,
                  scope = "ancestral", mode = "sets"), components = list(), observed_ids = character(),
    controls = list(crop = "fixture", water_model = "SoilWater", disambiguate = NULL,
                    simplify = FALSE, keep_common_causes = FALSE))
}
selection_record <- function(selection = "smallest") {
  values <- list(treatments = "T", outcomes = "Y", selection = selection)
  list(original_call = quote(get_dag()), resolved_call = as.call(c(list(as.name("get_dag")), values)),
       values = values, responses = data.frame(field = character(), prompt = character(), response = character()))
}

test_that("grouped selection preserves source and switches a saved alternative", {
  b <- selection_fixture("dag { T [exposure]; Y [outcome]; Z -> T; Z -> W -> Y; T -> Y }")
  before <- serialize(b$query, NULL)
  out <- NULL
  expect_output(out <- suppressWarnings(.gd_select_result(b, c(1, 1), selection_record())), "Causal query")
  expect_s3_class(out, "dagitty")
  expect_identical(serialize(original_query(out), NULL), before)
  expect_equal(dag_options(out)$group, c(1L, 1L))
  expect_equal(dag_call(out)$selection, c(1L, 1L))
  expect_equal(dag_checks(out)$direct$status, "changed")
  chosen <- Filter(function(z) z$chosen, dag_checks(out)$direct$records)
  expect_true(all(vapply(chosen, function(z) z$status == "preserved", logical(1))))
  other <- NULL
  expect_output(other <- suppressWarnings(select_dag(out, c(1, 2))), "selected")
  expect_false(identical(sort(names(out)), sort(names(other))))
  expect_identical(original_query(out), original_query(other))
  expect_equal(dag_call(other)$selection, c(1L, 2L))
  expect_error(select_dag(out, c(2, 1)), "absent")
  expect_error(select_dag(out, 2), "c\\(group, option\\)")
  # Explicit input means there is no notification claiming answers were collected.
  expect_message(suppressWarnings(.gd_select_result(b, c(1, 1), selection_record())), NA)
})

test_that("an empty adjustment set is not a failed search", {
  b <- selection_fixture("dag { T [exposure]; Y [outcome]; T -> Y }")
  expect_length(b$sets, 1L)
  expect_length(b$sets[[1L]], 0L)
  x <- NULL
  expect_warning(capture.output(x <- .gd_select_result(b, "smallest", selection_record())), NA)
  expect_s3_class(x, "dagitty")
  expect_equal(dag_options(x)$n_resolved, 0L)
  failed <- selection_fixture("dag { T [exposure]; Y [outcome]; U [latent]; U -> T; U -> Y; T -> Y }")
  expect_length(failed$sets, 0L)
  expect_warning(capture.output(y <- .gd_select_result(failed, "smallest", selection_record())), "No adjustment set found")
  expect_s3_class(y, "dagitty")
  expect_equal(dag_options(y)$status, "No adjustment set found")
})

test_that("new apparent adjustment results are checked against the reference", {
  source <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; U -> T; U -> Y; T -> Y }")
  wrong <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; T -> Y }")
  result <- .gd_check_family(source, wrong, list("U"), list(character()),
                            "total", TRUE, list("U"))
  expect_equal(result$status, "changed")
  expect_equal(result$records[[1]]$status, "not represented")
  expect_equal(result$introduced[[1]]$status, "not valid on reference")
})

test_that("validity and minimality and a capped prefix remain different", {
  source <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; U -> T; U -> Y; V; T -> Y }")
  result <- .gd_check_family(source, source, list(c("U", "V")), list("U"),
                            "total", TRUE, list(c("U", "V")))
  expect_equal(result$records[[1]]$status, "valid but not minimal")
  limited <- .gd_check_family(source, source, list("U"), list(), "total", FALSE, list("U"))
  expect_equal(limited$records[[1]]$status, "minimality not verified")
  expect_equal(limited$status, "not fully checked")
})

test_that("condensation does not invent a temporal lift", {
  source <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; U_t0 -> T; U_t1 -> Y; T -> Y }")
  display <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; U -> T; U -> Y; T -> Y }")
  expect_null(.gd_map_set("U", display, source))
  expect_equal(.gd_map_set(c("U_t0", "U_t1"), source, display), "U")
})

test_that("missing occasions cannot be silently renamed as later occasions", {
  from <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A_t2 -> T; A_t2 -> Y; T -> Y }")
  to <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; A_t3 -> T; A_t3 -> Y; T -> Y }")
  expect_null(.gd_map_set("A_t2", from, to))
  checked <- .gd_check_family(from, to, list("A_t2"), list("A_t3"), "total", TRUE, list("A_t2"))
  expect_equal(checked$records[[1]]$status, "not represented")
  expect_equal(checked$status, "changed")
})

test_that("an unobserved set is not available adjustment", {
  g <- dagitty::dagitty("dag { T [exposure]; Y [outcome]; U [latent]; U -> T; U -> Y; T -> Y }")
  expect_false(.gd_validate_sets(g, list("U"), "total"))
  expect_false(.gd_validate_sets(g, list("U"), "direct"))
})
