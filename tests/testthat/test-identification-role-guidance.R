################################################################################
# test-identification-role-guidance.R


# Check overlapping identification roles

# These roles belong to the private query record
test_that("separate strategies retain truthful overlapping roles", {
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] A -> B -> T",
    "A -> Y T -> Y }"
  ))
  result <- identification_result(dag)
  roles <- attr(result, "roles")
  instrument <- identification_strategy(result, "^instrumental variable$")

  expect_setequal(
    roles[variable == "B", role],
    c("adjustment variable", "instrument")
  )
  expect_true(dagitty::isAdjustmentSet(dag, Z = "B"))
  expect_true(instrument$feasible)
  expect_identical(instrument$variable_set, "{B} conditional on {A}")
  expect_identical(
    roles[variable == "B" & role == "instrument", conditioning_guidance],
    "used by the instrumental variable strategy; covariate adjustment is assessed separately"
  )
  expect_identical(
    roles[variable == "B" & role == "adjustment variable", conditioning_guidance],
    roles[variable == "A" & role == "adjustment variable", conditioning_guidance]
  )
})
