################################################################################
# test-identification-conditional-frontdoor.R


# Check complete conditional front door strategies and their identifying formula

conditional_frontdoor_result <- function(dag){ # Run the supported total effect check
  return(id_conditional_frontdoor(dag, id_admg(dag), "T", "Y"))
}

test_that("conditional front door finds the observed covariate and mediator", {
  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    C -> T
    C -> M
    U -> T
    U -> Y
    T -> M
    M -> Y
  }")
  result <- conditional_frontdoor_result(dag)
  expect_true(result$feasible)
  expect_identical(result$sets, list(list(mediators = "M", covariates = "C")))
  expect_length(dagitty::adjustmentSets(dag), 0L)
  expect_false(id_frontdoor(dag, id_admg(dag), "T", "Y")$feasible)
  expect_true(id_solve("Y", "T", id_admg(dag))$id)

  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    C -> T
    C -> M
    C -> Y
    U -> T
    U -> Y
    T -> M
    M -> Y
  }")
  expect_identical(conditional_frontdoor_result(dag)$sets, result$sets)
})

test_that("conditional front door retains complete alternative covariate sets", {
  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    C1 -> C2
    C2 -> T
    C1 -> M
    U -> T
    U -> Y
    T -> M
    M -> Y
  }")
  result <- conditional_frontdoor_result(dag)
  expect_identical(result$sets, list(
    list(mediators = "M", covariates = "C1"),
    list(mediators = "M", covariates = "C2")
  ))
})

test_that("conditional front door retains all smallest mediator alternatives", {
  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    C -> T
    C -> M1
    U -> T
    U -> Y
    T -> M1
    M1 -> M2
    M2 -> Y
  }")
  expect_identical(conditional_frontdoor_result(dag)$sets, list(
    list(mediators = "M1", covariates = "C"),
    list(mediators = "M2", covariates = "C")
  ))
})

test_that("conditional front door checks parallel mediators jointly", {
  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    C -> T
    C -> M1
    C -> M2
    U -> T
    U -> Y
    T -> M1
    T -> M2
    M1 -> Y
    M2 -> Y
  }")
  result <- conditional_frontdoor_result(dag)
  expect_identical(result$sets, list(
    list(mediators = c("M1", "M2"), covariates = "C")
  ))
  expect_true(id_solve("Y", "T", id_admg(dag))$id)
})

test_that("ordinary front door does not gain unnecessary conditioning variables", {
  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    C -> T
    C -> Y
    U -> T
    U -> Y
    T -> M
    M -> Y
  }")
  expect_true(id_frontdoor(dag, id_admg(dag), "T", "Y")$feasible)
  expect_identical(conditional_frontdoor_result(dag), list(feasible = FALSE, sets = list()))

  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    C -> T
    C -> M2
    U -> T
    U -> Y
    T -> M1
    M1 -> M2
    M2 -> Y
  }")
  expect_identical(id_frontdoor(dag, id_admg(dag), "T", "Y")$sets, list("M1"))
  expect_identical(conditional_frontdoor_result(dag)$sets, list(
    list(mediators = "M2", covariates = "C")
  ))
})

test_that("conditional front door excludes treatment descendants as covariates", {
  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    U -> T
    U -> Y
    T -> C
    C -> M
    M -> Y
  }")
  expect_identical(conditional_frontdoor_result(dag), list(feasible = FALSE, sets = list()))
  expect_true(id_frontdoor(dag, id_admg(dag), "T", "Y")$feasible)
})

test_that("conditional front door does not use latent covariates or mediators", {
  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    C [latent]
    C -> T
    C -> M
    U -> T
    U -> Y
    T -> M
    M -> Y
  }")
  expect_identical(conditional_frontdoor_result(dag), list(feasible = FALSE, sets = list()))

  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    M [latent]
    C -> T
    C -> M
    U -> T
    U -> Y
    T -> M
    M -> Y
  }")
  expect_identical(conditional_frontdoor_result(dag), list(feasible = FALSE, sets = list()))
})

test_that("conditional front door rejects a path bypassing the observed mediators", {
  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    C -> T
    C -> M
    U -> T
    U -> Y
    T -> M
    M -> Y
    T -> Y
  }")
  expect_identical(conditional_frontdoor_result(dag), list(feasible = FALSE, sets = list()))
})

test_that("conditional front door rejects unblocked mediator confounding", {
  dag <- dagitty::dagitty("dag {
    T [exposure]
    Y [outcome]
    U [latent]
    V [latent]
    C -> T
    C -> M
    U -> T
    U -> Y
    V -> M
    V -> Y
    T -> M
    M -> Y
  }")
  expect_identical(conditional_frontdoor_result(dag), list(feasible = FALSE, sets = list()))
  expect_false(id_solve("Y", "T", id_admg(dag))$id)
})

test_that("long ordinary or confounded chains do not search mediator subsets", {
  mediators <- paste0("M", seq_len(24L))
  chain <- paste(c("T", mediators), "->", c(mediators, "Y"), collapse = " ")
  ordinary <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] U [latent] C -> T C -> Y U -> T U -> Y",
    chain, "}"
  ))
  expect_identical(conditional_frontdoor_result(ordinary), list(feasible = FALSE, sets = list()))

  confounded <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] U [latent] C -> T C -> M1 U -> T U -> M1 U -> Y",
    chain, "}"
  ))
  expect_identical(conditional_frontdoor_result(confounded), list(feasible = FALSE, sets = list()))
})

test_that("conditional front door identifies the intervention probability exactly", {
  # Each probability comes from a binary structural model with independent errors.
  joint <- expand.grid(C = 0:1, U = 0:1, T = 0:1, M = 0:1, Y = 0:1)
  p_treatment <- 0.15 + 0.25 * joint$C + 0.30 * joint$U
  p_mediator <- 0.10 + 0.40 * joint$T + 0.20 * joint$C
  p_outcome <- 0.05 + 0.45 * joint$M + 0.25 * joint$U + 0.10 * joint$C
  joint$p <- ifelse(joint$C == 1, 0.4, 0.6) *
    ifelse(joint$U == 1, 0.3, 0.7) *
    ifelse(joint$T == 1, p_treatment, 1 - p_treatment) *
    ifelse(joint$M == 1, p_mediator, 1 - p_mediator) *
    ifelse(joint$Y == 1, p_outcome, 1 - p_outcome)
  expect_equal(sum(joint$p), 1)
  expect_true(all(joint$p > 0))

  identified <- numeric(2L)
  for(treatment in 0:1){
    for(covariate in 0:1){
      p_c <- sum(joint$p[joint$C == covariate])
      for(mediator in 0:1){
        given_t_c <- joint$T == treatment & joint$C == covariate
        p_m <- sum(joint$p[given_t_c & joint$M == mediator]) / sum(joint$p[given_t_c])
        outcome_mean <- 0
        for(reference in 0:1){
          given_reference_c <- joint$T == reference & joint$C == covariate
          p_reference <- sum(joint$p[given_reference_c]) / p_c
          given_m_reference_c <- given_reference_c & joint$M == mediator
          p_y <- sum(joint$p[given_m_reference_c & joint$Y == 1]) /
            sum(joint$p[given_m_reference_c])
          outcome_mean <- outcome_mean + p_y * p_reference
        }
        identified[[treatment + 1L]] <- identified[[treatment + 1L]] +
          p_c * p_m * outcome_mean
      }
    }
  }
  # Intervention removes the treatment mechanism. Averaging its outcome
  # equation over C, U and the mediator gives these exact probabilities.
  intervention <- 0.05 + 0.45 * (0.10 + 0.40 * (0:1) + 0.20 * 0.4) +
    0.25 * 0.3 + 0.10 * 0.4
  expect_equal(identified, intervention, tolerance = 1e-12)
  expect_equal(diff(identified), 0.18, tolerance = 1e-12)
})

test_that("conditional front door handles many irrelevant covariates without testing their combinations", {
  covariates <- paste0("C", seq_len(28L))
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] U [latent]",
    "U -> T U -> Y T -> M M -> Y C1 -> T C1 -> M C2 -> T C2 -> M",
    paste(covariates, "-> Y", collapse = " "), "}"
  ))
  calls <- 0L
  check <- dagitty::isAdjustmentSet
  testthat::local_mocked_bindings(
    isAdjustmentSet = function(...){
      calls <<- calls + 1L
      return(check(...))
    },
    .package = "dagitty"
  )
  expect_identical(conditional_frontdoor_result(dag)$sets, list(
    list(mediators = "M", covariates = c("C1", "C2"))
  ))
  expect_lt(calls, 20L)
})

test_that("conditional front door keeps alternatives among many covariates", {
  covariates <- paste0("C", seq.int(3L, 28L))
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] U [latent]",
    "U -> T U -> Y T -> M M -> Y C1 -> C2 C2 -> T C1 -> M",
    paste(covariates, "-> Y", collapse = " "), "}"
  ))
  expect_identical(conditional_frontdoor_result(dag)$sets, list(
    list(mediators = "M", covariates = "C1"),
    list(mediators = "M", covariates = "C2")
  ))
})

test_that("conditional front door rejects a long chain with unobserved outcome confounding", {
  mediators <- paste0("M", seq_len(24L))
  chain <- paste(c("T", mediators), "->", c(mediators, "Y"), collapse = " ")
  dag <- dagitty::dagitty(paste(
    "dag { T [exposure] Y [outcome] U [latent] V [latent]",
    "U -> T U -> Y V -> M1 V -> Y C -> T C -> M1", chain, "}"
  ))
  expect_identical(conditional_frontdoor_result(dag), list(feasible = FALSE, sets = list()))
})
