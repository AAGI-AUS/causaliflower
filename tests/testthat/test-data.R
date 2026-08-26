test_that("ESCDAGs matches its generation script and stays character-typed", {
  expect_s3_class(ESCDAGs, "data.frame")
  expect_identical(dim(ESCDAGs), c(3L, 5L))
  expect_identical(names(ESCDAGs),
                   c("name", "question", "description", "source", "required"))
  expect_identical(ESCDAGs$name,
                   c("Temporality", "Face-validity", "Recourse to theory"))
  expect_identical(ESCDAGs$required, c("yes", "yes", "no"))
  expect_false(any(vapply(ESCDAGs, is.factor, logical(1))))
  expect_match(ESCDAGs$description[1], "Glass et al., 2013", fixed = TRUE)
  expect_match(ESCDAGs$description[1],
               "Ferguson et al., 2020, p. 326", fixed = TRUE)
  expect_match(ESCDAGs$source[1], "10.1093/ije/dyz150", fixed = TRUE)
})

test_that("the shipped ESCDAGs drives the default assessment log", {
  dag <- string_to_dag("A -> B")
  log <- data.frame(v = "A", e = "->", w = "B",
                    Temporality = "y", `Face-validity` = "y",
                    `Recourse to theory` = "y",
                    check.names = FALSE, stringsAsFactors = FALSE)

  invisible(capture.output(
    replayed <- suppressMessages(assess_edges(
      dag,
      assess_causal_criteria = TRUE,
      save_answers = log
    ))
  ))
  expect_identical(paste(replayed$v, replayed$e, replayed$w), "A -> B")
})
