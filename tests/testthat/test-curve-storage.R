testthat::test_that("native controls retain exact zeros, precision and non-layout metadata", {
  text <- paste0('dag named_graph { bb="-2,-2,3,3"\n',
    '"T [->]" [exposure,pos="-1.123456789012345,0",label="Treatment",custom="keep"]\n',
    'Y [outcome,pos="2.123456789012345,0",label="Outcome"]\n',
    '"T [->]" -> Y [label="edge note", custom="alpha"]\n',
    '"T [->]" <-> Y [pos="1.123456789012345,0.555555555555555"]\n',
    '}\n')
  graph <- structure(text, class = "dagitty", labels = c("T [->]" = "long label", Y = "other"),
                     custom = list(version = 7))
  original <- serialize(graph, NULL)
  controls <- data.frame(v = "T [->]", e = "->", w = "Y", x = 0, y = 0.123456789012345)
  result <- causaliflower:::set_edge_controls(graph, controls)
  actual <- causaliflower:::edge_controls(result)
  testthat::expect_equal(nrow(actual), 2L)
  testthat::expect_identical(as.numeric(actual$x[actual$e == "->"]), 0)
  testthat::expect_identical(as.numeric(actual$y[actual$e == "->"]), controls$y)
  testthat::expect_identical(as.numeric(actual$x[actual$e == "<->"]), 1.123456789012345)
  testthat::expect_identical(attributes(result), attributes(graph))
  testthat::expect_identical(dagitty::coordinates(result), dagitty::coordinates(graph))
  for( preserved in c('bb="-2,-2,3,3"', 'label="Treatment"', 'custom="keep"',
                       'label="edge note"', 'custom="alpha"') ){
    testthat::expect_true(grepl(preserved, as.character(result), fixed = TRUE))
  }
  testthat::expect_identical(serialize(graph, NULL), original)
  second <- causaliflower:::set_edge_controls(result, controls)
  testthat::expect_identical(as.character(second), as.character(result))
  file <- tempfile(fileext = ".rds")
  saveRDS(result, file)
  testthat::expect_identical(readRDS(file), result)
  cleared <- causaliflower:::set_edge_controls(result, controls[FALSE, ], replace = controls[, c("v", "e", "w")])
  testthat::expect_identical(causaliflower:::edge_controls(cleared)$e, "<->")
})

testthat::test_that("native updates handle quoted escapes, chains, groups and comments", {
  fixtures <- c(
    'dag { "A \\"quote\\"" [pos="0,0"] B [pos="1,0"] "A \\"quote\\"" -> B [pos="0,0.7",note="a,b"] }',
    'dag { X [pos="0,0"] M [pos="1,0"] Y [pos="2,0"] X -> M -> Y [pos="0.5,1"] }',
    'dag { X [pos="0,0"] A [pos="1,1"] B [pos="1,-1"] X -> { A B } [pos="0.5,2"] }',
    'dag { X [pos="0,0"] Y [pos="1,0"] X -> Y [note="pos=\\"not metadata\\" // text", pos="0.5,0.4", other="a->b"] }'
  )
  for( text in fixtures ){
    graph <- structure(text, class = "dagitty")
    edges <- dagitty::edges(graph)
    controls <- edges[, c("v", "e", "w"), drop = FALSE]
    controls$x <- 0
    controls$y <- seq_len(nrow(controls)) / 7
    result <- causaliflower:::set_edge_controls(graph, controls)
    actual <- causaliflower:::edge_controls(result)
    match <- match(causaliflower:::edge_control_keys(controls$v, controls$e, controls$w),
                   causaliflower:::edge_control_keys(actual$v, actual$e, actual$w))
    testthat::expect_identical(as.numeric(actual$x[match]), as.numeric(controls$x))
    testthat::expect_identical(as.numeric(actual$y[match]), as.numeric(controls$y))
    testthat::expect_setequal(names(result), names(graph))
    testthat::expect_identical(dagitty::coordinates(result), dagitty::coordinates(graph))
    testthat::expect_identical(causaliflower:::set_edge_controls(result, controls), result)
  }
})

testthat::test_that("typed symmetric keys do not collapse directed or unusual names", {
  testthat::expect_identical(causaliflower:::edge_control_keys(character(), character(), character()), character())
  testthat::expect_identical(causaliflower:::edge_control_keys("a", "<->", "b"), causaliflower:::edge_control_keys("b", "<->", "a"))
  testthat::expect_false(causaliflower:::edge_control_keys("a", "->", "b") == causaliflower:::edge_control_keys("b", "->", "a"))
  testthat::expect_false(causaliflower:::edge_control_keys("a\rb", "->", "c") == causaliflower:::edge_control_keys("a", "->", "b\rc"))
  testthat::expect_false(causaliflower:::edge_control_keys("a", "->", "b") == causaliflower:::edge_control_keys("a", "<->", "b"))
})

testthat::test_that("empty graphs and invalid native comment syntax are handled explicitly", {
  graph <- structure('dag { X [pos="1,2"] }', class = "dagitty")
  controls <- causaliflower:::edge_controls(graph)
  testthat::expect_equal(nrow(controls), 0L)
  testthat::expect_identical(causaliflower:::set_edge_controls(graph, controls), graph)
  invalid <- structure('dag { X -> Y // unsupported comment\n }', class = "dagitty")
  testthat::expect_error(causaliflower:::edge_controls(invalid), "cannot read or update native edge control points", fixed = TRUE)
})

testthat::test_that("coordinate writes retain native text metadata and relative bends", {
  graph <- structure('dag { T [exposure,pos="0,0",label="t"] Y [outcome,pos="2,0"] T -> Y [pos="1,2",label="effect"] }',
                     class = "dagitty", custom = 9)
  result <- causaliflower:::set_coords_preserving_curves(graph, list(x = c(T = 3, Y = 3), y = c(T = 4, Y = 8)))
  control <- causaliflower:::edge_controls(result)
  testthat::expect_equal(control$x, -1)
  testthat::expect_equal(control$y, 6)
  testthat::expect_identical(attr(result, "custom"), 9)
  testthat::expect_true(grepl('label="t"', as.character(result), fixed = TRUE))
  testthat::expect_true(grepl('label="effect"', as.character(result), fixed = TRUE))
  testthat::expect_identical(dagitty::exposures(result), "T")
  testthat::expect_identical(dagitty::outcomes(result), "Y")
  testthat::expect_equal(causaliflower:::edge_controls(causaliflower:::set_coords_preserving_curves(result, dagitty::coordinates(result))), control)
  testthat::expect_error(causaliflower:::set_coords_preserving_curves(graph, list(x = c(T = 0, Y = 0), y = c(T = 0, Y = 0))),
                         "cannot preserve a curved edge's shape without finite, distinct endpoint coordinates", fixed = TRUE)
  missing <- structure('dag { T Y T -> Y [pos="0,1"] }', class = "dagitty")
  testthat::expect_error(causaliflower:::set_coords_preserving_curves(missing, list(x = c(T = 0, Y = 1), y = c(T = 0, Y = 0))),
                         "cannot preserve a curved edge's shape without finite, distinct endpoint coordinates", fixed = TRUE)
})

testthat::test_that("tiny and huge layout numbers are written without exponent notation", {
  graph <- structure('dag { T [pos="0,0"] M [pos="1,0"] Y [pos="2,0"] Z [pos="1,-1"] T -> M M -> Y T -> Y Z -> Y }',
                     class = "dagitty")
  result <- add_curves(graph, height = 1e-6, mode = "absolute")
  testthat::expect_identical(edge_keys(result), edge_keys(graph))
  testthat::expect_false(grepl("e-", as.character(result), fixed = TRUE))
  testthat::expect_false(grepl("e+", as.character(result), fixed = TRUE))
  control <- causaliflower:::edge_controls(result)
  testthat::expect_identical(as.numeric(control$y[control$v == "T" & control$w == "Y"]), -2e-06)
  coordinates <- list(x = c(T = 1e-6, M = 1e17, Y = -3.3e-7, Z = 5e18),
                      y = c(T = 4e-5, M = 0, Y = -1e-6, Z = 2.5))
  written <- causaliflower:::write_node_coordinates(graph, coordinates, "curve")
  updated <- dagitty::coordinates(written)
  testthat::expect_identical(as.numeric(updated$x[names(coordinates$x)]), as.numeric(coordinates$x))
  testthat::expect_identical(as.numeric(updated$y[names(coordinates$y)]), as.numeric(coordinates$y))
})
