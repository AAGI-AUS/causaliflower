testthat::test_that("graph editors retain surviving typed controls and exact node coordinates", {
  graph <- structure(paste0('dag { T [exposure,pos="-1.123456789012345,0"] ',
    'M [pos="0,1"] Y [outcome,pos="2.987654321098765,0"] ',
    'T -> M [pos="0,0.123456789012345"] M -> Y [pos="1.2,1.7"] ',
    'T -> Y [pos="1,2"] T <-> Y [pos="-0.0004,-1"] }'),
    class = "dagitty", labels = c(T = "Treatment", M = "Mediator", Y = "Outcome"),
    provenance = list(author = "fixture"))
  original <- serialize(graph, NULL)
  original_coords <- dagitty::coordinates(graph)
  original_controls <- causaliflower:::edge_controls(graph)
  keys <- function(x) causaliflower:::edge_control_keys(x$v, x$e, x$w)
  outputs <- list(
    remove_edges(graph, "M -> Y"),
    keep_edges(c("T -> M", "T -> Y", "Y <-> T"), graph),
    add_edges(graph, "M <-> T"),
    connect_nodes(graph, nodes = c("T", "M", "Y"), type = "ordered", print_edges = FALSE),
    add_nodes(graph, "A", parents = "M", children = "Y", type = "ordered", print_edges = FALSE)
  )
  for( result in outputs ){
    controls <- causaliflower:::edge_controls(result)
    matches <- match(keys(controls), keys(original_controls))
    testthat::expect_false(anyNA(matches))
    testthat::expect_identical(as.numeric(controls$x), as.numeric(original_controls$x[matches]))
    testthat::expect_identical(as.numeric(controls$y), as.numeric(original_controls$y[matches]))
    coords <- dagitty::coordinates(result)
    testthat::expect_identical(coords$x[names(original_coords$x)], original_coords$x)
    testthat::expect_identical(coords$y[names(original_coords$y)], original_coords$y)
    testthat::expect_identical(attr(result, "labels"), attr(graph, "labels"))
    testthat::expect_identical(attr(result, "provenance"), attr(graph, "provenance"))
    testthat::expect_identical(dagitty::exposures(result), "T")
    testthat::expect_identical(dagitty::outcomes(result), "Y")
  }
  added <- outputs[[3L]]
  added_controls <- causaliflower:::edge_controls(added)
  testthat::expect_false(any(
    causaliflower:::edge_control_keys(added_controls$v, added_controls$e,
                                      added_controls$w) ==
      causaliflower:::edge_control_keys("M", "<->", "T")
  ))
  testthat::expect_identical(serialize(graph, NULL), original)
  testthat::expect_warning(
    all_removed <- remove_edges(graph, dagitty::edges(graph)),
    "no remaining edges.*T.*M.*Y"
  )
  testthat::expect_equal(nrow(causaliflower:::edge_controls(all_removed)), 0L)
  testthat::expect_equal(nrow(dagitty::edges(all_removed)), 0L)
  testthat::expect_length(names(all_removed), 0L)
  testthat::expect_length(dagitty::coordinates(all_removed)$x, 0L)
  testthat::expect_length(dagitty::coordinates(all_removed)$y, 0L)
})

testthat::test_that("add_coords preserves full native metadata and moves controls without routing", {
  graph <- structure('dag named { bb="-2,-2,3,3" T [exposure,pos="0,0",label="Treatment",note="a"] M [pos="2,1"] Y [outcome,pos="4,0"] T -> M T -> Y [pos="0,1.25",label="direct",note="curve"] M -> Y }',
                     class = "dagitty", labels = c(T = "t", M = "m", Y = "y"))
  before <- dagitty::coordinates(graph)
  old <- causaliflower:::edge_controls(graph)
  result <- add_coords(graph, x_step = 3, y_step = 2)
  after <- dagitty::coordinates(result)
  actual <- causaliflower:::edge_controls(result)
  chord <- c(after$x["Y"] - after$x["T"], after$y["Y"] - after$y["T"])
  expected <- c(after$x["T"], after$y["T"]) + (1.25 / 4) * c(-chord[2L], chord[1L])
  testthat::expect_equal(as.numeric(actual$x), unname(expected[1L]))
  testthat::expect_equal(as.numeric(actual$y), unname(expected[2L]))
  testthat::expect_equal(nrow(actual), 1L)
  testthat::expect_identical(attr(result, "labels"), attr(graph, "labels"))
  for( piece in c('dag named', 'bb="-2,-2,3,3"', 'label="Treatment"', 'note="a"', 'label="direct"', 'note="curve"') ){
    testthat::expect_true(grepl(piece, as.character(result), fixed = TRUE))
  }
  testthat::expect_identical(dagitty::coordinates(graph), before)
  testthat::expect_identical(causaliflower:::edge_controls(graph), old)
})

testthat::test_that("join_graphs() keeps paths from dag and relative bends from new_dag", {
  first <- structure('dag { T [exposure,pos="-1,0"] Y [outcome,pos="1,0"] T -> Y }',
                      class = "dagitty", labels = c(T = "first", Y = "outcome"))
  donor <- structure('dag { T [exposure,pos="10,10"] Y [outcome,pos="20,10"] Z [pos="30,10"] T -> Y [pos="15,12"] Y -> Z [pos="25,15"] T <-> Y [pos="15,5"] }',
                      class = "dagitty", labels = c(Z = "donor"))
  result <- join_graphs(first, donor)
  controls <- causaliflower:::edge_controls(result)
  testthat::expect_false(any(controls$v == "T" & controls$w == "Y" & controls$e == "->"))
  testthat::expect_equal(nrow(controls), 2L)
  coords <- dagitty::coordinates(result)
  donor_curve <- controls[controls$e == "->", ]
  chord <- c(coords$x["Z"] - coords$x["Y"], coords$y["Z"] - coords$y["Y"])
  expected <- c(coords$x["Y"], coords$y["Y"]) + 0.5 * chord + 0.5 * c(-chord[2L], chord[1L])
  testthat::expect_equal(as.numeric(donor_curve$x), unname(expected[1L]))
  testthat::expect_equal(as.numeric(donor_curve$y), unname(expected[2L]))
  testthat::expect_identical(coords$x[c("T", "Y")], dagitty::coordinates(first)$x)
  testthat::expect_identical(coords$y[c("T", "Y")], dagitty::coordinates(first)$y)
  testthat::expect_identical(attr(result, "labels"), attr(first, "labels"))

  same_nodes <- structure('dag { T [exposure,pos="10,10"] Y [outcome,pos="20,10"] T -> Y [pos="15,12"] T <-> Y [pos="15,5"] }', class = "dagitty")
  result <- join_graphs(first, same_nodes)
  controls <- causaliflower:::edge_controls(result)
  testthat::expect_identical(controls$e, "<->")
  testthat::expect_equal(controls$x, 0)
  testthat::expect_equal(controls$y, -1)

  first_curved <- causaliflower:::set_edge_controls(first, data.frame(v = "T", e = "->", w = "Y", x = 0, y = 0.4))
  controls <- causaliflower:::edge_controls(join_graphs(first_curved, same_nodes))
  testthat::expect_equal(controls$y[controls$e == "->"], 0.4)
})

testthat::test_that("undefined donor transfers error without changing either input", {
  first <- structure('dag { X [pos="0,0"] Y [pos="2,0"] X -> Y }', class = "dagitty")
  donor <- structure('dag { X Y X <-> Y [pos="0,1"] }', class = "dagitty")
  before <- list(serialize(first, NULL), serialize(donor, NULL))
  testthat::expect_error(join_graphs(first, donor), "join_graphs(): cannot preserve a curved edge's shape without finite, distinct endpoint coordinates", fixed = TRUE)
  testthat::expect_identical(serialize(first, NULL), before[[1L]])
  testthat::expect_identical(serialize(donor, NULL), before[[2L]])
})
