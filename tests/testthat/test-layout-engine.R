node_places <- function(dag){
  coordinates <- dagitty::coordinates(dag)
  nodes <- intersect(names(coordinates$x), names(coordinates$y))
  paste(coordinates$x[nodes], coordinates$y[nodes])
}


test_that("levelling a symmetric edge never places a node on an occupied position", {
  symmetric <- dagitty::dagitty(
    'dag { "soil N" [outcome] Rain [exposure] S_2 "soil N" <-> Rain S_2 -> "soil N" }'
  )
  copied <- copy_nodes(symmetric, "Rain", copy_suffix = list(c("a", "b")),
                       print_edges = FALSE)
  expect_identical(anyDuplicated(node_places(copied)), 0L)
  expect_identical(anyDuplicated(node_places(add_coords(copied))), 0L)

  base <- add_coords(build_graph(
    variables = c("X1", "X2", "X3"), treatments = "T", outcomes = "Y",
    type = "saturated"
  ))
  base_coordinates <- dagitty::coordinates(base)
  added <- add_nodes(base, nodes = "W", role = "confounder", type = "full",
                     print_edges = FALSE)
  expect_identical(anyDuplicated(node_places(added)), 0L)
  added_coordinates <- dagitty::coordinates(added)
  expect_equal(added_coordinates$x[names(base_coordinates$x)], base_coordinates$x)
  expect_equal(added_coordinates$y[names(base_coordinates$y)], base_coordinates$y)
  expect_s3_class(add_curves(added), "dagitty")
})


test_that("add_nodes() and join_graphs() keep coordinates by default and regenerate them with renew_coords = TRUE", {
  spaced <- dagitty::dagitty(paste0(
    'dag { T [exposure,pos="0,0"] M [pos="2,1"] Y [outcome,pos="4,0"] ',
    'T -> M [pos="1,0.75"] M -> Y }'
  ))
  before <- dagitty::coordinates(spaced)
  existing <- c("T", "M", "Y")

  kept <- add_nodes(spaced, nodes = "W", role = "confounder", print_edges = FALSE)
  after <- dagitty::coordinates(kept)
  expect_equal(after$x[existing], before$x[existing])
  expect_equal(after$y[existing], before$y[existing])
  expect_true(is.finite(after$x[["W"]]) && is.finite(after$y[["W"]]))
  expect_identical(anyDuplicated(node_places(kept)), 0L)

  moved <- add_nodes(spaced, nodes = "W", role = "confounder", print_edges = FALSE,
                     renew_coords = TRUE)
  expect_equal(dagitty::coordinates(moved), dagitty::coordinates(add_coords(kept)))
  controls <- edge_controls(moved)
  expect_identical(paste(controls$v, controls$e, controls$w), "T -> M")
  expect_true(all(is.finite(c(controls$x, controls$y))))
  expect_error(
    add_nodes(spaced, nodes = "W", role = "confounder", print_edges = FALSE,
              renew_coords = NA),
    "add_nodes(): renew_coords must be TRUE or FALSE.", fixed = TRUE
  )

  donor <- dagitty::dagitty("dag { Y [outcome] Y -> Q }")
  joined <- join_graphs(spaced, donor)
  joined_coordinates <- dagitty::coordinates(joined)
  expect_equal(joined_coordinates$x[existing], before$x[existing])
  expect_equal(joined_coordinates$y[existing], before$y[existing])
  expect_true(is.finite(joined_coordinates$x[["Q"]]))
  expect_identical(anyDuplicated(node_places(joined)), 0L)

  relaid <- join_graphs(spaced, donor, renew_coords = TRUE)
  expect_equal(dagitty::coordinates(relaid), dagitty::coordinates(add_coords(joined)))
  relaid_controls <- edge_controls(relaid)
  expect_identical(paste(relaid_controls$v, relaid_controls$e, relaid_controls$w), "T -> M")
  expect_error(
    join_graphs(spaced, donor, renew_coords = "yes"),
    "join_graphs(): renew_coords must be TRUE or FALSE.", fixed = TRUE
  )
})


test_that("add_nodes() and join_graphs() regenerate every coordinate when a new node cannot be placed apart", {
  crowded <- dagitty::dagitty(
    'dag { A [exposure,pos="0,1"] B [pos="2,0"] C [outcome,pos="2,1"] A -> B }'
  )
  crowded_before <- dagitty::coordinates(crowded)
  expect_warning(
    added <- add_nodes(crowded, nodes = "W", role = "confounder",
                       print_edges = FALSE),
    paste0("add_nodes(): the new nodes could not be placed apart from existing ",
           "nodes; coordinates were regenerated for every node."),
    fixed = TRUE
  )
  expect_equal(dagitty::coordinates(added), dagitty::coordinates(add_coords(added)))
  expect_identical(anyDuplicated(node_places(added)), 0L)
  expect_false(isTRUE(all.equal(
    dagitty::coordinates(added)$y[names(crowded_before$y)], crowded_before$y
  )))

  host <- dagitty::dagitty(
    'dag { A [exposure,pos="4,1"] B [pos="1,2"] C [outcome,pos="1,0"] B -> C }'
  )
  host_before <- dagitty::coordinates(host)
  donor <- dagitty::dagitty("dag { B -> Q }")
  expect_warning(
    joined <- join_graphs(host, donor),
    paste0("join_graphs(): the new nodes could not be placed apart from existing ",
           "nodes; coordinates were regenerated for every node."),
    fixed = TRUE
  )
  expect_equal(dagitty::coordinates(joined), dagitty::coordinates(add_coords(joined)))
  expect_identical(anyDuplicated(node_places(joined)), 0L)
  expect_false(isTRUE(all.equal(
    dagitty::coordinates(joined)$y[names(host_before$y)], host_before$y
  )))
})
