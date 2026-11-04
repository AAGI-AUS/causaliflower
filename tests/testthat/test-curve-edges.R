curve_test_dag <- function(){
  dagitty::dagitty('dag {
    T [exposure,pos="0,0"]
    M [pos="2,0"]
    Y [outcome,pos="4,0"]
    E [pos="4,-4"]
    T -> M -> Y
    T -> Y
  }')
}

curve_test_control <- function(plan, from, to){
  Filter(function(control) control$from == from && control$to == to,
         plan$controls)[[1L]]
}

test_that("clear straight arrows are not curved at a shared endpoint", {
  dag <- dagitty::dagitty('dag {
    T [pos="0,0"] Y [pos="4,0"] U [pos="2,-2"]
    T -> Y
    T -> U
  }')
  plan <- causaliflower:::curve_edge_plan(causaliflower:::tidy_ggdagitty(dag))
  expect_true(all(vapply(plan$controls, function(control) !control$curved, logical(1L))))
  expect_true(all(table(plan$paths$edge_id) == 17L))
})

test_that("relative base height follows node and edge spans", {
  plan <- causaliflower:::curve_edge_plan(causaliflower:::tidy_ggdagitty(curve_test_dag()))
  direct <- curve_test_control(plan, "T", "Y")
  expect_equal(c(plan$X, plan$Y, plan$xref, plan$yref), c(4, 4, 4, 4))
  expect_equal(direct$base_height, 0.8)
  expect_equal(direct$height, 0.8)
  expect_true(direct$curved)
})

test_that("relative height responds more strongly to shorter spans", {
  dag <- dagitty::dagitty('dag {
    T [pos="0,0"] M [pos="2,0"] Y [pos="4,0"] H [pos="4,-4"]
    T -> M
    T -> Y
  }')
  plan <- causaliflower:::curve_edge_plan(causaliflower:::tidy_ggdagitty(dag))
  short <- curve_test_control(plan, "T", "M")
  full <- curve_test_control(plan, "T", "Y")
  expect_equal(short$base_height, 0.45)
  expect_gt(short$base_height, 0.4)
  expect_equal(full$base_height, 0.8)
})

test_that("vertical classification uses the edge endpoint scale", {
  dag <- dagitty::dagitty('dag {
    T [pos="0,0"] Y [pos="0.01,0"] F [pos="1000000000000,-1"]
    T -> Y
  }')
  plan <- causaliflower:::curve_edge_plan(causaliflower:::tidy_ggdagitty(dag))
  expect_false(curve_test_control(plan, "T", "Y")$lateral)
})

test_that("routing clears a crossing at the exact path midpoint", {
  dag <- dagitty::dagitty('dag {
    T [pos="0,0"] Y [pos="4,0"]
    U [pos="2,1"] V [pos="2,-1"]
    T -> Y
    U -> V
  }')
  plan <- causaliflower:::curve_edge_plan(causaliflower:::tidy_ggdagitty(dag))
  path <- plan$paths[plan$paths$from == "T" & plan$paths$to == "Y", , drop = FALSE]
  midpoint <- path[which.min(abs(path$index-0.5)), , drop = FALSE]
  expect_equal(midpoint$x, 2)
  expect_gt(abs(midpoint$y), 1)
  expect_false(curve_test_control(plan, "U", "V")$curved)
})

test_that("absolute curves have literal midpoint departures", {
  for( position in c("default", "flipped", "over", "under") ){
    plan <- causaliflower:::curve_edge_plan(causaliflower:::tidy_ggdagitty(curve_test_dag()), 1, "absolute", position)
    direct <- curve_test_control(plan, "T", "Y")
    expect_true(direct$curved)
    expect_equal(direct$height, 1)
    expect_equal(abs(direct$control_y), 2)
    path <- plan$paths[plan$paths$from == "T" & plan$paths$to == "Y", , drop = FALSE]
    expect_equal(abs(path$y[path$index == 0.5]), 1)
  }
})

test_that("sampled paths agree with independent quadratic interpolation", {
  plan <- causaliflower:::curve_edge_plan(causaliflower:::tidy_ggdagitty(curve_test_dag()))
  for( control in plan$controls ){
    path <- plan$paths[plan$paths$from == control$from & plan$paths$to == control$to, , drop = FALSE]
    a <- as.numeric(unlist(plan$nodes[match(control$from, plan$nodes$name), c("x", "y")], use.names = FALSE))
    b <- as.numeric(unlist(plan$nodes[match(control$to, plan$nodes$name), c("x", "y")], use.names = FALSE))
    cp <- c(control$control_x, control$control_y)
    t <- path$index
    first <- (1-t) %o% a+t %o% cp
    second <- (1-t) %o% cp+t %o% b
    expected <- first*(1-t)+second*t
    expect_equal(unname(as.matrix(path[, c("x", "y")])), expected,
                 tolerance = 1e-12, ignore_attr = TRUE)
  }
})

test_that("placement overrides reflect existing default choices", {
  data <- causaliflower:::tidy_ggdagitty(curve_test_dag())
  reference <- causaliflower:::curve_edge_plan(data)
  for( position in c("flipped", "over", "under") ){
    changed <- causaliflower:::curve_edge_plan(data, position = position)
    for( original in reference$controls ){
      control <- curve_test_control(changed, original$from, original$to)
      expect_identical(control$height, original$height)
      expect_identical(control$curved, original$curved)
      expected_sign <- if(position == "flipped") -original$sign else
        if(position == "over") 1 else -1
      expect_equal(control$sign, expected_sign)
    }
  }
})

test_that("routing ignores labels and does not mutate its input", {
  data <- causaliflower:::tidy_ggdagitty(curve_test_dag())
  before <- serialize(data, NULL)
  reference <- causaliflower:::curve_edge_plan(data)
  expect_identical(serialize(data, NULL), before)
  data$label <- paste("A long description for", data$name)
  data$label_width_mm <- 200
  data$label_height_mm <- 100
  attr(data, "seed") <- 98765
  changed <- causaliflower:::curve_edge_plan(data)
  expect_identical(changed$controls, reference$controls)
  expect_identical(changed$anchors, reference$anchors)
  expect_identical(changed$paths[, c("x", "y", "index")],
                   reference$paths[, c("x", "y", "index")])
})

test_that("relative routes follow independent axis rescaling and translation", {
  data <- causaliflower:::tidy_ggdagitty(curve_test_dag())
  reference <- causaliflower:::curve_edge_plan(data)
  changed <- data
  changed$x <- changed$x*10+100
  changed$xend <- changed$xend*10+100
  changed$y <- changed$y*5-50
  changed$yend <- changed$yend*5-50
  actual <- causaliflower:::curve_edge_plan(changed)
  expect_equal((actual$paths$x-100)/10, reference$paths$x, tolerance = 1e-12)
  expect_equal((actual$paths$y+50)/5, reference$paths$y, tolerance = 1e-12)
  expect_identical(actual$routing_metrics$changes, reference$routing_metrics$changes)
})

test_that("flat and vertical layouts retain fixed reference anchors", {
  flat <- dagitty::dagitty('dag {
    T [pos="0,0"] M [pos="2,0"] Y [pos="4,0"]
    T -> M -> Y
    T -> Y
  }')
  first <- causaliflower:::curve_edge_plan(causaliflower:::tidy_ggdagitty(flat), height = 0.1)
  second <- causaliflower:::curve_edge_plan(causaliflower:::tidy_ggdagitty(flat), height = 0.4)
  expect_equal(c(first$xref, first$yref), c(4, 1))
  expect_identical(first$anchors, second$anchors)
  expect_equal(first$anchors$y, c(-0.5, 0.5))

  vertical <- dagitty::dagitty('dag {
    T [pos="0,0"] M [pos="0,-2"] Y [pos="0,-4"]
    T -> M -> Y
    T -> Y
  }')
  plan <- causaliflower:::curve_edge_plan(causaliflower:::tidy_ggdagitty(vertical), height = 0.35,
                            mode = "absolute", position = "over")
  direct <- curve_test_control(plan, "T", "Y")
  expect_equal(c(plan$xref, plan$yref), c(4, 4))
  expect_true(direct$lateral)
  expect_equal(direct$control_x, -0.7)
  expect_equal(direct$height, 0.35)
})

test_that("add_curves writes native controls while preserving input and roles", {
  dag <- curve_test_dag()
  before <- as.character(dag)
  result <- add_curves(dag, height = 1, mode = "absolute", position = "over")
  expect_identical(as.character(dag), before)
  expect_equal(causaliflower:::dag_edges(result), causaliflower:::dag_edges(dag))
  expect_identical(dagitty::coordinates(result), dagitty::coordinates(dag))
  expect_identical(dagitty::exposures(result), dagitty::exposures(dag))
  expect_identical(dagitty::outcomes(result), dagitty::outcomes(dag))
  controls <- causaliflower:::edge_controls(result)
  direct <- controls[controls$v == "T" & controls$e == "->" & controls$w == "Y", ]
  expect_equal(direct$x, 2)
  expect_equal(direct$y, -2)
})

test_that("missing coordinates regenerate the complete layout", {
  dag <- dagitty::dagitty('dag {
    T [pos="100,50"] M Y [pos="200,60"]
    T -> M -> Y
    T -> Y [pos="150,40"]
  }')
  before <- as.character(dag)
  expected <- dagitty::coordinates(add_coords(causaliflower:::set_edge_controls(
    dag, data.frame(v = character(), e = character(), w = character(),
                    x = numeric(), y = numeric()),
    replace = causaliflower:::dag_edges(dag)
  )))
  result <- add_curves(dag)
  expect_identical(dagitty::coordinates(result), expected)
  expect_false(unname(dagitty::coordinates(result)$x["T"]) == 100)
  expect_identical(as.character(dag), before)
})

test_that("other typed-edge controls survive directed rerouting and straightening", {
  dag <- dagitty::dagitty('dag {
    T [pos="0,0"] M [pos="2,0"] Y [pos="4,0"]
    T -> M -> Y
    T -> Y [pos="2,-2"]
    T <-> Y [pos="1,2"]
  }')
  before <- causaliflower:::edge_controls(dag)
  before <- before[before$e == "<->", , drop = FALSE]
  result <- add_curves(dag)
  after <- causaliflower:::edge_controls(result)
  expect_equal(after[after$e == "<->", , drop = FALSE], before,
               ignore_attr = TRUE)
  straight <- add_curves(result, height = 0)
  controls <- causaliflower:::edge_controls(straight)
  expect_false(any(controls$e == "->"))
  expect_equal(controls, before, ignore_attr = TRUE)
})

test_that("zero height does not require a direction for coincident endpoints", {
  dag <- dagitty::dagitty('dag {
    T [pos="0,0"] M [pos="0,0"] Y [pos="4,0"]
    T -> M [pos="1,-1"]
    M -> Y
    T -> Y
  }')
  expect_error(add_curves(dag), "endpoints have identical coordinates", fixed = TRUE)
  result <- add_curves(dag, height = 0)
  expect_equal(nrow(causaliflower:::edge_controls(result)), 0L)
  expect_identical(dagitty::coordinates(result), dagitty::coordinates(dag))
})

test_that("add_curves input validation rejects invalid and unrepresentable settings", {
  dag <- curve_test_dag()
  expect_error(add_curves("dag { T -> Y }"), "dag must be a dagitty object", fixed = TRUE)
  for( height in list(NULL, numeric(), c(0.1, 0.2), NA_real_, NaN, Inf, -Inf,
                      -1, TRUE, "0.2", 1+1i) ){
    expect_error(add_curves(dag, height = height),
                  "height must be a single non-negative finite number", fixed = TRUE)
  }
  for( mode in list(NULL, NA_character_, c("relative", "absolute"), "other", TRUE) ){
    expect_error(add_curves(dag, mode = mode),
                  "mode must be 'relative' or 'absolute'", fixed = TRUE)
  }
  for( position in list(NULL, NA_character_, c("over", "under"), "auto", 1) ){
    expect_error(add_curves(dag, position = position),
                  "position must be 'default', 'flipped', 'over' or 'under'", fixed = TRUE)
  }
  expect_error(add_curves(dag, height = .Machine$double.xmax),
                "cannot be represented in this coordinate layout", fixed = TRUE)
  expect_error(add_curves(dag, height = c(named = 1L),
                            mode = c(named = "absolute"), position = c(named = "over")), NA)
})

test_that("empty graphs and graphs without edges do not enter curve planning", {
  empty <- dagitty::dagitty("dag {}")
  expect_identical(add_curves(empty), empty)
  isolated <- dagitty::dagitty("dag { A B }")
  result <- add_curves(isolated)
  expect_identical(names(result), names(isolated))
  expect_identical(dagitty::coordinates(result), dagitty::coordinates(add_coords(isolated)))
  expect_equal(nrow(causaliflower:::edge_controls(result)), 0L)
})


test_that("add_curves() accepts pdag and mag graphs and keeps their type", {
  pdag <- dagitty::dagitty("pdag { T [exposure] Y [outcome] Z T -> Y T -- Z Z -> Y }")
  curved <- add_curves(add_coords(pdag), height = 0.3)
  expect_s3_class(curved, "dagitty")
  expect_identical(dag_graph_type(curved), "pdag")
  expect_setequal(names(curved), c("T", "Y", "Z"))

  mag <- dagitty::dagitty("mag { T [exposure] Y [outcome] Z T -> Y Z <-> T Z -> Y }")
  expect_identical(dag_graph_type(add_curves(add_coords(mag), height = 0.3)), "mag")
})
