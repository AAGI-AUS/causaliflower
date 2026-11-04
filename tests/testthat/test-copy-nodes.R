copy_edge_keys <- function(dag){
  edges <- as.data.frame(dagitty::edges(dag))
  if( nrow(edges) == 0L ) return(character())
  sort(paste(edges$v, edges$e, edges$w))
}


copy_node_places <- function(dag){
  coordinates <- dagitty::coordinates(dag)
  nodes <- names(coordinates$x)
  paste(coordinates$x[nodes], coordinates$y[nodes])
}


test_that("copy_nodes() applies default and generated naming rules", {
  dag <- dagitty::dagitty("dag { P -> M; M -> Y }")

  default <- copy_nodes(dag, "M", print_edges = FALSE)
  expect_setequal(names(default), c("P", "M", "M_copy", "Y"))
  expect_setequal(
    copy_edge_keys(default),
    c("P -> M", "M -> Y", "P -> M_copy", "M_copy -> Y",
      "M -> M_copy")
  )

  ordered <- copy_nodes(dag, "M", suffix = "t", repeats = 2,
                        print_edges = FALSE)
  expect_setequal(names(ordered), c("P", "M_t1", "M_t2", "M_t3", "Y"))
  expect_true(all(c("M_t1 -> M_t2", "M_t2 -> M_t3") %in%
                    copy_edge_keys(ordered)))
  expect_false("M_t1 -> M_t3" %in% copy_edge_keys(ordered))

  parallel <- copy_nodes(dag, "M", suffix = list("t"), repeats = 2,
                         print_edges = FALSE)
  expect_true(all(c("M_t1 -> M_t2", "M_t1 -> M_t3") %in%
                    copy_edge_keys(parallel)))
  expect_false("M_t2 -> M_t3" %in% copy_edge_keys(parallel))
  expect_false("M_t3 -> M_t2" %in% copy_edge_keys(parallel))
})


test_that("copy suffix positions distinguish ordered and parallel copies", {
  dag <- dagitty::dagitty("dag { T -> M; M -> Y }")

  ordered <- copy_nodes(dag, "M", copy_suffix = c("a", "b"),
                        print_edges = FALSE)
  expect_true(all(c("M -> M_a", "M_a -> M_b") %in%
                    copy_edge_keys(ordered)))
  expect_false("M -> M_b" %in% copy_edge_keys(ordered))

  parallel <- copy_nodes(
    dag, "M", copy_suffix = list(c("survey", "sensor")),
    print_edges = FALSE
  )
  expect_true(all(c("M -> M_survey", "M -> M_sensor") %in%
                    copy_edge_keys(parallel)))
  expect_false(any(c("M_survey -> M_sensor", "M_sensor -> M_survey") %in%
                     copy_edge_keys(parallel)))

  tiered <- copy_nodes(
    dag, "M", copy_suffix = list(c("a", "b"), "c"),
    print_edges = FALSE
  )
  expect_true(all(c("M -> M_a", "M -> M_b", "M_a -> M_c",
                    "M_b -> M_c") %in% copy_edge_keys(tiered)))
  expect_false(any(c("M_a -> M_b", "M_b -> M_a", "M -> M_c") %in%
                     copy_edge_keys(tiered)))

  expect_warning(
    copy_nodes(dag, "M", copy_suffix = c("a", "b"), repeats = 2,
               print_edges = FALSE),
    "repeats is redundant"
  )
  expect_error(
    copy_nodes(dag, "M", copy_suffix = c("a", "b"), repeats = 3,
               print_edges = FALSE),
    "copy_suffix supplies 2 suffixes"
  )
})


test_that("descendant and unrelated anchors parent the first copy", {
  descendant <- dagitty::dagitty(
    "dag { P -> Z; Z -> Q; Q -> A; Z -> Y }"
  )
  copied <- copy_nodes(
    descendant, "Z", suffix = "pre", copy_suffix = c("post1", "post2"),
    anchor = "A", print_edges = FALSE
  )
  edges <- copy_edge_keys(copied)

  expect_true(all(c(
    "Z_pre -> Z_post1", "Z_post1 -> Z_post2",
    "A -> Z_post1",
    "P -> Z_post1", "P -> Z_post2",
    "Z_post1 -> Y", "Z_post2 -> Y"
  ) %in% edges))
  expect_false(any(c(
    "Z_post1 -> Q", "Z_post2 -> Q", "A -> Z_post2",
    "Z_post1 -> A", "Z_post2 -> A"
  ) %in% edges))
  expect_true(all(c("Z_pre -> Q", "Q -> A") %in% edges))

  unrelated <- dagitty::dagitty("dag { P -> Z; Z -> Y; A }")
  placed <- copy_nodes(unrelated, "Z", copy_suffix = "later", anchor = "A",
                       print_edges = FALSE)
  expect_true(all(c("A -> Z_later", "Z -> Z_later", "P -> Z_later",
                    "Z_later -> Y") %in% copy_edge_keys(placed)))
})


test_that("ancestor anchors child the last copy and reverse the boundary chain", {
  dag <- dagitty::dagitty("dag { A -> Q; Q -> M; P -> M; M -> Y }")
  copied <- copy_nodes(
    dag, "M", suffix = "post", copy_suffix = c("pre1", "pre2"),
    anchor = "A", print_edges = FALSE
  )
  edges <- copy_edge_keys(copied)

  expect_true(all(c(
    "M_pre1 -> M_pre2", "M_pre2 -> M_post",
    "M_pre2 -> A",
    "P -> M_pre1", "P -> M_pre2",
    "M_pre1 -> Y", "M_pre2 -> Y"
  ) %in% edges))
  expect_false(any(c("Q -> M_pre1", "Q -> M_pre2", "M_pre1 -> A",
                     "A -> M_pre1", "A -> M_pre2") %in% edges))
  expect_true(all(c("A -> Q", "Q -> M_post") %in% edges))
})


test_that("the anchor joins one end of the copies for each anchor kind", {
  descendant <- copy_nodes(
    dagitty::dagitty("dag { P -> Z; Z -> Q; Q -> A; Z -> Y }"), "Z",
    copy_suffix = c("c1", "c2", "c3"), anchor = "A", print_edges = FALSE
  )
  expect_setequal(
    copy_edge_keys(descendant),
    c("A -> Z_c1", "P -> Z", "P -> Z_c1", "P -> Z_c2", "P -> Z_c3",
      "Q -> A", "Z -> Q", "Z -> Y", "Z -> Z_c1", "Z_c1 -> Y",
      "Z_c1 -> Z_c2", "Z_c2 -> Y", "Z_c2 -> Z_c3", "Z_c3 -> Y")
  )

  unrelated <- copy_nodes(
    dagitty::dagitty("dag { P -> Z; Z -> Y; A }"), "Z",
    copy_suffix = c("c1", "c2", "c3"), anchor = "A", print_edges = FALSE
  )
  expect_setequal(
    copy_edge_keys(unrelated),
    c("A -> Z_c1", "P -> Z", "P -> Z_c1", "P -> Z_c2", "P -> Z_c3",
      "Z -> Y", "Z -> Z_c1", "Z_c1 -> Y", "Z_c1 -> Z_c2", "Z_c2 -> Y",
      "Z_c2 -> Z_c3", "Z_c3 -> Y")
  )

  ancestor <- copy_nodes(
    dagitty::dagitty("dag { A -> Q; Q -> M; P -> M; M -> Y }"), "M",
    copy_suffix = c("c1", "c2", "c3"), anchor = "A", print_edges = FALSE
  )
  expect_setequal(
    copy_edge_keys(ancestor),
    c("A -> Q", "M -> Y", "M_c1 -> M_c2", "M_c1 -> Y", "M_c2 -> M_c3",
      "M_c2 -> Y", "M_c3 -> A", "M_c3 -> M", "M_c3 -> Y", "P -> M",
      "P -> M_c1", "P -> M_c2", "P -> M_c3", "Q -> M")
  )

  tiered <- copy_nodes(
    dagitty::dagitty("dag { P -> Z; Z -> Q; Q -> A; Z -> Y }"), "Z",
    copy_suffix = list(c("s1", "s2"), "t"), anchor = "A", print_edges = FALSE
  )
  expect_true(all(c("A -> Z_s1", "A -> Z_s2") %in% copy_edge_keys(tiered)))
  expect_false("A -> Z_t" %in% copy_edge_keys(tiered))
})


test_that("chain places the copies in the copied node's own path", {
  worked <- copy_nodes(dagitty::dagitty("dag { T -> M -> Y }"), "M",
                       type = "chain", print_edges = FALSE)
  expect_setequal(
    copy_edge_keys(worked),
    c("T -> M", "M -> M_copy", "M_copy -> Y")
  )

  successive <- copy_nodes(
    dagitty::dagitty("dag { P -> M; M -> Y1; M -> Y2 }"), "M",
    copy_suffix = c("a", "b", "c"), type = "chain", print_edges = FALSE
  )
  expect_setequal(
    copy_edge_keys(successive),
    c("P -> M", "M -> M_a", "M_a -> M_b", "M_b -> M_c", "M_c -> Y1",
      "M_c -> Y2")
  )
})


test_that("chain places the copies after a descendant or unrelated anchor", {
  worked <- copy_nodes(
    dagitty::dagitty("dag { O -> A; A -> C1; A -> C2 }"), "O",
    anchor = "A", type = "chain", print_edges = FALSE
  )
  expect_setequal(
    copy_edge_keys(worked),
    c("O -> A", "A -> O_copy", "O_copy -> C1", "O_copy -> C2", "O -> O_copy")
  )

  tiered <- copy_nodes(
    dagitty::dagitty("dag { O -> A; A -> C1; A -> C2 }"), "O",
    copy_suffix = list(c("s1", "s2"), "t"), anchor = "A", type = "chain",
    print_edges = FALSE
  )
  expect_setequal(
    copy_edge_keys(tiered),
    c("O -> A", "A -> O_s1", "A -> O_s2", "O -> O_s1", "O -> O_s2",
      "O_s1 -> O_t", "O_s2 -> O_t", "O_t -> C1", "O_t -> C2")
  )

  unrelated <- copy_nodes(
    dagitty::dagitty("dag { P -> Z; Z -> Y; A -> B }"), "Z",
    anchor = "A", type = "chain", print_edges = FALSE
  )
  expect_setequal(
    copy_edge_keys(unrelated),
    c("P -> Z", "Z -> Y", "A -> Z_copy", "Z -> Z_copy", "Z_copy -> B")
  )
})


test_that("chain places the copies before an ancestor anchor", {
  worked <- copy_nodes(
    dagitty::dagitty("dag { P1 -> A; A -> X; X -> O }"), "O",
    anchor = "A", type = "chain", print_edges = FALSE
  )
  expect_setequal(
    copy_edge_keys(worked),
    c("P1 -> O_copy", "O_copy -> A", "O_copy -> O", "A -> X", "X -> O")
  )

  tiered <- copy_nodes(
    dagitty::dagitty("dag { P1 -> A; P2 -> A; A -> X; X -> O }"), "O",
    copy_suffix = list(c("p1", "p2")), anchor = "A", type = "chain",
    print_edges = FALSE
  )
  expect_setequal(
    copy_edge_keys(tiered),
    c("P1 -> O_p1", "P1 -> O_p2", "P2 -> O_p1", "P2 -> O_p2",
      "O_p1 -> A", "O_p2 -> A", "O_p1 -> O", "O_p2 -> O", "A -> X",
      "X -> O")
  )
})


test_that("chain keeps node statuses and reports the arrows it takes over", {
  dag <- build_graph("X", "T", "Y", mediators = "M",
                     latents = "U", type = "ordered")
  copied <- copy_nodes(dag, "M", type = "chain", print_edges = FALSE)

  expect_setequal(
    copy_edge_keys(copied),
    c("M -> M_copy", "M_copy -> Y", "U -> T", "U -> Y", "T -> M", "T -> Y",
      "X -> T", "X -> Y")
  )
  expect_identical(treatments(copied), "T")
  expect_identical(get_outcome_names(copied), "Y")
  expect_identical(unobserved(copied), "U")
  expect_true(dagitty::isAcyclic(copied))

  messages <- character()
  capture.output(
    withCallingHandlers(
      copy_nodes(dagitty::dagitty("dag { O -> A; A -> C1; A -> C2 }"), "O",
                 anchor = "A", type = "chain", print_edges = TRUE),
      message = function(condition){
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
  )
  expect_match(paste(messages, collapse = "\n"), "A -> C1' => 'O_copy -> C1")
  expect_match(paste(messages, collapse = "\n"), "A -> C2' => 'O_copy -> C2")
})


test_that("copy_nodes() accepts only chain as a type", {
  dag <- dagitty::dagitty("dag { T -> M -> Y }")
  expect_error(copy_nodes(dag, "M", type = "first", print_edges = FALSE),
               "type must be NULL")
  expect_error(copy_nodes(dag, "M", type = "saturated", print_edges = FALSE),
               "type must be NULL")
  expect_error(copy_nodes(dag, "M", type = c("chain", "chain"),
                          print_edges = FALSE),
               "type must be NULL")
  expect_error(copy_nodes(dag, "M", type = NA, print_edges = FALSE),
               "type must be NULL")
})


test_that("copy_nodes() retains declarations, coordinates, metadata, and curves", {
  dag <- dagitty::dagitty(
    paste0(
      'dag named {\n bb="-2,-2,5,3"\n',
      'T [exposure,adjusted,pos="0,0",label="Treatment",note="keep node"]\n',
      'M [latent,pos="2,1"]\n',
      'Y [outcome,selected,pos="4,0"]\n',
      'T -> M [pos="1,0.75",label="curved",note="keep edge"]\n',
      'M -> Y\n}'
    )
  )
  attr(dag, "labels") <- c(T = "t", M = "m", Y = "y")
  attr(dag, "provenance") <- list(author = "fixture")
  before <- serialize(dag, NULL)
  before_coordinates <- dagitty::coordinates(dag)
  before_controls <- edge_controls(dag)

  copied <- copy_nodes(
    dag, c("T", "M", "Y"), suffix = "original", copy_suffix = "copy",
    renew_coords = FALSE, print_edges = FALSE
  )

  expect_identical(serialize(dag, NULL), before)
  expect_setequal(treatments(copied), c("T_original", "T_copy"))
  expect_setequal(outcomes(copied), c("Y_original", "Y_copy"))
  expect_setequal(unobserved(copied), c("M_original", "M_copy"))
  expect_identical(dagitty::adjustedNodes(copied), "T_original")
  expect_identical(selected_nodes(copied), "Y_original")

  coordinates <- dagitty::coordinates(copied)
  rename <- c(T = "T_original", M = "M_original", Y = "Y_original")
  expect_equal(
    unname(coordinates$x[rename]), unname(before_coordinates$x[names(rename)])
  )
  expect_equal(
    unname(coordinates$y[rename]), unname(before_coordinates$y[names(rename)])
  )
  expect_false(anyNA(coordinates$x[c("T_copy", "M_copy", "Y_copy")]))
  expect_false(anyNA(coordinates$y[c("T_copy", "M_copy", "Y_copy")]))

  controls <- edge_controls(copied)
  expect_identical(controls$v, "T_original")
  expect_identical(controls$w, "M_original")
  expect_identical(controls$e, before_controls$e)
  expect_identical(controls$x, before_controls$x)
  expect_identical(controls$y, before_controls$y)
  expect_identical(attr(copied, "labels"),
                   c(T_original = "t", M_original = "m", Y_original = "y"))
  expect_identical(attr(copied, "provenance"), attr(dag, "provenance"))

  text <- as.character(copied)
  for( piece in c('dag named', 'bb="-2,-2,5,3"', 'label="Treatment"',
                  'note="keep node"', 'label="curved"',
                  'note="keep edge"') ){
    expect_true(grepl(piece, text, fixed = TRUE))
  }
  expect_true(all(c("T_copy -> M_copy", "M_copy -> Y_copy") %in%
                    copy_edge_keys(copied)))
})


test_that("copy_nodes() places copies in the source's role order", {
  dag <- build_graph(
    "X", "T", "Y",
    mediators = list(c(incoming("M1"), "M2"), outgoing("M3")),
    type = "saturated"
  )

  copied <- copy_nodes(
    dag, "M1", suffix = "original",
    copy_suffix = list(c("survey", "sensor"), "harvest"),
    print_edges = FALSE
  )

  expect_identical(
    test_role_positions(copied, "mediator"),
    list("M1_original", c("M1_survey", "M1_sensor"), "M1_harvest", "M3")
  )
  expect_true(all(c("T -> M1_original", "X -> M1_original",
                    "T -> M1_harvest", "X -> M1_harvest",
                    "M1_original -> M3", "M1_harvest -> M3") %in%
                    copy_edge_keys(copied)))
  expect_false("M1" %in% names(copied))

  before <- copy_nodes(
    dag, "M1", suffix = "original", copy_suffix = c("early", "late"),
    anchor = "T", print_edges = FALSE
  )
  expect_true(all(c("M1_early -> M1_late", "M1_late -> M1_original",
                    "M1_late -> T") %in% copy_edge_keys(before)))
  expect_false("M1_early -> T" %in% copy_edge_keys(before))
})


test_that("copy_nodes() copies each latent's own target mapping", {
  dag <- build_graph(
    "X", "T", "Y", mediators = "M",
    latents = list("V", list(outgoing("U"), c("T", "M"))),
    type = "saturated"
  )
  copied <- copy_nodes(
    dag, c("U", "V", "M"), suffix = "original",
    copy_suffix = list(c("survey", "sensor")), print_edges = FALSE
  )
  latent <- role_mapping_sources(copied, "latent")

  expect_setequal(latent, c("U_original", "U_survey", "U_sensor"))
  expect_false(any(c("V_original", "V_survey", "V_sensor") %in% latent))
  expect_true(all(c("U_original -> T", "U_original -> M_original",
                    "U_survey -> T", "U_sensor -> T",
                    "U_survey -> M_survey", "U_sensor -> M_sensor") %in%
                    copy_edge_keys(copied)))
  expect_false(any(c("U", "V", "M") %in% names(copied)))
})


test_that("copy_nodes() supports non-syntactic names", {
  dag <- dagitty::dagitty(
    'dag { "soil N" -> "grain yield"; "grain yield" -> outcome }'
  )
  copied <- copy_nodes(
    dag, "grain yield", suffix = "before", copy_suffix = "after",
    print_edges = FALSE
  )
  expect_setequal(
    names(copied),
    c("soil N", "grain yield_before", "grain yield_after", "outcome")
  )
  expect_true(all(c("soil N -> grain yield_before",
                    "soil N -> grain yield_after",
                    "grain yield_before -> outcome",
                    "grain yield_after -> outcome") %in%
                    copy_edge_keys(copied)))
})


test_that("copy_nodes() rejects invalid names, name collisions, and cycles", {
  collision <- dagitty::dagitty("dag { M; M_copy; M -> Y }")
  before_collision <- serialize(collision, NULL)
  expect_error(
    copy_nodes(collision, "M", print_edges = FALSE),
    "generated node names collide"
  )
  expect_identical(serialize(collision, NULL), before_collision)

  cyclic <- dagitty::dagitty("dag { S -> A; A -> T; S -> T }")
  before_cyclic <- serialize(cyclic, NULL)
  expect_error(
    copy_nodes(cyclic, c("S", "T"), copy_suffix = "copy", anchor = "A",
               print_edges = FALSE),
    "would create a directed cycle"
  )
  expect_identical(serialize(cyclic, NULL), before_cyclic)

  dag <- dagitty::dagitty("dag { T -> M -> Y }")
  expect_error(copy_nodes(dag, character(), print_edges = FALSE),
               "nodes must contain")
  expect_error(copy_nodes(dag, "absent", print_edges = FALSE),
               "not in the graph")
  expect_error(copy_nodes(dag, "M", anchor = "M", print_edges = FALSE),
               "anchor cannot")
  expect_error(copy_nodes(dag, "M", repeats = 0, print_edges = FALSE),
               "positive whole number")
  expect_error(copy_nodes(dag, "M", copy_suffix = list(list("a")),
                          print_edges = FALSE),
               "nested more than one level")
  expect_error(copy_nodes(dagitty::dagitty("pdag { T -- Y }"), "T",
                          print_edges = FALSE),
               "dag graph type")
})


test_that("copy_nodes() prints the edges an anchor reroutes", {
  dag <- dagitty::dagitty("dag { Z -> Q; Q -> A; Z -> Y }")
  messages <- character()
  printed <- capture.output(
    withCallingHandlers(
      copied <- copy_nodes(dag, "Z", copy_suffix = "post", anchor = "A",
                            print_edges = TRUE),
      message = function(condition){
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
  )
  expect_true(length(printed) > 0L)
  expect_match(paste(messages, collapse = "\n"), "rerouted edges")
  expect_match(paste(messages, collapse = "\n"), "Z_post -> Q")
  expect_match(paste(messages, collapse = "\n"), "A -> Z_post")
  expect_true("A -> Z_post" %in% copy_edge_keys(copied))
  expect_false("Z_post -> Q" %in% copy_edge_keys(copied))

  after_messages <- character()
  capture.output(
    withCallingHandlers(
      after <- copy_nodes(dagitty::dagitty("dag { Z -> Q; Q -> A; Z -> Y }"),
                          "Z", copy_suffix = c("post1", "post2"),
                          anchor = "A", print_edges = TRUE),
      message = function(condition){
      after_messages <<- c(after_messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
  )
  after_text <- paste(after_messages, collapse = "\n")
  expect_match(after_text, "Z_post1 -> Q' => 'A -> Z_post1")
  expect_false(grepl("A -> Z_post2", after_text, fixed = TRUE))
  expect_false(grepl("Z_post2 -> Q", after_text, fixed = TRUE))
  expect_true(all(c("A -> Z_post1", "Z -> Z_post1", "Z_post1 -> Z_post2",
                    "Z_post1 -> Y", "Z_post2 -> Y") %in%
                    copy_edge_keys(after)))
  expect_false(any(c("A -> Z_post2", "Z_post1 -> Q", "Z_post2 -> Q") %in%
                     copy_edge_keys(after)))
  after_lines <- strsplit(after_text, "\n", fixed = TRUE)[[1L]]
  after_targets <- regmatches(after_lines, regexpr("=> '[^']+'", after_lines))
  expect_true(all(substr(after_targets, 5L, nchar(after_targets) - 1L) %in%
                    copy_edge_keys(after)))

  before_messages <- character()
  capture.output(
    withCallingHandlers(
      before <- copy_nodes(dagitty::dagitty("dag { A -> P; P -> Z; Z -> Y }"),
                           "Z", copy_suffix = c("pre1", "pre2"),
                           anchor = "A", print_edges = TRUE),
      message = function(condition){
      before_messages <<- c(before_messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
  )
  before_text <- paste(before_messages, collapse = "\n")
  expect_match(before_text, "P -> Z_pre2' => 'Z_pre2 -> A")
  expect_false(grepl("Z_pre1 -> A", before_text, fixed = TRUE))
  expect_false(grepl("P -> Z_pre1", before_text, fixed = TRUE))
  expect_true(all(c("Z_pre2 -> A", "Z_pre1 -> Z_pre2", "Z_pre2 -> Z",
                    "Z_pre1 -> Y", "Z_pre2 -> Y") %in%
                    copy_edge_keys(before)))
  expect_false(any(c("Z_pre1 -> A", "P -> Z_pre1", "P -> Z_pre2") %in%
                     copy_edge_keys(before)))
  before_lines <- strsplit(before_text, "\n", fixed = TRUE)[[1L]]
  before_targets <- regmatches(before_lines,
                               regexpr("=> '[^']+'", before_lines))
  expect_true(all(substr(before_targets, 5L, nchar(before_targets) - 1L) %in%
                    copy_edge_keys(before)))
})


test_that("chain copies keep native attributes when an arrow is removed", {
  dag <- dagitty::dagitty(
    paste0(
      'dag named {\n bb="-2,-2,5,3"\n',
      'T [exposure,adjusted,pos="0,0",label="Treatment",note="keep node"]\n',
      'M [latent,pos="2,1"]\n',
      'Y [outcome,selected,pos="4,0"]\n',
      'T -> M [pos="1,0.75",label="curved",note="keep edge"]\n',
      'M -> Y\n}'
    )
  )

  copied <- copy_nodes(dag, "M", type = "chain", print_edges = FALSE)
  text <- as.character(copied)

  expect_match(text, "dag named", fixed = TRUE)
  expect_match(text, 'bb="-2,-2,5,3"', fixed = TRUE)
  expect_match(text, 'label="Treatment"', fixed = TRUE)
  expect_match(text, 'note="keep node"', fixed = TRUE)
  expect_match(text, 'label="curved"', fixed = TRUE)
  expect_match(text, 'note="keep edge"', fixed = TRUE)
  expect_true(all(c("T -> M", "M -> M_copy", "M_copy -> Y") %in%
                    copy_edge_keys(copied)))
  expect_false("M -> Y" %in% copy_edge_keys(copied))
  expect_identical(causaliflower:::selected_nodes(copied), "Y")
})


test_that("copy_nodes() regenerates coordinates by default and keeps them with renew_coords = FALSE", {
  crowded <- dagitty::dagitty(
    paste0('dag { A [pos="0,1"] B [pos="2,0"] C [pos="4,1"] ',
           'D [pos="6,0"] A -> D B -> D C -> D }')
  )

  renewed <- copy_nodes(crowded, "C", suffix = "v", anchor = "A",
                        print_edges = FALSE)
  renewed_coordinates <- dagitty::coordinates(renewed)
  expect_true(all(is.finite(renewed_coordinates$x)))
  expect_true(all(is.finite(renewed_coordinates$y)))
  expect_identical(anyDuplicated(copy_node_places(renewed)), 0L)
  expect_s3_class(add_curves(renewed), "dagitty")

  grDevices::pdf(NULL)
  device <- grDevices::dev.cur()
  on.exit({
    if( grDevices::dev.cur() == device ) grDevices::dev.off(device)
  }, add = TRUE)
  expect_s3_class(plot_dagitty(renewed), "ggplot")
  print(plot_dagitty(renewed))

  fallback <- paste0("copy_nodes(): the copies could not be placed apart ",
                     "from existing nodes; coordinates were regenerated ",
                     "for every node.")
  expect_warning(
    kept <- copy_nodes(crowded, "C", suffix = "v", anchor = "A",
                       renew_coords = FALSE, print_edges = FALSE),
    fallback, fixed = TRUE
  )
  expect_identical(anyDuplicated(copy_node_places(kept)), 0L)

  spaced <- dagitty::dagitty(
    'dag { T [pos="0,0"] M [pos="2,1"] Y [pos="4,0"] T -> M M -> Y }'
  )
  before <- dagitty::coordinates(spaced)
  placed <- expect_no_warning(
    copy_nodes(spaced, "M", copy_suffix = "copy", renew_coords = FALSE,
               print_edges = FALSE)
  )
  after <- dagitty::coordinates(placed)
  existing <- c("T", "M", "Y")
  expect_equal(unname(after$x[existing]), unname(before$x[existing]))
  expect_equal(unname(after$y[existing]), unname(before$y[existing]))

  curved <- dagitty::dagitty(
    paste0('dag { T [pos="0,0"] M [pos="2,1"] Y [pos="4,0"] ',
           'T -> M [pos="1,0.75"] M -> Y }')
  )
  moved <- copy_nodes(curved, "M", copy_suffix = "copy", print_edges = FALSE)
  controls <- edge_controls(moved)
  moved_coordinates <- dagitty::coordinates(moved)
  expect_identical(nrow(controls), 1L)
  expect_identical(paste(controls$v, controls$e, controls$w), "T -> M")
  expect_true(all(is.finite(c(controls$x, controls$y))))
  chord_x <- moved_coordinates$x[["M"]] - moved_coordinates$x[["T"]]
  chord_y <- moved_coordinates$y[["M"]] - moved_coordinates$y[["T"]]
  expect_gt(
    abs((controls$x - moved_coordinates$x[["T"]]) * chord_y -
          (controls$y - moved_coordinates$y[["T"]]) * chord_x),
    0
  )

  coincident <- dagitty::dagitty(
    paste0('dag { A [pos="0,0"] B [pos="0,0"] C [pos="4,0"] ',
           'A -> B [pos="1,1"] B -> C }')
  )
  expect_error(
    copy_nodes(coincident, "C", print_edges = FALSE),
    "cannot preserve a curved edge's shape", fixed = TRUE
  )

  unchanged <- paste0(
    'dag {\nA [pos="0,0"]\nB [pos="0,0"]\nC [pos="4,0"]\n',
    'A -> B [pos="1,1"]\nB -> C\n\n',
    '"C_copy" [pos="5.3330000000000002,0"]\n',
    '"B" -> "C_copy"\n"C" -> "C_copy"\n}\n'
  )
  original <- dagitty::coordinates(coincident)
  restored <- expect_no_warning(
    copy_nodes(coincident, "C", renew_coords = FALSE, print_edges = FALSE)
  )
  expect_identical(as.character(restored), unchanged)
  restored_coordinates <- dagitty::coordinates(restored)
  starting <- c("A", "B", "C")
  expect_equal(unname(restored_coordinates$x[starting]),
               unname(original$x[starting]))
  expect_equal(unname(restored_coordinates$y[starting]),
               unname(original$y[starting]))
  restored_controls <- edge_controls(restored)
  expect_identical(nrow(restored_controls), 1L)
  expect_identical(paste(restored_controls$v, restored_controls$e,
                         restored_controls$w), "A -> B")
  expect_equal(c(restored_controls$x, restored_controls$y), c(1, 1))

  expect_error(
    copy_nodes(crowded, "C", renew_coords = NA, print_edges = FALSE),
    "renew_coords must be TRUE or FALSE.", fixed = TRUE
  )
  expect_error(
    copy_nodes(crowded, "C", renew_coords = "yes", print_edges = FALSE),
    "renew_coords must be TRUE or FALSE.", fixed = TRUE
  )
})


test_that("copy_nodes() copies bidirected and undirected edges and drops them under chain", {
  full <- build_graph(variables = "X", treatments = "T", outcomes = "Y", type = "full")
  copied <- copy_nodes(full, "X", print_edges = FALSE)
  expect_true(all(c("T <-> X_copy", "X_copy <-> Y", "X -> X_copy") %in% edge_keys(copied)))

  undirected <- dagitty::dagitty("dag { A -- B B -> C }")
  copied_undirected <- copy_nodes(undirected, "B", print_edges = FALSE)
  expect_true(all(c("A -- B_copy", "B_copy -> C", "B -> B_copy") %in%
                    edge_keys(copied_undirected)))

  mixed <- dagitty::dagitty("dag { A -> B B -> C A <-> B }")
  chained <- copy_nodes(mixed, "B", type = "chain", print_edges = FALSE)
  chained_keys <- edge_keys(chained)
  expect_true(all(c("B -> B_copy", "B_copy -> C") %in% chained_keys))
  expect_false("A <-> B_copy" %in% chained_keys)
  expect_true("A <-> B" %in% chained_keys)
})
