test_that("remove_nodes removes named nodes and every incident edge", {
  dag <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] U [latent] T -> M -> Y; U -> M; U -> Y; M <-> N; M -- S; I }"
  )
  original <- serialize(dag, NULL)

  expect_identical(names(formals(remove_nodes)), c("dag", "nodes"))
  expect_no_warning(removed <- remove_nodes(dag, "M"))
  expect_setequal(names(removed), c("T", "Y", "U", "N", "S", "I"))
  expect_identical(edge_keys(removed), "U -> Y")
  expect_identical(treatments(removed), "T")
  expect_identical(outcomes(removed), "Y")
  expect_identical(unobserved(removed), "U")
  expect_identical(serialize(dag, NULL), original)
})


test_that("remove_nodes retains neighbours that become isolated", {
  dag <- dagitty::dagitty("dag { A -> B -> C; I }")

  removed <- remove_nodes(dag, "B")
  expect_setequal(names(removed), c("A", "C", "I"))
  expect_length(edge_keys(removed), 0L)

  isolated <- remove_nodes(removed, "I")
  expect_setequal(names(isolated), c("A", "C"))
  expect_length(edge_keys(isolated), 0L)

  empty <- remove_nodes(dag, names(dag))
  expect_length(names(empty), 0L)
  expect_length(edge_keys(empty), 0L)
  expect_identical(causaliflower:::dag_graph_type(empty), "dag")
})


test_that("remove_nodes validates inputs and handles duplicate or absent names", {
  dag <- dagitty::dagitty("dag { A -> B -> C; I }")
  original <- serialize(dag, NULL)

  expect_error(remove_nodes("dag { A -> B }", "A"), "only dagitty objects")
  expect_error(remove_nodes(dag), "argument.*missing")
  expect_error(remove_nodes(dag, NULL), "'nodes' must be a character vector",
               fixed = TRUE)
  expect_error(remove_nodes(dag, list("A")), "'nodes' must be a character vector",
               fixed = TRUE)
  expect_error(remove_nodes(dag, 1), "'nodes' must be a character vector",
               fixed = TRUE)
  expect_error(remove_nodes(dag, character()), "at least one node name")
  expect_error(remove_nodes(dag, c("A", NA_character_)),
               "non-missing, non-empty names")
  expect_error(remove_nodes(dag, c("A", "")), "non-missing, non-empty names")

  expect_no_warning(duplicates <- remove_nodes(dag, c("B", "B")))
  expect_identical(duplicates, remove_nodes(dag, "B"))
  expect_warning(
    mixed <- remove_nodes(dag, c("B", "absent", "absent")),
    "remove_nodes(): nodes not present in the graph were ignored: absent.",
    fixed = TRUE
  )
  expect_identical(mixed, duplicates)
  expect_warning(
    unchanged <- remove_nodes(dag, "absent"),
    "remove_nodes(): nodes not present in the graph were ignored: absent.",
    fixed = TRUE
  )
  expect_identical(serialize(unchanged, NULL), original)
  expect_identical(serialize(dag, NULL), original)
})


test_that("remove_nodes removes only the selected native status declarations", {
  dag <- dagitty::dagitty(paste0(
    "dag { T1 [exposure] T2 [exposure] Y1 [outcome] Y2 [outcome] ",
    "U1 [latent] U2 [latent] A1 [adjusted] A2 [adjusted] ",
    "S1 [selected] S2 [selected] T1 -> Y1; T2 -> Y2; U1 -> T1; U2 -> T2 }"
  ))

  removed <- remove_nodes(dag, c("T1", "Y1", "U1", "A1", "S1"))
  expect_setequal(names(removed), c("T2", "Y2", "U2", "A2", "S2"))
  expect_identical(edge_keys(removed), c("T2 -> Y2", "U2 -> T2"))
  expect_identical(treatments(removed), "T2")
  expect_identical(outcomes(removed), "Y2")
  expect_identical(unobserved(removed), "U2")
  expect_identical(dagitty::adjustedNodes(removed), "A2")
  expect_identical(causaliflower:::selected_nodes(removed), "S2")
})


test_that("remove_nodes preserves statuses in compact native declarations", {
  compact <- structure("dag {\nU [latent] S [selected]\nU -> T\n}",
                       class = "dagitty")
  original <- serialize(compact, NULL)

  expect_no_warning(removed <- remove_nodes(compact, "U"))
  expect_setequal(names(removed), c("S", "T"))
  expect_length(edge_keys(removed), 0L)
  expect_identical(causaliflower:::selected_nodes(removed), "S")
  expect_identical(serialize(compact, NULL), original)

  dag <- structure(paste0(
    'dag named {\n bb="-1,-1,4,3"\n',
    'U [latent] S [selected] T [exposure] Y [outcome] ',
    'X [adjusted] V [latent]\n',
    'T [pos="0.123456789012345,0",note="keep node"]\n',
    'Y [pos="3,0"]\n',
    'U -> T\nT -> Y [pos="1.5,1.123456789012345",note="keep edge"]\n}'),
    class = "dagitty", provenance = list(author = "fixture"))
  original <- serialize(dag, NULL)
  coordinates <- dagitty::coordinates(dag)
  controls <- causaliflower:::edge_controls(dag)

  expect_no_warning(removed <- remove_nodes(dag, "U"))
  status <- dagitty::dagitty(as.character(removed))
  expect_setequal(names(removed), c("S", "T", "Y", "X", "V"))
  expect_identical(edge_keys(removed), "T -> Y")
  expect_identical(dagitty::exposures(status), "T")
  expect_identical(dagitty::outcomes(status), "Y")
  expect_identical(dagitty::latents(status), "V")
  expect_identical(dagitty::adjustedNodes(status), "X")
  expect_identical(causaliflower:::selected_nodes(status), "S")
  expect_identical(dagitty::coordinates(removed)$x,
                   coordinates$x[names(removed)])
  expect_identical(dagitty::coordinates(removed)$y,
                   coordinates$y[names(removed)])
  actual_controls <- causaliflower:::edge_controls(removed)
  expect_identical(actual_controls$v, controls$v)
  expect_identical(actual_controls$e, controls$e)
  expect_identical(actual_controls$w, controls$w)
  expect_identical(as.numeric(actual_controls$x), as.numeric(controls$x))
  expect_identical(as.numeric(actual_controls$y), as.numeric(controls$y))
  for( piece in c('dag named', 'bb="-1,-1,4,3"', 'note="keep node"',
                  'note="keep edge"',
                  'S [selected] T [exposure] Y [outcome] X [adjusted] V [latent]') ){
    expect_match(as.character(removed), piece, fixed = TRUE)
  }
  expect_identical(attr(removed, "provenance"), attr(dag, "provenance"))
  expect_identical(serialize(dag, NULL), original)
})


test_that("remove_nodes preserves surviving native metadata, curves and coordinates", {
  dag <- structure(paste0(
    'dag named { bb="-2,-2,5,3"\n',
    'T [exposure,pos="-1.123456789012345,0",label="Treatment",note="keep node"]\n',
    'M [latent,pos="0,1",label="remove node"]\n',
    'Y [outcome,selected,pos="2.987654321098765,0"]\n',
    'I [adjusted,pos="4,2"]\n',
    'T -> M [pos="0,0.123456789012345",note="remove edge"]\n',
    'M -> Y [pos="1.2,1.7"]\n',
    'T -> Y [pos="1,2",label="keep directed"]\n',
    'T <-> Y [pos="-0.0004,-1",note="keep symmetric"] }'),
    class = "dagitty", labels = c(T = "Treatment", M = "Mediator", Y = "Outcome"),
    provenance = list(author = "fixture"), custom = list(nodes = c("M", "Y")))
  original <- serialize(dag, NULL)
  coordinates <- dagitty::coordinates(dag)
  controls <- causaliflower:::edge_controls(dag)
  retained <- controls[controls$v != "M" & controls$w != "M", , drop = FALSE]

  removed <- remove_nodes(dag, "M")
  text <- as.character(removed)
  expect_setequal(names(removed), c("T", "Y", "I"))
  expect_identical(edge_keys(removed), c("T -> Y", "T <-> Y"))
  for( piece in c('dag named', 'bb="-2,-2,5,3"', 'label="Treatment"',
                  'note="keep node"', 'label="keep directed"',
                  'note="keep symmetric"') ){
    expect_match(text, piece, fixed = TRUE)
  }
  expect_false(grepl('label="remove node"', text, fixed = TRUE))
  expect_false(grepl('note="remove edge"', text, fixed = TRUE))
  actual <- dagitty::coordinates(removed)
  expect_identical(actual$x, coordinates$x[names(removed)])
  expect_identical(actual$y, coordinates$y[names(removed)])
  actual_controls <- causaliflower:::edge_controls(removed)
  expect_identical(actual_controls$v, retained$v)
  expect_identical(actual_controls$e, retained$e)
  expect_identical(actual_controls$w, retained$w)
  expect_identical(as.numeric(actual_controls$x), as.numeric(retained$x))
  expect_identical(as.numeric(actual_controls$y), as.numeric(retained$y))
  expect_identical(attr(removed, "labels"), c(T = "Treatment", Y = "Outcome"))
  expect_identical(attr(removed, "provenance"), attr(dag, "provenance"))
  expect_identical(attr(removed, "custom"), attr(dag, "custom"))
  expect_identical(treatments(removed), "T")
  expect_identical(outcomes(removed), "Y")
  expect_identical(causaliflower:::selected_nodes(removed), "Y")
  expect_identical(dagitty::adjustedNodes(removed), "I")
  expect_length(unobserved(removed), 0L)
  expect_identical(serialize(dag, NULL), original)
})


test_that("remove_nodes matches complete quoted names and leaves text values alone", {
  dag <- structure(paste0(
    'dag named { "soil moisture" [latent,note="remove"] ',
    '"a [selected]" [selected] "a->b" [adjusted] "a;b" ',
    'Y [outcome,note="soil moisture a [selected] a->b a;b"] ',
    '"soil moisture" -> Y; "a [selected]" -> Y; ',
    '"a->b" -> Y; "a;b" -> Y; soil -> Y }'), class = "dagitty")

  removed <- remove_nodes(dag, c("soil moisture", "a [selected]", "a->b", "a;b"))
  expect_setequal(names(removed), c("soil", "Y"))
  expect_identical(edge_keys(removed), "soil -> Y")
  expect_identical(dagitty::outcomes(removed), "Y")
  expect_length(dagitty::latents(removed), 0L)
  expect_false(grepl("\"a [selected]\" [selected]", as.character(removed), fixed = TRUE))
  expect_length(dagitty::adjustedNodes(removed), 0L)
  expect_match(as.character(removed),
               'note="soil moisture a [selected] a->b a;b"', fixed = TRUE)
})


test_that("remove_nodes handles native chains and grouped edge statements", {
  chain <- structure("dag named { A -> B -> C -> D; E -> F; I }",
                     class = "dagitty")
  removed <- remove_nodes(chain, c("B", "E"))
  expect_setequal(names(removed), c("A", "C", "D", "F", "I"))
  expect_identical(edge_keys(removed), "C -> D")
  expect_match(as.character(removed), "dag named", fixed = TRUE)

  grouped <- structure("dag { A -> { B C }; B -> D; I }", class = "dagitty")
  removed <- remove_nodes(grouped, "B")
  expect_setequal(names(removed), c("A", "C", "D", "I"))
  expect_identical(edge_keys(removed), "A -> C")
})


test_that("remove_nodes preserves graph types and permits remaining cycles", {
  for( spec in c("dag { A -> B; B -> C; I }",
                 "pdag { A -- B; B -> C; I }",
                 "mag { A <-> B; B -> C; I }") ){
    dag <- dagitty::dagitty(spec)
    removed <- remove_nodes(dag, "B")
    expect_identical(causaliflower:::dag_graph_type(removed),
                     causaliflower:::dag_graph_type(dag))
    expect_setequal(names(removed), c("A", "C", "I"))
    expect_length(edge_keys(removed), 0L)

    empty <- remove_nodes(dag, names(dag))
    expect_identical(causaliflower:::dag_graph_type(empty),
                     causaliflower:::dag_graph_type(dag))
    expect_length(names(empty), 0L)
    expect_length(edge_keys(empty), 0L)
  }

  cyclic <- dagitty::dagitty("dag { A -> B; B -> A; B -> C }")
  expect_warning(remaining <- remove_nodes(cyclic, "C"), "contains cycles")
  expect_setequal(edge_keys(remaining), c("A -> B", "B -> A"))
  expect_no_warning(acyclic <- remove_nodes(cyclic, "B"))
  expect_true(dagitty::isAcyclic(acyclic))
  expect_setequal(names(acyclic), c("A", "C"))
})


test_that("remove_edges drops only nodes isolated by the current removal", {
  dag <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] U [latent] I [selected] T -> Y; U -> T; A <-> B; B -- C; D }"
  )
  original <- serialize(dag, NULL)

  expect_warning(
    directed <- remove_edges(dag, "U -> T"),
    "remove_edges(): removed nodes with no remaining edges: U.", fixed = TRUE
  )
  expect_setequal(names(directed), setdiff(names(dag), "U"))
  expect_identical(treatments(directed), "T")
  expect_identical(outcomes(directed), "Y")
  expect_length(unobserved(directed), 0L)
  expect_identical(causaliflower:::selected_nodes(directed), "I")

  expect_warning(
    bidirected <- remove_edges(dag, "B <-> A"),
    "remove_edges(): removed nodes with no remaining edges: A.", fixed = TRUE
  )
  expect_setequal(names(bidirected), setdiff(names(dag), "A"))
  expect_true("B -- C" %in% edge_keys(bidirected))

  expect_warning(
    undirected <- remove_edges(dag, "C -- B"),
    "remove_edges(): removed nodes with no remaining edges: C.", fixed = TRUE
  )
  expect_setequal(names(undirected), setdiff(names(dag), "C"))
  expect_true("A <-> B" %in% edge_keys(undirected))

  expect_warning(
    both <- remove_edges(dag, c("B <-> A", "C -- B")),
    "remove_edges(): removed nodes with no remaining edges: A, B, C.", fixed = TRUE
  )
  expect_setequal(names(both), c("T", "Y", "U", "I", "D"))
  expect_setequal(edge_keys(both), c("T -> Y", "U -> T"))

  expect_warning(
    empty <- remove_edges(dag, dagitty::edges(dag)),
    "remove_edges(): removed nodes with no remaining edges: A, B, C, T, U, Y.",
    fixed = TRUE
  )
  expect_setequal(names(empty), c("D", "I"))
  expect_length(edge_keys(empty), 0L)
  expect_identical(causaliflower:::selected_nodes(empty), "I")
  expect_identical(serialize(dag, NULL), original)
})


test_that("remove_edges preserves surviving curves and drops node label entries", {
  dag <- structure(paste0(
    'dag named { bb="-1,-1,5,3"\n',
    'T [exposure,pos="0,0",label="Treatment"]\n',
    'Y [outcome,pos="4,0"] U [latent,pos="0,2",note="remove node"]\n',
    'I [selected,pos="4,2"]\n',
    'U -> T [pos="-1,1",note="remove edge"]\n',
    'T -> Y [pos="2,0.123456789012345",note="keep edge"] }'),
    class = "dagitty", labels = c(T = "Treatment", U = "Unobserved", Y = "Outcome"),
    provenance = list(author = "fixture"))
  coordinates <- dagitty::coordinates(dag)
  controls <- causaliflower:::edge_controls(dag)
  retained <- controls[controls$v == "T", , drop = FALSE]

  expect_warning(removed <- remove_edges(dag, "U -> T"),
                 "no remaining edges.*U")
  expect_setequal(names(removed), c("I", "T", "Y"))
  actual_controls <- causaliflower:::edge_controls(removed)
  expect_identical(actual_controls$v, retained$v)
  expect_identical(actual_controls$e, retained$e)
  expect_identical(actual_controls$w, retained$w)
  expect_identical(as.numeric(actual_controls$x), as.numeric(retained$x))
  expect_identical(as.numeric(actual_controls$y), as.numeric(retained$y))
  actual <- dagitty::coordinates(removed)
  expect_identical(actual$x, coordinates$x[names(removed)])
  expect_identical(actual$y, coordinates$y[names(removed)])
  expect_identical(attr(removed, "labels"), c(T = "Treatment", Y = "Outcome"))
  expect_identical(attr(removed, "provenance"), attr(dag, "provenance"))
  expect_identical(treatments(removed), "T")
  expect_identical(outcomes(removed), "Y")
  expect_identical(causaliflower:::selected_nodes(removed), "I")
  for( piece in c('dag named', 'bb="-1,-1,5,3"', 'label="Treatment"',
                  'note="keep edge"') ){
    expect_match(as.character(removed), piece, fixed = TRUE)
  }
  expect_false(grepl('note="remove node"', as.character(removed), fixed = TRUE))
  expect_false(grepl('note="remove edge"', as.character(removed), fixed = TRUE))
})


test_that("keep_edges and native edge removal continue to retain isolated nodes", {
  dag <- dagitty::dagitty(
    "dag { T [exposure] Y [outcome] U [latent] I [selected] U -> T; T -> Y }"
  )
  kept <- keep_edges("T -> Y", dag)
  expect_setequal(names(kept), names(dag))
  expect_identical(unobserved(kept), "U")

  empty <- causaliflower:::native_remove_edges(dag, dagitty::edges(dag), "fixture")
  expect_setequal(names(empty), names(dag))
  expect_length(edge_keys(empty), 0L)
  expect_identical(treatments(empty), "T")
  expect_identical(outcomes(empty), "Y")
  expect_identical(unobserved(empty), "U")
  expect_identical(causaliflower:::selected_nodes(empty), "I")
})


test_that("remove_nodes can remove a latent node added with explicit children", {
  dag <- build_graph("X", "T", "Y", mediators = c("M1", "M2"),
                     type = "ordered")
  expect_warning(
    added <- add_nodes(dag, "U", role = "latent", children = c("M1", "M2"),
                       print_edges = FALSE),
    "parents/children input overrides role-based edges", fixed = TRUE
  )
  expect_setequal(setdiff(edge_keys(added), edge_keys(dag)),
                  c("U -> M1", "U -> M2"))

  removed <- remove_nodes(added, "U")
  expect_identical(edge_keys(removed), edge_keys(dag))
  expect_identical(names(removed), names(dag))
  expect_identical(treatments(removed), "T")
  expect_identical(outcomes(removed), "Y")
  expect_length(unobserved(removed), 0L)
})
