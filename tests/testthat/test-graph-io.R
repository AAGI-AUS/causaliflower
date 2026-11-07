test_that("string_to_dag parses chains and every supported operator", {
  chain <- string_to_dag("A -> B -> C")
  mixed <- string_to_dag("A <- B; C <-> D; E -- F")

  expect_identical(edge_keys(chain), c("A -> B", "B -> C"))
  expect_true("B -> A" %in% edge_keys(mixed))
  expect_true("C <-> D" %in% edge_keys(mixed))
  expect_true("E -- F" %in% edge_keys(mixed))
})

test_that("string_to_dag handles special and isolated names", {
  cases <- c(
    '"soil moisture" -> "air temp"',
    "soil-moisture -> air-temp",
    "rate% -> yield%"
  )
  expected <- list(
    c("soil moisture", "air temp"),
    c("soil-moisture", "air-temp"),
    c("rate%", "yield%")
  )

  for (i in seq_along(cases)) {
    expect_setequal(names(string_to_dag(cases[i])), expected[[i]])
  }
  expect_identical(names(string_to_dag("isolated node")), "isolated node")
  expect_length(names(string_to_dag("")), 0L)

  operator_name <- string_to_dag('"A->B" -> C')
  separator_name <- string_to_dag('"A;B" -> C')
  expect_setequal(names(operator_name), c("A->B", "C"))
  expect_identical(edge_keys(operator_name), "A->B -> C")
  expect_setequal(names(separator_name), c("A;B", "C"))
  expect_identical(edge_keys(separator_name), "A;B -> C")
  expect_identical(names(string_to_dag('"A->B"')), "A->B")
})

test_that("string parsing rejects malformed edge text", {
  expect_error(string_to_dag("A ->"), "node names cannot be empty")
  expect_error(string_to_dag("A => B"), "could not parse")
  expect_error(string_to_dag('"A -> B'), "unmatched quote")
})

test_that("keep_edges accepts strings, tables, dagitty objects, and reverse arrows", {
  dag <- build_graph(c("X1", "X2"), "T", "Y", type = "saturated")

  from_strings <- keep_edges(c("X1 -> T", "T -> Y"), dag)
  from_reverse <- keep_edges(c("T <- X1", "T -> Y"), dag)
  from_table <- keep_edges(
    data.frame(v = c("X1", "T"), e = "->", w = c("T", "Y")),
    dag
  )
  from_dag <- keep_edges(from_strings, dag)

  expected <- sort(c("X1 -> T", "T -> Y"))
  expect_identical(edge_keys(from_strings), expected)
  expect_identical(edge_keys(from_reverse), expected)
  expect_identical(edge_keys(from_table), expected)
  expect_identical(edge_keys(from_dag), expected)
})

test_that("edge editors use canonical edge argument names", {
  graph <- string_to_dag(c("A -> B", "B -> C"))

  expect_identical(names(formals(add_edges)), c("dag", "edges"))
  expect_identical(names(formals(keep_edges)),
                   c("edges", "dag", "edges_to_keep"))
  expect_identical(names(formals(remove_edges)), c("dag", "edges"))
  expect_identical(edge_keys(keep_edges(edges = "A -> B", dag = graph)),
                   "A -> B")
  expect_warning(
    removed <- remove_edges(dag = graph, edges = "A -> B"),
    "no remaining edges.*A"
  )
  expect_identical(edge_keys(removed), "B -> C")
  expect_warning(
    legacy <- keep_edges(edges_to_keep = "A -> B", dag = graph),
    paste0("keep_edges(): 'edges_to_keep' is deprecated as of 0.0.2. ",
           "Use 'edges' instead."),
    fixed = TRUE
  )
  expect_identical(edge_keys(legacy), "A -> B")
  expect_error(
    keep_edges("A -> B", graph, edges_to_keep = "B -> C"),
    "supply either 'edges' or its deprecated alias 'edges_to_keep'",
    fixed = TRUE
  )
})

test_that("add_edges accepts supported selectors and preserves the input", {
  graph <- string_to_dag(c("A -> B", "B -> C", "D"))
  before <- serialize(graph, NULL)
  from_chain <- add_edges(graph, "A -> D -> C")
  from_reverse <- add_edges(graph, "D <- A")
  from_table <- add_edges(graph, data.frame(v = "A", e = "->", w = "D"))
  from_dag <- add_edges(graph, string_to_dag("A -> D"))

  expect_setequal(edge_keys(from_chain),
                  c("A -> B", "A -> D", "B -> C", "D -> C"))
  for( result in list(from_reverse, from_table, from_dag) ){
    expect_true("A -> D" %in% edge_keys(result))
  }
  expect_identical(serialize(graph, NULL), before)
})

test_that("add_edges handles duplicates, invalid endpoints, and cycles atomically", {
  graph <- string_to_dag(c("A -> B", "B -> C", "D"))
  before <- serialize(graph, NULL)

  expect_warning(
    unchanged <- add_edges(graph, c("A -> B", "A -> B")),
    "already present"
  )
  expect_identical(serialize(unchanged, NULL), before)
  expect_warning(
    mixed <- add_edges(graph, c("A -> B", "A -> D")),
    "already present"
  )
  expect_true("A -> D" %in% edge_keys(mixed))
  expect_error(add_edges(graph, "A -> absent"), "not in the graph")
  expect_error(add_edges(graph, "A -> A"), "different nodes")
  expect_error(add_edges(graph, "C -> A"), "directed cycle")
  expect_error(add_edges(graph, character()), "at least one edge")
  expect_error(add_edges(graph, 1), "character vector")
  expect_identical(serialize(graph, NULL), before)
})

test_that("add_edges preserves graph type and native metadata", {
  pdag <- dagitty::dagitty("pdag { A -- B; B -> C; D }")
  mag <- dagitty::dagitty("mag { A <-> B; B -> C; D }")
  expect_identical(causaliflower:::dag_graph_type(add_edges(pdag, "C -> D")),
                   "pdag")
  expect_identical(causaliflower:::dag_graph_type(add_edges(mag, "C -> D")),
                   "mag")

  graph <- structure(
    'dag named { bb="-1,-1,3,3" A [label="first",note="node"] B C
       A -> B [label="old",note="edge"] }',
    class = "dagitty", provenance = "fixture"
  )
  result <- add_edges(graph, "B -> C")
  for( piece in c("dag named", 'bb="-1,-1,3,3"', 'label="first"',
                  'note="node"', 'label="old"', 'note="edge"') ){
    expect_true(grepl(piece, as.character(result), fixed = TRUE))
  }
  expect_identical(attr(result, "provenance"), "fixture")
})

test_that("keep_edges preserves nodes, roles, coordinates, and symmetric edges", {
  dag <- dagitty::dagitty(
    'dag { "soil moisture" [exposure] "yield%" [outcome]
           "air-temp" -> "soil moisture"
           "soil moisture" -> "yield%" A <-> B }'
  )
  dagitty::coordinates(dag) <- list(
    x = stats::setNames(seq_along(names(dag)), names(dag)),
    y = stats::setNames(rep(0, length(names(dag))), names(dag))
  )
  original_coordinates <- dagitty::coordinates(dag)

  kept <- keep_edges(
    c('"air-temp" -> "soil moisture"', '"soil moisture" -> "yield%"',
      "B <-> A"),
    dag
  )

  expect_setequal(names(kept), names(dag))
  expect_identical(treatments(kept), "soil moisture")
  expect_identical(causaliflower:::get_outcome_names(kept), "yield%")
  expect_equal(dagitty::coordinates(kept), original_coordinates)
  expect_true("A <-> B" %in% edge_keys(kept))
})

test_that("keep_edges has deliberate empty and missing-edge contracts", {
  dag <- build_graph(c("X1", "X2"), "T", "Y", type = "saturated")

  expect_error(keep_edges(character(), dag), "no edges were supplied")
  expect_error(keep_edges("A => B", dag), "could not parse")
  expect_error(keep_edges(data.frame(a = 1, b = 2), dag),
               "v, e, and w")

  expect_warning(empty <- keep_edges("absent -> edge", dag),
                 "not present")
  expect_length(edge_keys(empty), 0L)
  expect_setequal(names(empty), names(dag))
  expect_identical(treatments(empty), "T")
  expect_identical(causaliflower:::get_outcome_names(empty), "Y")

  expect_warning(
    expect_error(
      keep_edges(data.frame(v = NA, e = NA, w = NA), dag),
      "no complete edges"
    ),
    "containing missing values"
  )
})

test_that("remove_edges accepts strings, tables, dagitty objects, and reverse arrows", {
  dag <- string_to_dag(c("A -> B -> C", "A -> C"))

  expect_warning(from_chain <- remove_edges(dag, "A -> B -> C"),
                 "no remaining edges.*B")
  from_compact <- remove_edges(dag, "A->B")
  from_reverse <- remove_edges(dag, "B <- A")
  from_table <- remove_edges(
    dag,
    data.frame(v = "A", e = "->", w = "B")
  )
  selector <- string_to_dag("A -> B")
  from_dag <- remove_edges(dag, selector)
  from_duplicates <- remove_edges(dag, c("A -> B", "A -> B"))

  expect_identical(edge_keys(from_chain), "A -> C")
  expected <- sort(c("A -> C", "B -> C"))
  expect_identical(edge_keys(from_compact), expected)
  expect_identical(edge_keys(from_reverse), expected)
  expect_identical(edge_keys(from_table), expected)
  expect_identical(edge_keys(from_dag), expected)
  expect_identical(edge_keys(from_duplicates), expected)
  expect_setequal(edge_keys(dag), c("A -> B", "A -> C", "B -> C"))

  expect_warning(
    wrong_direction <- remove_edges(dag, "B -> A"),
    "not present"
  )
  expect_identical(as.character(wrong_direction), as.character(dag))
})

test_that("remove_edges matches symmetric edges without changing graph type", {
  pdag <- dagitty::dagitty("pdag { A -- B; B -> C }")
  mag <- dagitty::dagitty("mag { A <-> B; B -> C }")

  expect_warning(without_undirected <- remove_edges(pdag, "B -- A"),
                 "no remaining edges.*A")
  expect_warning(without_bidirected <- remove_edges(mag, "B <-> A"),
                 "no remaining edges.*A")

  expect_identical(causaliflower:::dag_graph_type(without_undirected), "pdag")
  expect_identical(causaliflower:::dag_graph_type(without_bidirected), "mag")
  expect_identical(edge_keys(without_undirected), "B -> C")
  expect_identical(edge_keys(without_bidirected), "B -> C")
})

test_that("remove_edges drops newly isolated nodes and their statuses", {
  dag <- dagitty::dagitty(
    "dag { T [exposure,adjusted] Y [outcome,selected] U [latent] U -> T; T -> Y }"
  )
  dagitty::coordinates(dag) <- list(
    x = c(T = 0, Y = 2),
    y = c(T = 0, Y = 0)
  )
  original <- serialize(dag, NULL)

  expect_warning(
    empty <- remove_edges(dag, c("U -> T", "T -> Y")),
    "no remaining edges.*T.*U.*Y"
  )

  expect_length(edge_keys(empty), 0L)
  expect_length(names(empty), 0L)
  expect_identical(causaliflower:::dag_graph_type(empty), "dag")
  expect_length(treatments(empty), 0L)
  expect_length(outcomes(empty), 0L)
  expect_length(unobserved(empty), 0L)
  expect_length(dagitty::adjustedNodes(empty), 0L)
  expect_length(causaliflower:::selected_nodes(empty), 0L)
  expect_length(dagitty::coordinates(empty)$x, 0L)
  expect_length(dagitty::coordinates(empty)$y, 0L)
  expect_identical(serialize(dag, NULL), original)
})

test_that("remove_edges reports invalid input and missing edges", {
  dag <- string_to_dag(c("A -> B", "B -> C"))

  expect_error(
    remove_edges(list(), "A -> B"),
    "only dagitty objects"
  )
  expect_error(remove_edges(dag), "argument.*missing")
  expect_error(remove_edges(dag, character()), "at least one edge")
  expect_error(remove_edges(dag, "  "), "at least one edge")
  expect_error(remove_edges(dag, "A => B"), "could not parse")
  expect_error(remove_edges(dag, 1), "character vector, an edge table")
  expect_error(
    remove_edges(dag, data.frame(a = 1, b = 2)),
    "v, e, and w"
  )
  expect_warning(
    expect_error(
      remove_edges(dag, data.frame(v = NA, e = NA, w = NA)),
      "no complete edges"
    ),
    "missing values"
  )

  expect_warning(
    unchanged <- remove_edges(dag, "absent -> edge"),
    "not present"
  )
  expect_identical(as.character(unchanged), as.character(dag))

  incomplete <- data.frame(
    v = c("A", NA), e = c("->", NA), w = c("B", NA)
  )
  expect_warning(
    expect_warning(with_complete <- remove_edges(dag, incomplete),
                   "no remaining edges.*A"),
    "missing values"
  )
  expect_identical(edge_keys(with_complete), "B -> C")
})

test_that("remove_edges warns when a pre-existing cycle remains", {
  cyclic <- string_to_dag(c("A -> B", "B -> A", "B -> C"))

  expect_warning(
    expect_warning(remaining_cycle <- remove_edges(cyclic, "B -> C"),
                   "no remaining edges.*C"),
    "contains cycles"
  )
  expect_setequal(edge_keys(remaining_cycle), c("A -> B", "B -> A"))
})


test_that("selected status ignores brackets inside quoted node names", {
  for( node in c("a[b", "a]b", "a[selected]", "a [selected]", "soil moisture",
                 "a;b", "a->b") ){
    selected <- dagitty::dagitty(paste0(
      'dag { "', node, '" [selected]; T -> Y; Y -> "', node, '" }'
    ))
    original <- selected
    expect_warning(removed <- remove_edges(selected, "T -> Y"),
                   "no remaining edges.*T")
    expect_identical(causaliflower:::selected_nodes(removed), node)
    expect_setequal(names(removed), setdiff(names(selected), "T"))
    expect_error(get_identification(selected, "T", "Y", quiet = TRUE),
                 "Adjusted/selected node")
    expect_identical(selected, original)

    ordinary <- dagitty::dagitty(paste0(
      'dag { "', node, '"; T -> Y; Y -> "', node, '" }'
    ))
    original <- ordinary
    expect_identical(causaliflower:::selected_nodes(ordinary), character())
    expect_true(attr(get_identification(ordinary, "T", "Y", quiet = TRUE),
                     "analyses")$total$identifiable)
    expect_identical(ordinary, original)
  }
})

test_that("selected status reads escaped quotes and literal backslashes", {
  for( node in c('a"b', 'a"[selected]b', 'a\\"[selected]b',
                 "a\\b", "a\\\\[selected]b") ){
    quoted <- paste0('"', gsub('"', '\\"', node, fixed = TRUE), '"')
    selected <- dagitty::dagitty(paste0(
      "dag { ", quoted, " [selected]; T -> Y }"
    ))
    original <- selected
    expect_identical(causaliflower:::selected_nodes(selected), node)
    expect_error(get_identification(selected, "T", "Y", quiet = TRUE),
                 "Adjusted/selected node")
    expect_error(remove_edges(selected, "T -> Y"),
                 "node names cannot contain quotes, backslashes, or line breaks")
    expect_identical(selected, original)

    ordinary <- dagitty::dagitty(paste0("dag { ", quoted, "; T -> Y }"))
    expect_identical(causaliflower:::selected_nodes(ordinary), character())
  }
})

test_that("edge editing preserves quoted selected nodes and other metadata", {
  dag <- dagitty::dagitty('dag {
    T [exposure,adjusted] Y [outcome] U [latent]
    "a [selected]" [selected] "isolated node"
    U -> T; T -> Y; Y -> "a [selected]"
  }')
  dagitty::coordinates(dag) <- list(
    x = c(T = 0, "a [selected]" = 3),
    y = c(T = 1, "a [selected]" = 2)
  )
  original <- dag
  selector <- data.table::as.data.table(dagitty::edges(dag))
  original_selector <- data.table::copy(selector)
  expect_warning(empty <- remove_edges(dag, selector), "no remaining edges")
  kept <- keep_edges('Y -> "a [selected]"', dag)

  expect_length(edge_keys(empty), 0L)
  expect_identical(names(empty), "isolated node")
  expect_length(treatments(empty), 0L)
  expect_length(outcomes(empty), 0L)
  expect_length(unobserved(empty), 0L)
  expect_length(dagitty::adjustedNodes(empty), 0L)
  expect_length(causaliflower:::selected_nodes(empty), 0L)
  expect_identical(edge_keys(kept), "Y -> a [selected]")
  expect_setequal(names(kept), names(dag))
  expect_identical(causaliflower:::dag_graph_type(kept), "dag")
  expect_identical(treatments(kept), "T")
  expect_identical(outcomes(kept), "Y")
  expect_identical(unobserved(kept), "U")
  expect_identical(dagitty::adjustedNodes(kept), "T")
  expect_identical(causaliflower:::selected_nodes(kept), "a [selected]")
  expect_equal(dagitty::coordinates(kept), dagitty::coordinates(dag))
  expect_identical(dag, original)
  expect_identical(selector, original_selector)
})

test_that("selected status ignores the word selected inside quoted values", {
  dag <- dagitty::dagitty('dag {
    T [label="a selected covariate",pos="0,0"]
    Y [note="not selected for adjustment",pos="1,0"]
    Z [pos="0,1"]
    T -> Y; Z -> Y
  }')

  expect_identical(causaliflower:::selected_nodes(dag), character())
  expect_true(attr(get_identification(dag, "T", "Y", quiet = TRUE),
                   "analyses")$total$identifiable)

  expect_warning(removed <- remove_edges(dag, "Z -> Y"),
                 "no remaining edges.*Z")
  expect_identical(edge_keys(removed), "T -> Y")
  expect_identical(causaliflower:::selected_nodes(removed), character())

  marked <- dagitty::dagitty(
    'dag { S [label="a selected covariate",selected] T -> S }'
  )
  expect_identical(causaliflower:::selected_nodes(marked), "S")
})

native_fixture <- function(){
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
  attr(dag, "provenance") <- list(author = "fixture")
  dag
}

test_that("remove_edges() keeps the graph's native attributes and coordinates", {
  dag <- native_fixture()
  coordinates <- dagitty::coordinates(dag)

  expect_warning(removed <- remove_edges(dag, "M -> Y"),
                 "no remaining edges.*Y")
  text <- as.character(removed)

  expect_match(text, "dag named", fixed = TRUE)
  expect_match(text, 'bb="-2,-2,5,3"', fixed = TRUE)
  expect_match(text, 'label="Treatment"', fixed = TRUE)
  expect_match(text, 'note="keep node"', fixed = TRUE)
  expect_match(text, 'label="curved"', fixed = TRUE)
  expect_match(text, 'note="keep edge"', fixed = TRUE)
  expect_identical(edge_keys(removed), "T -> M")
  expect_equal(dagitty::coordinates(removed),
               lapply(coordinates, function(axis) axis[c("M", "T")]))
  expect_identical(names(removed), c("M", "T"))
  expect_identical(treatments(removed), "T")
  expect_length(outcomes(removed), 0L)
  expect_identical(unobserved(removed), "M")
  expect_identical(dagitty::adjustedNodes(removed), "T")
  expect_length(causaliflower:::selected_nodes(removed), 0L)
  expect_identical(attr(removed, "provenance"), list(author = "fixture"))
  controls <- causaliflower:::edge_controls(removed)
  expect_identical(nrow(controls), 1L)
  expect_equal(c(controls$x, controls$y), c(1, 0.75))

  expect_warning(emptied <- remove_edges(removed, "T -> M"),
                 "no remaining edges.*M.*T")
  text <- as.character(emptied)

  expect_match(text, "dag named", fixed = TRUE)
  expect_match(text, 'bb="-2,-2,5,3"', fixed = TRUE)
  expect_false(grepl('label="Treatment"', text, fixed = TRUE))
  expect_false(grepl('note="keep node"', text, fixed = TRUE))
  expect_false(grepl('label="curved"', text, fixed = TRUE))
  expect_false(grepl('note="keep edge"', text, fixed = TRUE))
  expect_false(grepl('pos="1.000,0.750"', text, fixed = TRUE))
  expect_identical(edge_keys(emptied), character())
  expect_length(names(emptied), 0L)
  expect_length(dagitty::coordinates(emptied)$x, 0L)
  expect_length(dagitty::coordinates(emptied)$y, 0L)
  expect_identical(nrow(causaliflower:::edge_controls(emptied)), 0L)
})

test_that("keep_edges() keeps native attributes of the retained edges and nodes", {
  dag <- native_fixture()

  kept <- keep_edges("T -> M", dag)
  text <- as.character(kept)

  expect_match(text, "dag named", fixed = TRUE)
  expect_match(text, 'bb="-2,-2,5,3"', fixed = TRUE)
  expect_match(text, 'label="Treatment"', fixed = TRUE)
  expect_match(text, 'note="keep node"', fixed = TRUE)
  expect_match(text, 'label="curved"', fixed = TRUE)
  expect_match(text, 'note="keep edge"', fixed = TRUE)
  expect_identical(edge_keys(kept), "T -> M")
  expect_identical(names(kept), names(dag))
  expect_equal(dagitty::coordinates(kept), dagitty::coordinates(dag))
  expect_identical(causaliflower:::selected_nodes(kept), "Y")
  expect_identical(attr(kept, "provenance"), list(author = "fixture"))
})

test_that("native edge removal handles chains, reversed, symmetric, and quoted edges", {
  dag <- dagitty::dagitty(
    'dag { "a b" -> B -> C ; C <- D ; E <-> F ; G -- H ; H [label="h"] }'
  )

  chain <- remove_edges(dag, "B -> C")
  expect_identical(edge_keys(chain),
                   sort(c("a b -> B", "D -> C", "E <-> F", "G -- H")))
  expect_identical(names(chain), names(dag))

  expect_warning(symmetric <- remove_edges(dag, "F <-> E"),
                 "no remaining edges.*E.*F")
  expect_identical(edge_keys(symmetric),
                   sort(c("a b -> B", "B -> C", "D -> C", "G -- H")))
  expect_identical(names(symmetric), setdiff(names(dag), c("E", "F")))

  expect_warning(undirected <- remove_edges(dag, "G -- H"),
                 "no remaining edges.*G.*H")
  expect_identical(edge_keys(undirected),
                   sort(c("a b -> B", "B -> C", "D -> C", "E <-> F")))
  expect_identical(names(undirected), setdiff(names(dag), c("G", "H")))
  expect_false(grepl('label="h"', as.character(undirected), fixed = TRUE))

  expect_warning(isolated <- remove_edges(dag, c("a b -> B", "B -> C")),
                 "no remaining edges")
  expect_identical(edge_keys(isolated), sort(c("D -> C", "E <-> F", "G -- H")))
  expect_identical(names(isolated), setdiff(names(dag), c("a b", "B")))

  occurrences <- function(graph, node){
    text <- gsub("\\[[^]]*\\]", " ", as.character(graph))
    text <- gsub("<->|<-|->|--|[{}]", " ", text)
    sum(strsplit(text, "[[:space:];]+")[[1L]] == node)
  }

  expect_identical(occurrences(dag, "B"), 3L)
  expect_identical(occurrences(chain, "B"), 2L)
  expect_identical(occurrences(chain, "C"), 2L)
  expect_identical(occurrences(undirected, "G"), 0L)
  expect_identical(occurrences(isolated, "B"), 0L)

  chained <- structure("dag {\nA -> B -> C -> D\nE -> F\n}",
                       class = "dagitty")

  expect_warning(first <- remove_edges(chained, "A -> B"),
                 "no remaining edges.*A")
  expect_identical(edge_keys(first), sort(c("B -> C", "C -> D", "E -> F")))
  expect_identical(names(first), setdiff(names(chained), "A"))
  expect_identical(occurrences(first, "A"), 0L)

  expect_warning(last <- remove_edges(chained, "C -> D"),
                 "no remaining edges.*D")
  expect_identical(edge_keys(last), sort(c("A -> B", "B -> C", "E -> F")))
  expect_identical(names(last), setdiff(names(chained), "D"))
  expect_identical(occurrences(last, "D"), 0L)

  middle <- remove_edges(chained, "B -> C")
  expect_identical(edge_keys(middle), sort(c("A -> B", "C -> D", "E -> F")))
  expect_identical(names(middle), names(chained))
  expect_identical(occurrences(middle, "B"), 1L)
  expect_identical(occurrences(middle, "C"), 1L)
})

test_that("native edge removal falls back to rebuilding for grouped statements", {
  dag <- structure('dag { X -> { Y Z } }', class = "dagitty")
  expect_match(as.character(dag), "X -> { Y Z }", fixed = TRUE)

  expect_warning(removed <- remove_edges(dag, "X -> Y"),
                 "no remaining edges.*Y")

  expect_identical(edge_keys(removed), "X -> Z")
  expect_setequal(names(removed), c("X", "Z"))
  expect_false(grepl("{ Y Z }", as.character(removed), fixed = TRUE))
  expect_false(grepl("-> {", as.character(removed), fixed = TRUE))
})


test_that("add_nodes() and connect_nodes() keep the graph name, bounding box, and node attributes", {
  dag <- dagitty::dagitty(paste0(
    'dag named {\n bb="-2,-2,5,3"\n',
    'T [exposure,pos="0,0",label="Treatment",note="keep node"]\n',
    'Z [pos="1,1"]\n',
    'Y [outcome,pos="4,0"]\n',
    'T -> Y [pos="2,0.5"]\nZ -> T\n}'
  ))
  attr(dag, "provenance") <- list(author = "fixture")
  original <- as.character(dag)
  coordinates <- dagitty::coordinates(dag)

  added <- add_nodes(dag, nodes = "W", role = "confounder", print_edges = FALSE)
  text <- as.character(added)
  expect_true(startsWith(text, "dag named"))
  expect_true(grepl('bb="-2,-2,5,3"', text, fixed = TRUE))
  expect_true(grepl('label="Treatment"', text, fixed = TRUE))
  expect_true(grepl('note="keep node"', text, fixed = TRUE))
  expect_identical(attr(added, "provenance"), list(author = "fixture"))
  expect_setequal(names(added), c("T", "Z", "Y", "W"))
  expect_identical(treatments(added), "T")
  expect_identical(causaliflower:::get_outcome_names(added), "Y")
  expect_true(all(c("W -> T", "W -> Y") %in% edge_keys(added)))
  added_coordinates <- dagitty::coordinates(added)
  expect_equal(added_coordinates$x[names(coordinates$x)], coordinates$x)
  expect_equal(added_coordinates$y[names(coordinates$y)], coordinates$y)
  controls <- causaliflower:::edge_controls(added)
  expect_identical(paste(controls$v, controls$e, controls$w), "T -> Y")
  expect_equal(c(controls$x, controls$y), c(2, 0.5))

  expect_error(
    add_nodes(dag, nodes = 'weird"name', role = "confounder", print_edges = FALSE),
    "node names cannot contain quotes, backslashes, or line breaks", fixed = TRUE
  )

  latent_added <- add_nodes(dag, nodes = "L", role = "latent", print_edges = FALSE)
  expect_setequal(unobserved(latent_added), "L")
  expect_true(startsWith(as.character(latent_added), "dag named"))

  connected <- connect_nodes(dag, nodes = c("Z", "Y"), type = "full",
                             print_edges = FALSE)
  expect_true(startsWith(as.character(connected), sub("\\}\\s*$", "", original)))
  new_keys <- setdiff(edge_keys(connected), edge_keys(dag))
  expect_true(length(new_keys) > 0L)
  expect_true(all(grepl("<->", new_keys, fixed = TRUE)))
  expect_identical(attr(connected, "provenance"), list(author = "fixture"))
  expect_equal(dagitty::coordinates(connected), coordinates)
})


test_that("add_edges() warns instead of stopping on a graph that already has a cycle", {
  cyclic <- dagitty::dagitty("dag { A -> B B -> A C }")
  expect_warning(
    edited <- add_edges(cyclic, "A -> C"),
    "The outputted graph contains cycles", fixed = TRUE
  )
  expect_true("A -> C" %in% edge_keys(edited))
  expect_false(dagitty::isAcyclic(edited))

  acyclic <- string_to_dag(c("A -> B", "B -> C"))
  expect_error(
    add_edges(acyclic, "C -> A"),
    "add_edges(): the resulting graph would contain a directed cycle.",
    fixed = TRUE
  )
})


test_that("edge editors and inspection reject non-dagitty graph input", {
  dag <- string_to_dag("A -> B")
  expect_error(join_graphs("dag { A -> B }", dag),
               "join_graphs(): 'dag' must be a dagitty object.", fixed = TRUE)
  expect_error(join_graphs(dag),
               "join_graphs(): 'new_dag' must be a dagitty object.", fixed = TRUE)
  expect_error(join_graphs(dag, NULL),
               "join_graphs(): 'new_dag' must be a dagitty object.", fixed = TRUE)
  expect_error(keep_edges("A -> B", "dag { A -> B }"),
               "keep_edges(): 'dag' must be a dagitty object.", fixed = TRUE)
  expect_error(assess_edges(list()),
               "assess_edges(): 'dag' must be a dagitty object.", fixed = TRUE)

  shared <- paste0("Invalid input detected: only dagitty objects may be ",
                   "supplied as 'dag' input. Please check the input and try again.")
  expect_error(get_structure(NULL), shared, fixed = TRUE)
  expect_error(get_node_edges(NULL), shared, fixed = TRUE)
  expect_error(get_paths(character(0)), shared, fixed = TRUE)
})


test_that("remove_edges() and keep_edges() drop character NA entries with a warning", {
  dag <- string_to_dag(c("A -> B", "B -> C"))
  expect_warning(
    expect_warning(removed <- remove_edges(dag, c(NA, "A -> B")),
                   "no remaining edges.*A"),
    "remove_edges(): rows in 'edges' containing missing values were ignored.", fixed = TRUE
  )
  expect_identical(edge_keys(removed), "B -> C")

  expect_warning(
    kept <- keep_edges(c(NA, "A -> B"), dag),
    "keep_edges(): rows in 'edges' containing missing values were ignored.", fixed = TRUE
  )
  expect_identical(edge_keys(kept), "A -> B")

  expect_warning(
    expect_error(keep_edges(NA_character_, dag),
                 "keep_edges(): no complete edges were supplied.", fixed = TRUE),
    "containing missing values"
  )
  expect_warning(
    expect_error(remove_edges(dag, NA_character_), "must contain at least one edge"),
    "missing values were ignored"
  )
})
