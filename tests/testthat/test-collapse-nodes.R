test_that("collapse_nodes projects chains, common causes and collider paths", {
  cases <- list(
    list("dag { T -> M -> Y }", "M", "T -> Y"),
    list("dag { U -> A; U -> B }", "U", "A <-> B"),
    list("dag { T -> M; M -> A; M -> B }", "M",
         c("T -> A", "T -> B", "A <-> B")),
    list("dag { A -> C; B -> C }", "C", character()),
    list("dag { A <-> M; M -> B }", "M", "A <-> B"),
    list("dag { A -> M -> B; U -> A; U -> B }", c("M", "U"),
         c("A -> B", "A <-> B")),
    list("dag { T -> M1 -> M2 -> Y }", c("M1", "M2"), "T -> Y"),
    list("dag { T -> M1 -> Y; T -> M2 -> Y }", c("M1", "M2"), "T -> Y"),
    list("dag { T -> M1 -> K -> M2 -> Y }", c("M1", "M2"),
         c("T -> K", "K -> Y")),
    list("dag { U -> A; U -> B; U -> C }", "U",
         c("A <-> B", "A <-> C", "B <-> C")),
    list("dag { U0 -> U1 -> A; U0 -> U2 -> B }", c("U0", "U1", "U2"),
         "A <-> B"),
    list("dag { X -> U; U -> A; U -> B }", "U",
         c("X -> A", "X -> B", "A <-> B")),
    list("dag { A -> C; B -> C; C -> Y }", "C", c("A -> Y", "B -> Y")),
    list("dag { A <-> M; M <-> B }", "M", character()),
    list("dag { A -> M; M <-> B }", "M", character()),
    list("dag { A -> L; I }", "L", character()),
    list("dag { M -> A; M <-> B }", "M", "A <-> B")
  )
  for( case in cases ){
    dag <- dagitty::dagitty(case[[1L]])
    original <- serialize(dag, NULL)
    projected <- collapse_nodes(dag, case[[2L]])

    expect_identical(edge_keys(projected), sort(case[[3L]]))
    expect_setequal(names(projected), setdiff(names(dag), case[[2L]]))
    expect_identical(causaliflower:::dag_graph_type(projected), "dag")
    expect_true(dagitty::isAcyclic(projected))
    expect_identical(serialize(dag, NULL), original)
  }
})


test_that("collapse_nodes defaults to quiet output and prints only induced edges", {
  dag <- dagitty::dagitty(
    "dag { T -> M -> Y; U -> T; U -> Y; T -> Y }"
  )
  expect_identical(names(formals(collapse_nodes)), c("dag", "nodes", "print_edges"))
  expect_identical(formals(collapse_nodes)$print_edges, FALSE)
  expect_silent(default <- collapse_nodes(dag, c("M", "U")))
  expect_silent(quiet <- collapse_nodes(dag, c("M", "U"), print_edges = FALSE))
  output <- capture.output(
    expect_message(
      printed <- collapse_nodes(dag, c("M", "U"), print_edges = TRUE),
      "Printed new edges", fixed = TRUE
    )
  )
  expect_identical(default, quiet)
  expect_identical(printed, quiet)
  expect_match(paste(output, collapse = "\n"), "'T <-> Y'", fixed = TRUE)
  expect_false(grepl("'T -> Y'", paste(output, collapse = "\n"), fixed = TRUE))
  expect_silent(collapse_nodes(dag, "M", print_edges = TRUE))
})


test_that("collapse_nodes validates graph and node inputs before editing", {
  dag <- dagitty::dagitty("dag { T -> M -> Y; I }")
  original <- serialize(dag, NULL)

  expect_error(collapse_nodes("dag { T -> Y }", "T"), "only dagitty objects")
  expect_error(collapse_nodes(dag), "argument.*missing")
  for( nodes in list(NULL, list("M"), 1, factor("M")) ){
    expect_error(collapse_nodes(dag, nodes), "'nodes' must be a character vector",
                 fixed = TRUE)
  }
  expect_error(collapse_nodes(dag, character()), "at least one node name")
  expect_error(collapse_nodes(dag, c("M", NA_character_)),
               "non-missing, non-empty names")
  expect_error(collapse_nodes(dag, c("M", "")), "non-missing, non-empty names")
  for( value in list(NULL, NA, c(TRUE, FALSE), 1, "TRUE") ){
    expect_error(collapse_nodes(dag, "M", print_edges = value),
                 "'print_edges' must be TRUE or FALSE", fixed = TRUE)
  }
  for( spec in c("pdag { T -> M -> Y }", "mag { T -> M -> Y }") ){
    expect_error(collapse_nodes(dagitty::dagitty(spec), "M"),
                 "native graph type 'dag'", fixed = TRUE)
  }
  expect_error(collapse_nodes(dagitty::dagitty("dag { A -- B }"), "A"),
               "only directed and bidirected edges")
  expect_error(collapse_nodes(dagitty::dagitty("dag { A -> B -> A }"), "A"),
               "must be acyclic")
  self_loop <- dagitty::dagitty("dag { A <-> A; I }")
  expect_error(collapse_nodes(self_loop, "I"), "self-loops are not supported")
  expect_identical(edge_keys(self_loop), "A <-> A")
  expect_identical(serialize(dag, NULL), original)
})


test_that("collapse_nodes rejects conditional statuses in compact native text", {
  for( status in c("adjusted", "selected") ){
    dag <- structure(paste0("dag {\nU [latent] S [", status,
                             "]\nU -> T\n}"), class = "dagitty")
    original <- serialize(dag, NULL)
    expect_error(collapse_nodes(dag, "U"), "does not perform conditional projection")
    expect_error(collapse_nodes(dag, "S"), "does not perform conditional projection")
    expect_identical(serialize(dag, NULL), original)
  }
  ordinary <- dagitty::dagitty(paste0(
    'dag { "a [selected]" [note="selected"] T -> M -> Y }'
  ))
  expect_no_warning(projected <- collapse_nodes(ordinary, "M"))
  expect_true("a [selected]" %in% names(projected))
})


test_that("collapse_nodes ignores duplicate and absent names", {
  dag <- structure('dag named { bb="0,0,3,2" T -> M -> Y; I }',
                   class = "dagitty", provenance = "fixture")
  original <- serialize(dag, NULL)
  expect_silent(duplicates <- collapse_nodes(dag, c("M", "M")))
  expect_identical(duplicates, collapse_nodes(dag, "M"))
  expect_warning(
    mixed <- collapse_nodes(dag, c("M", "absent", "absent")),
    "collapse_nodes(): nodes not present in the graph were ignored: absent.",
    fixed = TRUE
  )
  expect_identical(mixed, duplicates)
  expect_warning(
    unchanged <- collapse_nodes(dag, "absent"),
    "collapse_nodes(): nodes not present in the graph were ignored: absent.",
    fixed = TRUE
  )
  expect_identical(serialize(unchanged, NULL), original)
  expect_identical(serialize(dag, NULL), original)
})


test_that("collapse_nodes removes only selected nodes and their native statuses", {
  dag <- dagitty::dagitty(paste0(
    "dag { T [exposure] Y [outcome] U [latent] V [latent] I ",
    "U -> T; U -> M; T -> M -> Y; V -> Y }"
  ))
  collapsed <- collapse_nodes(dag, c("M", "U"))
  expect_setequal(names(collapsed), c("T", "Y", "V", "I"))
  expect_identical(treatments(collapsed), "T")
  expect_identical(outcomes(collapsed), "Y")
  expect_identical(unobserved(collapsed), "V")
  expect_setequal(edge_keys(collapsed), c("T -> Y", "T <-> Y", "V -> Y"))

  no_treatment <- collapse_nodes(dag, "T")
  expect_length(treatments(no_treatment), 0L)
  expect_identical(outcomes(no_treatment), "Y")
  expect_setequal(unobserved(no_treatment), c("U", "V"))
  no_outcome <- collapse_nodes(dag, "Y")
  expect_length(outcomes(no_outcome), 0L)
  expect_identical(treatments(no_outcome), "T")
  expect_true("V" %in% names(no_outcome))

  latent_endpoint <- dagitty::dagitty("dag { A [latent] A <-> M; M -> B }")
  projected <- collapse_nodes(latent_endpoint, "M")
  expect_identical(edge_keys(projected), "A <-> B")
  expect_identical(unobserved(projected), "A")
})


test_that("collapse_nodes preserves survivor metadata and existing curved edges", {
  dag <- structure(paste0(
    'dag named {\n bb="-2,-2,5,3"\n',
    'T [exposure,pos="-1.123456789012345,0",label="Treatment",note="keep node"]\n',
    'M [pos="0,1",note="remove node"]\n',
    'Y [outcome,pos="2.987654321098765,0"]\n',
    'B [pos="3,2"] U [latent,pos="4,2"] I [pos="4,3"]\n',
    'T -> M [pos="0,0.123456789012345",note="remove edge"]\n',
    'M -> Y [pos="1.2,1.7"] M -> B\n',
    'T -> Y [pos="1,2",label="keep direct"]\n',
    'T <-> Y [pos="-0.0004,-1",note="keep confounding"]\n}'),
    class = "dagitty", labels = c(T = "Treatment", M = "Mediator", Y = "Outcome"),
    provenance = list(author = "fixture"), custom = list(nodes = c("M", "Y")))
  original <- serialize(dag, NULL)
  coordinates <- dagitty::coordinates(dag)
  controls <- causaliflower:::edge_controls(dag)
  retained <- controls[controls$v != "M" & controls$w != "M", , drop = FALSE]

  projected <- collapse_nodes(dag, "M")
  expect_setequal(edge_keys(projected), c("T -> Y", "T <-> Y", "T -> B", "B <-> Y"))
  expect_setequal(names(projected), setdiff(names(dag), "M"))
  actual_coordinates <- dagitty::coordinates(projected)
  expect_identical(actual_coordinates$x, coordinates$x[names(projected)])
  expect_identical(actual_coordinates$y, coordinates$y[names(projected)])
  actual_controls <- causaliflower:::edge_controls(projected)
  expect_identical(actual_controls$v, retained$v)
  expect_identical(actual_controls$e, retained$e)
  expect_identical(actual_controls$w, retained$w)
  expect_identical(as.numeric(actual_controls$x), as.numeric(retained$x))
  expect_identical(as.numeric(actual_controls$y), as.numeric(retained$y))
  expect_identical(attr(projected, "labels"), c(T = "Treatment", Y = "Outcome"))
  expect_identical(attr(projected, "provenance"), attr(dag, "provenance"))
  expect_identical(attr(projected, "custom"), attr(dag, "custom"))
  expect_identical(treatments(projected), "T")
  expect_identical(outcomes(projected), "Y")
  expect_identical(unobserved(projected), "U")
  for( piece in c('dag named', 'bb="-2,-2,5,3"', 'label="Treatment"',
                  'note="keep node"', 'label="keep direct"',
                  'note="keep confounding"') ){
    expect_match(as.character(projected), piece, fixed = TRUE)
  }
  expect_false(grepl('note="remove node"', as.character(projected), fixed = TRUE))
  expect_false(grepl('note="remove edge"', as.character(projected), fixed = TRUE))
  expect_identical(serialize(dag, NULL), original)
})


test_that("collapse_nodes keeps isolated nodes and can return an empty graph", {
  dag <- structure('dag named { bb="0,0,3,2" A -> L; I [latent] }',
                   class = "dagitty", provenance = "fixture")
  projected <- collapse_nodes(dag, "L")
  expect_setequal(names(projected), c("A", "I"))
  expect_length(edge_keys(projected), 0L)
  expect_identical(unobserved(projected), "I")

  empty <- collapse_nodes(dag, names(dag))
  expect_length(names(empty), 0L)
  expect_length(edge_keys(empty), 0L)
  expect_identical(causaliflower:::dag_graph_type(empty), "dag")
  expect_match(as.character(empty), "dag named", fixed = TRUE)
  expect_match(as.character(empty), 'bb="0,0,3,2"', fixed = TRUE)
  expect_identical(attr(empty, "provenance"), "fixture")

  one <- dagitty::dagitty("dag { U [latent] }")
  expect_length(names(collapse_nodes(one, "U")), 0L)
})


test_that("collapse_nodes handles quoted names and grouped native statements", {
  dag <- dagitty::dagitty(
    'dag { "a->b" -> "soil moisture" -> "a;b"; "a [selected]" }'
  )
  projected <- collapse_nodes(dag, "soil moisture")
  expect_setequal(names(projected), c("a->b", "a;b", "a [selected]"))
  expect_identical(edge_keys(projected), "a->b -> a;b")

  grouped <- structure("dag { A -> { M X }; M -> Y; I }", class = "dagitty")
  projected <- collapse_nodes(grouped, "M")
  expect_setequal(names(projected), c("A", "X", "Y", "I"))
  expect_setequal(edge_keys(projected), c("A -> X", "A -> Y"))
})


test_that("collapse_nodes agrees with hidden-parent DAG projection and separation", {
  set.seed(20260909)
  vertices <- c("A", "B", "C", "D", "E")
  pairs <- t(utils::combn(vertices, 2L))
  for( iteration in seq_len(24L) ){
    directed <- pairs[stats::runif(nrow(pairs)) < 0.3, , drop = FALSE]
    bidirected <- pairs[stats::runif(nrow(pairs)) < 0.3, , drop = FALSE]
    directed_text <- if( nrow(directed) > 0L ){
      paste(directed[, 1L], "->", directed[, 2L])
    }else{
      character()
    }
    bidirected_text <- if( nrow(bidirected) > 0L ){
      paste(bidirected[, 1L], "<->", bidirected[, 2L])
    }else{
      character()
    }
    dag <- dagitty::dagitty(paste("dag {", paste(c(vertices, directed_text,
                                                   bidirected_text), collapse = "; "), "}"))
    omitted <- sample(vertices, 2L)
    retained <- setdiff(vertices, omitted)
    hidden <- if( nrow(bidirected) > 0L ){
      paste0("hidden_", seq_len(nrow(bidirected)))
    }else{
      character()
    }
    expanded_text <- c(vertices, directed_text, paste(omitted, "[latent]"))
    if( length(hidden) > 0L ){
      expanded_text <- c(expanded_text, paste(hidden, "[latent]"),
                         paste(hidden, "->", bidirected[, 1L]),
                         paste(hidden, "->", bidirected[, 2L]))
    }
    expanded <- dagitty::dagitty(paste("dag {", paste(expanded_text, collapse = "; "), "}"))
    expected <- causaliflower:::id_as_dagitty(causaliflower:::id_admg(expanded))
    projected <- collapse_nodes(dag, omitted)
    successive <- collapse_nodes(collapse_nodes(dag, omitted[[1L]]), omitted[[2L]])
    reversed <- collapse_nodes(collapse_nodes(dag, omitted[[2L]]), omitted[[1L]])

    expect_identical(edge_keys(projected), edge_keys(expected))
    expect_setequal(names(projected), retained)
    expect_identical(edge_keys(successive), edge_keys(projected))
    expect_identical(edge_keys(reversed), edge_keys(projected))
    expect_identical(edge_keys(collapse_nodes(dag, rev(omitted))), edge_keys(projected))
    retained_pairs <- utils::combn(retained, 2L, simplify = FALSE)
    for( pair in retained_pairs ){
      other <- setdiff(retained, pair)
      for( given in list(character(), other) ){
        expect_identical(
          dagitty::dseparated(projected, pair[[1L]], pair[[2L]], given),
          dagitty::dseparated(dag, pair[[1L]], pair[[2L]], given)
        )
      }
    }
  }
})
