test_that("recombine_nodes() merges all original groups across tiers", {
  dag <- dagitty::dagitty(paste0(
    "dag { T [exposure]; Y [outcome]; ",
    "T -> A1; T -> A2; A1 -> B1; A1 -> B2; ",
    "A2 -> B1; A2 -> B2; B1 -> Y; B2 -> Y; I }"
  ))
  before <- dag
  groups <- recombine_nodes(dag, preview = TRUE)
  expect_identical(groups, list(A1_A2 = c("A1", "A2"),
                                B1_B2 = c("B1", "B2")))
  out <- recombine_nodes(dag)
  expected <- dagitty::dagitty(paste0(
    "dag { T [exposure]; Y [outcome]; T -> A1_A2 -> B1_B2 -> Y; I }"
  ))
  expect_setequal(names(out), names(expected))
  expect_setequal(edge_key(dag_edges(out)), edge_key(dag_edges(expected)))
  expect_setequal(dagitty::exposures(out), "T")
  expect_setequal(dagitty::outcomes(out), "Y")
  expect_identical(dag, before)

  manual <- dag
  for( name in names(groups) ){
    manual <- merge_nodes(manual, groups[[name]], name)
  }
  expect_identical(out, manual)
})

test_that("recombine_nodes() uses one pass and preserves bows between groups", {
  dag <- dagitty::dagitty("dag { A <-> B; B <-> C }")
  expect_identical(recombine_nodes(dag, preview = TRUE),
                   list(A_C = c("A", "C")))
  out <- recombine_nodes(dag)
  expected <- dagitty::dagitty("dag { A_C <-> B }")
  expect_setequal(names(out), names(expected))
  expect_setequal(edge_key(dag_edges(out)), edge_key(dag_edges(expected)))
  expect_identical(recombine_nodes(out, preview = TRUE),
                   list(A_C_B = c("A_C", "B")))

  dag <- dagitty::dagitty(paste0(
    "dag { A1 -> B1; A1 -> B2; A2 -> B1; A2 -> B2; ",
    "A1 <-> B1; A1 <-> B2; A2 <-> B1; A2 <-> B2 }"
  ))
  expect_identical(recombine_nodes(dag, preview = TRUE),
                   list(A1_A2 = c("A1", "A2"), B1_B2 = c("B1", "B2")))
  out <- recombine_nodes(dag)
  expected <- dagitty::dagitty("dag { A1_A2 -> B1_B2; A1_A2 <-> B1_B2 }")
  expect_setequal(edge_key(dag_edges(out)), edge_key(dag_edges(expected)))
})

test_that("recombine_nodes() keeps full graph context for scoped groups", {
  dag <- dagitty::dagitty(paste0(
    "dag { T [exposure]; Y [outcome]; T -> A -> Y; ",
    "T -> B -> Y; T -> C -> Y; X -> D -> Y; I }"
  ))
  groups <- recombine_nodes(dag, c("C", "A", "C"), preview = TRUE)
  expect_identical(groups, list(A_C = c("A", "C")))
  out <- recombine_nodes(dag, c("C", "A", "C"))
  expect_setequal(names(out), c("T", "Y", "A_C", "B", "X", "D", "I"))
  expect_identical(recombine_nodes(dag, c("A", "D"), TRUE), list())
  expect_identical(recombine_nodes(dag, c("T", "Y", "I"), TRUE), list())
  expect_identical(recombine_nodes(dag, "A"), dag)
})

test_that("recombine_nodes() reserves original names in radix order", {
  dag <- dagitty::dagitty(paste0(
    "dag { T [exposure]; Y [outcome]; T -> b -> Y; ",
    "T -> A -> Y; T -> a -> Y; A_a_b; A_a_b_2 }"
  ))
  groups <- recombine_nodes(dag, c("b", "a", "A"), preview = TRUE)
  expect_identical(groups, list(A_a_b_3 = c("A", "a", "b")))
  out <- recombine_nodes(dag, c("b", "a", "A"))
  expect_setequal(names(out), c("T", "Y", "A_a_b_3", "A_a_b", "A_a_b_2"))
  expected <- dagitty::dagitty(paste0(
    "dag { T [exposure]; Y [outcome]; T -> A_a_b_3 -> Y; A_a_b; A_a_b_2 }"
  ))
  expect_setequal(edge_key(dag_edges(out)), edge_key(dag_edges(expected)))
})

test_that("recombine_nodes() handles underscores in generated names", {
  dag <- dagitty::dagitty(paste0(
    "dag { T [exposure]; Y [outcome]; Z [outcome]; ",
    "T -> A -> Y; T -> B_C -> Y; T -> A_B -> Z; T -> C -> Z }"
  ))
  groups <- recombine_nodes(dag, preview = TRUE)
  expect_identical(groups, list(A_B_C = c("A", "B_C"),
                                A_B_C_2 = c("A_B", "C")))
  out <- recombine_nodes(dag)
  expected <- dagitty::dagitty(paste0(
    "dag { T [exposure]; Y [outcome]; Z [outcome]; ",
    "T -> A_B_C -> Y; T -> A_B_C_2 -> Z }"
  ))
  expect_setequal(names(out), names(expected))
  expect_setequal(edge_key(dag_edges(out)), edge_key(dag_edges(expected)))

  occupied <- dagitty::dagitty(paste0(
    "dag { T [exposure]; Y [outcome]; Z [outcome]; ",
    "T -> A -> Y; T -> B_C -> Y; T -> A_B -> Z; T -> C -> Z; ",
    "A_B_C; A_B_C_2 }"
  ))
  expect_identical(recombine_nodes(occupied, preview = TRUE),
                   list(A_B_C_3 = c("A", "B_C"), A_B_C_4 = c("A_B", "C")))
})

test_that("recombine_nodes() preserves first node metadata and curves", {
  dag <- structure(paste0(
    'dag named { bb="0,0,9,9"\n',
    'T [exposure,pos="0,0"]\nY [outcome,pos="9,0"]\n',
    'A1 [pos="2,0",label="alpha",note="first A"]\n',
    'A2 [pos="2,2",label="discard A",note="second A"]\n',
    'B1 [pos="5,0",label="beta",note="first B"]\n',
    'B2 [pos="5,2",label="discard B",note="second B"]\n',
    'I [pos="8,3",note="untouched"]\n',
    'T -> A1 [pos="1,1",note="incoming curve"]\nT -> A2\n',
    'A1 -> B1 [pos="3,1",note="between curve"]\nA1 -> B2\n',
    'A2 -> B1\nA2 -> B2\n',
    'B1 -> Y [pos="7,1",note="outgoing curve"]\nB2 -> Y\n',
    'I -> Y [note="untouched edge"]\n}'
  ), class = "dagitty")
  attr(dag, "labels") <- c(A1 = "alpha", A2 = "other A", B1 = "beta",
                            B2 = "other B", I = "isolate")
  attr(dag, "custom") <- list(nodes = c("A1", "B2"), data = 3L)
  before <- dag
  groups <- recombine_nodes(dag, c("B2", "A2", "B1", "A1"), TRUE)
  expect_identical(groups, list(A1_A2 = c("A1", "A2"),
                                B1_B2 = c("B1", "B2")))
  expect_identical(dag, before)
  out <- recombine_nodes(dag, c("B2", "A2", "B1", "A1"))
  text <- as.character(out)
  for( unchanged in c('dag named { bb="0,0,9,9"',
                      'T [exposure,pos="0,0"]', 'Y [outcome,pos="9,0"]',
                      'I [pos="8,3",note="untouched"]',
                      'I -> Y [note="untouched edge"]') ){
    expect_match(text, unchanged, fixed = TRUE)
  }
  expect_match(text, '"A1_A2" [pos="2,0",label="alpha",note="first A"]',
               fixed = TRUE)
  expect_match(text, '"B1_B2" [pos="5,0",label="beta",note="first B"]',
               fixed = TRUE)
  expect_match(text, 'T -> "A1_A2" [pos="1,1",note="incoming curve"]',
               fixed = TRUE)
  expect_match(text, '"A1_A2" -> "B1_B2" [pos="3,1",note="between curve"]',
               fixed = TRUE)
  expect_match(text, '"B1_B2" -> Y [pos="7,1",note="outgoing curve"]',
               fixed = TRUE)
  expect_false(grepl("discard|second", text))
  expect_equal(nrow(edge_controls(out, "test")), 3L)
  expect_identical(attr(out, "labels"), c(A1_A2 = "alpha", B1_B2 = "beta",
                                         I = "isolate"))
  expect_identical(attr(out, "custom"), attr(dag, "custom"))
  expect_identical(names(attributes(out)), names(attributes(dag)))
  expect_identical(dag, before)
})

test_that("recombine_nodes() executes protected and latent group rules", {
  dag <- dagitty::dagitty(paste0(
    "dag { T1 [exposure]; T2 [exposure]; Y1 [outcome]; Y2 [outcome]; ",
    "S1 [selected]; S2 [selected]; X1 [adjusted]; X2 [adjusted]; ",
    "U1 [latent]; U2 [latent]; I1; I2; ",
    "T1 -> K; T2 -> K; Y1 -> K; Y2 -> K; S1 -> K; S2 -> K; ",
    "X1 -> K; X2 -> K; U1 -> K; U2 -> K; A -> K; B -> K }"
  ))
  expect_identical(recombine_nodes(dag, preview = TRUE),
                   list(A_B = c("A", "B"), U1_U2 = c("U1", "U2")))
  out <- recombine_nodes(dag)
  expect_setequal(names(out), c("T1", "T2", "Y1", "Y2", "S1", "S2",
                               "X1", "X2", "U1_U2", "I1", "I2", "A_B", "K"))
  parsed <- dagitty::dagitty(as.character(out))
  expect_setequal(dagitty::exposures(parsed), c("T1", "T2"))
  expect_setequal(dagitty::outcomes(parsed), c("Y1", "Y2"))
  expect_setequal(selected_nodes(parsed), c("S1", "S2"))
  expect_setequal(dagitty::adjustedNodes(parsed), c("X1", "X2"))
  expect_setequal(dagitty::latents(parsed), "U1_U2")
})

test_that("recombine_nodes() leaves graphs with no matches unchanged", {
  expect_identical(names(formals(recombine_nodes)),
                   c("dag", "nodes", "preview"))
  expect_identical(formals(recombine_nodes)$preview, FALSE)
  dag <- dagitty::dagitty('dag { T -> "M 1" -> Y; T -> "M[2]" -> Y }')
  expect_silent(groups <- recombine_nodes(dag, preview = TRUE))
  expect_identical(groups, list(`M 1_M[2]` = c("M 1", "M[2]")))
  expect_silent(out <- recombine_nodes(dag))
  expect_identical(out, merge_nodes(dag, groups[[1L]], names(groups)[1L]))
  expect_identical(attributes(out), attributes(dag))
  expect_error(recombine_nodes(dag, print_edges = TRUE), "unused argument")

  for( text in c("dag {}", "dag named { A; B }", "dag { A -> B -> C }") ){
    dag <- dagitty::dagitty(text)
    attr(dag, "custom") <- list(value = 1L)
    expect_silent(groups <- recombine_nodes(dag, preview = TRUE))
    expect_identical(groups, list())
    expect_silent(out <- recombine_nodes(dag))
    expect_identical(out, dag)
  }
})

test_that("recombine_nodes() validates both modes without edits", {
  dag <- dagitty::dagitty("dag { T -> A -> Y; T -> B -> Y }")
  before <- dag
  for( preview in c(FALSE, TRUE) ){
    expect_error(recombine_nodes(dag, c("A", "absent"), preview), "not present")
    for( nodes in list(character(), NA_character_, "", 1L) ){
      expect_error(recombine_nodes(dag, nodes, preview), "node names")
    }
    expect_identical(dag, before)
    for( text in c("pdag { A -> B }", "mag { A -> B }",
                   "dag { A -- B }", "dag { A -> B -> A }",
                   "dag { A; B; I <-> I }") ){
      invalid <- dagitty::dagitty(text)
      original <- invalid
      expect_error(recombine_nodes(invalid, preview = preview),
                   "dag graph type|supported|acyclic")
      expect_identical(invalid, original)
    }
  }
  for( preview in list(NA, NULL, "yes", 1L, c(TRUE, FALSE)) ){
    expect_error(recombine_nodes(dag, preview = preview), "TRUE or FALSE")
    expect_identical(dag, before)
  }
})

test_that("recombine_nodes() validates generated names in both modes", {
  dag <- dagitty::dagitty('dag { T -> "A\\B" -> Y; T -> C -> Y }')
  before <- dag
  for( preview in c(FALSE, TRUE) ){
    expect_error(recombine_nodes(dag, preview = preview),
                 "cannot contain quotes, backslashes, or line breaks",
                 fixed = TRUE)
    expect_identical(dag, before)
  }
})
