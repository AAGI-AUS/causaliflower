# ==============================================================================
# identification.R


# get_identification() and the identification methods it calls
# ==============================================================================


#' Check if a causal effect is identifiable
#'
#' `get_identification()` checks whether a causal effect is identifiable under
#' the supplied `dagitty` DAG. Marginal effects are checked with the ID
#' algorithm. When `conditioning` contains node names, the conditional effect
#' is checked with IDC. When `experiments` contains node names,
#' z-identification is also checked. For marginal total effects, covariate
#' adjustment, front door, instrumental variable, and applicable surrogate
#' experiment strategies are reported separately from the overall verdict.
#'
#' With `effect = "direct"`, the function first uses latent projection to
#' obtain a graph containing only observed nodes. ID then checks the outcome
#' distribution when the treatment and the outcome's other parents are set by
#' intervention. The controlled direct effect is defined for this projected
#' graph. The other parents are printed as the `fixing set`, which can include
#' confounders and outcome predictors as well as mediators. They are set by
#' intervention, not conditioned on. The controlled direct effect can depend
#' on the values chosen for the fixing set. The function checks identification
#' without choosing those values or estimating the effect.
#'
#' Explicitly declared latent nodes are not held fixed. Projection can replace
#' paths through latent mediators with direct arrows, so the reported query can
#' include effects through latent mediators and need not be the controlled
#' direct effect in the supplied DAG.
#'
#' Direct queries also report a separate mediation check for natural direct
#' and indirect effects relative to selected mediators. This looks for a common
#' observed adjustment set for the treatment's effect on the mediators and the
#' joint effect of treatment and mediators on the outcome. Only nondescendants
#' of the treatment are considered as covariates. The search first tries the
#' empty set, then a single set containing all eligible observed ancestors of
#' the treatment, mediators and outcome. Each candidate must pass both checks.
#' If neither candidate passes, no covariate set of observed treatment
#' nondescendants passes both checks.
#' Other methods may still identify the natural direct and indirect effects.
#'
#' A successful mediation check requires a structural causal model with
#' independent background errors in an underlying DAG that includes the latent
#' common causes represented by bidirected edges. This entails cross world
#' assumptions. Consistency and positivity are also required. The check does
#' not verify these assumptions or estimate an effect. With several selected
#' mediators, the check concerns their joint natural effects, not an
#' effect decomposition for each mediator. The natural direct effect can still
#' include paths through mediators outside the selected set. The joint check
#' can fail even when a smaller mediator set defines the same natural effects
#' and passes the check. This can happen when the treatment affects one
#' selected mediator only through other selected mediators.
#'
#' @param dag A `dagitty` object representing a causal DAG. Bidirected edges
#'   are supported between observed nodes. All edges connected to explicitly
#'   declared latent nodes must be directed.
#' @param treatment A single character string naming an observed treatment
#'   node. If left `NULL`, the `dagitty` object must contain exactly one
#'   treatment marked with the `[exposure]` tag.
#' @param outcome A single character string naming an observed outcome node.
#'   If left `NULL`, the `dagitty` object must contain exactly one outcome
#'   marked with the `[outcome]` tag.
#' @param effect The effect to check: `"total"` (the default) or `"direct"`.
#' @param experiments An optional character vector naming observed nodes for
#'   which surrogate experimental data are available. It cannot include the
#'   treatment or outcome. When supplied, use `effect = "total"` and
#'   `conditioning = NULL`.
#' @param conditioning An optional character vector naming observed nodes to
#'   condition on in \eqn{P(y \mid do(x), z)}. It cannot include the treatment
#'   or outcome. When supplied, use `effect = "total"` and
#'   `experiments = NULL`.
#' @param quiet If `TRUE`, suppress the summary message.
#' @param mediator An optional nonempty character vector of observed mediators
#'   on directed paths from treatment to outcome. Requires `effect = "direct"`.
#'   Supply several mediator names to assess them jointly. With `NULL`, each
#'   observed mediator is checked separately, followed by all observed
#'   mediators together if there is more than one.
#'
#' @returns A `data.table` with columns `strategy`, `feasible`, `tier`,
#'   `variable_set`, and `note`. `feasible` records the graphical check for
#'   that row. The overall verdict is stored separately. For `effect =
#'   "direct"`, the general identification row also requires a direct arrow
#'   in the graph obtained by latent projection. If there is no such arrow,
#'   `FALSE` means the row for the direct effect is not applicable, rather than
#'   proving that the intervention is unidentifiable. The verdict attributes
#'   and hedge still report identification of the distribution under the
#'   intervention on the treatment and fixing set.
#'
#'   Attributes `identifiable`, `observationally_identifiable`,
#'   `z_identifiable`, and `data_availability` record the verdict and available
#'   data. `treatment`, `outcome`, `effect`, `conditional`, and `conditioning`
#'   record the query used by the function. For conditional effects,
#'   `conditioning_intervened` lists the conditioning variables that IDC moves
#'   into the intervention set. `conditioning_remaining` lists those that
#'   remain conditioned on. `hedge` is `NULL` when ID succeeds. When ID fails,
#'   it contains two `dagitty` graphs, `F` and `F_prime`, their common roots,
#'   and the query at the step where identification failed.
#'   `roles` is a table with columns `variable`, `role`, and
#'   `conditioning_guidance`. A variable may have different roles in separate
#'   strategies. An instrument role does not exclude membership in a valid
#'   adjustment set. `caveat` distinguishes graphical identification from
#'   estimation using data.
#'
#'   For a marginal total effect, `adjustment_sets` contains all minimal
#'   adjustment sets and `canonical_adjustment` contains dagitty's canonical
#'   set. Both attributes are lists of character vectors: an empty list means
#'   no valid set, while `list(character())` is a valid empty set. These
#'   attributes are absent for direct or conditional queries. The adjustment
#'   row displays at most five minimal sets and counts any remaining sets.
#'
#'   Direct queries additionally have a `mediation` attribute: a list of entries
#'   containing `mediators`, `feasible`, `adjustment_sets`, `note`, and
#'   `assumptions`. A successful entry contains one qualifying adjustment set,
#'   not necessarily a minimal set. `list(character())` denotes no adjustment
#'   required. An unsuccessful entry has an empty list of adjustment sets. If
#'   there are no observed mediators, `mediation` is an empty list and no
#'   mediation rows are added. Mediation results are reported separately. The
#'   verdict attributes, hedge, roles and summary describe the intervention on
#'   the treatment and fixing set.
#'
#' @references
#' Shpitser I, VanderWeele TJ (2011). "A Complete Graphical Criterion for the
#' Adjustment Formula in Mediation Analysis." *The International Journal of
#' Biostatistics*, **7**(1). \doi{10.2202/1557-4679.1297}.
#'
#' Shpitser I, Pearl J (2006). "Identification of Conditional Interventional
#' Distributions." arXiv:1206.6876.
#'
#' Shpitser I, Pearl J (2008). "Complete Identification Methods for the Causal
#' Hierarchy." *Journal of Machine Learning Research*, **9**, 1941--1979.
#'
#' Bareinboim E, Pearl J (2016). "Causal inference and the data-fusion
#' problem." *Proceedings of the National Academy of Sciences*, **113**(27),
#' 7345--7352. \doi{10.1073/pnas.1510507113}.
#'
#' Pearl J (2009). *Causality: Models, Reasoning, and Inference*. Second
#' edition. Cambridge University Press.
#'
#' @examples
#' dag <- dagitty::dagitty("dag {
#'   rainfall -> irrigation
#'   rainfall -> yield
#'   irrigation -> yield
#' }")
#' get_identification(dag, treatment = "irrigation", outcome = "yield")
#'
#' direct_dag <- dagitty::dagitty("dag {
#'   treatment -> growth -> yield
#'   treatment -> yield
#' }")
#' get_identification(
#'   direct_dag,
#'   treatment = "treatment",
#'   outcome = "yield",
#'   effect = "direct"
#' )
#' mediation <- get_identification(
#'   direct_dag, treatment = "treatment", outcome = "yield",
#'   effect = "direct", mediator = "growth", quiet = TRUE
#' )
#' attr(mediation, "mediation")
#'
#' @export
get_identification <- function(dag, treatment = NULL, outcome = NULL,
                               effect = c("total", "direct"),
                               experiments = NULL, conditioning = NULL,
                               quiet = FALSE, mediator = NULL){
  graph <- id_validate_graph(dag)

  if( missing(effect) ){
    effect <- "total"
  }
  if( !is.character(effect) || length(effect) != 1L || is.na(effect) ||
      !effect %in% c("total", "direct") ){
    stop("Invalid input detected: only 'total' or 'direct' may be supplied as 'effect' input. Please check the input and try again.",
         call. = FALSE)
  }
  effect <- as.character(effect)
  if( !is.logical(quiet) || length(quiet) != 1L || is.na(quiet) ){
    stop("Invalid input detected: only TRUE or FALSE may be supplied as 'quiet' input.",
         call. = FALSE)
  }

  query <- id_resolve_query(graph, treatment, outcome)
  treatment <- query$treatment
  outcome <- query$outcome

  has_conditioning <- !is.null(conditioning)
  has_experiments <- !is.null(experiments)

  if( effect != "total" && has_conditioning ){
    stop("Unsupported input combination detected: 'conditioning' must be 'NULL' when 'effect = \"direct\"'. Please remove the supplied 'conditioning' input or use 'effect = \"total\"' and try again.",
         call. = FALSE)
  }
  if( effect != "total" && has_experiments ){
    stop("Unsupported input combination detected: 'experiments' must be 'NULL' when 'effect = \"direct\"'. Please remove the supplied 'experiments' input or use 'effect = \"total\"' and try again.",
         call. = FALSE)
  }
  if( has_conditioning && has_experiments ){
    stop("Unsupported input combination detected: 'conditioning' and 'experiments' cannot both be supplied when 'effect = \"total\"'. Please remove either the supplied 'conditioning' input or the supplied 'experiments' input and try again.",
         call. = FALSE)
  }
  if( !is.null(mediator) && effect != "direct" ){
    stop("Unsupported input combination detected: 'mediator' requires 'effect = \"direct\"'.",
         call. = FALSE)
  }

  conditioning <- id_validate_nodes(
    conditioning,
    graph$observed,
    c(treatment, outcome),
    "conditioning"
  )
  experiments <- id_validate_nodes(
    experiments,
    graph$observed,
    c(treatment, outcome),
    "experiments"
  )

  G <- id_admg(graph$dag)
  strategy_dag <- id_as_dagitty(G)
  direct_path <- id_directed_reachable(G, treatment, outcome)
  direct_arrow <- any(G$D$from == treatment & G$D$to == outcome)
  mediation <- list()

  rows <- list()
  adjustment <- list(feasible = FALSE, sets = list())
  frontdoor <- list(feasible = FALSE, sets = list())
  instruments <- list(feasible = FALSE, entries = list())
  idc <- list(intervened = character(), remaining = character())
  z_result <- list(zid = NA, set = character())

  if( effect == "direct" ){
    outcome_parents <- sort(unique(G$D$from[G$D$to == outcome]))
    fixing_set <- setdiff(outcome_parents, treatment)
    intervention <- sort(unique(c(treatment, fixing_set)))
    solved <- id_solve(outcome, intervention, G)

    general_note <- if( isTRUE(solved$id) ){
      paste(
        "ID identifies the outcome distribution under intervention on",
        "the treatment and fixing set in the latent-projected graph"
      )
    }else{
      paste(
        "ID does not identify the outcome distribution under intervention on",
        "the treatment and fixing set in the latent-projected graph"
      )
    }
    if( !direct_arrow ){
      general_note <- paste(
        "no direct arrow in the latent-projected graph;", general_note
      )
    }
    rows[[length(rows) + 1L]] <- id_result_row(
      "general identification",
      direct_arrow && isTRUE(solved$id),
      "nonparametric",
      paste0("fixing set: ", id_render_set(fixing_set)),
      general_note
    )

    observational <- isTRUE(solved$id)
    identifiable <- observational
    hedge <- solved$hedge
    mediation <- id_mediation(strategy_dag, G, treatment, outcome, mediator)
    for(entry in mediation){
      rows[[length(rows) + 1L]] <- id_result_row(
        "mediation (natural effects)",
        entry$feasible,
        "assumption-dependent",
        paste0(
          "mediators: ", id_render_set(entry$mediators),
          "; adjustment: ", id_render_sets(entry$adjustment_sets)
        ),
        paste(entry$note, entry$assumptions, sep = ". ")
      )
    }
  }else if( length(conditioning) > 0L ){
    idc <- idc_solve(outcome, treatment, conditioning, G)
    solved <- idc
    detail <- id_idc_detail(idc$intervened, idc$remaining)
    general_note <- if( isTRUE(solved$id) ){
      "the IDC algorithm identifies the conditional effect"
    }else{
      "the IDC algorithm does not identify the conditional effect"
    }
    rows[[length(rows) + 1L]] <- id_result_row(
      "conditional identification",
      solved$id,
      "nonparametric",
      "not applicable",
      paste(general_note, detail, sep = "; ")
    )

    observational <- isTRUE(solved$id)
    identifiable <- observational
    hedge <- solved$hedge
  }else{
    solved <- id_solve(outcome, treatment, G)
    general_note <- if( !direct_path ){
      "the DAG contains no directed path from the treatment to the outcome"
    }else if( isTRUE(solved$id) ){
      "the ID algorithm identifies the total effect"
    }else{
      "the ID algorithm does not identify the total effect"
    }
    rows[[length(rows) + 1L]] <- id_result_row(
      "general identification",
      solved$id,
      "nonparametric",
      "not applicable",
      general_note
    )

    adjustment <- id_adjustment(strategy_dag, treatment, outcome)
    rows[[length(rows) + 1L]] <- id_result_row(
      "covariate adjustment",
      adjustment$feasible,
      "nonparametric",
      id_render_sets(adjustment$sets, max_sets = 5L),
      if( adjustment$feasible ){
        "minimal adjustment sets for the total effect"
      }else{
        "no valid adjustment set for the total effect"
      }
    )

    frontdoor <- id_frontdoor(strategy_dag, G, treatment, outcome)
    rows[[length(rows) + 1L]] <- id_result_row(
      "front-door",
      frontdoor$feasible,
      "nonparametric",
      id_render_sets(frontdoor$sets),
      if( frontdoor$feasible ){
        "valid front-door mediator set"
      }else{
        "no valid front-door mediator set"
      }
    )

    instruments <- id_instruments(strategy_dag, treatment, outcome)
    rows[[length(rows) + 1L]] <- id_result_row(
      "instrumental variable",
      instruments$feasible,
      "assumption-dependent",
      id_render_instruments(instruments$entries),
      if( instruments$feasible ){
        "graphical instrument; causal interpretation requires additional assumptions"
      }else{
        "no graphical instrument"
      }
    )

    observational <- isTRUE(solved$id)
    hedge <- solved$hedge
    if( length(experiments) > 0L ){
      z_result <- id_zidentifiable(
        G,
        treatment,
        outcome,
        experiments,
        observational
      )
      rows[[length(rows) + 1L]] <- id_result_row(
        "surrogate experiment (z-identification)",
        z_result$zid,
        "nonparametric",
        if( z_result$zid ){
          id_render_set(if( observational ) { experiments } else { z_result$set })
        }else{
          "none"
        },
        if( observational ){
          "the total effect is already identifiable without the declared surrogate experimental data"
        }else if( z_result$zid ){
          "identifies the total effect using the declared surrogate experimental data"
        }else{
          "does not identify the total effect using the declared surrogate experimental data"
        }
      )
    }
    identifiable <- observational || isTRUE(z_result$zid)
  }

  result <- data.table::rbindlist(rows, use.names = TRUE)
  roles <- id_roles(
    G,
    treatment,
    outcome,
    effect,
    adjustment,
    frontdoor,
    instruments
  )

  attr(result, "identifiable") <- isTRUE(identifiable)
  attr(result, "observationally_identifiable") <- isTRUE(observational)
  attr(result, "z_identifiable") <- if( length(experiments) > 0L ){
    isTRUE(z_result$zid)
  }else{
    NA
  }
  attr(result, "data_availability") <- if( length(experiments) > 0L ){
    paste0(
      "observational data + surrogate experimental data for ",
      id_render_set(experiments)
    )
  }else{
    "observational data"
  }
  attr(result, "hedge") <- hedge
  attr(result, "roles") <- roles
  if( effect == "direct" ){
    attr(result, "mediation") <- mediation
  }
  if( effect == "total" && length(conditioning) == 0L ){
    attr(result, "adjustment_sets") <- adjustment$sets
    attr(result, "canonical_adjustment") <- adjustment$canonical
  }
  attr(result, "treatment") <- treatment
  attr(result, "outcome") <- outcome
  attr(result, "effect") <- effect
  attr(result, "conditional") <- length(conditioning) > 0L
  attr(result, "conditioning") <- conditioning
  attr(result, "conditioning_intervened") <- sort(idc$intervened)
  attr(result, "conditioning_remaining") <- sort(idc$remaining)
  attr(result, "caveat") <- paste(
    "Identification is a graphical result. Estimation also requires",
    "positivity and data that support the comparisons used by the reported strategy."
  )

  if( !quiet ){
    id_summary_message(
      treatment = treatment,
      outcome = outcome,
      effect = effect,
      conditioning = conditioning,
      experiments = if( isTRUE(z_result$zid) ) { z_result$set } else { experiments },
      direct_path = direct_path,
      direct_arrow = direct_arrow,
      observational = observational,
      identifiable = identifiable
    )
  }

  return(result)
}


## Mediation checks ##
# The two adjustment checks follow Shpitser and VanderWeele (2011), Theorem 2.
# Search failures are deliberately not general natural-effect ID verdicts.
#' Check natural effects for selected mediators
#'
#' @param dag The projected graph as a dagitty object.
#' @param G The same graph as the internal list with vertices `V` and edges
#'   `D` and `B`.
#' @param treatment,outcome Observed treatment and outcome node names.
#' @param mediator Mediator names to check jointly, or `NULL` to check each
#'   observed mediator and then their joint set.
#' @return A list of checks with mediator names, feasibility, adjustment sets,
#'   notes and assumptions.
#' @noRd
id_mediation <- function(dag, G, treatment, outcome, mediator = NULL){
  on_path <- sort(Filter(
    function(node){
      matches <- id_directed_reachable(G, treatment, node) &&
        id_directed_reachable(G, node, outcome)

      return(matches)
    },
    setdiff(G$V, c(treatment, outcome))
  ))
  if( !is.null(mediator) ){
    valid <- is.character(mediator) && length(mediator) > 0L &&
      !anyNA(mediator) && all(nzchar(mediator)) && all(mediator %in% on_path)
    if( !valid ){
      stop("Invalid mediator input detected: supply observed node names on directed paths from the treatment to the outcome, excluding the treatment and outcome themselves.",
           call. = FALSE)
    }
    targets <- list(sort(unique(as.character(mediator))))
  }else{
    targets <- lapply(on_path, function(node){
      return(node)
    })
    if( length(on_path) > 1L ){
      targets[[length(targets) + 1L]] <- on_path
    }
  }

  assumptions <- paste(
    "Requires a structural causal model with independent background errors in",
    "the underlying DAG (including latent common causes), cross-world assumptions,",
    "consistency and positivity."
  )
  entries <- lapply(targets, function(mediators){
    sets <- id_mediation_adjustment(dag, G, treatment, outcome, mediators)
    feasible <- length(sets) > 0L
    entry <- list(
      mediators = mediators,
      feasible = feasible,
      adjustment_sets = sets,
      note = if( feasible ){
        paste(
          "the standard mediation formula identifies natural direct and indirect",
          "effects relative to the selected mediator set under the additional assumptions;",
          "the reported covariates satisfy both adjustment checks"
        )
      }else{
        paste(
          "no covariate set of observed treatment nondescendants satisfies both",
          "adjustment checks; natural direct and indirect effects are not established by",
          "this mediation-formula search, but other identification methods may apply"
        )
      },
      assumptions = assumptions
    )

    return(entry)
  })

  return(unname(entries))
}


#' Find covariates for the two mediation adjustment checks
#'
#' @param dag The projected graph as a dagitty object.
#' @param G The same graph as the internal list with vertices `V` and edges
#'   `D` and `B`.
#' @param treatment,outcome Observed treatment and outcome node names.
#' @param mediators Mediator node names to assess jointly.
#' @return A list containing one qualifying adjustment set, or an empty list
#'   if no set of eligible covariates passes both checks.
#' @noRd
id_mediation_adjustment <- function(dag, G, treatment, outcome, mediators){
  descendants <- Filter(
    function(node){
      reachable <- id_directed_reachable(G, treatment, node)

      return(reachable)
    }, G$V
  )
  eligible <- sort(setdiff(
    id_ancestors(G, c(treatment, mediators, outcome)),
    c(treatment, mediators, outcome, descendants)
  ))
  qualifies <- function(covariates){
    qualifies <- all(covariates %in% eligible) &&
      isTRUE(dagitty::isAdjustmentSet(
        dag, Z = covariates, exposure = treatment, outcome = mediators
      )) &&
      isTRUE(dagitty::isAdjustmentSet(
        dag, Z = covariates, exposure = c(treatment, mediators), outcome = outcome
      ))

    return(qualifies)
  }
  reduce_set <- function(covariates){
    for(node in covariates){
      smaller <- setdiff(covariates, node)
      if( qualifies(smaller) ){
        covariates <- smaller
      }
    }
    return(list(sort(covariates)))
  }

  # Among observed nondescendants of the treatment, the full eligible set
  # passes both checks whenever any eligible set does. If it fails, no
  # eligible candidate passes, so the search stops here.
  for(candidate in unique(list(character(), eligible))){
    if( qualifies(candidate) ){
      return(reduce_set(candidate))
    }
  }

  return(list())
}


## Query validation ##
#' Check graph inputs for causal identification
#'
#' @param dag A dagitty object.
#' @return A list containing the `dag`, `edges`, `vertices`, `latent`, and
#'   `observed` nodes.
#' @noRd
id_validate_graph <- function(dag){
  if( !dagitty::is.dagitty(dag) ){
    stop("Invalid input detected: only dagitty objects may be supplied as 'dag' input. Please check the input and try again.",
         call. = FALSE)
  }
  graph_edges <- dagitty::edges(dag)
  if( nrow(graph_edges) > 0L && any(!graph_edges$e %in% c("->", "<->")) ){
    stop("Unsupported edge type detected: only directed or bidirected edges may be supplied as 'dag' input. Please supply a dagitty object containing only directed and/or bidirected edges.",
         call. = FALSE)
  }
  if( !identical(dagitty::graphType(dag), "dag") ){
    stop("Unsupported graph type detected: 'get_identification()' requires a 'dagitty' DAG. Please supply one and try again.",
         call. = FALSE)
  }
  if( !isTRUE(dag_is_acyclic(dag)) ){
    stop("Cycles detected: only acyclic dagitty objects are supported as 'dag' input. Please check the input and try again.",
         call. = FALSE)
  }

  status_detected <- length(dagitty::adjustedNodes(dag)) > 0L ||
    length(selected_nodes(dag)) > 0L
  if( status_detected ){
    stop("Adjusted/selected node(s) detected: these are not supported by 'get_identification()' as 'dag' input. Please remove these from the supplied dagitty object and try again.",
         call. = FALSE)
  }

  latent <- as.character(dagitty::latents(dag))
  bidirected <- graph_edges$e == "<->"
  if( any(bidirected) && any(
    as.character(graph_edges$v[bidirected]) %in% latent |
      as.character(graph_edges$w[bidirected]) %in% latent
  ) ){
    stop("Unsupported bidirected edge detected: a bidirected edge in the DAG connects to a latent node. 'get_identification()' currently supports bidirected edges only between observed nodes. Please revise the supplied DAG and try again.",
         call. = FALSE)
  }

  vertices <- as.character(names(dag))
  graph_info <- list(
    dag = dag,
    edges = graph_edges,
    vertices = vertices,
    latent = latent,
    observed = setdiff(vertices, latent)
  )

  return(graph_info)
}


#' Resolve treatment and outcome names
#'
#' @param graph Validated graph list containing `dag` and `observed` node
#'   names.
#' @param treatment Optional treatment name; `NULL` uses the graph's exposure
#'   tag.
#' @param outcome Optional outcome name; `NULL` uses the graph's outcome tag.
#' @return A list containing one observed `treatment` and one observed
#'   `outcome`.
#' @noRd
id_resolve_query <- function(graph, treatment, outcome){
  if( is.null(treatment) ){
    treatment <- as.character(dagitty::exposures(graph$dag))
  }
  if( is.character(treatment) && length(treatment) == 0L ){
    stop("No treatment detected: 'get_identification()' requires exactly one named treatment, either tagged in the DAG or supplied through 'treatment'. Please ensure a single treatment is supplied and try again.",
         call. = FALSE)
  }
  if( is.character(treatment) && length(treatment) > 1L ){
    stop("Multiple treatments detected: 'get_identification()' requires exactly one named treatment, either tagged in the DAG or supplied through 'treatment'. Please ensure a single treatment is supplied and try again.",
         call. = FALSE)
  }
  if( !is.character(treatment) || length(treatment) != 1L ||
      is.na(treatment) || !nzchar(treatment) ||
      !treatment %in% graph$observed ){
    stop("Invalid treatment detected: 'get_identification()' requires exactly one named treatment, either tagged in the DAG or supplied through 'treatment'. Please ensure a single treatment is supplied and try again.",
         call. = FALSE)
  }

  treatment <- as.character(treatment)

  if( is.null(outcome) ){
    outcome <- as.character(dagitty::outcomes(graph$dag))
  }
  if( is.character(outcome) && length(outcome) == 0L ){
    stop("No outcome detected: 'get_identification()' requires exactly one named outcome, either tagged in the DAG or supplied through 'outcome'. Please ensure a single outcome is supplied and try again.",
         call. = FALSE)
  }
  if( is.character(outcome) && length(outcome) > 1L ){
    stop("Multiple outcomes detected: 'get_identification()' requires exactly one named outcome, either tagged in the DAG or supplied through 'outcome'. Please ensure a single outcome is supplied and try again.",
         call. = FALSE)
  }
  if( !is.character(outcome) || length(outcome) != 1L ||
      is.na(outcome) || !nzchar(outcome) || !outcome %in% graph$observed ){
    stop("Invalid outcome detected: 'get_identification()' requires exactly one named outcome, either tagged in the DAG or supplied through 'outcome'. Please ensure a single outcome is supplied and try again.",
         call. = FALSE)
  }
  outcome <- as.character(outcome)
  if( identical(treatment, outcome) ){
    stop("Identical treatment and outcome detected: 'get_identification()' requires the named treatment and outcome to refer to different nodes in the DAG. Please ensure different treatment and outcome nodes are supplied and try again.",
         call. = FALSE)
  }

  query <- list(treatment = treatment, outcome = outcome)

  return(query)
}


#' Check conditioning or experimental node names
#'
#' @param nodes Character vector of node names.
#' @param observed Observed node names allowed in the graph.
#' @param endpoints Treatment and outcome names to exclude.
#' @param argument Argument name for messages.
#' @return Sorted unique node names, or an empty character vector for `NULL`.
#' @noRd
id_validate_nodes <- function(nodes, observed, endpoints, argument){
  if( is.null(nodes) ){
    nodes_out <- character()
    return(nodes_out)
  }
  valid <- is.character(nodes) && !anyNA(nodes) && all(nzchar(nodes)) &&
    all(nodes %in% setdiff(observed, endpoints))
  if( !valid ){
    if( identical(argument, "conditioning") ){
      stop("Invalid conditioning input detected: 'conditioning' requires observed DAG node names as input. Please ensure each supplied name corresponds to an observed node in the supplied DAG other than the treatment and outcome, and try again.",
           call. = FALSE)
    }
    stop("Invalid experiments input detected: 'experiments' requires observed DAG node names as input. Please ensure each supplied name corresponds to an observed node in the supplied DAG other than the treatment and outcome, and try again.",
         call. = FALSE)
  }
  nodes_out <- sort(unique(nodes))

  return(nodes_out)
}


## Projected graphs ##
# G holds observed vertices in V, directed edges in D and bidirected edges in B


#' Create an empty directed-edge table
#'
#' @return A data frame with character `from` and `to` columns.
#' @noRd
id_empty_directed <- function(){
  edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)

  return(edges)
}


#' Create an empty bidirected-edge table
#'
#' @return A data frame with character `v` and `w` columns.
#' @noRd
id_empty_bidirected <- function(){
  edges <- data.frame(v = character(), w = character(), stringsAsFactors = FALSE)

  return(edges)
}


#' Remove repeated and self-directed edges
#'
#' @param arcs Directed-edge table with `from` and `to` columns, or `NULL`.
#' @return A sorted `from`, `to` data frame of unique non-self edges.
#' @noRd
id_normalize_directed <- function(arcs){
  if( is.null(arcs) || nrow(arcs) == 0L ){
    edges <- id_empty_directed()
    return(edges)
  }
  arcs <- arcs[arcs$from != arcs$to, c("from", "to"), drop = FALSE]
  arcs <- unique(arcs)
  edges <- arcs[order(arcs$from, arcs$to), , drop = FALSE]

  return(edges)
}


#' Order endpoints and remove repeated bidirected edges
#'
#' @param arcs Bidirected-edge table with `v` and `w` columns, or `NULL`.
#' @return A sorted `v`, `w` data frame of unique non-self edges.
#' @noRd
id_normalize_bidirected <- function(arcs){
  if( is.null(arcs) || nrow(arcs) == 0L ){
    edges <- id_empty_bidirected()
    return(edges)
  }
  ends <- t(apply(arcs[, c("v", "w"), drop = FALSE], 1L, sort))
  arcs <- data.frame(
    v = ends[, 1L],
    w = ends[, 2L],
    stringsAsFactors = FALSE
  )
  arcs <- arcs[arcs$v != arcs$w, , drop = FALSE]
  arcs <- unique(arcs)
  edges <- arcs[order(arcs$v, arcs$w), , drop = FALSE]

  return(edges)
}


#' Find observed nodes reached through latent intermediates
#'
#' @param start Starting node name.
#' @param directed_children Named list of directed children for each node.
#' @param is_latent Named logical vector identifying latent nodes.
#' @return Unique observed endpoint names reached through latent intermediate
#'   nodes.
#' @noRd
id_latent_reach <- function(start, directed_children, is_latent){
  reached <- character()
  visited <- character()
  frontier <- directed_children[[start]]
  if( is.null(frontier) ){
    frontier <- character()
  }

  while(length(frontier) > 0L){
    node <- frontier[[1L]]
    frontier <- frontier[-1L]
    if( node %in% visited ){
      next
    }
    visited <- c(visited, node)
    if( isTRUE(is_latent[[node]]) ){
      frontier <- c(frontier, directed_children[[node]])
    }else{
      reached <- c(reached, node)
    }
  }

  nodes <- unique(reached)

  return(nodes)
}


#' Project a graph onto its observed nodes
#'
#' @param dag A dagitty object.
#' @return A mixed-graph list with observed vertices `V`, directed edges `D`,
#'   and bidirected edges `B`.
#' @noRd
id_admg <- function(dag){
  graph_edges <- dagitty::edges(dag)
  vertices <- as.character(names(dag))
  latent <- as.character(dagitty::latents(dag))
  observed <- setdiff(vertices, latent)
  is_latent <- stats::setNames(vertices %in% latent, vertices)

  children <- stats::setNames(
    lapply(vertices, function(node) {
      child_nodes <- character()

      return(child_nodes)
    }),
    vertices
  )
  directed_edges <- graph_edges[graph_edges$e == "->", , drop = FALSE]
  if( nrow(directed_edges) > 0L ){
    for(index in seq_len(nrow(directed_edges))){
      from <- as.character(directed_edges$v[[index]])
      to <- as.character(directed_edges$w[[index]])
      children[[from]] <- unique(c(children[[from]], to))
    }
  }

  directed <- id_empty_directed()
  for(from in observed){
    reached <- id_latent_reach(from, children, is_latent)
    if( length(reached) > 0L ){
      directed <- rbind(
        directed,
        data.frame(from = from, to = reached, stringsAsFactors = FALSE)
      )
    }
  }

  bidirected <- id_empty_bidirected()
  supplied_bidirected <- graph_edges[graph_edges$e == "<->", , drop = FALSE]
  if( nrow(supplied_bidirected) > 0L ){
    bidirected <- rbind(
      bidirected,
      data.frame(
        v = as.character(supplied_bidirected$v),
        w = as.character(supplied_bidirected$w),
        stringsAsFactors = FALSE
      )
    )
  }
  for(node in latent){
    reached <- sort(intersect(
      id_latent_reach(node, children, is_latent),
      observed
    ))
    if( length(reached) >= 2L ){
      pairs <- utils::combn(reached, 2L)
      bidirected <- rbind(
        bidirected,
        data.frame(
          v = pairs[1L,],
          w = pairs[2L,],
          stringsAsFactors = FALSE
        )
      )
    }
  }

  graph <- list(
    V = sort(observed),
    D = id_normalize_directed(directed),
    B = id_normalize_bidirected(bidirected)
  )

  return(graph)
}


#' Restrict a mixed graph to selected vertices
#'
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param vertices Vertex names to retain.
#' @return A mixed-graph list containing the selected vertices and their
#'   connecting edges.
#' @noRd
id_induced <- function(G, vertices){
  vertices <- sort(intersect(vertices, G$V))
  graph <- list(
    V = vertices,
    D = id_normalize_directed(
      G$D[G$D$from %in% vertices & G$D$to %in% vertices, , drop = FALSE]
    ),
    B = id_normalize_bidirected(
      G$B[G$B$v %in% vertices & G$B$w %in% vertices, , drop = FALSE]
    )
  )

  return(graph)
}


#' Find directed ancestors in a mixed graph
#'
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param targets Target node names.
#' @param drop_into Nodes whose incoming directed edges are ignored.
#' @return Sorted ancestor names, including the target nodes.
#' @noRd
id_ancestors <- function(G, targets, drop_into = character()){
  targets <- intersect(targets, G$V)
  seen <- targets
  frontier <- targets

  while(length(frontier) > 0L){
    parents <- G$D$from[
      G$D$to %in% frontier & !G$D$to %in% drop_into
    ]
    parents <- setdiff(unique(parents), seen)
    seen <- c(seen, parents)
    frontier <- parents
  }

  nodes <- sort(unique(seen))

  return(nodes)
}


#' Check whether a directed path joins two node sets
#'
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param from Source node names.
#' @param to Target node names.
#' @return `TRUE` if a directed path joins the source and target sets,
#'   otherwise `FALSE`.
#' @noRd
id_directed_reachable <- function(G, from, to){
  if( length(from) == 0L || length(to) == 0L ||
      !all(c(from, to) %in% G$V) ){
    matches <- FALSE
    return(matches)
  }
  seen <- unique(from)
  frontier <- unique(from)
  while(length(frontier) > 0L){
    children <- unique(G$D$to[G$D$from %in% frontier])
    if( any(children %in% to) ){
      matches <- TRUE
      return(matches)
    }
    children <- setdiff(children, seen)
    seen <- c(seen, children)
    frontier <- children
  }
  matches <- FALSE

  return(matches)
}


#' Find bidirected connected components
#'
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @return A list of sorted vertex sets, one per bidirected component.
#' @noRd
id_ccomponents <- function(G){
  if( length(G$V) == 0L ){
    components_out <- list()
    return(components_out)
  }
  remaining <- G$V
  components <- list()

  while(length(remaining) > 0L){
    component <- remaining[[1L]]
    frontier <- component
    while(length(frontier) > 0L){
      neighbours <- unique(c(
        G$B$w[G$B$v %in% frontier],
        G$B$v[G$B$w %in% frontier]
      ))
      neighbours <- setdiff(neighbours, component)
      component <- c(component, neighbours)
      frontier <- neighbours
    }
    component <- sort(unique(component))
    components[[length(components) + 1L]] <- component
    remaining <- setdiff(remaining, component)
  }

  components_out <- components[order(vapply(components, paste, character(1L), collapse = "\r"))]

  return(components_out)
}


## Identification algorithms ##
#' Check an intervention query with the ID algorithm
#'
#' @param y Outcome node names.
#' @param x Intervention node names.
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @return A list containing the identification flag `id` and a `hedge` when
#'   identification fails.
#' @noRd
id_solve <- function(y, x, G){ # Return the ID verdict and recover a hedge on failure
  y <- sort(intersect(y, G$V))
  x <- sort(intersect(x, setdiff(G$V, y)))
  vertices <- G$V

  if( length(x) == 0L ){
    identification <- list(id = TRUE, hedge = NULL)
    return(identification)
  }

  ancestors_y <- id_ancestors(G, y)
  if( !setequal(ancestors_y, vertices) ){
    identification <- id_solve(
      y,
      intersect(x, ancestors_y),
      id_induced(G, ancestors_y)
    )
    return(identification)
  }

  ancestors_mutilated <- id_ancestors(G, y, drop_into = x)
  additional <- setdiff(setdiff(vertices, x), ancestors_mutilated)
  if( length(additional) > 0L ){
    identification <- id_solve(y, union(x, additional), G)
    return(identification)
  }

  remaining_graph <- id_induced(G, setdiff(vertices, x))
  components <- id_ccomponents(remaining_graph)
  if( length(components) > 1L ){
    for(component in components){
      solved <- id_solve(component, setdiff(vertices, component), G)
      if( !isTRUE(solved$id) ){
        return(solved)
      }
    }
    identification <- list(id = TRUE, hedge = NULL)
    return(identification)
  }
  if( length(components) != 1L ){
    stop("Internal error detected: the identification algorithm reached no valid case. Please report this issue.",
         call. = FALSE)
  }

  component <- components[[1L]]
  graph_components <- id_ccomponents(G)
  if( length(graph_components) == 1L &&
      setequal(graph_components[[1L]], vertices) ){
    hedge <- id_recover_hedge(G, component, y, x)
    identification <- list(id = FALSE, hedge = hedge)
    return(identification)
  }

  if( any(vapply(
    graph_components,
    function(candidate) {
      matches <- setequal(candidate, component)

      return(matches)
    },
    logical(1L)
  )) ){
    identification <- list(id = TRUE, hedge = NULL)
    return(identification)
  }

  larger <- Filter(
    function(candidate) {
      matches <- all(component %in% candidate)

      return(matches)
    },
    graph_components
  )
  if( length(larger) == 1L ){
    larger <- larger[[1L]]
    identification <- id_solve(
      y,
      intersect(x, larger),
      id_induced(G, larger)
    )
    return(identification)
  }

  stop("Internal error detected: the identification algorithm reached no valid case. Please report this issue.",
       call. = FALSE)
}


#' Find directed distances to a target set
#'
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param targets Target node names.
#' @return Named minimum directed distances, with `Inf` for nodes that cannot
#'   reach a target.
#' @noRd
id_distance_to <- function(G, targets){
  distance <- stats::setNames(rep(Inf, length(G$V)), G$V)
  distance[targets] <- 0
  changed <- TRUE
  while(changed){
    changed <- FALSE
    for(node in G$V){
      children <- G$D$to[G$D$from == node]
      if( length(children) == 0L ){
        next
      }
      candidate <- 1 + min(distance[children])
      if( is.finite(candidate) && candidate < distance[[node]] ){
        distance[[node]] <- candidate
        changed <- TRUE
      }
    }
  }

  return(distance)
}


#' Construct a hedge for an unidentified query
#'
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param S Bidirected component remaining after removing intervention nodes.
#' @param y Outcome node names.
#' @param x Intervention node names.
#' @return A list with dagitty graphs `F` and `F_prime`, their `roots`, and the
#'   reduced `query`.
#' @noRd
id_recover_hedge <- function(G, S, y, x){
  S <- sort(intersect(S, G$V))
  directed_s <- G$D[G$D$from %in% S & G$D$to %in% S, , drop = FALSE]
  roots <- sort(setdiff(S, unique(directed_s$from)))
  distance <- id_distance_to(G, roots)

  selected <- id_empty_directed()
  for(node in setdiff(G$V, roots)){
    children <- G$D$to[G$D$from == node]
    if( node %in% S ){
      children <- intersect(children, S)
    }
    children <- children[is.finite(distance[children])]
    if( length(children) == 0L ){
      stop("Internal error detected: the recovered hedge failed validation. Please report this issue.",
           call. = FALSE)
    }
    child_order <- order(distance[children], children)
    selected <- rbind(
      selected,
      data.frame(
        from = node,
        to = children[child_order[[1L]]],
        stringsAsFactors = FALSE
      )
    )
  }

  F <- list(
    V = G$V,
    D = id_normalize_directed(selected),
    B = G$B
  )
  F_prime <- id_induced(F, S)

  if( !id_validate_hedge(F, F_prime, G, roots, y, x) ){
    stop("Internal error detected: the recovered hedge failed validation. Please report this issue.",
         call. = FALSE)
  }

  hedge <- list(
    F = id_as_dagitty(F),
    F_prime = id_as_dagitty(F_prime),
    roots = roots,
    query = list(
      outcome = sort(y),
      intervention = sort(x),
      vertices = sort(G$V)
    )
  )

  return(hedge)
}


#' Create directed-edge comparison keys
#'
#' @param arcs Directed-edge table with `from` and `to` columns.
#' @return A character key for each directed edge.
#' @noRd
id_directed_keys <- function(arcs){
  if( nrow(arcs) == 0L ){
    keys <- character()
    return(keys)
  }
  keys <- paste(arcs$from, arcs$to, sep = "\r")

  return(keys)
}


#' Create bidirected-edge comparison keys
#'
#' @param arcs Bidirected-edge table with `v` and `w` columns.
#' @return A character key for each distinct, consistently ordered bidirected
#'   edge.
#' @noRd
id_bidirected_keys <- function(arcs){
  if( nrow(arcs) == 0L ){
    keys <- character()
    return(keys)
  }
  arcs <- id_normalize_bidirected(arcs)
  keys <- paste(arcs$v, arcs$w, sep = "\r")

  return(keys)
}


#' Find roots of a directed forest used in a hedge
#'
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @return Sorted node names with no outgoing directed edges.
#' @noRd
id_roots <- function(G){
  roots <- sort(setdiff(G$V, unique(G$D$from)))

  return(roots)
}


#' Check a proposed hedge against its graph and query
#'
#' @param F Larger mixed-graph forest.
#' @param F_prime Smaller mixed-graph forest.
#' @param source Mixed graph from which the forests were selected.
#' @param roots Expected root node names shared by both forests.
#' @param y Outcome node names.
#' @param x Intervention node names.
#' @return `TRUE` when the two forests satisfy the hedge checks, otherwise
#'   `FALSE`.
#' @noRd
id_validate_hedge <- function(F, F_prime, source, roots, y, x){
  if( length(F$V) <= length(F_prime$V) ||
      !all(F_prime$V %in% F$V) ||
      length(id_ccomponents(F)) != 1L ||
      length(id_ccomponents(F_prime)) != 1L ){
    valid <- FALSE
    return(valid)
  }

  children_F <- table(factor(F$D$from, levels = F$V))
  children_F_prime <- table(factor(F_prime$D$from, levels = F_prime$V))
  if( any(children_F > 1L) || any(children_F_prime > 1L) ||
      !setequal(id_roots(F), roots) ||
      !setequal(id_roots(F_prime), roots) ){
    valid <- FALSE
    return(valid)
  }

  if( length(intersect(F$V, x)) == 0L ||
      length(intersect(F_prime$V, x)) > 0L ||
      !all(roots %in% id_ancestors(source, y, drop_into = x)) ){
    valid <- FALSE
    return(valid)
  }

  source_d <- id_directed_keys(source$D)
  source_b <- id_bidirected_keys(source$B)
  if( !all(id_directed_keys(F$D) %in% source_d) ||
      !all(id_bidirected_keys(F$B) %in% source_b) ){
    valid <- FALSE
    return(valid)
  }

  induced <- id_induced(F, F_prime$V)
  valid <- setequal(id_directed_keys(induced$D), id_directed_keys(F_prime$D)) &&
    setequal(id_bidirected_keys(induced$B), id_bidirected_keys(F_prime$B))

  return(valid)
}


#' Quote a node name for dagitty text
#'
#' @param node Node name to format.
#' @return A double-quoted node name with embedded double quotes escaped.
#' @noRd
id_quote_node <- function(node){
  text <- paste0("\"", gsub("\"", "\\\"", node, fixed = TRUE), "\"")

  return(text)
}


#' Convert a mixed graph to a dagitty object
#'
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @return A dagitty object containing the vertices and directed and bidirected
#'   edges.
#' @noRd
id_as_dagitty <- function(G){
  lines <- vapply(G$V, id_quote_node, character(1L))
  if( nrow(G$D) > 0L ){
    lines <- c(
      lines,
      paste(
        vapply(G$D$from, id_quote_node, character(1L)),
        "->",
        vapply(G$D$to, id_quote_node, character(1L))
      )
    )
  }
  if( nrow(G$B) > 0L ){
    lines <- c(
      lines,
      paste(
        vapply(G$B$v, id_quote_node, character(1L)),
        "<->",
        vapply(G$B$w, id_quote_node, character(1L))
      )
    )
  }
  dag <- dagitty::dagitty(paste0("dag {\n", paste(lines, collapse = "\n"), "\n}"))

  return(dag)
}


#' Remove edges for an IDC separation check
#'
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param intervention Nodes whose incoming directed and incident bidirected
#'   edges are removed.
#' @param remove_outgoing Nodes whose outgoing directed edges are removed.
#' @return A dagitty object after removing the specified incoming and outgoing
#'   edges.
#' @noRd
idc_mutilate <- function(G, intervention, remove_outgoing){
  directed <- G$D[
    !G$D$to %in% intervention & !G$D$from %in% remove_outgoing,
    ,
    drop = FALSE
  ]
  bidirected <- G$B[
    !G$B$v %in% intervention & !G$B$w %in% intervention,
    ,
    drop = FALSE
  ]
  dag <- id_as_dagitty(list(V = G$V, D = directed, B = bidirected))

  return(dag)
}


#' Check a conditional intervention query with IDC
#'
#' @param y Outcome node names.
#' @param x Intervention node names.
#' @param z Conditioning node names.
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @return A list containing `id`, `hedge`, and the conditioning sets
#'   `intervened` and `remaining`.
#' @noRd
idc_solve <- function(y, x, z, G){
  intervention <- x
  remaining <- sort(z)
  intervened <- character()

  repeat {
    shifted <- FALSE
    for(candidate in remaining){
      mutilated <- idc_mutilate(G, intervention, candidate)
      given <- union(intervention, setdiff(remaining, candidate))
      movable <- dagitty::dseparated(
        mutilated,
        X = y,
        Y = candidate,
        Z = given
      )
      if( isTRUE(movable) ){
        intervention <- union(intervention, candidate)
        intervened <- c(intervened, candidate)
        remaining <- setdiff(remaining, candidate)
        shifted <- TRUE
        break
      }
    }
    if( !shifted ){
      break
    }
  }

  solved <- id_solve(union(y, remaining), intervention, G)
  identification <- list(
    id = solved$id,
    hedge = solved$hedge,
    intervened = sort(intervened),
    remaining = sort(remaining)
  )

  return(identification)
}


## Individual identification strategies ##
#' Check covariate adjustment for the total effect
#'
#' @param dag A dagitty object.
#' @param treatment Observed treatment node name.
#' @param outcome Observed outcome node name.
#' @return A list with a `feasible` flag, all minimal adjustment `sets`, and
#'   the `canonical` adjustment set.
#' @noRd
id_adjustment <- function(dag, treatment, outcome){
  sets <- dagitty::adjustmentSets(
    dag,
    exposure = treatment,
    outcome = outcome,
    type = "minimal",
    effect = "total"
  )
  sets <- lapply(unclass(sets), function(set) {
    adjustment_set <- sort(as.character(set))

    return(adjustment_set)
  })
  canonical <- dagitty::adjustmentSets(
    dag,
    exposure = treatment,
    outcome = outcome,
    type = "canonical",
    effect = "total"
  )
  canonical <- lapply(unclass(canonical), function(set){
    adjustment_set <- sort(as.character(set))

    return(adjustment_set)
  })
  adjustment <- list(
    feasible = length(sets) > 0L,
    sets = unname(sets),
    canonical = unname(canonical)
  )

  return(adjustment)
}


#' Find mediator sets satisfying the front door checks
#'
#' @param dag A dagitty object.
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param treatment Observed treatment node name.
#' @param outcome Observed outcome node name.
#' @return A list with a `feasible` flag and the smallest qualifying mediator
#'   `sets`.
#' @noRd
id_frontdoor <- function(dag, G, treatment, outcome){
  if( any(G$D$from == treatment & G$D$to == outcome) ){
    frontdoor <- list(feasible = FALSE, sets = list())
    return(frontdoor)
  }

  candidates <- Filter(
    function(node){
      matches <- id_directed_reachable(G, treatment, node) &&
        id_directed_reachable(G, node, outcome)

      return(matches)
    },
    setdiff(G$V, c(treatment, outcome))
  )
  # A set has no open back-door path from the treatment exactly when each of
  # its members has none.
  candidates <- Filter(
    function(node){
      matches <- isTRUE(dagitty::isAdjustmentSet(
        dag,
        Z = character(),
        exposure = treatment,
        outcome = node
      ))

      return(matches)
    },
    candidates
  )
  # A member with a back-door path to the outcome that is open given the
  # treatment and avoids the other candidates fails the last check in every
  # candidate set containing it. Removing it can expose further such paths.
  repeat {
    blocked <- Filter(
      function(node){
        others <- setdiff(candidates, node)
        H <- id_induced(G, setdiff(G$V, others))
        H$D <- H$D[H$D$from != node, , drop = FALSE]
        matches <- !isTRUE(dagitty::dseparated(
          id_as_dagitty(H),
          X = node,
          Y = outcome,
          Z = treatment
        ))

        return(matches)
      },
      candidates
    )
    if( length(blocked) == 0L ){
      break
    }
    candidates <- setdiff(candidates, blocked)
  }
  # Every front-door set lies within the remaining candidates, so if these
  # leave a directed path open, no subset intercepts it.
  if( length(candidates) == 0L || id_directed_reachable(
    id_induced(G, setdiff(G$V, candidates)),
    treatment,
    outcome
  ) ){
    frontdoor <- list(feasible = FALSE, sets = list())
    return(frontdoor)
  }
  sets <- list()
  if( length(candidates) > 0L ){
    for(size in seq_along(candidates)){
      candidate_sets <- utils::combn(candidates, size, simplify = FALSE)
      for(mediators in candidate_sets){
        without_mediators <- id_induced(G, setdiff(G$V, mediators))
        intercepts <- !id_directed_reachable(
          without_mediators,
          treatment,
          outcome
        )
        if( !intercepts ){
          next
        }

        no_treatment_mediator_confounding <- isTRUE(
          dagitty::isAdjustmentSet(
            dag,
            Z = character(),
            exposure = treatment,
            outcome = mediators
          )
        )
        treatment_blocks_mediator_outcome <- isTRUE(
          dagitty::isAdjustmentSet(
            dag,
            Z = treatment,
            exposure = mediators,
            outcome = outcome
          )
        )
        if( no_treatment_mediator_confounding &&
            treatment_blocks_mediator_outcome ){
          sets[[length(sets) + 1L]] <- sort(mediators)
        }
      }
      if( length(sets) > 0L ){
        break
      }
    }
  }

  frontdoor <- list(feasible = length(sets) > 0L, sets = sets)

  return(frontdoor)
}


#' Find instruments and their conditioning sets
#'
#' @param dag A dagitty object.
#' @param treatment Observed treatment node name.
#' @param outcome Observed outcome node name.
#' @return A list with a `feasible` flag and instrument `entries`.
#' @noRd
id_instruments <- function(dag, treatment, outcome){
  entries <- dagitty::instrumentalVariables(
    dag,
    exposure = treatment,
    outcome = outcome
  )
  entries <- lapply(unclass(entries), function(entry){
    instrument <- list(
      instrument = as.character(entry$I),
      conditioning = sort(as.character(entry$Z))
    )

    return(instrument)
  })
  instruments <- list(feasible = length(entries) > 0L, entries = unname(entries))

  return(instruments)
}


#' Check identification using surrogate experiments
#'
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param treatment Observed treatment node name.
#' @param outcome Observed outcome node name.
#' @param experiments Observed nodes with available surrogate experimental
#'   data.
#' @param observational Whether the query is identifiable from observational
#'   data.
#' @return A list with a `zid` flag and the qualifying experimental node `set`.
#'   The set is empty when observational identification succeeds or no
#'   experimental set qualifies.
#' @noRd
id_zidentifiable <- function(G, treatment, outcome, experiments,
                              observational){
  if( isTRUE(observational) ){
    identification <- list(zid = TRUE, set = character())
    return(identification)
  }

  candidates <- sort(intersect(
    experiments,
    setdiff(G$V, c(treatment, outcome))
  ))
  graph_without_treatment <- id_induced(G, setdiff(G$V, treatment))

  for(size in seq_along(candidates)){
    subsets <- utils::combn(candidates, size, simplify = FALSE)
    for(experiment_set in subsets){
      bypasses_treatment <- any(vapply(
        experiment_set,
        function(node){
          reaches_outcome <- id_directed_reachable(graph_without_treatment, node, outcome)

          return(reaches_outcome)
        },
        logical(1L)
      ))
      if( bypasses_treatment ){
        next
      }

      randomized <- list(
        V = G$V,
        D = G$D[!G$D$to %in% experiment_set, , drop = FALSE],
        B = G$B[
          !G$B$v %in% experiment_set & !G$B$w %in% experiment_set,
          ,
          drop = FALSE
        ]
      )
      if( isTRUE(id_solve(outcome, treatment, randomized)$id) ){
        identification <- list(zid = TRUE, set = sort(experiment_set))
        return(identification)
      }
    }
  }

  identification <- list(zid = FALSE, set = character())

  return(identification)
}


## Result formatting ##
#' Format one node name within a printed set
#'
#' @param node Node name to format.
#' @return The node name, quoted when it contains set punctuation.
#' @noRd
id_render_set_node <- function(node){
  if( grepl("[{},;]", node) ){
    text <- id_quote_node(node)
    return(text)
  }

  return(node)
}


#' Format a set of node names
#'
#' @param set Node names to display.
#' @return One character string enclosing sorted unique names in braces.
#' @noRd
id_render_set <- function(set){
  set <- sort(unique(as.character(set)))
  if( length(set) == 0L ){
    text <- "{}"
    return(text)
  }
  rendered <- vapply(set, id_render_set_node, character(1L))
  text <- paste0("{", paste(rendered, collapse = ", "), "}")

  return(text)
}


#' Format several sets of node names
#'
#' @param sets List of node-name sets.
#' @param max_sets Maximum number of sets to display; `Inf` displays all sets.
#' @return One character string separating the sets with semicolons, or
#'   `"none"`. If sets are omitted, their count follows the displayed sets.
#' @noRd
id_render_sets <- function(sets, max_sets = Inf){
  if( length(sets) == 0L ){
    text <- "none"
    return(text)
  }
  shown <- min(length(sets), max_sets)
  text <- paste(
    vapply(sets[seq_len(shown)], id_render_set, character(1L)),
    collapse = "; "
  )
  if( shown < length(sets) ){
    text <- paste0(text, "; ... (", length(sets) - shown, " more)")
  }

  return(text)
}


#' Format instruments and their conditioning sets
#'
#' @param entries List of entries containing `instrument` and `conditioning`
#'   node names.
#' @return One character string describing the instruments, or `"none"`.
#' @noRd
id_render_instruments <- function(entries){
  if( length(entries) == 0L ){
    text <- "none"
    return(text)
  }
  rendered <- vapply(entries, function(entry){
    instrument <- id_render_set(entry$instrument)
    if( length(entry$conditioning) == 0L ){

      return(instrument)
    }else{
      text <- paste0(
        instrument,
        " conditional on ",
        id_render_set(entry$conditioning)
      )

      return(text)
    }
  }, character(1L))
  text <- paste(rendered, collapse = "; ")

  return(text)
}


#' Construct one identification strategy row
#'
#' @param strategy Name of the identification strategy.
#' @param feasible Logical result of the strategy check.
#' @param tier Evidence or assumption category for the strategy.
#' @param variable_set Printed variable-set description.
#' @param note Explanation of the strategy result.
#' @return A data.table with `strategy`, `feasible`, `tier`, `variable_set`,
#'   and `note` columns.
#' @noRd
id_result_row <- function(strategy, feasible, tier, variable_set, note){
  row <- data.table::data.table(
    strategy = strategy,
    feasible = isTRUE(feasible),
    tier = tier,
    variable_set = variable_set,
    note = note
  )

  return(row)
}


#' Describe the conditioning split produced by IDC
#'
#' @param intervened Conditioning nodes moved into the intervention set.
#' @param remaining Nodes remaining in the conditioning set.
#' @return A character description of the moved and remaining conditioning
#'   sets.
#' @noRd
id_idc_detail <- function(intervened, remaining){
  if( length(intervened) > 0L && length(remaining) > 0L ){
    text <- paste0(
      "IDC moves ", id_render_set(intervened),
      " into the intervention set; ", id_render_set(remaining),
      " remain in the conditioning set"
    )
    return(text)
  }
  if( length(intervened) > 0L ){
    text <- paste0(
      "IDC moves ", id_render_set(intervened),
      " into the intervention set; no variables remain in the conditioning set"
    )
    return(text)
  }
  text <- paste0(
    "IDC moves no conditioning variables; ",
    id_render_set(remaining),
    " remain in the conditioning set"
  )

  return(text)
}


#' Summarise variable roles for an identification query
#'
#' @param G Observed-node mixed graph as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param treatment Observed treatment node name.
#' @param outcome Observed outcome node name.
#' @param effect Effect type, `"total"` or `"direct"`.
#' @param adjustment Covariate-adjustment result containing `sets`.
#' @param frontdoor Front-door result containing mediator `sets`.
#' @param instruments Instrument result containing `entries`.
#' @return A data.table with `variable`, `role`, and `conditioning_guidance`
#'   columns.
#' @noRd
id_roles <- function(G, treatment, outcome, effect, adjustment,
                      frontdoor, instruments){
  variables <- sort(setdiff(G$V, c(treatment, outcome)))
  role_variable <- character()
  role_name <- character()
  role_guidance <- character()

  add_role <- function(nodes, role, guidance){ # Append to the vectors in this id_roles() call
    nodes <- sort(intersect(unique(nodes), variables))
    if( length(nodes) == 0L ){
      out <- NULL
      return(invisible(out))
    }
    role_variable <<- c(role_variable, nodes)
    role_name <<- c(role_name, rep(role, length(nodes)))
    role_guidance <<- c(role_guidance, rep(guidance, length(nodes)))
    out <- NULL

    return(invisible(out))
  }

  adjustment_nodes <- unique(unlist(adjustment$sets, use.names = FALSE))
  add_role(
    adjustment_nodes,
    "adjustment variable",
    "member of at least one minimal adjustment set for the total effect; use a complete set"
  )

  if( identical(effect, "total") ){
    mediators <- Filter(
      function(node){
        matches <- id_directed_reachable(G, treatment, node) &&
          id_directed_reachable(G, node, outcome)

        return(matches)
      },
      variables
    )
    add_role(
      mediators,
      "mediator",
      "do not include in an adjustment set for the total effect"
    )
    add_role(
      unique(unlist(frontdoor$sets, use.names = FALSE)),
      "front-door mediator",
      "used by the front-door strategy; this is not ordinary covariate adjustment"
    )
    add_role(
      vapply(
        instruments$entries,
        function(entry) {
          instrument <- entry$instrument

          return(instrument)
        },
        character(1L)
      ),
      "instrument",
      "used by the instrumental variable strategy; covariate adjustment is assessed separately"
    )
  }else{
    outcome_parents <- intersect(
      unique(G$D$from[G$D$to == outcome]),
      variables
    )
    add_role(
      outcome_parents,
      "outcome parent",
      "held fixed in the direct effect query"
    )
  }

  add_role(
    setdiff(variables, role_variable),
    "other",
    "no variable-level conclusion; assess as part of a set"
  )

  roles <- data.table::data.table(
    variable = role_variable,
    role = role_name,
    conditioning_guidance = role_guidance
  )
  if( nrow(roles) > 0L ){
    roles <- unique(roles)
    roles <- roles[order(roles$variable, roles$role),]
  }

  return(roles)
}


#' Print the causal identification verdict
#'
#' @param treatment Observed treatment node name.
#' @param outcome Observed outcome node name.
#' @param effect Effect type, `"total"` or `"direct"`.
#' @param conditioning Observed conditioning node names.
#' @param experiments Observed nodes with available surrogate experimental
#'   data.
#' @param direct_path Whether a directed treatment-to-outcome path exists.
#' @param direct_arrow Whether a directed treatment-to-outcome edge exists.
#' @param observational Whether the query is identifiable from observational
#'   data.
#' @param identifiable Whether the query is identifiable from the available
#'   data.
#' @return Invisible `NULL`.
#' @noRd
id_summary_message <- function(treatment, outcome, effect, conditioning,
                                experiments, direct_path, direct_arrow,
                                observational, identifiable){
  if( identical(effect, "direct") ){
    if( !direct_arrow ){
      message(
        "The latent-projected graph contains no direct arrow from '",
        treatment, "' to '", outcome, "'. The verdict attributes refer to ",
        "the intervention on the treatment and fixing set."
      )
    }else if( observational ){
      message(
        "The outcome distribution for the controlled direct effect of '",
        treatment, "' on '", outcome, "' is identifiable from observational ",
        "data in the latent-projected graph (with the outcome's other ",
        "parents held fixed)."
      )
    }else{
      message(
        "The outcome distribution for the controlled direct effect of '",
        treatment, "' on '", outcome, "' is not identifiable from ",
        "observational data in the latent-projected graph (with the ",
        "outcome's other parents held fixed)."
      )
    }
    out <- NULL
    return(invisible(out))
  }

  if( length(conditioning) > 0L ){
    if( observational ){
      message(
        "The effect of '", treatment, "' on '", outcome, "' given ",
        id_render_set(conditioning),
        " is identifiable from observational data under this DAG."
      )
    }else{
      message(
        "The effect of '", treatment, "' on '", outcome, "' given ",
        id_render_set(conditioning),
        " is not identifiable from observational data under this DAG."
      )
    }
    out <- NULL
    return(invisible(out))
  }

  if( !direct_path ){
    message(
      "The DAG contains no directed path from '", treatment, "' to '",
      outcome, "', so it encodes no total causal effect of '", treatment,
      "' on '", outcome, "'."
    )
  }else if( observational ){
    message(
      "The total effect of '", treatment, "' on '", outcome,
      "' is identifiable from observational data under this DAG."
    )
  }else if( identifiable ){
    message(
      "The total effect of '", treatment, "' on '", outcome,
      "' is not identifiable from observational data alone, but is ",
      "identifiable using surrogate experimental data for ",
      id_render_set(experiments), " declared through 'experiments'."
    )
  }else{
    message(
      "The total effect of '", treatment, "' on '", outcome,
      "' is not identifiable from the supplied data under this DAG."
    )
  }
  out <- NULL

  return(invisible(out))
}
