# ==============================================================================
# identification.R


# get_identification() and the identification methods it calls
# ==============================================================================


#' Check if a causal effect is identifiable
#'
#' `get_identification()` checks identification strategies for observed treatment
#' and outcome nodes in the supplied DAG. Results are grouped by effect, with a
#' separate list for each complete strategy. Identification verdicts are stored
#' and displayed separately. With one treatment and one outcome, and without
#' `conditioning` or `experiments`, the default call checks total, direct and
#' indirect effects. With `conditioning` or `experiments`, only total effects
#' are supported. With several treatments or outcomes, it checks one joint
#' total effect.
#'
#' The default strategy list includes successful checks. Use `id_status = "all"`
#' to include unsuccessful checks, candidates and checks not assessed. This
#' selection does not hide the separate identification verdicts. Each strategy
#' contains its method and one complete set of variables. Their order does not
#' rank the methods. See [identification_methods] for status definitions,
#' assumptions and sources.
#'
#' @section Inspecting the results:
#' Printing the result shows identification verdicts and numbered strategies
#' within each effect group. HTML documents show the same results and fields.
#' Access the nested lists normally, for example `result$direct[[1]]`.
#' A subset made with `[`, such as `result["direct"]`, is an independent plain
#' list without the report's class or analysis attributes.
#'
#' `summary(result)` returns all strategy statuses within the same effect and
#' method selection. It uses stored analyses without performing new
#' identification calculations. `print(result, details = TRUE)` prints this
#' summary. If the nested contents have been edited, both use the current
#' contents without restoring removed strategies or displaying stored verdicts.
#'
#' `as.data.frame(result)` and `data.table::as.data.table(result)` give one row
#' per complete strategy, with columns `effect`, `method`, `assumptions`,
#' `variables`, `id_status` and `notes`. The `variables` column retains a named
#' list of character vectors for each row. The other columns are character
#' vectors, with assumptions and notes combined into text.
#'
#' Identification verdicts are available in `attr(result, "verdicts")`, grouped
#' by effect and named by the check. The table conversions retain this attribute
#' while the report contents are unchanged. Verdicts are not strategy rows.
#'
#' The function does not fit a model, return an identifying formula or prescribe
#' an estimator. Results depend on the supplied graph and each method's
#' assumptions. There is no single yes/no verdict for the whole result.
#'
#' @param dag A `dagitty` object representing a causal DAG. Bidirected
#'   edges are supported between observed nodes. All edges connected to explicitly
#'   declared latent nodes must be directed.
#' @param treatment A nonempty character vector naming observed treatment nodes.
#'   If `NULL`, uses the DAG's `[exposure]` nodes. Repeated names are used once.
#' @param outcome A nonempty character vector naming observed outcome nodes.
#'   If `NULL`, uses the DAG's `[outcome]` nodes. Treatment and outcome must not
#'   share a node.
#' @param effect Omit this argument or use `NULL` for an overview of total,
#'   direct and indirect effects when there is one treatment and one outcome.
#'   With several treatments or outcomes, or with `conditioning` or
#'   `experiments`, the function checks only the total effect. Supply `"total"`,
#'   `"direct"` or `"indirect"` to request only that effect group.
#'   Supply one value at a time.
#' @param experiments Optional observed node names with experimental data from
#'   interventions on those nodes. Cannot include treatment or outcome.
#'   Requires a total effect query and `conditioning = NULL`.
#' @param conditioning Optional observed node names to condition on in
#'   \eqn{P(y \mid do(x), z)}. Cannot include treatment or outcome.
#'   Requires a total effect query and `experiments = NULL`.
#' @param quiet If `TRUE`, suppress the message emitted during the call.
#'   The returned object can still be printed, including by R's automatic display.
#' @param mediator Optional nonempty character vector of observed mediators on
#'   directed paths from treatment to outcome. Defines the mediator set for
#'   natural direct or indirect effects. With `NULL`, each observed mediator
#'   is checked separately, followed by their joint set if there is more than
#'   one. Latent mediators cannot be selected. Requires a selection that includes
#'   natural direct or indirect effects. This does not change either front door
#'   search, the single door search or the stored parent intervention diagnostic.
#' @param method Optional character vector selecting strategies by their exact
#'   names: `"adjustment"`, `"front_door"`, `"conditional_front_door"`,
#'   `"instrumental_variable"`, `"single_door"` or `"mediation_formula"`.
#'   With `NULL`, includes all compatible strategies. Several values select
#'   alternatives. Only combinations compatible with `effect` are retained.
#'   An error is returned if none remain. This selection does not remove the
#'   identification verdicts for the selected effects.
#' @param id_status Optional character vector selecting strategy statuses.
#'   The default `"identifiable"` includes successful strategies. Use `NULL`
#'   or `"all"` to include all statuses. Other strategy statuses are
#'   `"not_established"`, `"not_assessed"` and `"candidate"`.
#'   The values `"not_identifiable"` and `"zero_by_graph"` are accepted but
#'   do not match strategy records. They occur only in the separate verdicts,
#'   which this argument does not filter. Several values select alternatives,
#'   with `"all"` taking precedence.
#'
#' @returns A `causaliflower_identification` object containing nested lists.
#'   The first level names effect groups (`total`, `direct` or `indirect`).
#'   Each group contains complete strategies. Each strategy has five fields:
#'   \itemize{
#'     \item `method`: the identification strategy.
#'     \item `assumptions`: a character vector of required assumptions.
#'     \item `variables`: a named list of character vectors, labelled by role.
#'       Node names are preserved exactly. An empty covariate vector means no
#'       adjustment is required for that strategy.
#'     \item `id_status`: the result of the strategy check.
#'       Definitions and their scope are in [identification_methods].
#'     \item `notes`: a character vector describing restrictions or the
#'       reason for the result.
#'   }
#'   Selected effect groups remain present when no strategy matches. Their lists
#'   are empty, while the relevant identification verdicts remain available.
#'
#' @section Stored analyses:
#'   The `verdicts` attribute contains the assessed identification checks,
#'   grouped by effect. Each record has `method`, `query`, `id_status`,
#'   `assumptions`, `variables` and `notes`. ID, IDC and z-ID concern the
#'   requested total effect query. Graph structure records describe an effect
#'   or coefficient that is zero by graph. An unassessed algorithm has no
#'   verdict record. These verdicts are independent of the strategy filters.
#'
#'   Attributes `treatment`, `outcome`, `scope` and `effect` describe the
#'   treatment and outcome sets and the effects checked. The `selection`
#'   attribute records the selectors. The `analyses` attribute contains named
#'   records for the included checks:
#'   \itemize{
#'     \item `total`: the query and its `identifiable`,
#'       `observationally_identifiable`, `z_identifiable`,
#'       `data_availability`, `hedge`, `roles`, `adjustment_sets`,
#'       `canonical_adjustment`, `frontdoor_sets`, `instruments`,
#'       `experiments`, `experiment_set` and `caveat`.
#'     \item `conditional_frontdoor`: a `feasible` flag and paired `sets` of
#'       `mediators` and `covariates` for the conditional front door strategies.
#'     \item `direct`: the supplied DAG's `direct_arrow`, single door
#'       `feasible`, `status` and complete `adjustment_sets`, with
#'       `treatment` and `outcome`.
#'     \item `controlled`: the same fields as `total`, for the separate
#'       intervention on treatment and the outcome's other parents in the graph
#'       obtained by latent projection. This diagnostic is not displayed as a
#'       direct effect or included in the nested results. See
#'       [identification_methods].
#'     \item `mediation`: entries with `mediators`, `feasible`,
#'       `adjustment_sets`, `note` and `assumptions`. Empty if there are no
#'       observed mediators. A passing entry returns one qualifying set,
#'       not necessarily a minimal one.
#'   }
#'   Query records specify `intervention`, `fixing_set` and `conditioning`.
#'   The fields `conditioning_intervened` and `conditioning_remaining` record
#'   which conditioning variables IDC moved into the intervention set and which
#'   remain conditioned on. The fields `direct_arrow` and `direct_path` refer to
#'   the graph obtained by latent projection. When ID fails, `hedge` holds the
#'   two graphs (`F` and `F_prime`) that demonstrate why the distribution is not
#'   identifiable, their common roots and the query at the step where ID failed.
#'   Roles have `variable`, `role` and `conditioning_guidance` columns.
#'   A variable may have different roles in different strategies. Instrument
#'   entries hold `instrument` and `conditioning`.
#'
#'   Adjustment and ordinary front door sets are lists of character vectors.
#'   Conditional front door sets pair mediator and covariate vectors as described
#'   above. `list()` contains no stored set. Consult the corresponding result
#'   to distinguish an unsuccessful check from a skipped one.
#'   `list(character())` is a valid empty set. These records and the nested
#'   strategies preserve variable names without parsing displayed text.
#'
#' @section Joint total effects:
#' With more than one distinct treatment or outcome, the query is the joint
#' distribution \eqn{P(Y \mid do(X))}: intervene simultaneously on every treatment
#' and consider every outcome together. It is not a collection of pairwise
#' results. Supply at least one observed treatment and one observed outcome,
#' with no node in both sets. Ordering and repeated names do not change the sets.
#'
#' Joint support is limited to unconditional total effect ID and covariate
#' adjustment. `effect = "direct"` and `effect = "indirect"` are rejected,
#' as are `conditioning`, `experiments` or `mediator` inputs other than
#' `NULL`. Front door, conditional front door, instrumental variable, direct,
#' natural, conditional (IDC) and experimental (z-ID) checks are not assessed
#' when there is more than one treatment or outcome. `summary(result)` includes
#' only strategy statuses within the selected total effect methods. A `not_assessed`
#' status is not a failure of identification. The total analysis records
#' `method_support` (`assessed`, `not_assessed`, `reason`) and
#' `query$joint = TRUE`. These extra fields are absent when there is only one
#' treatment and one outcome.
#' Variable roles record membership of complete adjustment sets and otherwise
#' give neutral guidance. Joint mediator and instrument roles are not inferred.
#'
#' @section Compatibility:
#' The earlier default checked only the total effect. Omitting `effect` or
#' setting it to `NULL` now gives an overview. Use `effect = "total"` to select
#' only total effects. The earlier table columns `strategy`, `feasible`,
#' `tier`, `variable_set` and `note` are replaced by the nested lists
#' described above. Use the conversion methods for a table.
#' Verdicts and structured variable sets remain available in the `analyses`
#' attribute. For example, use `attr(result, "analyses")$total$identifiable`
#' for the total effect identification result from the stated data, or
#' `attr(result, "analyses")$direct$adjustment_sets` for single door sets.
#' The stored `controlled` record is available only through this attribute.
#'
#' @seealso [identification_methods], [minimal_sets()]
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
#'   mediator = "growth", quiet = TRUE
#' )
#' mediation$direct
#' mediation$indirect
#' as.data.frame(mediation)
#'
#' @export
get_identification <- function(dag, treatment = NULL, outcome = NULL,
                               effect = c("total", "direct", "indirect"),
                               experiments = NULL, conditioning = NULL,
                               quiet = FALSE, mediator = NULL,
                               method = NULL,
                               id_status = c("identifiable", "all", "not_identifiable",
                                 "not_established", "not_assessed", "candidate",
                                 "zero_by_graph")){ # Assess identification and return the selected strategies
  if( missing(effect) ){
    effect <- NULL
  }
  if( !is.null(effect) && (!is.character(effect) || length(effect) != 1L ||
      is.na(effect) || !effect %in% c("total", "direct", "indirect")) ){
    stop("'effect' must be NULL for an overview, or one of 'total', 'direct' and 'indirect'.",
         call. = FALSE)
  }
  if( missing(id_status) ){
    id_status <- "identifiable"
  }
  id_status <- id_report_statuses(id_status)
  if( !is.logical(quiet) || length(quiet) != 1L || is.na(quiet) ){
    stop("Invalid input detected: only TRUE or FALSE may be supplied as 'quiet' input.",
         call. = FALSE)
  }
  if( !is.null(effect) ){
    effect <- as.character(effect)
  }
  joint <- length(treatment) > 1L || length(outcome) > 1L ||
    (is.null(treatment) && dagitty::is.dagitty(dag) &&
       length(dagitty::exposures(dag)) > 1L) ||
    (is.null(outcome) && dagitty::is.dagitty(dag) &&
       length(dagitty::outcomes(dag)) > 1L)
  if( joint ){
    query <- id_resolve_query(id_validate_graph(dag), treatment, outcome)
    treatment <- query$treatment
    outcome <- query$outcome
    joint <- length(treatment) > 1L || length(outcome) > 1L
    if( joint ){
      id_validate_joint_scope(effect, conditioning, experiments, mediator)
    }
  }
  scope <- if( !is.null(effect) ){
    effect
  }else if( joint || !is.null(conditioning) || !is.null(experiments) ){
    "total"
  }else{
    c("total", "direct", "indirect")
  }
  if( (joint || !is.null(conditioning) || !is.null(experiments)) &&
      any(scope != "total") ){
    stop("Joint endpoints, conditioning and experiments support only total effect queries.",
         call. = FALSE)
  }
  selection <- id_report_selection(scope, method, mediator)
  scope <- selection$effect
  selection$id_status <- id_status

  analyses <- list()
  if( "total" %in% scope ){
    total <- id_identification_query(
      dag, treatment, outcome, effect = "total", experiments = experiments,
      conditioning = conditioning, quiet = TRUE
    )
    analyses$total <- id_analysis_record(total)
    query <- analyses$total$query
    if( !isTRUE(query$joint) && !length(query$conditioning) ){
      analyses$conditional_frontdoor <- id_conditional_frontdoor(
        dag, id_admg(dag), query$treatment, query$outcome
      )
    }
  }
  if( any(c("direct", "indirect") %in% scope) ){
    controlled <- id_identification_query(
      dag, treatment, outcome, effect = "direct", experiments = experiments,
      conditioning = conditioning, quiet = TRUE, mediator = mediator
    )
    analyses$controlled <- id_analysis_record(controlled)
    analyses$mediation <- attr(controlled, "mediation", exact = TRUE)
    query <- analyses$controlled$query
    if( "direct" %in% scope ){
      analyses$direct <- id_direct_adjustment(dag, query$treatment, query$outcome)
    }
  }

  first_analysis <- if( !is.null(analyses$total) ){
    analyses$total
  }else{
    analyses$controlled
  }
  result <- id_identification_result(
    analyses, first_analysis$query$treatment, first_analysis$query$outcome,
    selection, dag
  )
  if( !quiet ){
    treatment_label <- first_analysis$query$treatment
    outcome_label <- first_analysis$query$outcome
    if( joint ){
      treatment_label <- id_display_nodes(treatment_label)
      outcome_label <- id_display_nodes(outcome_label)
    }
    message("Identification checks for ", treatment_label, " -> ", outcome_label,
            ": ", paste(scope, collapse = " and "),
            ". Verdicts and strategies are reported separately.")
  }
  return(result)
}


#' Understanding identification strategies and results
#'
#' A causal effect is identifiable if its value is uniquely determined by the
#' available observational or experimental distributions under the stated causal
#' assumptions. All causal models consistent with those assumptions and
#' distributions must agree on the effect.
#'
#' [get_identification()] checks identification strategies using the supplied
#' graph. It reports identification verdicts separately from the strategies
#' and their variable sets. It does not derive an identifying formula or
#' estimate an effect from data. This page describes the results, assumptions
#' and sources for those checks. The results do not rank estimators. The
#' supplied graph is assumed to represent the relevant causal structure,
#' including unobserved common causes.
#'
#' Total, direct and indirect effects form separate groups. Each strategy's
#' notes describe the effect or coefficient it concerns. A single door check
#' concerns a direct path coefficient under a linear model. A mediation formula
#' check concerns a natural direct or indirect effect relative to the selected
#' mediators. These checks have different assumptions and need not give the
#' same result.
#'
#' @section Result definitions:
#' \describe{
#'   \item{\code{identifiable}}{The stated method establishes identification of this
#'     effect or coefficient under the graph and required assumptions. For natural effects
#'     these include the structural and cross world assumptions described below.}
#'   \item{\code{not_identifiable}}{The supported identification algorithm proves that
#'     the requested distribution cannot be determined uniquely from the stated
#'     data under its assumptions. This does not rule out a
#'     different target or stronger model restrictions.}
#'   \item{\code{not_established}}{This restricted criterion or method does not
#'     establish identification. It is not a claim that every method fails.}
#'   \item{\code{not_assessed}}{This check was not requested or is not supported
#'     for the query. No identification verdict was calculated for this method.}
#'   \item{\code{candidate}}{An instrument passes the graphical search. The linear
#'     model and instrument relevance requirements described below still apply.
#'     This function does not test them in data.}
#'   \item{\code{zero_by_graph}}{The absence of the relevant causal path or original
#'     direct arrow forces the stated effect or coefficient to be zero.
#'     This is a structural implication, not an estimate from data. It must
#'     not be confused with an unidentified intervention distribution.}
#' }
#'
#' @section Variables and assumptions:
#' Each strategy stores variables in a named list. Roles include `covariates`,
#' `mediators`, `instrument` and `conditioning`. Verdicts can also record
#' `intervention`, `experiments`, `experiment_set`, `treatment` and `outcome`.
#' Each role contains a character vector of node names. Names are preserved
#' exactly. Each strategy contains one complete set of variables. Strategies
#' are not instructions to adjust for the union of their sets. An empty
#' `covariates` vector means no covariate adjustment is required for that
#' strategy. An unsuccessful search has no covariate set and a note explaining
#' the result. Natural effects are defined relative to the selected mediators.
#'
#' For conditional front door adjustment, `mediators` and `covariates` form
#' one complete strategy. The covariates are averaged over in the identifying
#' formula. They do not change the query into a conditional effect.
#'
#' Consistency means the observed response agrees with the response under its
#' corresponding well defined intervention. Positivity requires support for
#' the conditional distributions used by the identifying method. Appropriate
#' forms of these assumptions apply in addition to the graph. The `assumptions`
#' and `notes` fields highlight requirements and restrictions. They are not an
#' exhaustive checklist and an omitted restriction is not evidence that it holds.
#' No data are inspected.
#'
#' @section Identification verdicts:
#' The function uses the ID algorithm to check whether the outcome distribution
#' under treatment intervention is identifiable in a nonparametric causal model.
#' Nonparametric means no particular functional form, such as linear equations,
#' is imposed beyond the stated causal model. Conditional effects are checked
#' with IDC when `conditioning` is supplied. Those variables are conditioned on,
#' not automatically treated as an adjustment set or set by intervention.
#' Query records retain the requested conditioning, any variables IDC moves
#' into the intervention set and those that remain conditioned on. These
#' verdicts are separate from the strategy list. A positive ID verdict can
#' occur when none of the listed strategies succeeds. The current function
#' returns that verdict without deriving the identifying formula.
#' Sources: Shpitser and Pearl (2006), the joint intervention
#' paper for ID and the conditional intervention paper, Figure 3 and Theorems
#' 7--8, for IDC.
#'
#' @section Covariate adjustment:
#' The listed set satisfies the graphical adjustment criterion for the total
#' effect. An empty set means no covariate adjustment is required by this
#' criterion. Failure of adjustment does not imply failure of general
#' identification. Source: Perkovic et al. (2015), Definition 3.3 and Theorem 3.4.
#'
#' @section Front door adjustment:
#' The listed mediators satisfy the ordinary front door criterion for the total
#' effect. The search returns qualifying sets of the smallest size found for
#' this criterion. It does not implement every generalised front door method.
#' Source: Pearl (2009), Definition 3.3.3 and Theorem 3.3.4, pp. 82--83.
#'
#' @section Conditional front door adjustment:
#' Each strategy pairs observed mediators with observed covariates that are
#' ancestors of treatment or outcome and are not descendants of treatment.
#' The mediators intercept every directed path from treatment to outcome.
#' The covariates must be an adjustment set for treatment's effect on the
#' mediators. Treatment and those covariates must be an adjustment set for
#' the joint effect of the mediators on the outcome. Several mediators are
#' checked jointly. Failure of this check does not rule out other identifying
#' formulas.
#'
#' The search returns qualifying strategies with the smallest combined number
#' of mediators and covariates. A mediator set that passes without covariates
#' is reported under front door adjustment. The covariates are averaged over
#' in the identifying formula, so this remains a total effect query.
#' An unsuccessful search does not rule out other identification strategies.
#' Source: Xu et al. (2024), Definition 4, Theorem 3 and equation 3, with the
#' covariates unaffected by treatment as in Section 4.1.
#'
#' @section Single door adjustment:
#' This uses `dagitty::adjustmentSets(effect = "direct")` on the supplied DAG,
#' matching [minimal_sets()]. The target is the coefficient of the original
#' arrow from treatment to outcome in a linear structural equation model.
#' Linearity restricts this result. It is not required for every definition of
#' a direct effect and is not tested by the graph. With no original direct
#' arrow the coefficient is zero, and no single door search is needed.
#' An unsuccessful search does not exclude other methods for direct effects.
#' Source: Pearl (2009), Theorem 5.3.1, pp. 150--151.
#'
#' @section Instrumental variables:
#' Each strategy gives one candidate instrument and the covariates to
#' condition on when using it. The target is the linear total effect, under a
#' linear structural equation model. The instrument must also have nonzero
#' correlation with treatment after conditioning on those covariates.
#' The function does not test that correlation
#' or assess problems caused by weak instruments.
#' A candidate does not override a negative nonparametric identification
#' verdict: those results concern different assumptions. No candidate found
#' means this search did not establish a route.
#' Source: van der Zander et al. (2015), Definitions 3.2--3.3.
#'
#' @section Surrogate experiments:
#' Supplying `experiments` declares available data from interventions on those
#' variables. The z-identification check asks whether these data, together with
#' observational data, identify the total effect distribution. The record
#' distinguishes observational identification, declared experiments and a
#' qualifying experimental set. If observational identification already
#' succeeds, the declared experiments are not needed for that result.
#' Source: Bareinboim and Pearl (2012), Theorem 3.
#'
#' @section Natural direct and indirect effects:
#' A natural direct effect compares outcomes under different treatment levels
#' while keeping the selected mediators at the values they would have taken
#' under a reference treatment.
#' A natural indirect effect compares outcomes as the mediators change between
#' the values they would take under different treatment levels, while treatment
#' is held fixed. The mediation formula concerns these effects relative to the
#' selected mediator set. The same observed covariate set must pass adjustment
#' checks for the treatment's effect on the mediators and the joint effect of
#' treatment and mediators on the outcome. Eligible covariates cannot be reached
#' by following arrows from treatment. An unsuccessful check rules out this
#' criterion for the target, not all identification methods.
#'
#' A positive result requires a structural causal model with independent
#' background errors in an underlying DAG, including the latent common causes
#' represented by bidirected edges. This entails cross world assumptions,
#' which connect responses under different hypothetical treatments.
#' Consistency and positivity are also required. The graph cannot verify them.
#' Several selected mediators are checked jointly, not as separate contributions
#' to a decomposition. A natural direct effect can include paths through
#' mediators outside the selected set. The joint criterion can be conservative
#' when treatment affects a selected mediator only through others in the set.
#' A smaller set can pass when the joint set fails.
#' Source: Shpitser and VanderWeele (2011), Theorem 2.
#'
#' @section Stored parent intervention diagnostic:
#' The `analyses` attribute retains a `controlled` record when direct or indirect
#' effects are checked. It is not included in the nested strategies or table
#' conversions, and neither printing nor `summary()` displays it.
#' This is a diagnostic for the existing
#' parent intervention query, not a reported controlled direct effect.
#' Latent projection first produces a graph containing only observed nodes.
#' ID checks whether the outcome distribution is identifiable when treatment
#' and the outcome's other parents
#' in that graph are set by intervention. A parent is a node with an arrow
#' directly into the outcome. The other parents can include confounders and
#' predictors as well as mediators. Intervention is different from conditioning.
#' A controlled contrast also requires chosen treatment levels and fixed
#' parent values. The function checks identifiability without choosing those
#' values.
#' A zero contrast does not establish that the distribution itself is identified.
#'
#' Explicit latent nodes are not fixed. Latent projection can replace paths
#' through latent mediators by arrows, so this query need not represent a
#' controlled direct effect in the original supplied DAG. It remains distinct
#' from the direct path coefficient in the supplied graph and from natural
#' effects.
#' Source for distribution identification: Shpitser and Pearl (2006), ID.
#'
#' @references
#' ID: Shpitser I, Pearl J (2006). "Identification of Joint Interventional
#' Distributions in Recursive Semi-Markovian Causal Models." *AAAI*,
#' 1219--1226. <https://cdn.aaai.org/AAAI/2006/AAAI06-191.pdf>.
#'
#' IDC: Shpitser I, Pearl J (2006). "Identification of Conditional
#' Interventional Distributions." *UAI*, 437--444, Figure 3 and Theorems 7--8.
#' <https://arxiv.org/abs/1206.6876>.
#'
#' Adjustment: Perkovic E, Textor J, Kalisch M, Maathuis MH (2015).
#' "A Complete Generalized Adjustment Criterion." *UAI*, Definition 3.3
#' and Theorem 3.4. <https://arxiv.org/abs/1507.01524>.
#'
#' Front door and single door: Pearl J (2009). *Causality: Models,
#' Reasoning, and Inference*. Second edition. Cambridge University Press.
#' Front door: Definition 3.3.3 and Theorem 3.3.4, pp. 82--83.
#' Single door: Theorem 5.3.1, pp. 150--151.
#' \doi{10.1017/CBO9780511803161}.
#'
#' Conditional front door: Xu Z, Cheng D, Li J, Liu J, Liu L, Yu K (2024).
#' "Causal Inference with Conditional Front-Door Adjustment and Identifiable
#' Variational Autoencoder." ICLR, Definition 4, Theorem 3 and equation 3.
#' <https://proceedings.iclr.cc/paper_files/paper/2024/file/5186db2bfac87765a744547e82ac8848-Paper-Conference.pdf>.
#'
#' Conditional instruments: van der Zander B, Textor J, Liskiewicz M (2015).
#' "Efficiently Finding Conditional Instruments for Causal Inference."
#' *IJCAI*, 3243--3249, Definitions 3.2--3.3.
#' <https://www.ijcai.org/Proceedings/15/Papers/457.pdf>.
#'
#' z-ID: Bareinboim E, Pearl J (2012). "Causal Inference by Surrogate
#' Experiments: z-Identifiability." *UAI*, 113--120, Theorem 3.
#' <https://arxiv.org/abs/1210.4842>.
#'
#' Natural effects: Shpitser I, VanderWeele TJ (2011). "A Complete Graphical
#' Criterion for the Adjustment Formula in Mediation Analysis."
#' *The International Journal of Biostatistics*, **7**(1), Article 16,
#' Theorem 2. \doi{10.2202/1557-4679.1297}.
#'
#' @seealso [get_identification()], [minimal_sets()]
#' @name identification_methods
#' @keywords documentation
NULL


## Identification result methods ##

#' @export
print.causaliflower_identification <- function(x, ..., details = FALSE){ # Print separate verdicts and numbered strategies for each effect
  id_display_validate_details(details)
  report <- if( details && id_display_is_current(x) ){
    summary(x)
  }else{
    x
  }
  verdicts <- if( id_display_is_current(report) ){
    attr(report, "verdicts", exact = TRUE)
  }else{
    list()
  }
  cat(id_display_heading(report), "\n", sep = "")
  if( !length(report) ){
    cat("No effect groups remain in this report.\n")
  }
  for(effect in names(report)){
    cat("\n", id_report_label(effect), "\n", sep = "")
    for(verdict in verdicts[[effect]]){
      fields <- id_verdict_fields(verdict)
      for(field in names(fields)){
        lines <- strwrap(paste0(field, ": ", fields[[field]]),
                         width = max(40L, getOption("width") - 2L))
        cat(paste0("  ", lines), sep = "\n")
      }
      cat("\n")
    }
    cat("  Strategies\n")
    alternatives <- report[[effect]]
    if( !length(alternatives) ){
      cat("    No strategies match the selected methods and identification statuses.\n")
    }
    for(index in seq_along(alternatives)){
      cat("    [", index, "]\n", sep = "")
      fields <- id_report_fields(alternatives[[index]])
      for(field in names(fields)){
        lines <- strwrap(paste0(field, ": ", fields[[field]]),
                         width = max(40L, getOption("width") - 4L))
        cat(paste0("    ", lines), sep = "\n")
      }
      if( index < length(alternatives) ){
        cat("\n")
      }
    }
  }
  cat("\nMethods and assumptions: ?identification_methods\n")
  return(invisible(x))
}


#' @export
summary.causaliflower_identification <- function(object, ...){ # Return all stored strategy statuses without recalculating identification
  if( !id_display_is_current(object) ){
    return(data.table::copy(object))
  }
  selection <- attr(object, "selection", exact = TRUE)
  selection$id_status <- NULL
  result <- id_identification_result(
    attr(object, "analyses", exact = TRUE),
    attr(object, "treatment", exact = TRUE),
    attr(object, "outcome", exact = TRUE),
    selection, attr(object, "dag", exact = TRUE)
  )
  return(result)
}


#' @exportS3Method knitr::knit_print
knit_print.causaliflower_identification <- function(x, ..., details = FALSE){ # Format the same verdicts and strategies as labelled HTML blocks
  id_display_validate_details(details)
  if( !knitr::is_html_output() ){
    return(print(x, ..., details = details))
  }
  report <- if( details && id_display_is_current(x) ){
    summary(x)
  }else{
    x
  }
  verdicts <- if( id_display_is_current(report) ){
    attr(report, "verdicts", exact = TRUE)
  }else{
    list()
  }
  output <- paste0("<div class=\"causaliflower-identification\">\n<p><strong>",
                   id_display_escape(id_display_heading(report)), "</strong></p>")
  if( !length(report) ){
    output <- paste0(output, "\n<p>No effect groups remain in this report.</p>")
  }
  for(effect in names(report)){
    output <- paste0(output, "\n<h3>", id_display_escape(id_report_label(effect)), "</h3>")
    for(verdict in verdicts[[effect]]){
      fields <- id_verdict_fields(verdict)
      output <- paste0(output, "\n<dl class=\"causaliflower-verdict\">")
      for(field in names(fields)){
        output <- paste0(output, "\n<dt>", id_display_escape(field), "</dt><dd>",
                         id_display_escape(fields[[field]]), "</dd>")
      }
      output <- paste0(output, "\n</dl>")
    }
    output <- paste0(output, "\n<p><strong>Strategies</strong></p>")
    alternatives <- report[[effect]]
    if( !length(alternatives) ){
      output <- paste0(output, "\n<p>No strategies match the selected methods and ",
                       "identification statuses.</p>")
    }
    for(index in seq_along(alternatives)){
      output <- paste0(output, "\n<p>[", index, "]</p>\n<dl>")
      fields <- id_report_fields(alternatives[[index]])
      for(field in names(fields)){
        output <- paste0(output, "\n<dt>", id_display_escape(field), "</dt><dd>",
                         id_display_escape(fields[[field]]), "</dd>")
      }
      output <- paste0(output, "\n</dl>")
    }
  }
  output <- paste0(output, "\n<p>Methods and assumptions: ",
                   "<code>?identification_methods</code>.</p>\n</div>\n")
  return(knitr::asis_output(output))
}


#' @export
`[.causaliflower_identification` <- function(x, ...){ # Return independent selections without unrelated report metadata
  contents <- lapply(x, function(value){ return(value) })
  result <- data.table::copy(contents[...])
  return(result)
}


#' @export
as.data.frame.causaliflower_identification <- function(x, row.names = NULL,
                                                      optional = FALSE, ...){ # Flatten complete strategies without losing variable roles
  rows <- list()
  for(effect in names(x)){
    for(entry in x[[effect]]){
      row <- data.frame(effect = effect, method = entry$method,
        assumptions = paste(entry$assumptions, collapse = ". "),
        id_status = entry$id_status, notes = paste(entry$notes, collapse = ". "),
        stringsAsFactors = FALSE)
      row$variables <- I(list(data.table::copy(entry$variables)))
      row <- row[c("effect", "method", "assumptions", "variables", "id_status", "notes")]
      rows[[length(rows) + 1L]] <- row
    }
  }
  if( !length(rows) ){
    result <- data.frame(effect = character(), method = character(),
      assumptions = character(), id_status = character(), notes = character(),
      stringsAsFactors = FALSE)
    result$variables <- I(list())
    result <- result[c("effect", "method", "assumptions", "variables", "id_status", "notes")]
  }else{
    result <- do.call(rbind, rows)
  }
  if( !is.null(row.names) ){
    rownames(result) <- row.names
  }
  if( id_display_is_current(x) ){
    attr(result, "verdicts") <- data.table::copy(attr(x, "verdicts", exact = TRUE))
  }
  return(result)
}


#' @exportS3Method data.table::as.data.table
as.data.table.causaliflower_identification <- function(x, ...){ # Return a strategy table with separate verdicts in its attributes
  frame <- as.data.frame(x, ...)
  result <- data.table::as.data.table(frame)
  attr(result, "verdicts") <- data.table::copy(attr(frame, "verdicts", exact = TRUE))
  return(result)
}


## Identification for one query ##

id_identification_query <- function(dag, treatment = NULL, outcome = NULL,
                                    effect = c("total", "direct"),
                                    experiments = NULL, conditioning = NULL,
                                    quiet = FALSE, mediator = NULL){ # Retain the ID, IDC and mediation calculations for one target
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
  joint <- length(treatment) > 1L || length(outcome) > 1L
  if( joint ){
    id_validate_joint_scope(effect, conditioning, experiments, mediator)
  }

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
  direct_arrow <- any(G$D$from %in% treatment & G$D$to %in% outcome)
  mediation <- list()

  rows <- list()
  adjustment <- list(feasible = FALSE, sets = list())
  frontdoor <- list(feasible = FALSE, sets = list())
  instruments <- list(feasible = FALSE, entries = list())
  idc <- list(intervened = character(), remaining = character())
  z_result <- list(zid = NA, set = character())
  fixing_set <- character()
  intervention <- treatment

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
      isTRUE(solved$id),
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

    if( !joint ){
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
          "graphical instrument for a linear effect; assess linearity and instrument strength before estimation"
        }else{
          "no graphical instrument"
        }
      )
    }

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
    instruments,
    conditioning,
    idc$intervened
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
  attr(result, "query") <- list(
    treatment = treatment, outcome = outcome, effect = effect,
    intervention = intervention, fixing_set = fixing_set,
    conditioning = conditioning,
    conditioning_intervened = sort(idc$intervened),
    conditioning_remaining = sort(idc$remaining),
    direct_arrow = direct_arrow, direct_path = direct_path
  )
  if( joint ){
    attr(result, "query")$joint <- TRUE
    attr(result, "method_support") <- list(
      assessed = c("ID", "covariate adjustment"),
      not_assessed = c("front-door", "instrumental variable", "direct effects",
                       "mediation", "IDC", "z-ID"),
      reason = "Joint endpoint support is limited to unconditional total effect ID and covariate adjustment."
    )
  }
  attr(result, "frontdoor_sets") <- frontdoor$sets
  attr(result, "instruments") <- instruments$entries
  attr(result, "experiments") <- experiments
  attr(result, "experiment_set") <- z_result$set
  attr(result, "caveat") <- paste(
    "Identification is a graphical result. Estimation also requires",
    "positivity and data that support the comparisons used by the reported strategy."
  )

  attr(result, "display_rows") <- lapply(result, data.table::copy)
  if( !quiet ){
    id_summary_message(result)
  }

  return(result)
}


## Mediation checks ##
# The two adjustment checks follow Shpitser and VanderWeele (2011), Theorem 2.
# Failure of these checks does not establish failure of general identification.
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

id_validate_joint_scope <- function(effect, conditioning, experiments, mediator){ # Refuse unsupported joint methods instead of running pairwise checks
  if( identical(effect, "direct") || !is.null(conditioning) ||
      !is.null(experiments) || !is.null(mediator) ){
    stop("Joint endpoint queries support only unconditional total effect ID and covariate adjustment. Use effect = 'total' (or NULL), with conditioning, experiments and mediator all NULL. Joint direct, mediation, conditional and experimental methods are not assessed.",
         call. = FALSE)
  }
  return(invisible(NULL))
}
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
#' @param treatment Optional treatment name. `NULL` uses the graph's exposure
#'   tag.
#' @param outcome Optional outcome name. `NULL` uses the graph's outcome tag.
#' @return A list containing one observed `treatment` and one observed
#'   `outcome`.
#' @noRd
id_resolve_query <- function(graph, treatment, outcome){ # Resolve and validate the observed query endpoints
  if( is.null(treatment) ){
    treatment <- as.character(dagitty::exposures(graph$dag))
  }
  if( is.character(treatment) && length(treatment) == 0L ){
    stop("No treatment detected: 'get_identification()' requires exactly one named treatment, either tagged in the DAG or supplied through 'treatment'. Please ensure a single treatment is supplied and try again.",
         call. = FALSE)
  }
  if( is.character(treatment) && length(treatment) > 1L &&
      (anyNA(treatment) || any(!nzchar(treatment)) ||
       any(!treatment %in% graph$observed)) ){
    stop("Invalid treatment set detected: supply a nonempty character vector of observed node names in the DAG.",
         call. = FALSE)
  }
  if( !is.character(treatment) || length(treatment) == 0L ||
      anyNA(treatment) || any(!nzchar(treatment)) ||
      any(!treatment %in% graph$observed) ){
    stop("Invalid treatment detected: 'get_identification()' requires exactly one named treatment, either tagged in the DAG or supplied through 'treatment'. Please ensure a single treatment is supplied and try again.",
         call. = FALSE)
  }

  treatment <- sort(unique(as.character(treatment)))

  if( is.null(outcome) ){
    outcome <- as.character(dagitty::outcomes(graph$dag))
  }
  if( is.character(outcome) && length(outcome) == 0L ){
    stop("No outcome detected: 'get_identification()' requires exactly one named outcome, either tagged in the DAG or supplied through 'outcome'. Please ensure a single outcome is supplied and try again.",
         call. = FALSE)
  }
  if( is.character(outcome) && length(outcome) > 1L &&
      (anyNA(outcome) || any(!nzchar(outcome)) ||
       any(!outcome %in% graph$observed)) ){
    stop("Invalid outcome set detected: supply a nonempty character vector of observed node names in the DAG.",
         call. = FALSE)
  }
  if( !is.character(outcome) || length(outcome) == 0L ||
      anyNA(outcome) || any(!nzchar(outcome)) || any(!outcome %in% graph$observed) ){
    stop("Invalid outcome detected: 'get_identification()' requires exactly one named outcome, either tagged in the DAG or supplied through 'outcome'. Please ensure a single outcome is supplied and try again.",
         call. = FALSE)
  }
  outcome <- sort(unique(as.character(outcome)))
  if( length(treatment) == 1L && identical(treatment, outcome) ){
    stop("Identical treatment and outcome detected: 'get_identification()' requires the named treatment and outcome to refer to different nodes in the DAG. Please ensure different treatment and outcome nodes are supplied and try again.",
         call. = FALSE)
  }
  if( length(intersect(treatment, outcome)) ){
    stop("Overlapping treatment and outcome sets detected: joint endpoint queries require disjoint observed treatment and outcome sets.",
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


#' Create an empty table of directed edges
#'
#' @return A data frame with character `from` and `to` columns.
#' @noRd
id_empty_directed <- function(){
  edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)

  return(edges)
}


#' Create an empty table of bidirected edges
#'
#' @return A data frame with character `v` and `w` columns.
#' @noRd
id_empty_bidirected <- function(){
  edges <- data.frame(v = character(), w = character(), stringsAsFactors = FALSE)

  return(edges)
}


#' Remove repeated edges and edges from a node to itself
#'
#' @param arcs Table of directed edges with `from` and `to` columns, or `NULL`.
#' @return A sorted `from`, `to` data frame of unique edges between different nodes.
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
#' @param arcs Table of bidirected edges with `v` and `w` columns, or `NULL`.
#' @return A sorted `v`, `w` data frame of unique edges between different nodes.
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
#' @return A mixed graph list with observed vertices `V`, directed edges `D`,
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param vertices Vertex names to retain.
#' @return A mixed graph list containing the selected vertices and their
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
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


#' Create keys for comparing directed edges
#'
#' @param arcs Table of directed edges with `from` and `to` columns.
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


#' Create keys for comparing bidirected edges
#'
#' @param arcs Table of bidirected edges with `v` and `w` columns.
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @return Sorted node names with no outgoing directed edges.
#' @noRd
id_roots <- function(G){
  roots <- sort(setdiff(G$V, unique(G$D$from)))

  return(roots)
}


#' Check a proposed hedge against its graph and query
#'
#' @param F Larger forest in a mixed graph.
#' @param F_prime Smaller forest in a mixed graph.
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
#' @return A node name in double quotes with embedded double quotes escaped.
#' @noRd
id_quote_node <- function(node){
  text <- paste0("\"", gsub("\"", "\\\"", node, fixed = TRUE), "\"")

  return(text)
}


#' Convert a mixed graph to a dagitty object
#'
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
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

id_direct_adjustment <- function(dag, treatment, outcome){ # Check single door adjustment for the supplied graph's direct arrow
  graph <- id_validate_graph(dag)
  edges <- dagitty::edges(graph$dag)
  direct_arrow <- any(edges$v == treatment & edges$w == outcome & edges$e == "->")
  sets <- list()
  if( direct_arrow ){
    sets <- dagitty::adjustmentSets(
      graph$dag, exposure = treatment, outcome = outcome,
      type = "minimal", effect = "direct"
    )
    sets <- lapply(sets, function(nodes){
      return(sort(as.character(nodes)))
    })
    sets <- unname(sets[order(lengths(sets))])
  }
  feasible <- NA
  status <- "structural_zero"
  if( direct_arrow ){
    feasible <- length(sets) > 0L
    status <- if( length(sets) ){
      "available"
    }else{
      "not_found"
    }
  }
  result <- list(
    treatment = treatment, outcome = outcome, direct_arrow = direct_arrow,
    feasible = feasible, status = status, adjustment_sets = sets
  )
  return(result)
}
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
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
  # A set has no open back door path from the treatment exactly when each of
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
  # A member with a back door path to the outcome that is open given the
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
  # Every front door set lies within the remaining candidates, so if these
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


#' Find conditional front door mediator and covariate sets
#'
#' @param dag A dagitty object.
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param treatment Observed treatment node name.
#' @param outcome Observed outcome node name.
#' @return A list with a `feasible` flag and complete paired `sets` of mediators
#'   and covariates. The search returns every pair with the smallest total
#'   number of variables satisfying these checks. Covariates must be observed
#'   ancestors of treatment or outcome and cannot be descendants of treatment.
#'   Mediator sets that qualify without covariates belong to ordinary front door
#'   adjustment and are not repeated here.
#' @noRd
id_conditional_frontdoor <- function(dag, G, treatment, outcome){
  result <- list(feasible = FALSE, sets = list())
  if( any(G$D$from == treatment & G$D$to == outcome) ){
    return(result)
  }

  descendants <- Filter(
    function(node){
      return(id_directed_reachable(G, treatment, node))
    }, setdiff(G$V, treatment)
  )
  mediator_candidates <- sort(intersect(
    descendants,
    setdiff(id_ancestors(G, outcome), c(treatment, outcome))
  ))
  covariate_candidates <- sort(setdiff(
    id_ancestors(G, c(treatment, outcome)),
    c(treatment, outcome, descendants)
  ))
  if( !length(mediator_candidates) || !length(covariate_candidates) ){
    return(result)
  }

  # A joint mediator set cannot qualify unless treatment's effect on each
  # member can be identified by covariate adjustment.
  mediator_candidates <- Filter(
    function(node){
      sets <- dagitty::adjustmentSets(
        dag, exposure = treatment, outcome = node, type = "canonical"
      )
      return(length(sets) > 0L)
    }, mediator_candidates
  )
  if( !length(mediator_candidates) ){
    return(result)
  }

  # A back door path that avoids all other candidate mediators must be
  # blocked by pretreatment covariates. Remove a node when no allowed set
  # can do this, then repeat because its removal can rule out other nodes.
  repeat {
    blocked <- Filter(
      function(node){
        others <- setdiff(mediator_candidates, node)
        H <- id_induced(G, setdiff(G$V, others))
        H$D <- H$D[H$D$from != node, , drop = FALSE]
        ancestors <- id_ancestors(H, c(node, treatment, outcome))
        covariates <- intersect(covariate_candidates, ancestors)
        return(!isTRUE(dagitty::dseparated(
          id_as_dagitty(H), X = node, Y = outcome,
          Z = c(treatment, covariates)
        )))
      }, mediator_candidates
    )
    if( !length(blocked) ){
      break
    }
    mediator_candidates <- setdiff(mediator_candidates, blocked)
  }
  if( !length(mediator_candidates) ){
    return(result)
  }

  # Removing every eligible mediator must interrupt all directed causal paths.
  if( id_directed_reachable(
    id_induced(G, setdiff(G$V, mediator_candidates)), treatment, outcome
  ) ){
    return(result)
  }

  # A set made entirely of these nodes needs no covariates. In particular,
  # ordinary front door chains should not trigger a search for extra variables.
  ordinary_candidates <- Filter(
    function(node){
      return(isTRUE(dagitty::isAdjustmentSet(
        dag, Z = character(), exposure = treatment, outcome = node
      )) && isTRUE(dagitty::isAdjustmentSet(
        dag, Z = treatment, exposure = node, outcome = outcome
      )))
    }, mediator_candidates
  )
  if( setequal(ordinary_candidates, mediator_candidates) ){
    return(result)
  }

  # Search by the total number of measured variables. Prepare each mediator
  # set once and reuse its covariate search at successive sizes.
  searches <- vector("list", length(mediator_candidates))
  largest_size <- length(mediator_candidates) + length(covariate_candidates)
  for(size in seq.int(2L, largest_size)){
    for(mediator_size in seq_len(min(length(mediator_candidates), size - 1L))){
      covariate_size <- size - mediator_size
      if( covariate_size > length(covariate_candidates) ){
        next
      }
      if( is.null(searches[[mediator_size]]) ){
        entries <- list()
        mediator_sets <- utils::combn(mediator_candidates, mediator_size, simplify = FALSE)
        for(mediators in mediator_sets){
          if( all(mediators %in% ordinary_candidates) ){
            next
          }
          if( id_directed_reachable(
            id_induced(G, setdiff(G$V, mediators)), treatment, outcome
          ) ){
            next
          }
          ordinary_frontdoor <- isTRUE(dagitty::isAdjustmentSet(
            dag, Z = character(), exposure = treatment, outcome = mediators
          )) && isTRUE(dagitty::isAdjustmentSet(
            dag, Z = treatment, exposure = mediators, outcome = outcome
          ))
          if( ordinary_frontdoor ){
            next
          }
          covariates <- id_frontdoor_covariate_search(
            dag, G, treatment, outcome, mediators, covariate_candidates
          )
          if( !is.null(covariates) ){
            entries[[length(entries) + 1L]] <- list(
              mediators = mediators, covariates = covariates
            )
          }
        }
        searches[[mediator_size]] <- entries
      }
      for(entry in searches[[mediator_size]]){
        for(covariates in entry$covariates(covariate_size)){
          result$sets[[length(result$sets) + 1L]] <- list(
            mediators = entry$mediators, covariates = covariates
          )
        }
      }
    }
    if( length(result$sets) ){
      result$feasible <- TRUE
      break
    }
  }
  return(result)
}


#' Prepare a covariate search for fixed front door mediators
#'
#' @param dag A dagitty object.
#' @param G Mixed graph of observed nodes.
#' @param treatment Observed treatment node name.
#' @param outcome Observed outcome node name.
#' @param mediators A candidate mediator set.
#' @param candidates Allowed covariates, none descendants of treatment.
#' @return A function returning all qualifying sets of a requested size, or
#'   NULL when an adjustment constraint cannot be satisfied.
#' @noRd
id_frontdoor_covariate_search <- function(dag, G, treatment, outcome, mediators, candidates){
  force(mediators)
  treatment_graph <- G
  treatment_graph$D <- G$D[G$D$from != treatment, , drop = FALSE]
  mediator_graph <- G
  mediator_graph$D <- G$D[!G$D$from %in% mediators, , drop = FALSE]
  treatment_dag <- id_as_dagitty(treatment_graph)
  mediator_dag <- id_as_dagitty(mediator_graph)

  can_complete <- function(required, allowed){ # Reject a branch with no possible separator
    treatment_ancestors <- id_ancestors(treatment_graph, c(treatment, mediators, required))
    treatment_covariates <- intersect(allowed, treatment_ancestors)
    if( !isTRUE(dagitty::dseparated(
      treatment_dag, X = treatment, Y = mediators, Z = treatment_covariates
    )) ){
      return(FALSE)
    }
    mediator_ancestors <- id_ancestors(mediator_graph, c(mediators, outcome, treatment, required))
    mediator_covariates <- intersect(allowed, mediator_ancestors)
    return(isTRUE(dagitty::dseparated(
      mediator_dag, X = mediators, Y = outcome,
      Z = c(treatment, mediator_covariates)
    )))
  }

  # A separator, if one exists within the allowed nodes, exists within the
  # ancestors of its endpoints and required nodes. Adding the other allowed
  # ancestors preserves separation in the corresponding moral graph.
  if( !can_complete(character(), candidates) ){
    return(NULL)
  }
  required <- Filter(function(node){
    return(!can_complete(character(), setdiff(candidates, node)))
  }, candidates)
  optional <- setdiff(candidates, required)

  search <- function(included, remaining, size){ # Enumerate only branches that can still qualify
    if( length(included) > size || length(included) + length(remaining) < size ){
      return(list())
    }
    if( length(included) == size ){
      covariates <- sort(included)
      if( isTRUE(dagitty::isAdjustmentSet(
        dag, Z = covariates, exposure = treatment, outcome = mediators
      )) && isTRUE(dagitty::isAdjustmentSet(
        dag, Z = c(treatment, covariates), exposure = mediators, outcome = outcome
      )) ){
        return(list(covariates))
      }
      return(list())
    }
    if( !can_complete(included, c(included, remaining)) ){
      return(list())
    }
    node <- remaining[[1L]]
    remaining <- remaining[-1L]
    with_node <- search(c(included, node), remaining, size)
    without_node <- search(included, remaining, size)
    return(c(with_node, without_node))
  }

  find_sets <- function(size){ # Reuse this mediator set's constraints at each requested size
    if( size < max(1L, length(required)) || size > length(candidates) ){
      return(list())
    }
    return(search(required, optional, size))
  }
  return(find_sets)
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
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

id_analysis_record <- function(result){ # Keep the structured calculation records for a reported target
  fields <- c("query", "identifiable", "observationally_identifiable", "z_identifiable",
    "data_availability", "hedge", "roles", "adjustment_sets", "canonical_adjustment",
    "frontdoor_sets", "instruments", "experiments", "experiment_set", "caveat") # Stored calculation attributes
  record <- lapply(fields, function(field){
    return(attr(result, field, exact = TRUE))
  })
  names(record) <- fields
  support <- attr(result, "method_support", exact = TRUE)
  if( !is.null(support) ){
    record$method_support <- support
  }
  return(record)
}


id_query_label <- function(query){ # Describe the intervention and conditioning sets in a distribution query
  intervention <- query$intervention
  if( is.null(intervention) ){
    intervention <- query$treatment
  }
  intervention <- c(intersect(query$treatment, intervention),
                    setdiff(intervention, query$treatment))
  given <- character()
  if( length(intervention) ){
    given <- paste0("do(", id_display_nodes(intervention), ")")
  }
  if( length(query$conditioning) ){
    given <- c(given, id_display_nodes(query$conditioning))
  }
  given_label <- ""
  if( length(given) ){
    given_label <- paste0(" | ", paste(given, collapse = ", "))
  }
  label <- paste0("P(", id_display_nodes(query$outcome), given_label, ")")
  return(label)
}


id_display_nodes <- function(nodes){ # Quote node names only when they need delimiters
  rendered <- vapply(as.character(nodes), function(node){
    if( grepl("^[A-Za-z][A-Za-z0-9_.]*$", node) ){
      return(node)
    }
    return(encodeString(node, quote = '"'))
  }, character(1L))
  label <- paste(rendered, collapse = ", ")
  return(label)
}


id_report_statuses <- function(id_status){ # Validate status filters and retain structural zeros with identified results
  choices <- c("identifiable", "not_identifiable", "not_established",
               "not_assessed", "candidate", "zero_by_graph") # Stored identification statuses
  if( is.null(id_status) ){
    return(NULL)
  }
  if( !is.character(id_status) || !length(id_status) || anyNA(id_status) ||
      any(!id_status %in% c(choices, "all")) ){
    stop("'id_status' must be NULL, 'all', or one or more identification statuses.",
         call. = FALSE)
  }
  if( "all" %in% id_status ){
    return(NULL)
  }
  id_status <- unique(unname(id_status))
  if( "identifiable" %in% id_status ){
    id_status <- union(id_status, "zero_by_graph")
  }
  return(id_status)
}


id_report_selection <- function(scope, method, mediator){ # Resolve compatible effect and strategy selections
  method_choices <- c("adjustment", "front_door", "conditional_front_door",
    "instrumental_variable", "single_door", "mediation_formula") # Supported strategies
  if( !is.null(method) && (!is.character(method) || !length(method) ||
      anyNA(method) || any(!method %in% method_choices)) ){
    stop("'method' must be NULL or one or more supported strategy names.", call. = FALSE)
  }
  if( is.null(method) ){
    method <- method_choices
  }else{
    method <- unique(unname(method))
  }
  compatible <- vapply(scope, function(effect){
    strategies <- switch(effect,
      total = c("adjustment", "front_door", "conditional_front_door", "instrumental_variable"),
      direct = c("single_door", "mediation_formula"),
      indirect = "mediation_formula"
    )
    return(any(method %in% strategies))
  }, logical(1L))
  scope <- scope[compatible]
  if( !length(scope) ){
    stop("The selected methods do not apply to the requested effects.", call. = FALSE)
  }
  if( !is.null(mediator) &&
      (!any(c("direct", "indirect") %in% scope) || !"mediation_formula" %in% method) ){
    stop("'mediator' requires a natural direct or natural indirect effect selection.",
         call. = FALSE)
  }
  return(list(effect = scope, method = method))
}


id_identification_result <- function(analyses, treatment, outcome, selection, dag){ # Attach the calculations and separate verdicts to the strategy report
  result <- id_identification_report(analyses, selection)
  attr(result, "analyses") <- data.table::copy(analyses)
  attr(result, "verdicts") <- id_identification_verdicts(analyses, selection)
  attr(result, "treatment") <- treatment
  attr(result, "outcome") <- outcome
  attr(result, "scope") <- selection$effect
  attr(result, "effect") <- if( length(selection$effect) > 1L ){
    "overview"
  }else{
    selection$effect
  }
  attr(result, "selection") <- data.table::copy(selection)
  attr(result, "dag") <- data.table::copy(dag)
  attr(result, "display_rows") <- data.table::copy(lapply(result, function(value){
    return(value)
  }))
  class(result) <- "causaliflower_identification"
  return(result)
}


id_identification_verdicts <- function(analyses, selection){ # Keep algorithm verdicts and structural implications outside the strategies
  result <- stats::setNames(vector("list", length(selection$effect)), selection$effect)
  distribution_assumptions <- c("Nonparametric causal model", "Consistency", "Positivity") # Requirements beyond the supplied graph
  add <- function(effect, method, query, id_status, assumptions, variables, notes){
    result[[effect]][[method]] <<- list(method = method, query = query,
      id_status = id_status, assumptions = assumptions, variables = variables, notes = notes)
    return(invisible(NULL))
  }
  status <- function(value){
    return(if( isTRUE(value) ){ "identifiable" }else{ "not_identifiable" })
  }
  total <- analyses$total
  if( !is.null(total) ){
    query <- total$query
    variables <- list(intervention = query$treatment)
    if( length(query$conditioning) ){
      variables$conditioning <- query$conditioning
      add("total", "idc", query, status(total$identifiable), distribution_assumptions,
          variables, "Uses observational data. No identifying formula is returned")
    }else{
      add("total", "id", query, status(total$observationally_identifiable),
          distribution_assumptions, variables,
          "Uses observational data. No identifying formula is returned")
      if( identical(query$direct_path, FALSE) ){
        add("total", "graph_structure", NULL, "zero_by_graph", character(), variables,
            "No directed path from treatment to outcome. Treatment does not change the outcome distribution")
      }
      if( length(total$experiments) ){
        experimental_variables <- list(experiments = total$experiments)
        if( isTRUE(total$observationally_identifiable) ){
          experimental_variables$experiment_set <- character()
          note <- "Observational data already identify this distribution. The declared experiments are not needed"
        }else if( isTRUE(total$z_identifiable) ){
          experimental_variables$experiment_set <- total$experiment_set
          note <- "Uses observational data and the qualifying experimental data. No identifying formula is returned"
        }else{
          note <- "The stated observational and experimental data do not identify this distribution"
        }
        add("total", "z_id", query, status(total$z_identifiable), distribution_assumptions,
            experimental_variables, note)
      }
    }
  }
  direct <- analyses$direct
  if( !is.null(direct) &&
      (identical(direct$status, "structural_zero") || identical(direct$direct_arrow, FALSE)) ){
    add("direct", "graph_structure", NULL, "zero_by_graph", "Linear structural equations",
        list(treatment = direct$treatment, outcome = direct$outcome),
        "No direct arrow from treatment to outcome in the supplied graph. Its direct path coefficient is zero")
  }
  return(result)
}


id_identification_report <- function(analyses, selection){ # Keep complete strategy records under their effect
  result <- stats::setNames(rep(list(list()), length(selection$effect)), selection$effect)
  add <- function(effect, method, id_status, variables, assumptions,
                  notes = character()){ # Append one matching strategy without combining alternatives
    if( !effect %in% selection$effect || !method %in% selection$method ||
        (!is.null(selection$id_status) && !id_status %in% selection$id_status) ){
      return(invisible(NULL))
    }
    entry <- list(method = method, assumptions = assumptions, variables = variables,
                  id_status = id_status, notes = notes)
    result[[effect]][[length(result[[effect]]) + 1L]] <<- entry
    return(invisible(NULL))
  }
  base_assumptions <- c("Consistency", "Positivity") # Shared requirements for distribution checks
  direct_assumptions <- "Linear structural equations" # Restriction of single door identification
  natural_assumptions <- c("Cross world structural assumptions", base_assumptions) # Natural effect requirements
  iv_assumptions <- c("Linear structural equations", "Nonzero instrument relevance") # Not verified from data
  total <- analyses$total
  if( !is.null(total) ){
    query <- total$query
    joint <- isTRUE(query$joint)
    if( length(query$conditioning) ){
      for(method in c("adjustment", "front_door", "conditional_front_door")){
        add("total", method, "not_assessed", list(), base_assumptions,
            "Only IDC is assessed for a conditional query")
      }
      add("total", "instrumental_variable", "not_assessed", list(), iv_assumptions,
          "Instrumental variables are not assessed for a conditional query")
    }else{
      if( length(total$adjustment_sets) ){
        for(nodes in total$adjustment_sets){
          add("total", "adjustment", "identifiable", list(covariates = nodes),
              base_assumptions, "A sufficient adjustment set for the total effect")
        }
      }else{
        add("total", "adjustment", "not_established", list(), base_assumptions,
            "No qualifying adjustment set. Other identification methods may apply")
      }
      if( length(total$frontdoor_sets) ){
        for(nodes in total$frontdoor_sets){
          add("total", "front_door", "identifiable", list(mediators = nodes),
              base_assumptions, "The mediators satisfy the ordinary front door criterion for the total effect")
        }
      }else if( joint ){
        add("total", "front_door", "not_assessed", list(), base_assumptions,
            "Front door adjustment is not assessed for joint endpoints")
      }else{
        add("total", "front_door", "not_established", list(), base_assumptions,
            "No qualifying front door set. Other identification methods may apply")
      }
      conditional <- analyses$conditional_frontdoor
      if( !is.null(conditional) && length(conditional$sets) ){
        for(pair in conditional$sets){
          add("total", "conditional_front_door", "identifiable", pair, base_assumptions,
              "Covariates are averaged over in the identifying formula")
        }
      }else if( joint ){
        add("total", "conditional_front_door", "not_assessed", list(), base_assumptions,
            "Conditional front door adjustment is not assessed for joint endpoints")
      }else{
        add("total", "conditional_front_door", "not_established", list(), base_assumptions,
            "No qualifying pair requiring covariates. Sets needing no covariates are reported under front door adjustment")
      }
      if( length(total$instruments) ){
        for(entry in total$instruments){
          add("total", "instrumental_variable", "candidate",
              list(instrument = entry$instrument, conditioning = entry$conditioning),
              iv_assumptions, "The graphical search found this instrument. Its relevance has not been checked in data")
        }
      }else if( joint ){
        add("total", "instrumental_variable", "not_assessed", list(), iv_assumptions,
            "Instrumental variables are not assessed for joint endpoints")
      }else{
        add("total", "instrumental_variable", "not_established", list(), iv_assumptions,
            "No qualifying instrument found by this search")
      }
    }
  }
  direct <- analyses$direct
  if( !is.null(direct) && isTRUE(direct$direct_arrow) ){
    if( isTRUE(direct$feasible) ){
      for(nodes in direct$adjustment_sets){
        add("direct", "single_door", "identifiable", list(covariates = nodes),
            direct_assumptions, "Concerns the coefficient on the original treatment to outcome arrow")
      }
    }else{
      add("direct", "single_door", "not_established", list(), direct_assumptions,
          "No qualifying single door set for the original arrow coefficient. Other identification methods may apply")
    }
  }
  for(entry in analyses$mediation){
    variables <- list(mediators = entry$mediators)
    if( isTRUE(entry$feasible) ){
      for(nodes in entry$adjustment_sets){
        variables$covariates <- nodes
        add("direct", "mediation_formula", "identifiable", variables, natural_assumptions,
            "Natural direct effect relative to the selected mediator set")
        add("indirect", "mediation_formula", "identifiable", variables, natural_assumptions,
            "Natural indirect effect relative to the selected mediator set")
      }
    }else{
      add("direct", "mediation_formula", "not_established", variables, natural_assumptions,
          "No eligible covariate set passes both adjustment checks for the natural direct effect. Other identification methods may apply")
      add("indirect", "mediation_formula", "not_established", variables, natural_assumptions,
          "No eligible covariate set passes both adjustment checks for the natural indirect effect. Other identification methods may apply")
    }
  }
  if( !is.null(analyses$controlled) && !length(analyses$mediation) ){
    note <- "No observed mediator is available on a directed path from treatment to outcome"
    add("direct", "mediation_formula", "not_assessed", list(), natural_assumptions, note)
    add("indirect", "mediation_formula", "not_assessed", list(), natural_assumptions, note)
  }
  return(result)
}


id_display_is_current <- function(x){ # Reuse stored analyses only while the selected strategies remain unchanged
  snapshot <- attr(x, "display_rows", exact = TRUE)
  contents <- lapply(x, function(value){ return(value) })
  current <- inherits(x, "causaliflower_identification") &&
    is.list(attr(x, "analyses", exact = TRUE)) &&
    is.list(attr(x, "selection", exact = TRUE)) &&
    is.list(snapshot) && identical(contents, snapshot)
  return(current)
}


id_report_label <- function(value){ # Render stable identifiers as readable labels
  label <- switch(value,
    total = "Total effect",
    direct = "Direct effect",
    indirect = "Indirect effect",
    graph_structure = "Graph structure",
    id = "ID",
    idc = "IDC",
    adjustment = "Covariate adjustment",
    front_door = "Front door adjustment",
    conditional_front_door = "Conditional front door adjustment",
    instrumental_variable = "Instrumental variable",
    z_id = "z-ID",
    single_door = "Single door adjustment",
    mediation_formula = "Mediation formula",
    identifiable = "Identifiable",
    not_identifiable = "Not identifiable",
    not_established = "Not established",
    not_assessed = "Not assessed",
    candidate = "Candidate",
    zero_by_graph = "Zero by graph",
    covariates = "Covariates",
    mediators = "Mediators",
    instrument = "Instrument",
    conditioning = "Conditioning",
    intervention = "Intervention",
    experiments = "Declared experiments",
    experiment_set = "Experiments needed",
    treatment = "Treatment",
    outcome = "Outcome",
    value
  )
  return(label)
}


id_report_fields <- function(entry){ # Format one strategy without changing its stored values
  variables <- character()
  for(role in names(entry$variables)){
    variables <- c(variables, paste0(id_report_label(role), ": ",
                                     id_display_variables(entry$variables[[role]])))
  }
  if( !length(variables) ){
    variables <- "None reported"
  }
  notes <- if( length(entry$notes) ){
    paste(entry$notes, collapse = ". ")
  }else{
    "None"
  }
  fields <- list(
    Method = id_report_label(entry$method),
    Assumptions = paste(entry$assumptions, collapse = ". "),
    Variables = paste(variables, collapse = ". "),
    "Identification status" = id_report_label(entry$id_status),
    Notes = notes
  )
  return(fields)
}


id_verdict_fields <- function(entry){ # Describe one query verdict without presenting it as a strategy
  label <- if( identical(entry$method, "graph_structure") ){
    "Graph structure"
  }else{
    paste0("Identification check (", id_report_label(entry$method), ")")
  }
  fields <- stats::setNames(list(id_report_label(entry$id_status)), label)
  if( !is.null(entry$query) ){
    fields$Query <- id_query_label(entry$query)
  }
  details <- id_report_fields(entry)
  if( length(entry$assumptions) ){
    fields$Assumptions <- details$Assumptions
  }
  if( identical(entry$method, "z_id") ){
    fields$Variables <- details$Variables
  }
  fields$Notes <- details$Notes
  return(fields)
}


id_display_validate_details <- function(details){ # Validate the option to print the full summary
  if( !is.logical(details) || length(details) != 1L || is.na(details) ){
    stop("'details' must be TRUE or FALSE.", call. = FALSE)
  }
  return(invisible(NULL))
}


id_display_heading <- function(x){ # Name the treatment and outcome sets
  heading <- paste0("Identification: ", id_display_nodes(attr(x, "treatment", exact = TRUE)),
    " -> ", id_display_nodes(attr(x, "outcome", exact = TRUE)))
  return(heading)
}


id_display_variables <- function(nodes){ # Distinguish a valid empty set from named variables
  if( !length(nodes) ){
    return("none needed")
  }
  return(id_display_nodes(nodes))
}


id_display_escape <- function(text){ # Escape each label and value once for HTML output
  text <- gsub("&", "&amp;", text, fixed = TRUE)
  text <- gsub("<", "&lt;", text, fixed = TRUE)
  text <- gsub(">", "&gt;", text, fixed = TRUE)
  text <- gsub("\"", "&quot;", text, fixed = TRUE)
  escaped <- gsub("'", "&#39;", text, fixed = TRUE)
  return(escaped)
}


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
#' @param sets List of sets of node names.
#' @param max_sets Maximum number of sets to display. `Inf` displays all sets.
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
#' @param variable_set Printed description of the variable set.
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
#' @param G Mixed graph of observed nodes as a list of vertices `V`, directed edges
#'   `D`, and bidirected edges `B`.
#' @param treatment Observed treatment node name.
#' @param outcome Observed outcome node name.
#' @param effect Effect type, `"total"` or `"direct"`.
#' @param adjustment Covariate adjustment result containing `sets`.
#' @param frontdoor Front door result containing mediator `sets`.
#' @param instruments Instrument result containing `entries`.
#' @return A data.table with `variable`, `role`, and `conditioning_guidance`
#'   columns.
#' @noRd
id_roles <- function(G, treatment, outcome, effect, adjustment,
                      frontdoor, instruments, conditioning = character(),
                      intervened = character()){ # Describe each observed variable's role in the requested analysis
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

  if( length(conditioning) > 0L ){
    add_role(
      intersect(conditioning, intervened),
      "conditioning variable",
      "conditioned on in the requested query. IDC permits an equivalent intervention on this variable"
    )
    add_role(
      setdiff(conditioning, intervened),
      "conditioning variable",
      "remains conditioned on in the IDC query. This is not a recommendation to adjust for it for the total effect"
    )
  }else if( identical(effect, "total") &&
            length(treatment) == 1L && length(outcome) == 1L ){
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
  }else if( !identical(effect, "total") ){
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
#' @param result The identification result with its resolved query.
#' @return Invisible `NULL`.
#' @noRd
id_summary_message <- function(result){ # Describe the supported distribution query and its verdict
  query <- attr(result, "query")
  verdict <- if( isTRUE(attr(result, "observationally_identifiable")) ){
    "identified from observational data under this graph."
  }else if( isTRUE(attr(result, "identifiable")) ){
    paste0("identified using surrogate experimental data for ",
           id_render_set(attr(result, "experiment_set")),
           ". Not identified from observational data alone.")
  }else{
    "not identified from the declared data under this graph alone."
  }
  message(id_query_label(query), ": nonparametric distribution ", verdict)
  if( identical(query$effect, "direct") && !query$direct_arrow ){
    message("No direct arrow in the graph after latent projection: the controlled contrast is zero.")
  }else if( identical(query$effect, "total") && !query$direct_path &&
            length(query$conditioning) == 0L ){
    message("No directed path: the total causal contrast is zero under this graph.")
  }
  return(invisible(NULL))
}
