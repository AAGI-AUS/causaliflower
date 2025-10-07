#' LM with dagitty
#'
#' causalify_lm()
#'
#' @importFrom data.table data.table fcase
#' @param dag dagitty object
#' @param data A data frame or data table, with variables corresponding to nodes in the supplied dag.
#' @return Model object of either
#' @export
causalify_lm <- function( dag, data ){
  .datatable.aware <- TRUE

      setClass("definition",
               slots = list(outcome = "character",
                            treatment = "character",
                            vars = "character",
                            data = "data.table",
                            model = "function",
                            lm = "list"),
               where = causali.envr
      )

      setMethod("initialize", "definition",
                function(.Object, outcome, treatment, vars, data, model) {
                  fo <- as.formula(paste(outcome, "~", treatment, paste(" + ", vars, collapse = " "), collapse = " "))
                  model_lm <- model(fo,
                                    data  = data)
                  .Object@lm <- list(model_lm)
                  .Object
                },
                where = causali.envr
      )


      causali_lm <- with(causali.envr, function(dag, data){

        feature_map_list <- getFeatureMap(dag)

        lm <- new("definition",
                  outcome = as.character(feature_map_list$outcome$ancestor),
                  treatment = as.character(feature_map_list$treatment$ancestor),
                  vars = as.character(feature_map_list$confounder$ancestor),
                  data = tab_work,
                  model = lm)

        return(lm)

      } )

      lm_out <- causali_lm(dag, tab_work)

      return(lm_out)
}


#' GLM with dagitty
#'
#' causalify_glm()
#'
#' @importFrom data.table data.table fcase
#' @importFrom mgcv gam gamm
#' @importFrom nlme corAR1
#' @param dag dagitty object
#' @param model Current options for model include "lm", "glm", "lme", or "gamm". The last two options are taken from the 'mgcv' package, with "lme" returning the linear mixed-effect model output from mgcv::gamm().
#' @param data A data frame or data table, with variables corresponding to nodes in the supplied dag.
#' @param family Used for glm and gamm models, defaults to gaussion().
#' @return Model object of either
#' @export
causalify_glm <- function( dag, data, family = gaussion() ){
  .datatable.aware <- TRUE

  # glm
  setClass("definition",
           slots = list(outcome = "character",
                        treatment = "character",
                        vars = "character",
                        data = "data.table",
                        family = "function",
                        model = "function",
                        glm = "list"),
           where = causali.envr
  )

  setMethod("initialize", "definition",
            function(.Object, outcome, treatment, vars, data, family, model) {
              fo <- as.formula(paste(outcome, "~", treatment, paste(" + ", vars, collapse = " "), collapse = " "))
              model_glm <- model(fo,
                                 data  = data,
                                 family = family)
              .Object@glm <- list(model_glm)
              .Object
            },
            where = causali.envr
  )


  causali_glm <- with(causali.envr, function(dag, data, family){

    feature_map_list <- getFeatureMap(dag)

    glm <- new("definition",
               outcome = as.character(feature_map_list$outcome$ancestor),
               treatment = as.character(feature_map_list$treatment$ancestor),
               vars = as.character(feature_map_list$confounder$ancestor),
               data = tab_work,
               family = gaussian(),
               model = glm)

    return(glm)

  } )

  glm_out <- causali_glm(dag, tab_work)

  return(glm_out)

}


#' GAMM with dagitty
#'
#' causalify_gamm()
#'
#' @importFrom data.table data.table fcase
#' @importFrom mgcv gam gamm
#' @importFrom nlme corAR1
#' @param dag dagitty object
#' @param data A data frame or data table, with variables corresponding to nodes in the supplied dag.
#' @param random_effects A vector or list of variables to be treated as random effects. Only used when "gamm" or "gam" is specified, the  are not required to be included in the dag.
#' @param method Defaults to method = "REML", only used for gamm models.
#' @param family Used for glm and gamm models, defaults to gaussion().
#' @param corr_structure Limited functionality. Currently only used when model = "lme" or "gamm" to specify nlme::corAR1(0, form=~1|year) if a random_effects input contains the string 'year'.
#' @return Model object of either
#' @export
causalify_lme <- function(dag, data, random_effects = NULL, method = "REML", family = gaussion(), corr_structure = "corAR1"){
  .datatable.aware <- TRUE

  # gamm / lme
  setClass("definition",
           slots = list(outcome = "character",
                        treatment = "character",
                        vars = "character",
                        data = "data.table",
                        random = "list",
                        correlation = "function",
                        family = "function",
                        method = "character",
                        model = "function",
                        gamm = "list"),
           where = causali.envr
  )

  setMethod("initialize", "definition",
            function(.Object, outcome, treatment, vars, data, random, correlation, family, method, model) {
              fo <- as.formula(paste(outcome, "~", treatment, paste(" + ", vars, collapse = " "), collapse = " "))
              model_gamm <- model(fo,
                                  data  = data,
                                  random = random,
                                  family = family,
                                  method = method)
              .Object@gamm <- list(model_gamm)
              .Object
            },
            where = causali.envr
  )


  {
    random_effects <- list("year", "loc", "rep")

    random_list <- lapply(seq_along( random_effects ), function(x){

      random_list <- list()

      random_list[x] <- c(~1)

    })

    random_list <- unlist(random_list)

    names(random_list) <- random_effects


    if( "year" %in% random_effects ){

      cor_var <- random_effects[ random_effects %in% "year" ]

      year <- function (x)
      {
        UseMethod(as.character(cor_var))
      }

      corr_structure <- nlme::corAR1(0, form=~1|year)

    }else{

      corr_structure <- NULL

    }

  }


  causalify_gamm <- with(causali.envr, function(dag, data, random, family, method, corr_structure){

    feature_map_list <- getFeatureMap(dag)

    gamm <- new("definition",
                outcome = as.character(feature_map_list$outcome$ancestor),
                treatment = as.character(feature_map_list$treatment$ancestor),
                vars = as.character(feature_map_list$confounder$ancestor),
                data = tab_work,
                random = random_list,
                family = family,
                method = method,
                correlation = cor_AR1,
                model = mgcv::gamm)

    return(gamm)

  } )

  lme_out <- causalify_gamm(dag, tab_work, random, family, method, corr_structure)$lme

  return(lme_out)

}


#' GAMM with dagitty
#'
#' causalify_gamm()
#'
#' @importFrom data.table data.table fcase
#' @importFrom mgcv gam gamm
#' @importFrom nlme corAR1
#' @param dag dagitty object
#' @param data A data frame or data table, with variables corresponding to nodes in the supplied dag.
#' @param random_effects A vector or list of variables to be treated as random effects. Only used when "gamm" or "gam" is specified, the  are not required to be included in the dag.
#' @param method Defaults to method = "REML", only used for gamm models.
#' @param family Used for glm and gamm models, defaults to gaussion().
#' @param corr_structure Limited functionality. Currently only used when model = "lme" or "gamm" to specify nlme::corAR1(0, form=~1|year) if a random_effects input contains the string 'year'.
#' @return Model object of either
#' @export
causalify_gamm <- function(dag, data, random_effects = NULL, method = "REML", family = gaussion(), corr_structure = "corAR1"){
  .datatable.aware <- TRUE

  # gamm / lme
  setClass("definition",
           slots = list(outcome = "character",
                        treatment = "character",
                        vars = "character",
                        data = "data.table",
                        random = "list",
                        correlation = "function",
                        family = "function",
                        method = "character",
                        model = "function",
                        gamm = "list"),
           where = causali.envr
  )

  setMethod("initialize", "definition",
            function(.Object, outcome, treatment, vars, data, random, correlation, family, method, model) {
              fo <- as.formula(paste(outcome, "~", treatment, paste(" + ", vars, collapse = " "), collapse = " "))
              model_gamm <- model(fo,
                                  data  = data,
                                  random = random,
                                  family = family,
                                  method = method)
              .Object@gamm <- list(model_gamm)
              .Object
            },
            where = causali.envr
  )


  {
    random_effects <- list("year", "loc", "rep")

    random_list <- lapply(seq_along( random_effects ), function(x){

      random_list <- list()

      random_list[x] <- c(~1)

    })

    random_list <- unlist(random_list)

    names(random_list) <- random_effects


    if( "year" %in% random_effects ){

      cor_var <- random_effects[ random_effects %in% "year" ]

      year <- function (x)
      {
        UseMethod(as.character(cor_var))
      }

      corr_structure <- nlme::corAR1(0, form=~1|year)

    }else{

      corr_structure <- NULL

    }

  }


  causalify_gamm <- with(causali.envr, function(dag, data, random, family, method, corr_structure){

    feature_map_list <- getFeatureMap(dag)

    gamm <- new("definition",
                outcome = as.character(feature_map_list$outcome$ancestor),
                treatment = as.character(feature_map_list$treatment$ancestor),
                vars = as.character(feature_map_list$confounder$ancestor),
                data = tab_work,
                random = random_list,
                family = family,
                method = method,
                correlation = cor_AR1,
                model = mgcv::gamm)

    return(gamm)

  } )

  gamm_out <- causalify_gamm(dag, tab_work, random, family, method, corr_structure)

  return(gamm_out)

}
