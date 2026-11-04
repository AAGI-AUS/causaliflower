# ==============================================================================
# package.R


# Package imports and names registered for R's code checks
# ==============================================================================


#' Package imports
#'
#' @importFrom stats complete.cases na.omit reshape
#' @importFrom utils globalVariables tail
#' @noRd
NULL

globalVariables(c( # Register names used in table expressions for R's code checks
  "..cols_required", ".edge_group", "ESCDAGs", "display_label", "index",
  "instrument",
  "label", "name", "outcome", "role", "treatment", "x", "xend", "y",
  "yend"
))

#' @keywords internal
"_PACKAGE"
