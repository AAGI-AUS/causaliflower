#' Package imports
#'
#' @importFrom stats complete.cases na.omit reshape
#' @importFrom utils globalVariables
#' @noRd
NULL

globalVariables(c(
  "..cols_required", "ESCDAGs", "instrument", "label", "outcome",
  "role", "treatment", "x", "xend", "y", "yend"
))

#' @keywords internal
"_PACKAGE"
