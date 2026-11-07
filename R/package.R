# ==============================================================================
# package.R


# Package imports and names registered for R's code checks
# ==============================================================================


#' Package imports
#'
#' @importFrom data.table :=
#' @importFrom stats ave complete.cases na.omit reshape
#' @importFrom utils flush.console globalVariables head tail
#' @noRd
NULL

globalVariables(c( # Register names used in table expressions and bundled data for R's code checks
  "..cols_required", ".edge_group", "ESCDAGs", "display_label", "from", "index",
  "instrument", "is_lagged", "kind", "lag", "management_defaults",
  "label", "name", "node_aliases", "outcome", "role", "to", "treatment", "x", "xend", "y",
  "yend"
))

#' @useDynLib causaliflower, .registration = TRUE, .fixes = "C_"
#' @keywords internal
"_PACKAGE"
