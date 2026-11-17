# ==============================================================================
# zzz.R


# Bind prepared recovery metadata when the package loads
# ==============================================================================


#' Bind effective recovery metadata when the namespace loads
#'
#' @param libname,pkgname Library path and package name supplied by R's namespace
#'   loader. This hook does not use these arguments.
#' @return `NULL`, invisibly, after binding the prepared table in the namespace
#'   and lazy-data environments.
#' @noRd
.onLoad <- function(libname, pkgname){
  return(.cf_bind_recovery_policy(environment(.onLoad)))
}
