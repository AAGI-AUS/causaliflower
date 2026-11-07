################################################################################
# zzz.R


# Bind prepared recovery metadata when the package loads

.onLoad <- function(libname, pkgname){ # Bind the reviewed recovery interpretation when prepared SCM data are loaded
  return(.cf_bind_recovery_policy(environment(.onLoad)))
}
