#################################################################
# creates a new environment on package load for needed variables
causali.envr <- NULL
.onLoad <- function(...) {

  causali.envr <<- new.env()

}
#################################################################
