# Date: 2017-Aug-03
#'
#' Is an Object of Class "cif"?
#'
#' Checks whether an object is of Class "cif".
#'
#' @param x An R object
#'

is.cif <- function(x)
  inherits(x, "cif")
