#' Get kernel distance from maximum
#' 
#' Calculates distance d from kernel dimensions.
#' 
#' @param kdim kernel dimensions.
#' 
#' @return 
#' A matrix.
#'

getDistKernelFromMax <- function(kdim) {
  kdim <- ceiling(kdim)
  wDim <- kdim * 2 + 1
  locSeq <- seq(-kdim, kdim)
  y <- matrix(locSeq, wDim, wDim)
  xx <- t(y)
  d <- (xx^2 + y^2)^0.5
  return(d)
}
