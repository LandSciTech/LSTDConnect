#' Convolution kernels
#' 
#' Creates a convolution kernel
#' 
#' @param r The radius of the kernel, in physical units.
#' @param mu Alternative to r, so that r = 3 * mu / 2 (equation for the 
#'    average distance from the center of a circle).
#' @param cellDim Default to 1, cell dimensions used as a scaling factor 
#'    for the kernel size
#' @param negligible Negligible value when rounding edges of exponential kernels.
#' 
#' @return 
#' A kernel (`matrix`).
#' 
#' @rdname kernels
#' @export
uniformKernel <- function(r=NULL, mu=NULL, cellDim = 1) {
  if (!is.null(mu)) {
    # https://math.stackexchange.com/questions/3019165/average-distance-from-center-of-circle
    r <- 3 * mu / 2
  }
  hdim <- ceiling(r / cellDim)
  weights <- getDistKernelFromMax(hdim)
  weights <- weights <= r / cellDim
  weights <- weights / sum(weights)
  return(weights)
}

#' @rdname kernels
#' @export
exponentialKernel <- function(mu, cellDim = 1, negligible = 10^-10) {
  # Exponential kernel from Hughes et al 2015 American Naturalist
  muCell <- mu / cellDim
  
  dmax <- -0.5 * muCell * log(pi * muCell^2 * negligible / 2)
  
  if (dmax < 0) {
    stop("Set negligible so that pi*mu^2*negligible/2 <=1")
  }
  
  kdim <- floor(dmax / 2) * 2
  d <- getDistKernelFromMax(kdim)
  m <- (2 / (pi * muCell^2)) * exp(-2 * d / muCell)
  m[m < negligible] <- 0
  
  k <- m / sum(m)
  return(k)
}