#' Convolution kernels
#' 
#' Creates a convolution kernel
#' 
#' @param dmax TODO
#' @param cellDim TODO
#' @param useAveDist TODO
#' @param dbar TODO
#' @param negligible TODO
#' @param returnScale TODO
#' @param dmax TODO
#' 
#' @return 
#' A kernel (`matrix`).
#' 
#' @rdname kernels
#' @export
uniformKernel <- function(dmax, cellDim = 1, useAveDist = F) {
  if (useAveDist) {
    # https://math.stackexchange.com/questions/3019165/average-distance-from-center-of-circle
    dmax <- 3 * dmax / 2
  }
  hdim <- ceiling(dmax / cellDim)
  weights <- getDistKernelFromMax(hdim)
  weights <- weights <= dmax / cellDim
  weights <- weights / sum(weights)
  return(weights)
}

#' @rdname kernels
#' @export
exponentialKernel <- function(dbar, cellDim = 1, negligible = 10^-10, 
                              returnScale = F, dmax = NULL) {
  # Exponential kernel from Hughes et al 2015 American Naturalist
  dbarCell <- dbar / cellDim
  if (is.null(dmax)) {
    dmax <- -0.5 * dbarCell * log(pi * dbarCell^2 * negligible / 2)
    if (dmax < 0) {
      stop("Set negligible so that pi*dbar^2*negligible/2 <=1")
    }
  }
  kdim <- floor(dmax / 2) * 2
  d <- getDistKernelFromMax(kdim)
  m <- (2 / (pi * dbarCell^2)) * exp(-2 * d / dbarCell)
  m[m < negligible] <- 0
  if (returnScale) {
    return(sum(m))
  }
  k <- m / sum(m)
  return(k)
}