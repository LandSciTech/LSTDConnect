#' TODO
#' 
#' TODO
#' 
#' @param kdim TODO
#' 
#' @return 
#' TODO
#' 
#' @export
getDistKernelFromMax <- function(kdim) {
  # kdim=2
  kdim <- ceiling(kdim)
  wDim <- kdim * 2 + 1
  locSeq <- seq(-kdim, kdim)
  y <- matrix(locSeq, wDim, wDim)
  xx <- t(y)
  d <- (xx^2 + y^2)^0.5
  return(d)
}

# Helper functions for kernel work ----------------------------------------

# this is the function called by focal for intactness metric
qprime <- function(x, centrePt, dvec, z) {
  # x=exQuality[1:nrow(foc.w),1:ncol(foc.w)]
  if (is.na(x[centrePt])) {
    return(NA)
  }
  recs <- which(!is.na(x))
  # return(sum((x[centrePt] * x[recs])^z * dvec[recs]) / sum(dvec[recs]))
  return(sum((x[centrePt] * x[recs])^z * dvec[recs])) # assume dvec has already been standardized for efficiency
}

# this is the function called by focal for DWS convolution
focalConvolve <- function(x, dvec) {
  # x=exQuality[1:nrow(foc.w),1:ncol(foc.w)]
  recs <- which(!is.na(x))
  return(sum(x[recs] * dvec[recs]))
}