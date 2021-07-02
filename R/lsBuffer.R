#' @include AAAClassDefinitions.R
#' @include helperFns.R
NULL
# Note this is copied from lstools to avoid dependence on a private repository - 
# reintegrate once that library is public.

#' Apply buffer.
#'
#' @details
#' For a Raster* object use velox package.
#'
#' @param x RasterLayer.
#' @param width Numeric. Width of buffer to apply, in the units of x.
#' @param subset Logical expression. Optional. Passed as an argument to the subset 
#'   function in order to select a subset of raster values or polygons. Default 
#'   for rasters is expression(x!=0) - 0s and NAs are ignored. Default for 
#'   polygons is NULL (all polygons selected).
#'   
#' @return 
#' A binary RasterLayer object showing buffered area.
#'   
#' @examples
#' # TODO: examples
#' 
#' @export
setGeneric("lsBuffer", function(x, width, subset = NULL) standardGeneric("lsBuffer"))

#' @rdname lsBuffer
#' @export
setMethod("lsBuffer", signature(x = "RasterLayer"), 
          function(x, width, subset) {
  # x=patches;subset=NULL;width=buffDist#x=x[[kk]]

  if (is.null(subset)) {
    subset <- expression(x != 0)
  }

  xtemp <- x # distinguish between NA values and 0 values. Note velox does not handle NA values well.
  x[is.na(x)] <- 0

  x[!eval(subset)] <- 0 # TO DO test this.
  if (length(raster::unique(x)) > 0) {
    if (width > 0) {
      weights <- uniformKernel(width, raster::res(x)[1])

      vx <- velox::velox(x)
      vx$sumFocal(weights = weights, bands = c(1))

      x <- vx$as.RasterLayer()
      x[x != 0] <- 1
    } else {
      x[x > 1] <- 1
    }
  }
  x[is.na(x)] <- 0
  x[is.na(xtemp)] <- NA
  return(x)
})
