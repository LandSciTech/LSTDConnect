#' Distance weighted sum convolution or intactness.
#'
#' @param quality RasterLayer.
#' @param d Numeric. Distance parameter. For euclidean kernel this is average 
#'   dispersal distances. For uniform kernel this is width of buffer or average 
#'   distance if useAveDist=T. In the units of x.
#' @param kernel String or matrix. "Exponential" or "Uniform". If a matrix is 
#'   supplied it will be used as the convolution kernel, and d/negligible/useAveDist 
#'   arguments will be ignored.
#' @param patches RasterLayer. Optional. Used to mask return layer if provided.
#' @param negligible Numeric. Truncation value for exponential kernel.
#' @param convolutionMethod Character. focal, velox or spatialFil. focal default. 
#'   Use spatialFil for wrapped boundaries. velox and focal return identical 
#'   results, but velox is much faster.
#' @param useAveDist Logical. If TRUE d is average dispersal distance for 
#'   uniform kernels. Otherwise (default) d is buffer width.
#' @param z Numeric. z parameter for intactness metric. If z=NA the output 
#'   will be simple distance weighted sum.
#' 
#' @details
#' Exponential kernel as in Hughes et al 2015.
#' Uniform kernel gives buffered sum.
#' Intactness as in Beyer et al 2019.
#' 
#' @examples
#' # TODO: examples
#' 
#' @return 
#' A RasterLayer - value of each pixel is kernel weighted neighbourhood average 
#' of quality, or intactness if z!=NA.
#' If quality is a rasterbrick, A convolved RasterBrick - value of each pixel is 
#' kernel weighted neighbourhood average of quality.
#' 
#' @export
setGeneric("applyKernel", 
           function(quality, d, kernel, patches = NULL, negligible = 10^-10, 
                    convolutionMethod = "focal", useAveDist = F, z = NA) {
             standardGeneric("applyKernel")}
)

#' @rdname applyKernel
#' @export
setMethod("applyKernel", signature(quality = "RasterLayer"), 
          function(quality, d, kernel, patches, negligible, 
                   convolutionMethod, useAveDist, z) {
  # quality=exQuality;kernelShape="Exponential";negligible=10^-10
  
  if (class(kernel) == "matrix") {
    if (sum(kernel, na.rm = T) != 1) {
      warning("Sum of the convolution kernel is ", sum(matrix, na.rm = T), 
              ". To be interpreted as a dispersal kernel the sum should be 1.")
    }
    k <- kernel
  } else {
    kernelTypes <- c("Exponential", "Uniform")
    if (!is.element(kernel, kernelTypes)) {
      stop("kernel shape ", kernel, " not recognized. Options are:", 
           paste(kernelTypes, collapse = ","))
    }
    cellDim <- raster::res(quality)[1]
    dbar <- d / cellDim
    
    if (kernel == "Exponential") {
      k <- exponentialKernel(dbar, negligible = negligible)
    }
    if (kernel == "Uniform") {
      k <- uniformKernel(dbar, useAveDist = useAveDist)
    }
  }
  
  if (!is.na(z)) {
    kShape <- k
    kShape[kShape != 0] <- 1
    dvec <- as.vector(k)
    centrePt <- ceiling(length(dvec) / 2) # NOTE: this only works for symmetric kernels
    trMap <- raster::focal(quality, kShape, fun = function(x) {
      qprime(x, centrePt, dvec, z)
    }, pad = FALSE)
    # TO DO: testing qprime
  } else {
    convMethods <- c("velox", "spatialFil", "focal")
    if (!is.element(convolutionMethod, convMethods)) {
      stop("convolutionMethod not recognized. Options are: ", 
           paste(convMethods, sep = ","))
    }
    if (convolutionMethod == "spatialFil") {
      # TO DO: stop if spatialFil package not loaded.
      temp <- spatialfil::convKernel(sigma = 1, k = c("gaussian"))
      temp$matrix <- k
      temp$kernel <- as.character(kernel)
      temp$sigma <- dbar
      mTemp <- raster::as.matrix(quality)
      trTemp <- spatialfil::applyFilter(mTemp, temp)
      trTemp <- raster::raster(trTemp)
      trMap <- quality
      raster::values(trMap) <- raster::values(trTemp)
    }
    if (convolutionMethod == "velox") {
      # TO DO: stop if velox not loaded.
      vx <- velox::velox(quality)
      vx$sumFocal(weights = k, bands = c(1))
      trMap <- vx$as.RasterLayer()
    }
    if (convolutionMethod == "focal") {
      # kShape=k
      # kShape[kShape!=0] = 1
      # dvec <- as.vector(k)
      # trMap <- raster::focal(quality, kShape, fun= function(x){focalConvolve(x,dvec)}, pad=FALSE)
      quality[is.na(quality)] <- 0
      trMap <- raster::focal(quality, k, fun = sum, pad = FALSE)
      # TO DO: figure out what to do about PAD argument.
    }
  }
  
  if (!is.null(patches)) {
    trMap[(patches == 0) | is.na(patches)] <- NA
  }
  return(trMap)
})

#' @rdname applyKernel
#' @export
setMethod("applyKernel", signature(quality = "RasterBrick"), 
          function(quality, d, kernel, patches, negligible, 
                   convolutionMethod, useAveDist, z) {
  # patches=exPatch
  outBrick <- quality
  for (nn in names(quality)) {
    outBrick[[nn]] <- applyKernel(quality[[nn]], d, kernel, patches[[nn]], 
                                  negligible, convolutionMethod, useAveDist, z)
  }
  return(outBrick)
})
