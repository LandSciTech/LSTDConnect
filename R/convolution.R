#' @include AAAClassDefinitions.R
#' @include helperFns.R
NULL

#' Distance weighted sum convolution or intactness. 
#'
#' @details
#' Exponential kernel as in Hughes et al 2015.
#'
#' Uniform kernel gives buffered sum.
#' 
#' Intactness as in Beyer et al 2019.
#'
#' @param quality RasterLayer.
#' @param d Numeric. Distance parameter. For euclidean kernel this is average dispersal distances. For uniform kernel this is width of buffer or average distance if useAveDist=T. In the units of x.
#' @param kernel String or matrix. "Exponential" or "Uniform". If a matrix is supplied it will be used as the convolution kernel, and d/negligible/useAveDist arguments will be ignored.
#' @param patches RasterLayer. Optional. Used to mask return layer if provided.
#' @param negligible Numeric. Truncation value for exponential kernel.
#' @param convolutionMethod Character. focal, velox or spatialFil. focal default. Use spatialFil for wrapped boundaries. velox and focal return identical results, but velox is much faster.
#' @param useAveDist Logical. If TRUE d is average dispersal distance for uniform kernels. Otherwise (default) d is buffer width.
#' @param z Numeric. z parameter for intactness metric. If z=NA the output will be simple distance weighted sum.
#' @examples
#' # TODO: examples
#' @export
setGeneric('applyKernel',function(quality,d,kernel,patches=NULL,negligible=10^-10,convolutionMethod="focal",useAveDist=F,z=NA) standardGeneric('applyKernel'))

#' @return A RasterLayer - value of each pixel is kernel weighted neighbourhood average of quality, or intactness if z!=NA.
#' @rdname applyKernel
#' @export
setMethod('applyKernel', signature(quality="RasterLayer"), function(quality,d,kernel,patches,negligible,convolutionMethod,useAveDist,z) {
  #quality=exQuality;kernelShape="Exponential";negligible=10^-10

  if(class(kernel)=="matrix"){
    if(sum(kernel,na.rm=T)!=1){
      warning("Sum of the convolution kernel is ",sum(matrix,na.rm=T),". To be interpreted as a dispersal kernel the sum should be 1.")
    }
    k=kernel
  }else{
    kernelTypes = c("Exponential","Uniform")
    if(!is.element(kernel,kernelTypes)){
      stop("kernel shape ",kernel," not recognized. Options are:",paste(kernelTypes,collapse=","))
    }
    cellDim = res(quality)[1]
    dbar=d/cellDim
    
    if(kernelShape=="Exponential"){
      k = exponentialKernel(dbar,negligible=negligible)
    }
    if(kernelShape=="Uniform"){
      k = uniformKernel(dbar,useAveDist=useAveDist)
    }
  }
  
  if(!is.na(z)){
    kShape=k
    kShape[kShape!=0] = 1
    dvec <- as.vector(k)
    centrePt = ceiling(length(dvec)/2) #NOTE: this only works for symmetric kernels
    trMap <- focal(quality, kShape, fun= function(x){qprime(x,centrePt,dvec,z)}, pad=FALSE)
    #TO DO: testing qprime
  }else{
    convMethods = c("velox","spatialFil","focal")
    if(!is.element(convolutionMethod,convMethods)){
      stop("convolutionMethod not recognized. Options are: ",paste(convMethods,sep=","))
    }
    if(convolutionMethod=="spatialFil"){
      #TO DO: stop if spatialFil package not loaded.
      temp = spatialfil::convKernel(sigma = 1, k = c("gaussian"))
      temp$matrix=k
      temp$kernel=as.character(kernel)
      temp$sigma = dbar
      mTemp = raster::as.matrix(quality)
      trTemp = spatialfil::applyFilter(mTemp,temp)
      trTemp=raster::raster(trTemp)
      trMap = quality
      values(trMap)=values(trTemp)
    }
    if(convolutionMethod=="velox"){
      #TO DO: stop if velox not loaded.
      vx = velox::velox(quality)
      vx$sumFocal(weights=k, bands=c(1))
      trMap = vx$as.RasterLayer()
    }
    if(convolutionMethod=="focal"){
      #kShape=k
      #kShape[kShape!=0] = 1
      #dvec <- as.vector(k)
      #trMap <- focal(quality, kShape, fun= function(x){focalConvolve(x,dvec)}, pad=FALSE)
      trMap <- focal(quality, k, fun= sum, pad=FALSE)
      #TO DO: figure out what to do about PAD argument. 
    }
  }
  
  if(!is.null(patches)){
    trMap[(patches==0)|is.na(patches)]=NA
  }
  return(trMap)
})

#' @return A convolved RasterBrick - value of each pixel is kernel weighted neighbourhood average of quality.
#' @rdname applyKernel
#' @export
setMethod('applyKernel', signature(quality="RasterBrick"), function(quality,d,kernel,patches,negligible,convolutionMethod,useAveDist,z) {
   #patches=exPatch
   outBrick=quality
   for(nn in names(quality)){
     outBrick[[nn]]=applyKernel(quality[[nn]],d,kernel,patches[[nn]],negligible,convolutionMethod,useAveDist,z)
   }
   return(outBrick)
})


# this is the function called by focal for intactness metric
#' @export
qprime <- function(x,centrePt,dvec,z){
  #x=exQuality[1:nrow(foc.w),1:ncol(foc.w)]  
  if (is.na(x[centrePt])) return(NA)
  recs <- which(!is.na(x))
  #return(sum((x[centrePt] * x[recs])^z * dvec[recs]) / sum(dvec[recs]))
  return(sum((x[centrePt] * x[recs])^z * dvec[recs])) #assume dvec has already been standardized for efficiency
}

# this is the function called by focal for DWS convolution
#' @export
focalConvolve <- function(x,dvec){
  #x=exQuality[1:nrow(foc.w),1:ncol(foc.w)]  
  recs <- which(!is.na(x))
  return(sum( x[recs] * dvec[recs]))
}

#' @export
getDistKernelFromMax<-function(kdim){
  #kdim=2
  kdim=ceiling(kdim)
  wDim = kdim*2+1
  locSeq = seq(-kdim,kdim)
  y = matrix(locSeq, wDim, wDim)
  xx=t(y)
  d = (xx^2+y^2)^0.5
  return(d)
}

#' @export
uniformKernel<-function(dmax,cellDim=1,useAveDist=F){
  if(useAveDist){
    #https://math.stackexchange.com/questions/3019165/average-distance-from-center-of-circle
    dmax=3*dmax/2
  }
  hdim = ceiling(dmax/cellDim)
  weights =getDistKernelFromMax(hdim)
  weights=weights <= dmax/cellDim
  weights=weights/sum(weights)
  return(weights)
}

#' @export
exponentialKernel<-function(dbar,cellDim=1,negligible=10^-10,returnScale=F,dmax=NULL){
  #Exponential kernel from Hughes et al 2015 American Naturalist
  dbarCell = dbar/cellDim
  if(is.null(dmax)){
    dmax = -0.5*dbarCell*log(pi*dbarCell^2*negligible/2)
    if(dmax<0){
      stop("Set negligible so that pi*dbar^2*negligible/2 <=1")
    }
  }
  kdim = floor(dmax/2)*2
  d=getDistKernelFromMax(kdim)
  m = (2/(pi*dbarCell^2))*exp(-2*d/dbarCell)
  m[m<negligible]=0
  if(returnScale){
    return(sum(m))
  }
  k = m/sum(m)
  return(k)
}
