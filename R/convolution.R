#' @include AAAClassDefinitions.R
#' @include helperFns.R
NULL

#' Convolution of a raster with a kernel.
#'
#' @details
#' Exponential kernel as in Hughes2015.
#'
#' Uniform kernel gives buffered sum.
#'
#' @param quality RasterLayer.
#' @param d Numeric. Distance parameter. For euclidean kernel this is average dispersal distances. For uniform kernel this is width of buffer or average distance if useAveDist=T. In the units of x.
#' @param kernelShape String. Exponential or Uniform
#' @param patches RasterLayer. Optional. Used to mask return layer if provided.
#' @param negligible Numeric. Truncation value for exponential kernel.
#' @param useSpatialfil Logical. Use spatialfil instead of velox. Slower (?), wrapped boundaries.
#' @param useAveDist Logical. If TRUE d is average dispersal distance for uniform kernels. Otherwise (default) d is buffer width.
#' @examples
#' # TODO: examples
#' @export
setGeneric('applyKernel',function(quality,d,kernelShape,patches=NULL,negligible=10^-10,useSpatialfil=F,useAveDist=F) standardGeneric('applyKernel'))

#' @return A convolved RasterLayer - value of each pixel is kernel weighted neighbourhood average of quality.
#' @rdname applyKernel
#' @export
setMethod('applyKernel', signature(quality="RasterLayer"), function(quality,d,kernelShape,patches,negligible,useSpatialfil,useAveDist) {
  #quality=exQuality;kernelShape="Exponential";negligible=10^-10

  kernelTypes = c("Exponential","Uniform")
  if(!is.element(kernelShape,kernelTypes)){
    stop("kernelShape ",kernelShape," not recognized. Options are:",paste(kernelTypes,collapse=","))
  }
  cellDim = res(quality)[1]
  dbar=d/cellDim

  if(kernelShape=="Exponential"){
    k = exponentialKernel(dbar,negligible=negligible)
  }
  if(kernelShape=="Uniform"){
    k = uniformKernel(dbar,useAveDist=useAveDist)
  }

  if(useSpatialfil){
    temp = spatialfil::convKernel(sigma = 1, k = c("gaussian"))
    temp$matrix=k
    temp$kernel=kernelShape
    temp$sigma = dbar

    mTemp = raster::as.matrix(quality)

    trTemp = spatialfil::applyFilter(mTemp,temp)
    trTemp=raster::raster(trTemp)
    trMap = quality
    values(trMap)=values(trTemp)
  }else{
    vx = velox::velox(quality)
    vx$sumFocal(weights=k, bands=c(1))

    trMap = vx$as.RasterLayer()
  }
  if(!is.null(patches)){
    trMap[(patches==0)|is.na(patches)]=NA
  }
  return(trMap)
})

#' @return A convolved RasterBrick - value of each pixel is kernel weighted neighbourhood average of quality.
#' @rdname applyKernel
#' @export
setMethod('applyKernel', signature(quality="RasterBrick"), function(quality,d,kernelShape,patches,negligible,useSpatialfil,useAveDist) {
   #patches=exPatch
   outBrick=quality
   for(nn in names(quality)){
     outBrick[[nn]]=applyKernel(quality[[nn]],d,kernelShape,patches[[nn]],negligible,useSpatialfil,useAveDist)
   }
   return(outBrick)
})


#' @export
getDistKernelFromMax<-function(kdim){
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
