#Install packages code - this only needs to happen once
#install.packages(c("NLMR","landscapetools","devtools"))
#devtools::install_github("LandSciTech/LSTDConnect")
#devtools::document();devtools::load_all()

library(NLMR)
library(landscapetools)
library(raster)
library(LSTDConnect)

#Example landcover map
exQuality=nlm_mpd(ncol=200,nrow=200,resolution=1,roughness=0.9,rand_dev=1)
plot(exQuality)

#Redistribution kernel examples
exponentialKernel(0.3,negligible=10^-6) #exponential kernel with mean dispersal distance of 1. Note the tail of the kernel (density < 10^-6) is truncated.
uniformKernel(1,useAveDist=T) #uniform kernel with mean dispersal distance of 1
uniformKernel(1,useAveDist=F) #uniform kernel with max dispersal distance of 1
sum(exponentialKernel(1)) #sum of discretized redistribution kernel is always 1

#I haven't yet written documentation for exponentialKernel and uniformKernel functions.
#The functions are small, though. You can see how they work by looking at the code.
exponentialKernel
uniformKernel

#See documentation for applyKernel function. Such as it is...
?applyKernel

#Expample, exponential convolution kernel, euclidean distance
dbar= 5
#absorbing boundaries - note edge effects
trMap = applyKernel(exQuality,dbar,kernelShape="Exponential",useSpatialfil = F)
sum(values(exQuality))
sum(values(trMap)) # Note loss from absorbing boundaries
plot(trMap)
#The value at each pixel of trMap is the distance weighted sum of the input map in the vicinity of that pixel.

#wrapped boundaries - wrapped boundary conditions are useful avoiding edge effects on small simulated landscapes. Don't use with real data.
dwMean5 = applyKernel(exQuality,dbar,kernelShape="Exponential",useSpatialfil = T)
sum(values(exQuality))
sum(values(dwMean5)) # The kernel is a redistribution function so the sum of the input and output maps should match more or less - there is a little numerical error.
plot(stack(exQuality,dwMean5))

#uniform convolution kernel, euclidean distance. This is a threshold based analysis or buffered sum.
tbMax5 = applyKernel(exQuality,5,kernelShape="Uniform",useSpatialfil = T) #wrapped boundaries, buffer width = 5
tbMean5 = applyKernel(exQuality,5,kernelShape="Uniform",useSpatialfil = T,useAveDist=T) #wrapped boundaries, mean dispersal distance = 5

viewStack = stack(dwMean5,tbMean5,tbMax5);names(viewStack)=c("exp mean 5","unif mean 5","unif max 5")
plot(viewStack)
#Note exponential and uniform kernels with mean distance of 5 are more similar to one another than either is to the uniform kernel with a 5 unit radius.

