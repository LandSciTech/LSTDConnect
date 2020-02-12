
#install.packages(c("NLMR","landscapetools","devtools"))
#devtools::install_github("LandSciTech/LSTDConnect")
#devtools::document();devtools::load_all()

library(NLMR)
library(landscapetools)
library(raster)

exQuality=nlm_mpd(ncol=200,nrow=200,resolution=1,roughness=0.9,rand_dev=1)

plot(exQuality)

exPatch = util_binarize(exQuality,0.17)-1
plot(exPatch)

#exponential convolution kernel, euclidean distance
dbar= 5
exponentialKernel(0.2)
#wrapped boundaries
trMap = applyKernel(exQuality,dbar,kernelShape="Exponential",useSpatialfil = T)
sum(values(exQuality))
sum(values(trMap)) # kernel is a redistribution function. Total sum should match (more or less).
plot(trMap)

#absorbing boundaries
trMap = applyKernel(exQuality,dbar,kernelShape="Exponential",useSpatialfil = F)
sum(values(exQuality))
sum(values(trMap)) # Loss from absorbing boundaries
plot(trMap)

#Use patches as mask
trMap = applyKernel(exQuality,dbar,patches=exPatch,kernelShape="Exponential",useSpatialfil = T)
plot(trMap)

#uniform convolution kernel, euclidean distance. This is a buffered sum.
dmax = 10
uniformKernel(2) #buffer width of 2
uniformKernel(2,useAveDist=T) #mean dispersal distance of 2
trMap = applyKernel(exQuality,dmax,kernelShape="Uniform",useSpatialfil = T) #wrapped boundaries
sum(values(exQuality))
sum(values(trMap)) # kernel is a redistribution function. Total sum should match (more or less).
plot(trMap)

#buffered sum of a binary map gives the density of that feature in the vicinity of each pixel.
#e.g. patch density
trMap = applyKernel(exPatch,dmax,kernelShape="Uniform",useSpatialfil = T) #wrapped boundaries
plot(trMap)
patchMask=exPatch;patchMask[patchMask==0]=NA
plot(patchMask,add=T)

#TO DO: implement uniform kernel (i.e. buffer sums). Send to Scott.
#TO DO: apply exponential and uniform variants to simulated landscape sets. Correlate with PARC given same benefit map.
