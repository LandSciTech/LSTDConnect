library(raster)
library(NLMR)
#Example landcover map
exQuality=nlm_mpd(ncol=500,nrow=500,resolution=1,roughness=0.9,rand_dev=1)
exQuality = exQuality*140;exQuality[exQuality>100]=100
#exQuality[!is.na(exQuality)]=80
#exQuality[100,100]=100
plot(exQuality)

foc.w = exponentialKernel(10,negligible=10^-6) #exponential kernel with mean dispersal distance of 1. Note the tail of the kernel (density < 10^-6) is truncated.

cc = applyKernel(exQuality,kernel=foc.w,convolutionMethod = "focal") #simple distance weighted sum using focal
dd = applyKernel(exQuality,kernel=foc.w,convolutionMethod = "velox") #velox implementation
ee = applyKernel(exQuality,kernel=foc.w,z=0.5) #qprime intactness
plot(ee)
plot(stack(exQuality,ee))

#TO DO: figure out transformation from our quality metric to Beyer's exp transformed metric.
#We use an exponential transformation to scale HFI (range [0, 50]) to quality (range [0, 1]) such that quality is 1 when HFI is 0, and quality is 20% when HFI is 4.
gamma <- -log(c(0.2))/4

dd = data.frame(quality=seq(0,100))
dd$hfi = (100-dd$quality)*0.5
dd$q = exp(-gamma*dd$hfi)
plot(q~quality,data=dd)

r.hfi = (100-exQuality)*0.5
hist(exQuality@data@values)

hist(r.hfi@data@values)
#r.hfi <- raster("hfp2013_meris.tif")
r.q <- exp(-gamma * r.hfi)
hist(r.q@data@values)

plot(stack(exQuality,r.hfi,r.q))
#TO DO: figure out what to do about pad argument.
#TO DO: proper testing of qprime intactness.


