#Install packages code - this only needs to happen once
#install.packages(c("NLMR","landscapetools","devtools"))
#devtools::install_github("LandSciTech/LSTDConnect")
#devtools::document();devtools::load_all()

library(NLMR)
library(landscapetools)
library(raster)
library(LSTDConnect)
library(samc)

#Example landcover map
exQuality=nlm_mpd(ncol=200,nrow=200,resolution=1,roughness=0.9,rand_dev=1)
plot(exQuality)

oneMap = exQuality;oneMap[exQuality!=1]=1

#Example - variable mortality & occurence
backgroundMortality = 10^-5
t = 91
resistance <- oneMap
absorption <- oneMap*(backgroundMortality+(1/15)*(1-backgroundMortality)*(1-exQuality))
plot(absorption)
occurrence <- exQuality

#return RasterLayer
trMap = applySAMC(occurrence,absorption,resistance,t) 
plot(trMap)

#no need for raster package if all inputs are matrices. output will also be a matrix
trMap = applySAMC(as.matrix(occurrence),as.matrix(absorption),as.matrix(resistance),t) 

#if doing multiple calculations for the same resistance/absorption surfaces, will be faster to construct samcObj only once
samcObj <- samc::samc(as.matrix(resistance), as.matrix(absorption), tr_fun = function(x) 1/mean(x))
times = c(1,20,40,60,80,100)
mStack = stack(occurrence)
for(t in times){
  mStack[[paste("time",t)]] = applySAMC(occurrence,t=t,samcObj=samcObj)
}
plot(mStack)

#TO DO:
# - create & use lookup table to find t that most closely corresponds to selected d.

