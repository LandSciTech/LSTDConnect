#' Least cost path distance matrix from rasters.
#'
#' @details
#' TO DO: add description
#'
#' @param patches RasterLayer. Binary patch map.
#' @param cost RasterLayer. Cost. Must be >0. Set cost=1 for euclidean distance.
#' @param maxDist Numeric. Maximum distance between i and j in map units.
#' @param bufferedPatches RasterLayer. Optional. If provided maxDist will be ignored. bufferedPatches defines the calculation region - use to restrict calculations to a portion of the map area.
#' @param neighbourhood Character. 'rook','queen', or 'octagon'. 'octagon' option is a modified version of the queen's 8 cell neighbourhood in which diagonals weights are 2^0.5x higher than horizontal/vertical weights.
#' @examples
#' # TO DO: examples
#' @export
setGeneric('leastCostPathDistances',function(patches,cost,maxDist,bufferedPatches,neighbourhood='octagon') standardGeneric('leastCostPathDistances'))

#' @return List of effective distance d_ij matrices for each subgraph.
#' @rdname leastCostPathDistances
setMethod('leastCostPathDistances', signature(patches="RasterLayer"), function(patches,cost,maxDist,bufferedPatches,neighbourhood) {
  #patches=cPatches; cost=cCost;bufferedPatches=!is.na(cPatches);maxDist=maxDist
  checkAllign = raster::compareRaster(patches,cost)
  if(!checkAllign){
    stop("All input rasters must have the same same extent, number of rows and columns,projection, resolution, and origin.")
  }
  if(is.null(bufferedPatches)){
    #Buffer patches
    costMin = raster::cellStats(cost,min)
    if(costMin==0){stop("Cost must be >0.")}
    buffDist = maxDist/costMin
    bufferedPatches = lsBuffer(patches,maxDist)
  }

  bufferedPatches[bufferedPatches==0]=NA
  #plot(bufferedPatches)
  patches[is.na(bufferedPatches)]=NA
  cost[is.na(bufferedPatches)]=NA
  #scale cost to mapunits
  cRes = raster::res(patches)
  if(isFALSE(all.equal(cRes[1],cRes[2]))){stop("Method assumes square pixels. x and y map resolution must be equal.")}
  cost = cost*cRes[1]
  #sort( sapply(ls(),function(x){object.size(get(x))}))

  #cost[!is.na(cost)]=cRes[1]
  #library(tidyverse);library(data.table)
  sim = roadCLUS.getGraph(sim=list(costSurface=cost),neighbourhood = neighbourhood)

  #sim = roadCLUS.getGraph(sim=list(costSurface=cost))
  #id clusters and loop over each cluster
  buffClumps = raster::clump(bufferedPatches)
  #plot(buffClumps)
  clumps = raster::unique(buffClumps)
  igList=list()
  for (ic in clumps){
    #ic=clumps[1]
    cPatches=patches;cCost=cost
    cPatches[buffClumps!=ic]=NA
    cCost[buffClumps!=ic]=NA

    #get IDs of from vertices and to vertices.
    fromV = getVs(cPatches>0)
    toV = getVs(cCost,omit0=F)


    #first try to allocation matrices - fail early if too big
    checkMemory=NULL;gc()
    checkMemory = try(matrix(0.0000001,nrow=length(fromV),ncol=length(toV)*1.5),silent=T)

    if(inherits(checkMemory,"try-error")){
      stop(paste("Expecting memory problems. Try reducing focal patch size or the spatial scale parameter.",checkMemory)) 
    }
    checkMemory=NULL;gc()

    #Goal here is to fail gracefully if memory requirements are excessive, before igraph crashes unceremoniously.
    #The problem is I don't know how much memory igraph really needs. Clearly larger than output matrix size, but how big?
    
    #sort( sapply(ls(),function(x){object.size(get(x))}))

    igDists = igraph::distances(sim$g,v=fromV,to=toV)

    igDists[is.infinite(igDists)]=NA

    cPatches=NULL;cCost=NULL

    rownames(igDists)=fromV
    colnames(igDists)=toV

    igDists[igDists>maxDist]=NA
    igList[[paste0("clump",ic)]]=igDists
  }

  return(igList)

})
