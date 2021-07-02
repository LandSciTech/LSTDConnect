#' @include AAAClassDefinitions.R
#' @include helperFns.R
NULL

# @name ParcConnectedness
# @rdname ParcConnectedness-class
setMethod(f = "initialize", signature = "ParcConnectedness", 
          definition = function(.Object, patches, d_ij, H_j, M_i = list(), 
                                Mbar_p = data.frame(), Mbar = data.frame(), 
                                metric = "") {
            .Object@patches <- patches
            .Object@d_ij <- d_ij
            .Object@H_j <- H_j
            .Object@M_i <- M_i
            .Object@Mbar_p <- Mbar_p
            .Object@Mbar <- Mbar
            .Object@metric <- metric
            return(.Object)
          })

#' Calculate PARC-Connectedness
#'
#'
#' @param x RasterStack, matrix or ParcConnectedness object. RasterStack of 
#'   patches, cost, and habitatValue will be passed to distanceMethod for 
#'   calculation of d_ij effective distance matrix. Or pass in precalculated 
#'   effective distance matrix d_ij or ParcConnectedness object for increased 
#'   computational efficiency.
#' @param habitatValue vector. Vector of H_j values. Indices must correspond 
#'   to rows(?) of d_ij matrix x. Ignored if x is a RasterStack or 
#'   ParcConnectedness object.
#' @param distanceMethod Character string. Options are leastCostPaths or ... 
#'   Argument is ignored if x is a matrix or ParcConnectedness object.
#' @param maxDist Number. Maximum distance between i and j in units of 
#'   RasterStack x. Ignored if x is a matrix or ParcConnectedness object.
#' @param alpha Number or vector of these. 1/alpha is in units of effective 
#'   distance. See Drielsma et al 2007 for details.
#' @param metric Character string. Options are "connectivity" or "colonization 
#'   potential" (equations 3 or 4 of Drielsma)
#' @param region RasterLayer or character string. Optional. Areas to treat as 
#'   distinct in calculations. If x is a ParcConnectedness object region should 
#'   be one or more names of the d_ij list.
#' @param memoryLimit Number. Total number of GB to allocate to d_ij list. When 
#'   size of d_ij list exceeds this limit calculations will be done, but matrices 
#'   will not be retained and returned.
#' @param stopOnMemoryLimit Boolean. If true algorithm will return error if 
#'   expected memory requirement is greater than available memory.
#' @param neighbourhood Character. 'rook','queen', or 'octagon'. 'octagon' 
#'   option is a modified version of the queen's 8 cell neighbourhood in which 
#'   diagonals weights are 2^0.5x higher than horizontal/vertical weights.
#' @return A ParcConnectedness object including patches, effective distance 
#'   d_ij matrix, habitat value H_j vector, M_i maps, average metric value for 
#'   each patch Mbar_p, data frame average metric of all patches Mbar.
#' 
#' @details
#' TO DO: add description
#' Average colonization potential will be calculated for each stratum (or patch id) 
#' in patches layer.
#'
#' cost must be >0. Set cost = 1 for euclidean distance.
#' maxDist is in map units (metres etc). See leastCostPathDistances() for details.
#'
#' Algorithm notes:
#'   The limiting step is calculation and storage of d_ij matrix the set of pixels 
#'   within a buffered clump. Given N protected pixels in a clump, and B pixels 
#'   within maxDist of the protected pixels, the size of d_ij will be N*B.
#'   For example, d_ij for a circular park with 10km diameter, maxDist=500km, 
#'   1km resolution, cost>=1: B~800000, N~80, N*B~63 million. Obviously we will 
#'   need to do something about this for maps with big park clusters - it works 
#'   fine for small examples.
#' 
#' @examples
#' # TO DO: examples
#' # TO DO: define parcConnectedness for other input signatures.
#' 
#' @export
setGeneric("parcConnectedness", 
           function(x, habitatValue = NULL, distanceMethod = "leastCostPaths", 
                    maxDist = 500, alpha = c(2, 25, 100), metric = "connectivity", 
                    region = NULL, memoryLimit = 1, stopOnMemoryLimit = T, 
                    neighbourhood = "octagon") {
             standardGeneric("parcConnectedness")}
)

#' @rdname parcConnectedness
setMethod("parcConnectedness", signature(x = "ParcConnectedness"), 
          function(x, habitatValue, distanceMethod, maxDist, alpha, metric, 
                   region, memoryLimit, stopOnMemoryLimit, neighbourhood) {
  # x=initPC
  metricOptions <- c("connectivity", "colonization potential")
  if (!is.element(metric, metricOptions)) {
    stop(paste0("Metric ", metric, " not recognized. Options are: ", paste(metricOptions, collapse = ",")))
  }
  if (class(x@d_ij) == "matrix") {
    x$d_ij <- list(clump1 = x@d_ij)
  }
  if (class(x@H_j) == "matrix") {
    x@H_j <- list(clump1 = x@H_j)
  }
  
  setMismatch <- setdiff(names(x@d_ij), names(x@H_j))
  if (length(setMismatch) > 0) {
    stop("Elements of d_ij and H_j lists do not match:", paste(setMismatch, collape = T))
  }
  
  # Make blank R_i maps for merging.
  blankMap <- x@patches
  blankMap[] <- NA
  R_i <- raster::brick(blankMap)
  names(R_i) <- paste0("alpha", alpha[[1]])
  
  if (length(alpha) > 1) {
    for (a in 2:length(alpha)) {
      R_i <- raster::addLayer(R_i, blankMap)
      names(R_i)[length(names(R_i))] <- paste0("alpha", alpha[[a]])
    }
  }
  
  # loop over clumps
  clumpSet <- setdiff(names(x@d_ij), c(NA))
  for (cl in clumpSet) {
    # cl="clump1"
    
    x@H_j[[cl]] <- x@H_j[[cl]][is.element(colnames(x@d_ij[[cl]]), colnames(x@d_ij[[cl]]))]
    
    pId <- getVs(x@patches, id = F, locs = rownames(x@d_ij[[cl]]))
    R_all <- data.frame(alpha = NA, patchId = NA, R_i = NA)
    
    for (ca in alpha) {
      # ca = alpha[1]
      cName <- paste0("alpha", ca)
      cName <- gsub("-", ".", cName, fixed = T)
      
      # scaleKernel = exponentialKernel(2/ca,cellDim=res(x@patches)[1],
      #                                 returnScale=T,negligible=10^-6)#put in the vicinity of 1
      # Exponential kernel from Hughes et al 2015 American Naturalist
      # w_ij =(2/(pi*(1/ca)^2))* exp(-ca*x@d_ij[[cl]])
      w_ij <- exp(-ca * x@d_ij[[cl]])
      
      # Drop d_ij to save memory if it is not being used.
      if ((memoryLimit == 0) && (length(clumpSet) == 1) && (length(alpha) == 1)) {
        x@d_ij <- list()
      }
      print("Memory use before bottleneck")
      print(sort(sapply(ls(), function(aa) {
        object.size(get(aa))
      })))
      w_ij <- t(w_ij)
      
      gamma_ij <- w_ij * x@H_j[[cl]]
      
      w_ij <- NULL
      # gamma_ij = t(t(w_ij)*x@H_j[[cl]])
      
      gamma_i <- colSums(gamma_ij, na.rm = T)
      gamma_ij <- NULL
      if (metric == "colonization potential") {
        R_iv <- x@H_j[[cl]] * gamma_i
      } else {
        R_iv <- gamma_i
      }
      
      #########
      # Convert R_iv vector to map
      R_ivd <- data.frame(id = as.numeric(names(R_iv)), R = R_iv)
      R_ivd <- merge(data.frame(id = 1:ncell(x@patches)), R_ivd, all.x = T)
      tempM <- blankMap
      tempM[] <- R_ivd$R
      
      R_i[[cName]][!is.na(tempM)] <- tempM[!is.na(tempM)]
      #########
      # Summarize by patches
      Rbar_pbit <- data.frame(alpha = ca, patchId = pId, R_i = R_iv)
      R_all <- rbind(R_all, Rbar_pbit)
    }
  }
  Rbar_p <- plyr::ddply(subset(R_all, !is.na(alpha)), plyr::.(alpha, patchId), 
                        plyr::summarize, Mbar_p = mean(R_i), 
                        patchArea = sum(!is.na(R_i)))
  Rbar <- plyr::ddply(Rbar_p, plyr::.(alpha), plyr::summarize, 
                      Mbar = sum(Mbar_p * patchArea) / sum(patchArea))
  
  x@M_i <- R_i
  x@Mbar_p <- Rbar_p
  x@Mbar <- Rbar
  x@metric <- metric
  if (memoryLimit == 0) {
    x@d_ij <- list()
  }
  return(x)
})

#' @rdname parcConnectedness
setMethod("parcConnectedness", signature(x = "RasterStack"), 
          function(x, habitatValue, distanceMethod, maxDist, alpha, metric, 
                   region, memoryLimit, stopOnMemoryLimit, neighbourhood) {
  # x=stack(fPatch,exCost,qSurface);distanceMethod="leastCostPaths";maxDist=maxDistHold*res(qSurface)[1];alpha=2/(dbar*res(qSurface)[1]);region=NULL;memoryLimit=0;metric="connectivity";stopOnMemoryLimit=F;neighbourhood="octagon"
  # x=exDat;distanceMethod="leastCostPaths";maxDist=25;alpha=5;region=NULL;memoryLimit=0;metric="colonization potential"
  # x=stack(testPatch,oneMap,exQuality);neighbourhood="octagon";distanceMethod="leastCostPaths";maxDist=nrow(exponentialKernel(dbar,negligible=10^-6));alpha=2/dbar;memoryLimit=0;stopOnMemoryLimit=T
  names(x) <- c("patches", "cost", "habitatValue")
  
  distMethods <- c("leastCostPaths")
  if (!is.element(distanceMethod, distMethods)) {
    stop("Distance method", distanceMethod, "not recognized. Options are:", paste(distMethods, collapse = ","))
  }
  
  if (is.null(region)) {
    region <- !is.na(x$cost)
  }
  patches <- x$patches
  x <- raster::dropLayer(x, "patches")
  patches[!region] <- NA
  
  isPatches <- cellStats(patches, "sum")
  if (isPatches == 0) {
    stop("There are no patches for parc calculation A.")
  }
  
  # To get buffer width need to scale maxDist by minimum cost -
  costMin <- raster::cellStats(x$cost, min)
  if (costMin == 0) {
    stop("Cost must be >0.")
  }
  buffDist <- maxDist / costMin
  
  bufferedPatches <- lsBuffer(patches, buffDist)
  buffClumps <- raster::clump(bufferedPatches)
  # plot(patches)
  # plot(bufferedPatches)
  bufferedPatches <- NULL
  clumps <- setdiff(raster::unique(buffClumps), c(NA, 0))
  outPC <- NULL
  # plot(x$habitatValue)
  # x$habitatValue[!is.na(buffClumps)]=100
  # plot(buffClumps)
  #  cellStats(buffClumps,stat="sum")
  #  cellStats(patches,stat="sum")
  # cellStats(patches&!is.na(x$habitatValue),stat="sum")
  for (ic in clumps) {
    # ic=clumps[1]
    
    cPatches <- patches
    cPatches[buffClumps != ic] <- NA
    cPatches[is.na(buffClumps)] <- NA
    isPatches <- cellStats(cPatches, "sum")
    if (isPatches == 0) {
      stop("There are no patches for parc calculation B.")
    }
    cCost <- x$cost
    cCost[buffClumps != ic] <- NA
    cCost[is.na(buffClumps)] <- NA
    
    # sort( sapply(ls(),function(x){object.size(get(x))}))
    if (distanceMethod == "leastCostPaths") {
      d_ij <- leastCostPathDistances(cPatches, cCost, maxDist, 
                                     bufferedPatches = !is.na(cPatches), 
                                     neighbourhood = neighbourhood)
    }
    
    
    # leastCostPathDistances includes check for exceeding memory requirements. 
    # Figure out how to reduce problem at this point.
    
    cCost <- NULL
    # sort( sapply(ls(),function(x){object.size(get(x))}))
    
    # str(d_ij[[1]])
    # plot(image(1:nrow(d_ij[[1]]), 1:ncol(d_ij[[1]]), d_ij[[1]], axes = FALSE, xlab="", ylab=""))
    H_j <- list()
    for (nn in names(d_ij)) {
      # nn="clump1"
      hv <- x$habitatValue
      hv[!is.na(cPatches) & is.na(hv)] <- 0
      colnames(d_ij[[nn]])
      H_j[[nn]] <- getVs(hv, id = F, locs = colnames(d_ij[[nn]]), omit0 = F)
      hv <- NULL
    }
    
    names(d_ij) <- paste0("clump", ic)
    names(H_j) <- paste0("clump", ic)
    initPC <- new("ParcConnectedness", patches = cPatches, d_ij = d_ij, H_j = H_j)
    d_ij <- NULL
    H_j <- NULL
    cPatches <- NULL # clear memory
    if (length(clumps) == 1) {
      x <- NULL
      patches <- NULL
      buffClumps <- NULL
      region <- NULL
    }
    # sort( sapply(ls(),function(x){object.size(get(x))}))
    # alpha=2/c(5000,1000);stopOnMemoryLimit=F
    
    
    cPC <- parcConnectedness(initPC, alpha = alpha, memoryLimit = memoryLimit, 
                             metric = metric, stopOnMemoryLimit = stopOnMemoryLimit, 
                             neighbourhood = neighbourhood)
    
    outPC <- mergeParcConnectedness(cPC, outPC, memoryLimit = memoryLimit, 
                                    clumpName = ic)
  }
  
  outPC@Mbar <- plyr::ddply(outPC@Mbar_p, plyr::.(alpha), plyr::summarize, 
                            Mbar = sum(Mbar_p * patchArea) / sum(patchArea))
  
  return(outPC)
})
