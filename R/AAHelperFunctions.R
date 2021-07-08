#' TODO
#' 
#' TODO
#' 
#' @param sim TODO
#' @param neighbourhood TODO
#' 
#' @return 
#' TODO
#' 
#' @import data.table
#' @export
roadCLUS.getGraph <- function(sim, neighbourhood) {
  
  # Assign pronouns to NULL to avoid CRAN notes
  id <- NULL
  
  ### Set the grpah which determines least cost paths
  # Creates a graph (sim$g) in inititation phase which can be updated and 
  # solved for paths
  sim$paths.v <- NULL
  #------prepare the cost surface raster
  ras.matrix <- raster::as.matrix(sim$costSurface) # get the cost surface as a matrix using the raster package
  
  weight <- c(t(ras.matrix)) # transpose then vectorize which matches the same order as adj
  weight <- data.table::data.table(weight) # convert to a data.table - faster for large objects than data.frame
  # weight<-data.table(getValues(sim$costSurface)) #Try
  weight[, id := seq_len(.N)] # get the id for ther verticies which is used to merge with the edge list from adj
  
  #------get the adjacency using SpaDES function adj
  # rooks case
  if (!is.element(neighbourhood, c("rook", "octagon", "queen"))) {
    stop("neighbourhood type not recognized")
  }
  
  edges <- SpaDES.tools::adj(returnDT = TRUE, numCol = ncol(ras.matrix), 
                             numCell = ncol(ras.matrix) * nrow(ras.matrix), 
                             directions = 4, 
                             cells = 1:as.integer(ncol(ras.matrix) * nrow(ras.matrix)))
  edges <- data.table::data.table(edges)
  # edges[from < to, c("from", "to") := .(to, from)]
  edges[edges$from < edges$to, ] <- edges[edges$from < edges$to, c("to", "from")]
  edges <- unique(edges)
  edges.w1 <- merge(x = edges, y = weight, by.x = "from", by.y = "id") # merge in the weights from a cost surface
  data.table::setnames(edges.w1, c("from", "to", "w1")) # reformat
  edges.w2 <- data.table::setDT(merge(x = edges.w1, y = weight, by.x = "to", by.y = "id")) # merge in the weights to a cost surface
  data.table::setnames(edges.w2, c("from", "to", "w1", "w2")) # reformat
  edges.w2$weight <- (edges.w2$w1 + edges.w2$w2) / 2 # take the average cost between the two pixels
  
  if (neighbourhood == "rook") {
    edges.weight <- edges.w2
  } else {
    # bishop's case - multiply weights by 2^0.5
    if (neighbourhood == "octagon") {
      mW <- 2^0.5
    } else {
      mW <- 1
    }
    weight$weight <- weight$weight * mW
    edges <- SpaDES.tools::adj(returnDT = TRUE, numCol = ncol(ras.matrix), 
                               numCell = ncol(ras.matrix) * nrow(ras.matrix), 
                               directions = "bishop", 
                               cells = 1:as.integer(ncol(ras.matrix) * nrow(ras.matrix)))
    edges <- data.table::data.table(edges)
    
    # edges[from < to, c("from", "to") := .(to, from)]
    edges[edges$from < edges$to, ] <- edges[edges$from < edges$to, c("to", "from")]
    edges <- unique(edges)
    
    edges.w1 <- merge(x = edges, y = weight, by.x = "from", by.y = "id") # merge in the weights from a cost surface
    data.table::setnames(edges.w1, c("from", "to", "w1")) # reformat
    
    edges.w3 <- data.table::setDT(merge(x = edges.w1, y = weight, by.x = "to", by.y = "id")) # merge in the weights to a cost surface
    data.table::setnames(edges.w3, c("from", "to", "w1", "w2")) # reformat
    edges.w3$weight <- (edges.w3$w1 + edges.w3$w2) / 2 # take the average cost between the two pixels
    
    #------get the edges list
    edges.weight <- rbind(edges.w2, edges.w3)
  }
  
  edges.weight <- edges.weight[stats::complete.cases(edges.weight), c(1:2, 5)] # get rid of NAs caused by barriers. Drop the w1 and w2 costs.
  edges.weight[, id := seq_len(.N)] # set the ids of the edge list. Faster than using as.integer(row.names())
  
  #------make the graph
  sim$g <- igraph::graph.edgelist(as.matrix(edges.weight)[, 1:2], dir = FALSE) # create the graph using to and from columns. Requires a matrix input
  igraph::E(sim$g)$weight <- as.matrix(edges.weight)[, 3] # assign weights to the graph. Requires a matrix input
  
  #------clean up
  rm(edges.w1, edges.w2, edges.w3, edges, weight, ras.matrix) # remove unused objects
  gc() # garbage collection
  return(invisible(sim))
}

# Helpers -----------------------------------------------------------------

stripX <- function(x, parts = 1, sep = "_") {
  # x=testSet$landscape[1]
  xparts <- strsplit(x, sep)[[1]]
  xparts <- xparts[(parts + 1):length(xparts)]
  return(paste(xparts, collapse = sep))
}

mergeParcConnectedness <- function(cPC, outPC = NULL, memoryLimit, clumpName = "c") {
  # outPC=NULL;memoryLimit=1;clumpName=ic
  # outPC=cPC
  cPC@Mbar_p$clump <- clumpName
  cPC@Mbar$clump <- clumpName
  
  if (is.null(outPC)) {
    return(cPC)
  }
  
  outPC@patches <- raster::merge(cPC@patches, outPC@patches)
  if (as.numeric(object.size(outPC@d_ij)) / 10^9 < memoryLimit) {
    # names(cPC@d_ij)="clump2"
    outPC@d_ij <- c(outPC@d_ij, cPC@d_ij)
    outPC@H_j <- c(outPC@H_j, cPC@H_j)
  } else {
    addNames <- names(cPC@d_ij)
    addBit <- list("Memory limit exceeded.")
    for (adn in addNames) {
      names(addBit) <- adn
      outPC@d_ij <- c(outPC@d_ij, adn)
      outPC@H_j <- c(outPC@H_j, adn)
    }
  }
  
  addAlphas <- names(cPC@M_i)
  for (aName in addAlphas) {
    # aName = addAlphas[1]
    if (!is.element(aName, names(outPC@M_i))) {
      outPC@M_i <- raster::addLayer(cPC@M_i[[aName]])
    } else {
      outPC@M_i[[aName]] <- raster::merge(outPC@M_i[[aName]], cPC@M_i[[aName]])
    }
  }
  
  outPC@Mbar_p <- rbind(outPC@Mbar_p, cPC@Mbar_p)
  outPC@Mbar <- rbind(outPC@Mbar, cPC@Mbar)
  return(outPC)
}

getVs <- function(inMap, id = T, locs = NULL, omit0 = T) {
  # TO DO: likely there is a more efficient straightforward way to do this. Fix if it matters.
  # inMap = x$patches; id = F; locs = names(R_iv);omit0=T
  
  ras.matrix <- raster::as.matrix(inMap)
  weight <- c(t(ras.matrix))
  weight <- data.frame(weight = weight)
  weight$id <- seq_len(nrow(weight))
  subset(weight, weight > 0)
  weight <- subset(weight, !is.na(weight))
  if (omit0) {
    weight <- subset(weight, weight > 0)
  }
  if (!is.null(locs)) {
    # locs = c(3,7)
    locs <- as.numeric(locs)
    missingLocs <- setdiff(locs, weight$id)
    if (length(missingLocs) > 0) {
      stop("inMap does not include these selected locations:", paste(missingLocs, collapse = ","))
    }
    weight <- merge(data.frame(id = locs), weight)
  }
  if (id) {
    return(weight$id)
  } else {
    return(weight$weight)
  }
}
