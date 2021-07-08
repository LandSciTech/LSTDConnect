#' Manipulate costs information
#' 
#' Two functions to manipulate the information related to costs
#' 
#' @param composites TODO
#' @param landclasses TODO
#' @param outDir TODO
#' @param outName TODO
#' @param costOptions TODO
#' @param costScns TODO
#' @param extraProtectionValue TODO
#' @param exemptionClasses TODO
#' @param costFile TODO
#' 
#' @return 
#' TODO
#' 
#' 
#' @rdname costs
#' @export
setCosts <- function(composites, landclasses, outDir, outName, costOptions, 
                     costScns = NULL, extraProtectionValue = c(), 
                     exemptionClasses = c("urban", "cropland", "cutblock")) {
  
  # costScns=c("humanFootprintVaried.waterOmit");exemptionClasses=c("urban","cropland","cutblock");extraProtectionValue=c()
  if (is.null(costScns)) {
    costScns <- setdiff(names(costOptions), c("Landclass", "speciesType"))
    print(costScns)
  }
  costOptions$speciesType[is.na(costOptions$speciesType)] <- ""
  
  outPaths <- list()
  lcids <- unique(subset(landclasses, select = c("Landclass", "LCID")))
  outCostValues <- subset(costOptions, select = c("Landclass", "speciesType"))
  for (cc in costScns) {
    # cc=costScns[1]
    if (!is.element(cc, names(costOptions))) {
      stop("Cost scenario ", cc, " not recognized. Options are: ", 
           paste(setdiff(names(costList$costOptions), c("Landclass", "speciesType")), collapse = ","))
    }
    
    if (is.element(cc, names(extraProtectionValue))) {
      protectionLevels <- extraProtectionValue[[cc]]
    } else {
      protectionLevels <- c(0)
    }
    
    aCosts <- subset(costOptions, select = c("Landclass", "speciesType", cc))
    names(aCosts)[names(aCosts) == cc] <- "cost"
    aCosts <- subset(aCosts, !is.na(aCosts$cost))
    
    spTypes <- unique(aCosts$speciesType)
    for (ee in protectionLevels) {
      for (sp in spTypes) {
        # ee=protectionLevels[1];sp=spTypes[1]
        if (ee != 0) {
          outNameC <- paste0(c(cc, ee), collapse = ".")
        } else {
          outNameC <- paste0(c(cc), collapse = ".")
        }
        
        if (!is.na(sp) && (sp != "")) {
          outNameC <- paste0(c(outNameC, sp), collapse = ".")
        }
        print(paste("Calculating", outNameC, "costs."))
        
        cCosts <- subset(aCosts, aCosts$speciesType == sp, select = c("Landclass", "cost"))
        
        # apply extra degradation to non-protected natural areas
        if (ee != 0) {
          ee <- strsplit(ee, split = "_")[[1]]
          eeMax <- cCosts$cost[cCosts$Landclass == ee[1]]
          cCosts$NPCost <- cCosts$cost + (eeMax - cCosts$cost) * as.numeric(ee[2]) / 100
          cCosts$cost[!grepl(" P", cCosts$Landclass, fixed = T) & 
                        (!is.element(cCosts$Landclass, exemptionClasses))] <- 
            cCosts$NPCost[!grepl(" P", cCosts$Landclass, fixed = T) & 
                            (!is.element(cCosts$Landclass, exemptionClasses))]
          cCosts$NPCost <- NULL
        }
        # get LCIDs
        oddLCs <- setdiff(cCosts$Landclass, lcids$Landclass)
        if (length(oddLCs) > 0) {
          stop("Landcover classes in costOptions not recognized: ", paste(oddLCs, sep = ","))
        }
        cCosts <- merge(cCosts, lcids)
        rcl <- subset(cCosts, select = c("LCID", "cost"))
        
        outCosts <- composites
        outCosts <- raster::subs(outCosts, rcl)
        outPaths[[outNameC]] <- paste0(outDir, "/", outName, "_", outNameC, ".RData")
        save(outCosts, file = outPaths[[outNameC]])
        
        # add costs to output table
        names(cCosts)[names(cCosts) == "cost"] <- outNameC
        cCosts$speciesType <- sp
        outCostValues <- merge(outCostValues, cCosts, all.x = T)
      }
    }
  }
  
  write.csv(outCostValues, paste0(outDir, "/", outName, ".csv"))
  return(outPaths)
}

#' @rdname costs
#' @export
getCostOptions <- function(costFile) {
  costOptions <- read.csv(costFile, stringsAsFactors = F)
  costOptionDescriptions <- costOptions[1, ]
  costOptions <- costOptions[2:nrow(costOptions), ]
  
  optionCols <- setdiff(names(costOptions), c("Landclass", "speciesType"))
  
  for (o in optionCols) {
    costOptions[[o]] <- as.numeric(costOptions[[o]])
  }
  
  return(list(costOptions = costOptions, costOptionDescriptions = costOptionDescriptions))
}
