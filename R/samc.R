#' @useDynLib LSTDConnect
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

#' Spatial Absorbing Markov Chain (SAMC) in parallel
#' 
#' Implements the SAMC algorithm from [samc][samc::samc()].
#' 
#' @param directions Either 4 (rook) or 8 (queen), the directions of movement.
#' @param kernel Alternative to "directions", a square kernel.
#' @param resistance_na_mask Default to 0, the value to give to NAs and NaNs values 
#'     in resistance.
#' @param absorption_na_mask Default to 0, the value to give to NAs and NaNs 
#'     values in absorption.
#' @param fidelity_na_mask Default to 0, the value to give to NAs and NaNs 
#'     values in fidelity.
#' @param symmetric Default to true, whether movement is symmetrical 
#'     between cells.
#' 
#' @param occ The initial state of the model (abundance matrix).
#' @param dead Optionnal matrix containing information on deaths.
#' @param time A positive integer or a vector of positive integers representing
#'     time steps.
#' 
#' @inheritParams samc::samc
#' 
#' 
#' @param samc A list resulting from calling samc.
#' @inheritParams samc::distribution
#' 
#' @return 
#' `samc` will return a list.
#' `distribution` will also return a list.
#' 
#' @rdname samc
#' @export
samc <- function(resistance, absorption = NULL, fidelity = NULL, 
                 directions = 8, kernel = NULL, 
                 resistance_na_mask = 0, absorption_na_mask = 0, 
                 fidelity_na_mask = 0, symmetric = TRUE) {
  
  if (is.matrix(resistance)){
    if (!is.numeric(resistance)) {
      stop("'resistance' must be a numeric matrix")
    }
  } else if (class(resistance) == "RasterLayer"){
    resistance <- as.matrix(resistance)
  }
  
  if(is.null(directions) & is.null(kernel)){
    stop(paste0("You must provide one of 'directions' (num value of 4 or 8) OR",
                " 'kernel' (an odd sided matrix) "))
  }
  
  if (is.null(kernel)){
    
    if(!is.numeric(directions)){
      stop("kernel must be numeric")
    }
    
    if (directions == 4) {
      
      kernel <- matrix(
        c(
          0, 1, 0,
          1, 0, 1,
          0, 1, 0
        ), 3, 3
      )
      
    } else if (directions == 8) {
      
      kernel <- matrix(
        c(
          1 / sqrt(2), 1, 1 / sqrt(2),
          1.0000000, 0, 1.0000000,
          1 / sqrt(2), 1, 1 / sqrt(2)
        ), 3, 3
      )
      
    } else {
      
      stop("'directions' must be equal to either 4 or 8")
      
    }
    
  }
  
  if (!is.numeric(kernel)) {
    stop("kernel must be numeric")
    
  } else if (is.matrix(kernel)) {
    
    if (!(nrow(kernel) %% 2) || !(ncol(kernel) %% 2)) {
      stop("If the kernel is a matrix, it must be an n*m matrix where both n and m are odd")
    }
  } else {
    stop("kernel should be a matrix")
  }
  
  # -------------------------------------------------------------------------
  
  if (is.null(absorption)) {
    absorption <- matrix(0, nrow(resistance), ncol(resistance))
  } else if (is.numeric(absorption) && !is.matrix(absorption)) {
    absorption <- matrix(absorption, nrow(resistance), ncol(resistance))
  } else if (!is.numeric(absorption) || !is.matrix(absorption)) {
    stop("absorption must be a numeric matrix of the same dimentions as resistance, or a numeric to be used to create such a matrix, or null to default to no absorption")
  }
  
  if (!all(dim(resistance) == dim(absorption))) {
    stop("If absorption is a matrix, it must be of the same dimentions as resistance")
  }
  
  if (is.null(fidelity)) {
    fidelity <- matrix(0, nrow(resistance), ncol(resistance))
  } else if (is.numeric(fidelity) && !is.matrix(fidelity)) {
    fidelity <- matrix(fidelity, nrow(resistance), ncol(resistance))
  } else if (!is.numeric(fidelity) || !is.matrix(fidelity)) {
    stop("fidelity must be a numeric matrix of the same dimentions as resistance, or a numeric to be used to create such a matrix, or null to default to no fidelity")
  }
  
  if (!all(dim(resistance) == dim(fidelity))) {
    stop("If fidelity is a matrix, it must be of the same dimentions as resistance")
  }
  
  resistance[!is.finite(resistance)] <- resistance_na_mask
  absorption[!is.finite(absorption)] <- absorption_na_mask
  fidelity[!is.finite(fidelity)] <- fidelity_na_mask
  
  if (any(0 > resistance)) {
    stop("resistance must > 0")
  }
  
  if (any(0 > absorption || absorption > 1)) {
    stop("absorption must be in the range [0,1]")
  }
  
  if (any(0 > fidelity || fidelity > 1)) {
    stop("fidelity must be in the range [0,1]")
  }
  
  if (any(fidelity + absorption > 1)) {
    stop("fidelity+absorption must be <= 1")
  }
  
  return(cache_samc_cpp(kernel, resistance, fidelity, absorption, symmetric))
}

#' @rdname samc
#' @export
distribution <- function(samc, occ, time = 1, dead = NULL) {
  
  if (!is.numeric(time)) {
    stop("time must be numeric")
  } else if(length(time) <= 0){
    stop("time must have a length of at least 1")
  }
  
  if (is.matrix(occ)){
    if (!is.numeric(occ)) {
      stop("'occ' must be a numeric matrix")
    }
  } else if (class(occ) == "RasterLayer"){
    occ <- as.matrix(occ)
  }
  
  warned_about_rounding <- FALSE
  for (i in 1:length(time)) {
    if (time[i] %% 1 && !warned_about_rounding) {
      warned_about_rounding <- TRUE
      warning("time values should be whole numbers. The fractional part will be truncated off.")
    }
  }
  time <- as.integer(floor(time))
  
  sizes <- samc_cache_sizes_cpp(samc)
  
  if (nrow(occ) != sizes[1]) {
    stop("occ has the wrong height")
  } else if (ncol(occ) != sizes[2]) {
    stop("occ has the wrong width")
  }
  
  if (is.null(dead)) {
    dead <- matrix(0, nrow = sizes[1], ncol = sizes[2])
  } else if (nrow(dead) != sizes[1]) {
    stop("Dead has the wrong height")
  } else if (ncol(dead) != sizes[2]) {
    stop("Dead has the wrong width")
  }
  
  return(samc_step_cpp(time, samc, occ, dead))
}

# Helpers -----------------------------------------------------------------

.compare_implementations <- function(t, resistance, absorption, fidelity,
                                     occ, directions = 8) {
  
  samc_obj <- samc::samc(resistance, absorption, fidelity, 
                         tr_fun = function(x) 1 / mean(x), override = TRUE, directions = directions)
  
  acc <- matrix(0, nrow(occ), ncol(occ))
  
  for (i in 1:length(occ)) {
    acc <- acc + occ[i] * matrix(samc::distribution(samc_obj, origin = i, time = t), 
                                 nrow(occ), ncol(occ), byrow = TRUE)
  }
  
  return(acc - LSTDConnect::distribution(
    LSTDConnect::samc(resistance, fidelity, absorption, directions), 
    occ, time = 1)[["occ"]][[1]])
}
