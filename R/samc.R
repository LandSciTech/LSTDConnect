#' @useDynLib LSTDConnect
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp


#' @export
samc_cache <- function(permiability, absorbtion=NULL, fidelity=NULL, kernel=8, permiability_na_mask = 0, absorbtion_na_mask = 0, fidelity_na_mask = 0){
  
  
  if(!is.numeric(permiability) || !is.matrix(permiability)){
    stop("permiability must be a numeric matrix")
  }
  
  if(!is.numeric(kernel)){
    stop("kernel must be numeric")
  }else if(is.matrix(kernel)){
    if(!(nrow(kernel) %% 2) || !(ncol(kernel) %% 2)){
      stop("If the kernel is a matrix, it must be an n*m matrix where both n and m are odd")
    }
  }else if(kernel == 4){
    kernel <- matrix(
      c(0,1,0,
        1,0,1,
        0,1,0), 3, 3)
  }else if(kernel == 8){
    kernel <- matrix(
      c(1/sqrt(2), 1, 1/sqrt(2),
        1.0000000, 0, 1.0000000,
        1/sqrt(2), 1, 1/sqrt(2)), 3, 3)
  }else{
    stop("If kernel is not a matrix, it must be equal to 4 or 8")
  }
  
  if(is.null(absorbtion)){
    absorbtion <- matrix(0, nrow(permiability), ncol(permiability))
  }else if(is.numeric(absorbtion) && !is.matrix(absorbtion)){
    absorbtion <- matrix(absorbtion, nrow(permiability), ncol(permiability))
  }else if(!is.numeric(absorbtion) || !is.matrix(absorbtion)){
    stop("absorbtion must be a numeric matrix of the same dimentions as permiability, or a numeric to be used to create such a matrix, or null to default to no absorbtion")
  }
  
  if(!all(dim(permiability) == dim(absorbtion))){
    stop("If absorbtion is a matrix, it must be of the same dimentions as permiability.")
  }
  
  if(is.null(fidelity)){
    fidelity <- matrix(0, nrow(permiability), ncol(permiability))
  
  }else if(is.numeric(fidelity) && !is.matrix(fidelity)){
    fidelity <- matrix(fidelity, nrow(permiability), ncol(permiability))
    
  }else if(!is.numeric(fidelity) || !is.matrix(fidelity)){
    stop("fidelity must be a numeric matrix of the same dimentions as permiability, or a numeric to be used to create such a matrix, or null to default to no fidelity")
  }
  
  if(!all(dim(permiability) == dim(fidelity))){
    stop("If fidelity is a matrix, it must be of the same dimentions as permiability.")
  }
  
  permiability[!is.finite(permiability)] <- permiability_na_mask
  absorbtion[  !is.finite(absorbtion)]   <- absorbtion_na_mask
  fidelity[    !is.finite(fidelity)]     <- fidelity_na_mask
  
  if(any(0 > absorbtion || absorbtion > 1)){
    stop("absorbtion must be in the range [0,1]")
  }
  
  if(any(0 > fidelity || fidelity > 1)){
    stop("fidelity must be in the range [0,1]")
  }
  
  if(any(fidelity+absorbtion > 1)){
    stop("fidelity+absorbtion must be <= 1")
  }
  
  return(cache_samc_cpp(kernel, permiability, fidelity, absorbtion))
}


#' @export
samc_step <- function(steps = 1, cache, population, dead=NULL){
  if(!is.numeric(steps)){
    stop("steps must be numeric")
  # }else if(length(steps) <= 0){
    stop("steps must have a length of at least 1")
  }
  
  warned_about_rounding <- FALSE
  for(i in 1:length(steps)){
    if(steps[i]%%1 && !warned_about_rounding){
      warned_about_rounding <- TRUE
      warning("steps values should be whole numbers. The fractional part will be truncated off.")
    }
  }
  steps <- as.integer(floor(steps))
  
  sizes <- samc_cache_sizes_cpp(cache)
  #print(sizes)
  # {ca->nrow, ca->ncol, ca->left_extra_cols, ca->right_extra_cols};
  
  if(nrow(population) != sizes[1]){
    stop("Population has the wrong height")
  }else if(ncol(population) != sizes[2]){
    stop("Population has the wrong width")
  }
  
  if(is.null(dead)){
    dead <- matrix(0, nrow=sizes[1], ncol=sizes[2])
  }else if(nrow(dead) != sizes[1]){
    stop("Dead has the wrong height")
  }else if(ncol(dead) != sizes[2]){
    stop("Dead has the wrong width")
  }
  
  return(samc_step_cpp(steps, cache, population, dead))
}