#' @useDynLib LSTDConnect
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp


#' @export
samc_cache <- function(permiability, lethality=NULL, kernel=matrix(c(0,1,0,1,1,1,0,1,0),3,3)){
  if(!is.numeric(permiability) || !is.matrix(permiability)){
    stop("permiability must be a numeric matrix")
  }
  
  if(!is.numeric(kernel) || !is.matrix(kernel)){
    stop("kernel must be a numeric matrix, and should have odd size")
  }
  
  if(!(nrow(kernel) %% 2) || !(ncol(kernel) %% 2)){
    warning("kernel should be an n*m matrix where both n and m are odd")
  }
  
  if(is.null(lethality)){
    lethality <- matrix(0, nrow(permiability), ncol(permiability))
  }else if(!is.numeric(lethality) || !is.matrix(lethality)){
    stop("lethality must be a numeric matrix, or null to default to no lethality")
  }
  
  return(cache_samc_cpp(kernel, permiability, lethality))
}


#' @export
samc_step <- function(steps = 1, cache, population, dead=NULL){
  if(!is.numeric(steps)){
    stop("steps must be numeric")
  }else if(length(steps) <= 0){
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