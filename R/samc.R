#' @useDynLib LSTDConnect
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp

#' @export
hello_world_rcpp <- function(){
  rcpp_hello_world_cpp()
}