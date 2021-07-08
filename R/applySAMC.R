#' SAMC disperser distribution
#'
#' @param occurrence RasterLayer or matrix. Initial disperser distribution.
#' @param absorption RasterLayer or matrix. Mortality should be between 0 and 1. 
#'   Optional if samcObj is specified.
#' @param resistance RasterLayer or matrix. Resistance should be greater than 0. 
#'   Optional if samcObj is specified.
#' @param t Numeric. Number of timesteps. Only one of t and d should be specified.
#' @param d Numeric. Mean displacement given resistance of 1. In the units of 
#'   occurrence.
#' @param patches RasterLayer. Optional. Used to mask return layer if provided.
#' @param samcObj samc object. Optional. If supplied occurrence and absorption 
#'   arguments will be ingored. Use to increase computational efficiency if 
#'   interested in results for more than one occurrence layer or t/d with same 
#'   resistance and absorption parameters.
#'   
#' @details
#' A wrapper around distribution function of samc R package. Allows users to 
#' specify a mean displacement distance rather than a time.
#'   
#' @examples
#' # TODO: examples
#' 
#' @return 
#' A disperser distribution RasterLayer or matrix. If occurrence is a RasterLayer, 
#' the returned disperser distribution will also be a RasterLayer.
#' 
#' @export
setGeneric("applySAMC", function(occurrence, absorption = NULL, resistance = NULL, 
                                 t = NULL, d = NULL, patches = NULL, samcObj = NULL) {
  standardGeneric("applySAMC")}
)

#' @rdname applySAMC
#' @export
setMethod("applySAMC", signature(occurrence = "matrix"), 
          function(occurrence, absorption, resistance, t, d, patches, samcObj) {
            # d=NULL;patches=NULL;samcObj=NULL;occurrence=as.matrix(occurrence)
            if (is.null(t) & is.null(d)) {
              stop("Specify t or d.")
            }
            if (!is.null(t) & !is.null(d)) {
              stop("Specify only one of t or d.")
            }
            
            if (is.null(samcObj)) {
              if (class(absorption) == "RasterLayer") {
                absorption <- raster::as.matrix(absorption)
              }
              if (class(resistance) == "RasterLayer") {
                resistance <- raster::as.matrix(resistance)
              }
              if (!is.matrix(absorption)) {
                stop("absorption should be a RasterLayer or matrix.")
              }
              if (!is.matrix(resistance)) {
                stop("resistance should be a RasterLayer or matrix.")
              }
              
              if ((min(absorption) <= 0) | (max(absorption) > 1)) {
                stop("absorption should be between 0 and 1.")
              }
              if (min(resistance) <= 0) {
                stop("resistance should be greater than 0.")
              }
              
              samcObj <- samc::samc(resistance, absorption, tr_fun = function(x) 1 / mean(x))
            }
            if (!is.null(d)) {
              stop("TO DO: create & use lookup table to find t that most closely corraster::responds to selected dbar.")
              cellDim <- raster::res(quality)[1]
              dbar <- d / cellDim
              # TO DO: create & use lookup table to find t that most closely corresponds to selected dbar.
            }
            
            trMap <- samc::distribution(samcObj, occurrence, time = t)
            if (!is.null(patches)) {
              patches <- as.matrix(patches)
              trMap[(patches == 0) | is.na(patches)] <- NA
            }
            
            return(trMap)
          })

#' @rdname applySAMC
#' @export
setMethod("applySAMC", signature(occurrence = "RasterLayer"), 
          function(occurrence, absorption, resistance, t, d, patches, samcObj) {
            # d=NULL;patches=NULL;samcObj=NULL
            if (is.null(samcObj)) {
              if (class(absorption) == "RasterLayer") {
                absorption <- raster::as.matrix(absorption)
              }
              if (class(resistance) == "RasterLayer") {
                resistance <- raster::as.matrix(resistance)
              }
              if (!is.matrix(absorption)) {
                stop("absorption should be a RasterLayer or matrix.")
              }
              if (!is.matrix(resistance)) {
                stop("resistance should be a RasterLayer or matrix.")
              }
              
              if ((min(absorption) <= 0) | (max(absorption) > 1)) {
                stop("absorption should be between 0 and 1.")
              }
              if (min(resistance) <= 0) {
                stop("resistance should be greater than 0.")
              }
              
              samcObj <- samc::samc(resistance, absorption, tr_fun = function(x) 1 / mean(x))
            }
            
            trMap <- applySAMC(raster::as.matrix(occurrence), t = t, d = d, patches = patches, samcObj = samcObj)
            
            trMap <- samc::map(samcObj, trMap)
            
            oMap <- occurrence
            oMap@data@values <- trMap@data@values
            return(oMap)
          })
