#'@importClassesFrom raster RasterBrick RasterStack
setClassUnion("missingOrNULL", c("missing", "NULL"))
setClassUnion("missingOrNULLOrChar", c("missing", "NULL","character"))
setClassUnion("listOrBrickOrStack", c("list", "RasterBrick","RasterStack"))
setClassUnion("listOrMatrix", c("list", "matrix"))
setClassUnion("listOrNumeric", c("list", "numeric"))

#NOTE: Constructors for each class are defined in the R file bearing the name of the class (lower case).

#' ParcConnectedness
#'
#' PARC-Connectedness analysis information including distance d_ij matrix, habitat value H_j vector, metric M_i RasterBrick or list of vectors, data frame average metric of each patch Mbar_p, data frame average metric of all patches Mbar.
#'
#' @examples
#' #TO DO - examples
#'
#' @slot patches RasterLayer. Patch id map.
#' @slot d_ij matrix or list of these. Effective distances among locations.
#' @slot H_j vector or list of these. Vector of H_j values. Indices must correspond to rows(?) of d_ij matrix.
#' @slot M_i RasterBrick or list of vectors. Element names are alpha values.
#' @slot Mbar_p data.frame. Average metric value for each alpha and each stratum (or patch id) in patches layer.
#' @slot Mbar data.frame. Average metric for each alpha across all patches.
#' @slot metric character. Metric type.Options are "connectivity" or "colonization potential" (equations 3 or 4 of Drielsma)
#' @name ParcConnectedness-class
#' @rdname ParcConnectedness-class
#' @export ParcConnectedness
#' @importFrom methods new slot slotNames
ParcConnectedness <- setClass("ParcConnectedness", representation(patches="RasterLayer",d_ij="listOrMatrix",H_j="listOrNumeric",M_i="listOrBrickOrStack",Mbar_p="data.frame",Mbar="data.frame",metric="character"))
