#' @importClassesFrom raster RasterBrick RasterStack
#' @importFrom utils object.size read.csv write.csv
setClassUnion("missingOrNULL", c("missing", "NULL"))
setClassUnion("missingOrNULLOrChar", c("missing", "NULL", "character"))
setClassUnion("listOrBrickOrStack", c("list", "RasterBrick", "RasterStack"))
setClassUnion("listOrMatrix", c("list", "matrix"))
setClassUnion("listOrNumeric", c("list", "numeric"))
