#' Generic for writing an R object to one or more csv files
#'
#' @param object [`object`] An R object to dump to csv
#' @param ... [`pairlist`] Allow definition of new arguments in setMethod
#'
#' @export
setGeneric("writeToCsv", function(object, ...) standardGeneric('writeToCsv'))


#'
#'
#'
#'
#'
setMethod("writeToCsv", signature(object="PharmacoSet", path="character"), function(object, path){

    .writeAnnotationsToCsv

    .writeSensitivityToCsv

    .writeSensitivityToCsv

    .writePerturbationToCsv

    .writeDrugToCsv

    .writeDataSetTypeToCsv

    .writeCellToCSV

})