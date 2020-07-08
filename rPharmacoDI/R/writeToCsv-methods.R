#'
#'
#'
#'
#'
setMethod("writeToCsv")


#'
#'
#'
#' @import
#' @expprt
setMethod("writeToCsv", signature(object="PharmacoSet", path="character"), function(object, path){

    .writeAnnotationsToCsv(object, path)

    .writeSensitivityToCsv <- function(object, path)
        fwrite(data.table(object, rownames="row.names"), file=paste)

    .writeSensitivityToCsv()

    .writePerturbationToCsv()

    .writeDrugToCsv()

    .writeDataSetTypeToCsv(object, path)

    .writeCellToCSV(object, path)

})

.writeAnnotationToCsv <- function(object, path) {

}

