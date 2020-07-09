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
setMethod("writeToCsv", signature(object="PharmacoSet", filePath="character"), function(object, path) {

    objectName <- name(object)

    .writeMolecularProfilesToCsv(molecularProfilesSlot(object), filePath, objectName)

    # annotation and datasetType
    .writeAnnotationToTxt(annotation(object), datasetType(object), filePath, objectName)

    .writeSensitivityToCsv <- function(object, filePath)
        fwrite(data.table(object, rownames="row.names"), file=paste)

    .writeSensitivityToCsv(sensitivitySlot(object), filePath, objectName)

    .writePerturbationToCsv(slot(object, 'perturbation'), filePath, objectName)

    .writeDrugToCsv(slot(object, 'drug'), filePath, objectName)



    .writeCellToCSV(object, filePath)

})

#'
#'
#'
#'
#'
#'
#'
#' @export
.writeMolecularProfilesToCsv <- function(SElist, path, objectName) {
    SEnames <- names(SElist)
    for (idx in seq_along(SElist)) {
        fileName <- paste0(objectName, "_SElong_", SEnames[[idx]], '.csv')
        message(paste0("Saving: ", fileName))
        fwrite(as(SElist[[idx]], 'data.table'), file=file.path(path, fileName))
    }
}

#' Write the annotation listt in a PharmacoSet to a text file
#'
#' Each list item is one line, lines in each item are separated with `sepChar` so they
#'    can be split into to their original stucture.
#'
#' @param annotations [`list`] Annotations as returned by the `annotaiton` function
#' @param path [`character`] The path to save the output file to.
#' @param objectName [`character`] The name of the PSet the annotations are from.
#' @param sepChar [`character`] Separator to use when pasting together multiple lines of a list item.
#'
#' @return Writes to disk, does not return.
#'
#' @keywords internal
#' @export
.writeAnnotationToTxt <- function(annotations, dsType, filePath, objectName, sepChar="|||") {
    file <- file.path(filePath, paste0(objectName, '_annotations.txt'))

    # Date Created
    dateCreated <- annotations$dateCreated

    # Session Info
    sessInfo <- paste0(capture.output(annotations$sessionInfo), collapse=sepChar)

    # Call
    creationCall <- paste0(capture.output(annotations$call), collapse=sepChar)

    # Version
    version <- annotations$version

    annots <- rbind(objectName, dateCreated, sessInfo, creationCall, version, dsType)

    write.table(annots, file=file, sep="\n", row.names=FALSE, col.names=FALSE)
}

