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
.writeMolecularProfilesToCsv <- function(SElist, filePath, objectName) {
    SEnames <- names(SElist)
    for (idx in seq_along(SElist)) {
        fileName <- paste0(objectName, "_SElong_", SEnames[[idx]], '.csv')
        message(paste0("Saving: ", fileName))
        fwrite(as(SElist[[idx]], 'data.table'), file=file.path(filePath, fileName))
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


#' Convert a PSets' sensitivty slot into a long data.table and save to disk as a .csv
#'
#' @param sensSlot [`list`] PharmacoSet sensitivity slot, as returned by `sensitivitySlot`.
#'
#'
#' @export
.writeSensitivityToCsv <- function(sensSlot, filePath, objectName) {
    longSensDT <- .sensSlotToLong(sensSlot)
    fwrite(longSensDT, file=file.path(filePath, paste0(objectName, '_sensitivity.csv')))
}

## TODO:: Can this be generalized to perturbation as well?
#' Merge all the data in the `sensitivity` slot into a single long format `data.table`
#'
#' @param sensSlot [`list`] PharmacoSet sensitivity slot, as returned by `sensitivitySlot`.
#'
#' @return [`data.table`] Long format of data in sensitivity slot
#'
#' @export
.sensSlotToLong <- function(sensSlot) {

    # -- sensitivityRaw
    rawDTdose <- data.table(sensSlot$raw[,,1], keep.rownames='.exp_id')
    rawDTviab <- data.table(sensSlot$raw[,,2], keep.rownames='.exp_id')

    moltenRawDose <- melt.data.table(rawDTdose, id.vars='.exp_id', value.vars=colnames(sensSlot$raw[,,1]),
                                     variable.name='.sample', value.name='concentration')
    moltenRawViab <- melt.data.table(rawDTviab, id.vars='.exp_id', value.vars=colnames(sensSlot[,,2]),
                                     variable.name='.sample', value.name='viability')

    # Memory manage - not triggering gc() because it is slow, could if we run out of system memory
    rm(rawDTdose, rawDTviab)

    # Set keys then join on them
    setkey(moltenRawDose, .exp_id, .sample)
    setkey(moltenRawViab, .exp_id, .sample)
    rawDT <- merge.data.table(moltenRawDose, moltenRawViab)

    # Memory manage
    rm(moltenRawDose, moltenRawViab)

    # -- sensitivityInfo
    infoDT <- data.table(sensSlot$info, keep.rownames=".exp_id")
    colnames(infoDT)[2:ncol(infoDT)] <- paste0('info_', colnames(infoDT)[2:ncol(infoDT)])

    # -- sensitivityProfiles
    profDT <- data.table(sensSlot$profiles, keep.rownames=".exp_id")
    colnames(profDT)[2:ncol(profDT)] <- paste0('prof_', colnames(profDT)[2:ncol(profDT)])


    # -- sensNumber
    numDT <- data.table(sensSlot$n, keep.rownames='.cancer_type')
    moltenNumDT <- melt(numDT, id.vars='.cancer_type', measure.vars=colnames(sensSlot$n),
                        variable.name=".drug_id", value.name="n")
    rm(numDT)

    # -- longSensDT

    # Join sensitivityInfo with sensitivityProfiles
    setkey(infoDT, .exp_id)
    setkey(profDT, .exp_id)
    annotDT <- merge.data.table(infoDT, profDT)

    # Memory manage
    rm(infoDT, profDT)

    # Join annotDT with numDT
    setkey(annotDT, info_cellid, info_drugid)
    setkey(moltenNumDT, .cancer_type, .drug_id)
    metaDT <- merge.data.table(annotDT, moltenNumDT,
                                 by.x=c('info_cellid', 'info_drugid'),
                                 by.y=c('.cancer_type', '.drug_id'))

    # Memory manage
    rm(annotDT, moltenNumDT)

    # Join annotsDT with rawDT
    setkey(metaDT, .exp_id)
    setkey(rawDT, .exp_id)
    longSensDT <- merge.data.table(rawDT, metaDT)

    return(longSensDT)
}