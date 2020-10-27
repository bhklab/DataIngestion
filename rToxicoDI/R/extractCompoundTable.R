#' Construct the compound table from a `ToxicoSet` object
#'
#' @param tSet `ToxicoSet` object
#' @param outDir `character` path to save the table to. Defaults to `tempdir()`.
#' @param fileName `character` name of the output file. Defaults to `name(tSet)`.
#'
#' @return Nothing; writes to disk
#'
#' @import data.table
#' @md
#' @export
extractCompoundTable <- function(tSet, outDir=tempdir(), fileName=name(tSet)) {

    if (!is(tSet, 'ToxicoSet'))
        stop('[rPharmacoDI::extractCellTable] tSet must be a ToxicoSet object!')

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    compoundInfo <- as.data.table(drugInfo(tSet))[, 'drugid']
    setnames(compoundInfo, 'drugid', 'name')

    # process the file name
    fileName <- split(fileName, ' ')
    if (length(fileName) > 1)
        fileName <- paste(fileName[-length(fileName)], collapse='_')
    else
        fileName <- unlist(fileName)

    fwrite(compoundInfo, file=file.path(outDir, paste0(fileName, '.csv')))
}


#' Contruct the compound for for each `ToxicoSet` from a list of `ToxicoSet`s
#'
#' @param tSets `list` of `ToxicoSet` objects
#' @param outDir `character` path to save the table to. Defaults to `tempdir()`.
#'
#' @return Nothing; writes to disk
#'
#' @md
#' @export
extractAllCompoundTables <- function(tSets, outDir=tempdir()) {

    ## FIXME:: Error handling

    for (tSet in tSets) extractCompoundTable(tSet, outDir=outDir)
}