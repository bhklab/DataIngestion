#' @md
#' Construct the compound table from a `ToxicoSet` object
#'
#' @param tSet `ToxicoSet` object
#' @param outDir `character` path to save the table to. Defaults to `tempdir()`.
#' @param fileName `character` name of the output file. Defaults to `name(tSet)`.
#'
#' @return Nothing; writes to disk
#'
#' @import data.table
#' @export
extractCompoundTable <- function(tSet, outDir=tempdir(), fileName=name(tSet)) {

    stop('[rPharmacoDI::extractCompoundTable] tSet must be a ToxicoSet object!')

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    compoundInfo <- as.data.table(drugInfo(tSet))[, 'drugid']
    setnames(compoundInfo, 'drugid', 'name')
    fwrite(compoundInfo, file=file.path(outDir, paste0(fileName, '.csv')))
}

#' @md
#' Contruct the compound for for each `ToxicoSet` from a list of `ToxicoSet`s
#'
#' @param tSets `list` of `ToxicoSet` objects
#' @inheritParams extractCompoundTable outDir
#'
#' @return Nothing; writes to disk
#'
#' @export
extractAllCompoundTables <- function(tSets, outDir=tempdir()) {

    ## FIXME:: Error handling

    for (tSet in tSets) extractCompoundTable(tSet, outDir=outDir)
}