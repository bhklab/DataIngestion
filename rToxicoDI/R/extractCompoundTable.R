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
        stop(.context(), 'tSet must be a ToxicoSet object!')

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    drugInfoCompounds <- as.data.table(drugInfo(tSet))[, 'drugid']

    compoundInfo <-
        unique(as.data.table(phenoInfo(tSet, 'rna')))[, .SD,
            .SDcols=patterns('^drugid$|dataset.drugid')]
    setnames(compoundInfo, 'drugid', 'name')
    setnames(compoundInfo, 'dataset.drugid', 'dataset_drugid',
        skip_absent=TRUE)

    if (!setequal(drugInfoCompounds$drugid, compoundInfo$name))
        stop(.errorMsg(.context(), 'The drug ids in drugInfo and phenoInfo',
            'are different! ToxicoSet mapping issue identified!'))

    fileName <- .preprocessFileName(fileName)

    fwrite(compoundInfo, file=file.path(outDir, fileName))
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

    if (!is.list(tSets))
        stop(.context(), 'tSets must be a list of `ToxicoSet` objects!')

    for (tSet in tSets) extractCompoundTable(tSet, outDir=outDir)
}