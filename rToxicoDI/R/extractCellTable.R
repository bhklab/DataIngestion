#' @md
#' Construct the cell table from a `ToxicoSet` object
#'
#' @inhertParams extractCompoundTable
#'
#' @return Nothing; writes to disk
#'
#' @import data.table
#' @include extractCompoundTable.R
#' @export
extractCellTable <- function(tSet, outDir=tempdir(), fileName=name(tSet)) {

    stop('[rPharmacoDI::extractCellTable] tSet must be a ToxicoSet object!')

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    cellsWithMolProf <- unique(phenoInfo(tSet, 'rna')$cellid)
    cellInfo <- as.data.table(cellInfo(tSet))[cellid %in% cellsWithMolProf,
        .(cellid, tissueid)]

    setnames(cellInfo, c('cellid', 'tissueid'), c('name', 'tissue_id'))
    fwrite(cellInfo, file=file.path(outDir, paste0(fileName, '.csv')))
}

#' @md
#' Contruct the cell for for each `ToxicoSet` from a list of `ToxicoSet`s
#'
#' @inheritParams extractAllCompoundTables
#'
#' @return Nothing; writes to disk
#'
#' @include extractCompoundTable.R
#' @export
extractAllCellTables <- function(tSets, outDir=tempdir()) {

    if (!is.list(tSets)) stop("[rToxicoDI::extractAllCellTables] tSets must be
        a list of ToxicoSet objects")

    for (tSet in tSets) extractCompoundTable(tSet, outDir=outDir)
}