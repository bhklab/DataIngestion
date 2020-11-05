#' Construct the cell table from a `ToxicoSet` object
#'
#' @inheritParams extractCompoundTable
#'
#' @return Nothing; writes to disk
#'
#' @import data.table
#' @include extractCompoundTable.R
#' @md
#' @export
extractCellTable <- function(tSet, outDir=tempdir(), fileName=name(tSet)) {

    # handle errors
    if (!is(tSet, 'ToxicoSet'))
        stop(.context(), 'tSet must be a ToxicoSet object!')

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    cellsWithMolProf <- unique(phenoInfo(tSet, 'rna')$cellid)
    cellInfo <- as.data.table(cellInfo(tSet))[cellid %in% cellsWithMolProf,
        .(cellid, tissueid)]
    ## FIXME:: Need to add species to cellInfo in all tSets
    cellInfo$species <- if (length(cellInfo$species) > 0) cellInfo$species else
        'R.norvegicus'

    # rename columns by reference
    setnames(cellInfo, c('cellid', 'tissueid'), c('name', 'tissue_id'))

    # process the file name
    fileName <- .preprocessFileName(fileName)

    # save to csv
    fwrite(cellInfo, file=file.path(outDir, fileName))
}

#' Contruct the cell for for each `ToxicoSet` from a list of `ToxicoSet`s
#'
#' @inheritParams extractAllCompoundTables
#'
#' @return Nothing; writes to disk
#'
#' @include extractCompoundTable.R
#' @md
#' @export
extractAllCellTables <- function(tSets, outDir=tempdir()) {

    if (!is.list(tSets))
        stop(.context(), 'tSets must be a list of ToxicoSet objects!')

    for (tSet in tSets) extractCellTable(tSet, outDir=outDir)
}