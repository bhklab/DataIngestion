#' Construct the sample table from a `ToxicoSet` object
#'
#' @inheritParams extractCompoundTable
#'
#' @return Nothing; writes to disk
#'
#' @import data.table
#' @include extractCompoundTable.R
#' @md
#' @export
extractSampleTable <- function(tSet, outDir=tempdir(), fileName=name(tSet)) {

    # handle errors
    if (!is(tSet, 'ToxicoSet')) stop('[rPharmacoDI::extractSampleTable]
        tSet must be a ToxicoSet object!')

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # get the required data from the tSet
    sampleInfo <- as.data.table(phenoInfo(tSet, 'rna'))[, c('samplename',
        'drugid', 'cellid', 'dose_level', 'duration', 'individual_id',
        'concentration')]

    # rename columns by reference
    setnames(sampleInfo,
        old=c('samplename', 'drugid', 'cellid', 'dose_level', 'duration',
            'individual_id'),
        new=c('name', 'compound_id', 'cell_id', 'dose', 'time', 'replicate'))

    fileName <- .preprocessFileName(fileName)

    # write to disk
    fwrite(sampleInfo, file=file.path(outDir, fileName))
}

#' Contruct the sample table for for each `ToxicoSet` from a list of `ToxicoSet`s
#'
#' @inheritParams extractAllCompoundTables
#'
#' @return Nothing; writes to disk
#'
#' @include extractCompoundTable.R
#' @md
#' @export
extractAllSampleTables <- function(tSets, outDir=tempdir()) {

    if (!is.list(tSets)) stop('\n[rToxicoDI::extractAllSampleTables] tSets must
        be a list of `ToxicoSet` objects!')

    for (tSet in tSets) extractSampleTable(tSet, outDir=outDir)
}