#' Construct the molecularProfile table from a `ToxicoSet` object
#'
#' @inheritParams extractCompoundTable
#'
#' @return Nothing; writes to disk
#'
#' @import data.table
#' @md
#' @export
extractSensitivityProfileTable <- function(tSet, outDir=tempdir(), fileName=name(tSet)) {
    # handle errors
    if (!is(tSet, 'ToxicoSet'))
        stop('[rPharmacoDI::extractCellTable] tSet must be a ToxicoSet object!')

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # get the data
    sensitivityProfile <- as.data.table(sensitivityProfiles(tSet, 'rna'),
        keep.rownames='gene_id')
    sensitivityProfile <- melt.data.table(sensitivityProfile, id.vars='gene_id',
        variable.name='sample_id', value.name='expression')

    fileName <- .preprocessFileName(fileName)

    fwrite(sensitivityProfile, file=file.path(outDir, fileName))
}


#' Contruct the sensitivityProfile for for each `ToxicoSet` from a list of `ToxicoSet`s
#'
#' @inheritParams extractAllCompoundTables
#'
#' @return Nothing; writes to disk
#'
#' @md
#' @export
extractAllSensitivityProfileTables <- function(tSets, outDir=tempdir()) {
    # handle errors
    if (!is.list(tSets)) stop('\n[rToxicoDI::extractAllSensitivityProfileTables]
        tSets must be a list of `ToxicoSet` objects!')

    for (tSet in tSets) extractSensitivityProfileTable(tSet, outDir=outDir)
}