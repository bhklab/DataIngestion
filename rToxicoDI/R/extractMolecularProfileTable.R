#' Construct the molecularProfile table from a `ToxicoSet` object
#'
#' @inheritParams extractCompoundTable
#'
#' @return Nothing; writes to disk
#'
#' @import data.table
#' @md
#' @export
extractMolecularProfileTable <- function(tSet, outDir=tempdir(), fileName=name(tSet)) {
    # handle errors
    if (!is(tSet, 'ToxicoSet'))
        stop('[rPharmacoDI::extractCellTable] tSet must be a ToxicoSet object!')

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # get the data
    molecularProfile <- as.data.table(molecularProfiles(tSet, 'rna'),
        keep.rownames='gene_id')
    molecularProfile <- melt.data.table(molecularProfile, id.vars='gene_id',
        variable.name='sample_id', value.name='expression')

    fileName <- .preprocessFileName(fileName)

    fwrite(molecularProfile, file=file.path(outDir, fileName))
}


#' Contruct the molecularProfile for for each `ToxicoSet` from a list of `ToxicoSet`s
#'
#' @inheritParams extractAllCompoundTables
#'
#' @return Nothing; writes to disk
#'
#' @md
#' @export
extractAllMolecularProfileTables <- function(tSets, outDir=tempdir()) {
    # handle errors
    if (!is.list(tSets)) stop('\n[rToxicoDI::extractAllMolecularProfileTables]
        tSets must be a list of `ToxicoSet` objects!')

    for (tSet in tSets) extractMolecularProfileTable(tSet, outDir=outDir)
}