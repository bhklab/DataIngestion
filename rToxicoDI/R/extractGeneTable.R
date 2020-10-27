#' Construct the gene table from a `ToxicoSet` object
#'
#' @inheritParams extractCompoundTable
#'
#' @return Nothing; writes to disk
#'
#' @import data.table
#' @include extractCompoundTable.R
#' @md
#' @export
extractGeneTable <- function(tSet, outDir=tempdir(), fileName=name(tSet)) {

    # handle errors
    if (!is(tSet, 'ToxicoSet'))
        stop('[rPharmacoDI::extractGeneTable] tSet must be a ToxicoSet object!')

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # get the required data from the tSet
    geneInfo <- as.data.table(featureInfo(tSet, 'rna'))[,
        c('gene_id', "Symbol", "EntrezGene.ID", "transcript_name",
        "transcript_id")]

    # rename columns by reference
    setnames(geneInfo,
        old= c("Symbol", "EntrezGene.ID", "transcript_name", "transcript_id"),
        new=c('symbol', 'entrez_gid', 'transcript_name', 'transcript_id'))

    # process the file name
    fileName <- split(fileName, ' ')
    if (length(fileName) > 1)
        fileName <- paste(fileName[-length(fileName)], collapse='_')
    else
        fileName <- unlist(fileName)

    # write to disk
    fwrite(geneInfo, file=file.path(outDir, paste0(fileName, '.csv')))
}

#' Contruct the gene table for for each `ToxicoSet` from a list of `ToxicoSet`s
#'
#' @inheritParams extractAllCompoundTables
#'
#' @return Nothing; writes to disk
#'
#' @include extractCompoundTable.R
#' @export
extractAllGeneTables <- function(tSets, outDir=tempdir()) {

    if (!is.list(tSets)) stop('\n[rToxicoDI::extractAllGeneTables] tSets must be
        a list of ToxicoSet objects')

    for (tSet in tSets) extractGeneTable(tSet, outDir=outDir)
}