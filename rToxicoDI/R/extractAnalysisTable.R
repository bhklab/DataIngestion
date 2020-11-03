#' Construct the analysis table from a `ToxicoSet` object
#'
#' @inheritParams extractCompoundTable
#'
#' @return Nothing; writes to disk
#'
#' @import data.table
#' @md
#' @export
extractAnalysisTable <- function(tSet, outDir=tempdir(), fileName=name(tSet)) {

    if (!is(tSet, 'ToxicoSet'))
        stop(.context(), ' tSet must be a ToxicoSet object!')

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # get the data
    analysis <- computeLimmaDiffExpr(tSet)

    # rename columns
    setnames(analysis,
        old=c('gene', 'compound', 'cell', 'duration'),
        new=c('gene_id', 'compound_id', 'cell_id', 'time'),
        skip_absent=TRUE)

    fileName <- .preprocessFileName(fileName)

    fwrite(analysis, file=file.path(outDir, fileName))
}


#' Contruct the analysis table for each `ToxicoSet` from a list of `ToxicoSet`s
#'
#' @inheritParams extractAllCompoundTables
#'
#' @return Nothing; writes to disk
#'
#' @md
#' @export
extractAllAnalysisTables <- function(tSets, outDir=tempdir()) {

    if (!is.list(tSets))
        stop(.context(), 'tSets must be a list of `ToxicoSet` objects!')

    # TODO:: Can I parallelize this without crashing due to memory error?
    for (tSet in tSets) extractAnalysisTable(tSet, outDir=outDir)
}
