#' Extract all data needed to build the ToxicoDB database tables from a set of
#'   `ToxicoSet` objects
#'
#' @param path `character`
#'
#'
#' @md
#' @import data.table
#' @export
extractAllTables <- function(path='tSets',
    pattern='EM.*rds|.*ldh.*rds|.*drugMatrix.*rds', outDir='procdata')
{
    # get paths
    message('Configuring file paths...')
    tSetPaths <- list.files(path, pattern, full.names=TRUE)
    tables <- c('compound', 'gene', 'sample', 'molecProf', 'analysis',
        'sensitivity', 'cell', 'pathway')
    outDirs <- as.list(file.path(outDir, tables))
    names(outDirs) <- tables

    # load tSets
    message('Loading TSets...')
    tSets <- lapply(tSetPaths, readRDS)

    ## TODO:: refactor to call lapply over tables

    # compound
    message("Extracting compound tables to ", outDirs$compound)
    extractAllCompoundTables(tSets, outDirs$compound)

    # gene
    message("Extracting gene tables to ", outDirs$gene)
    extractAllGeneTables(tSets, outDirs$gene)

    # sample
    message("Extracting sample tables to ", outDirs$sample)
    extractAllSampleTables(tSets, outDirs$sample)

    # molecProf
    message("Extracting molecular profile tables to ", outDirs$molecProf)
    extractAllMolecProfTables(tSets, outDirs$molecProf)

    # analysis
    message("Extracting analysis tables to ", outDirs$analysis)
    extractAllAnalysisTables(tSets, outDirs$analysis)

    # sensitivity
    message("Extracting sensitivity tables to ", outDirs$sensitivity)
    extractAllSensitivityTables(tSets, outDirs$sensitivity)

    # cell
    message("Extracting cell tables to ", outDirs$cell)
    extractAllCellTables(tSets, outDirs$cell)

    # pathway
    message("Extracting pathway table to ", outDirs$pathway)
    extractPathwayTable(outDir=outDirs$pathway)
}

if (sys.nframe() == 0) {
    library(ToxicoDI)

    path='tSets'
    pattern='EM.*rds|.*ldh.*rds|.*drugMatrix.*rds'
    outDir='procdata'


}