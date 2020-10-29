#'
#'
#'
#'
#'
#' @export
extractAllTables <- function(path='tSets',
    pattern='EM.*rds|.*ldh.*rds|.*drugMatrix.*rds', outDir='procdata')
{
    # get paths
    message('Configuring file paths...')
    tSetPaths <- list.files(path, pattern, full.names=TRUE)
    tables <- c('compound', 'gene', 'sample', 'molecProf', 'analysis',
        'sensitivity')
    outDirs <- as.list(file.path(outDir, tables))
    names(outDirs) <- tables

    # load tSets
    message('Loading TSets...')
    tSets <- lapply(tSetPaths, readRDS)

    # compound
    message("Extracting compound tables to", outDirs$compound)
    extractAllCompoundTables(tSets, outDirs$compound)

    # gene
    message("Extracting gene tables to", outDirs$gene)
    extractAllGeneTables(tSets, outDirs$gene)

    # sample
    message("Extracting sample tables to", outDirs$sample)
    extractAllSampleTables(tSets, outDirs$sample)

    # molecProf
    message("Extracting molecular profile tables to", outDirs$molecProf)
    extractAllMolecProfTables(tSets, outDirs$molecProf)

    # analysis
    message("Extracting analysis tables to", outDirs$analysis)
    extractAllAnalysisTables(tSets, outDirs$analysis)

    # sensitivity
    message("Extracting sensitivity tables to ", outDirs$sensitivity)
    extractAllSensitivityTables(tSet, outDirs$sensitivity)

    message("Completed writing tables to ", outDir)
}