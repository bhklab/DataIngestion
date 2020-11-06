#' Build all tables related to samples for ToxicoDB and write to disk
#'
#' @param path `character` path to the directory where the ToxicoSet data was
#'  extracted
#' @param outDir `character` path to save table csvs to
#' @param ... `pairlist` force further arguments to be named
#'
#' @export
buildSensitivityTables <- function(path='procdata', outDir='latest', ...)
{
    # -- ensure the output directory exists
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # -- load the compound table for each tSet
    files <- list.files(file.path(path, 'sensitivity'), pattern='csv', full.names=TRUE)
    sensitivityTables <- lapply(files, fread)
    names(sensitivityTables) <-
        trimws(gsub('^.*/|.csv$', '', files))

    # -- load dependent tables
    ## TODO:: Parameterize these paths
    sample <- fread(file.path(outDir, 'sample.csv'))

    # -- build viability table
    ## FIXME:: Why are there NA sample_id in this table?
    viability <- na.omit(rbindlist(sensitivityTables)[, .(sample_id, viability)])
    ## TODO:: Start using .csvy to specify column types
    # fix column types
    viability[, sample_id := as.character(sample_id)]

    # map sampe_id to viability table
    setkeyv(sample, 'name')
    setkeyv(viability, 'sample_id')
    viability[sample, sample_id := i.id]

    # -- write to disk
    for (table in 'viability') {
        fwrite(get(table), file.path(outDir, paste0(table, '.csv')))
    }
}

if (sys.nframe() == 0) {
    library(data.table)
    path='procdata'
    outDir='latest'
    sample='sample'

}