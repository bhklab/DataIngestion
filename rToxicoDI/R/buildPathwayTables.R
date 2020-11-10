#' Build all tables related to samples for ToxicoDB and write to disk
#'
#' @param path `character` path to the directory where the ToxicoSet data was
#'  extracted
#' @param outDir `character` path to save table csvs to
#' @param ... `pairlist` force further arguments to be named
#'
#' @export
buildPathwayTables <- function(path='procdata', outDir='latest', ...)
{
    # -- ensure the output directory exists
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # -- load the compound table for each tSet
    files <- list.files(file.path(path, 'pathway'), pattern='csv',
        full.names=TRUE)

    # --
    whichPathwayPath <- grepl('pathway.csv', files)
    pathway <- fread(files[whichPathwayPath])
    files <- files[!whichPathwayPath]

    # --
    pathwayStats <- lapply(files, fread)
    names(pathwayStats) <- gsub( '.*/|.csv', '', files)


}

if (sys.nframe() == 0) {
    library(data.table)
    library()
    path='procdata'
}