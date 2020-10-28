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

    if (!is(tSet, 'ToxicoSet'))
        stop(.context(), ' tSet must be a ToxicoSet object!')

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # get the data
    sensInfo <- as.data.table(sensitivityInfo(tSet), keep.rownames='rownames')
    if (nrow(sensInfo) < 1) {  # stop if there are no profiles
        message(.context(), name(tSet), ' has no sensitivty profiles.')
        return()
    }
    sensitivityProfile <- as.data.table(sensitivityRaw(tSet)[,, 2],
        keep.rownames='rownames')

    # wide -> long format
    sensInfo <-

    # join the sensitivity annotations to the data
    setkeyv(sensitivityProfile, 'rownames')
    setkeyv(sensInfo, 'rownames')
    sensitivityProfile <- sensInfo[sensitivityProfile]


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
    if (!is.list(tSets))
        stop(.context(), 'tSets must be a list ToxicoSet objects!')

    for (tSet in tSets) extractSensitivityProfileTable(tSet, outDir=outDir)
}