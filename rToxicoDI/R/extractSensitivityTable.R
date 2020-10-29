#' Construct the molecularProfile table from a `ToxicoSet` object
#'
#' @inheritParams extractCompoundTable
#'
#' @return Nothing; writes to disk
#'
#' @import data.table
#' @md
#' @export
extractSensitivityTable <- function(tSet, outDir=tempdir(), fileName=name(tSet)) {

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
    sensProfile <- as.data.table(sensitivityRaw(tSet)[,, 2],
        keep.rownames='rownames')

    # wide -> long format
    idCols <- setdiff(colnames(sensInfo), c("Control", "Low", "Middle", "High"))
    sensInfo <- melt(sensInfo, id.vars=idCols, variable.name='dose',
        value.name='sample_id')
    ## TODO:: Infer the dose levels based on the values of dose, instead of assuming they are in order
    setnames(sensProfile, c('doses1', 'doses2', 'doses3'),
        c('Low', 'Middle', 'High'), skip_absent=TRUE)
    sensProfile <- melt(sensProfile, id.vars='rownames',
        variable.name='dose', value.name='viability')

    # join the sensitivity annotations to the data
    setkeyv(sensProfile, c('rownames', 'dose'))
    setkeyv(sensInfo, c('rownames', 'dose'))
    sensProfile <- sensInfo[sensProfile, !'rownames']
    setnames(sensProfile,
        old=c('cellid', 'drugid', 'duration_h'),
        new=c('cell_id', 'compund_id', 'time'),
        skip_absent=TRUE)

    fileName <- .preprocessFileName(fileName)

    fwrite(sensProfile, file=file.path(outDir, fileName))
}


#' Contruct the sensitivityProfile for for each `ToxicoSet` from a list of
#'   `ToxicoSet`s
#'
#' @inheritParams extractAllCompoundTables
#'
#' @return Nothing; writes to disk
#'
#' @md
#' @export
extractAllSensitivityTables <- function(tSets, outDir=tempdir()) {
    # handle errors
    if (!is.list(tSets))
        stop(.context(), 'tSets must be a list ToxicoSet objects!')

    for (tSet in tSets) extractSensitivityProfileTable(tSet, outDir=outDir)
}