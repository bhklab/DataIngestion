#' Download all ToxicoSet objects using the ToxicoGx R package.
#' 
#' @param saveDir `character` Path to save the tSets to. Defaults to 'tSets'.
#' 
#' @importFrom ToxicoGx availableTSets downloadTSet
#' @export
downloadAllTSets <- function(saveDir='tSets') {
    tSetNames <- availableTSets(saveDir)$ToxicoSet.Name

    for (name in tSetNames) x <- downloadTSet(name, saveDir)
}