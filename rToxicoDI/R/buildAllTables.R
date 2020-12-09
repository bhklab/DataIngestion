#' Perform joins on all of the individual database tables between all ToxicoSets
#' 
#' Outputs the final database tables to the 'latest' directory
#'
#' @md
#' @export
buildAllTables <- function() {

    buildCompoundTables()

    buildCellTables()

    buildSampleTables()

    buildMolecProfTables()

    buildSensitivityTables()

    # Current this function needs to be run last, because it remaps names
    #   for both datasets and compounds
    buildPathwayTables()
}