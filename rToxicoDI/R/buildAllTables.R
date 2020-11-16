#'
#'
#'
#'
#'
buildAllTables <- function() {

    buildCompoundTables()

    buildCellTables()

    buildSampleTables()

    buildMolecProfTables()

    buildSampleTables()

    # Current this function needs to be run last, because it remaps names
    #   for both datasets and compounds
    buildPathwayTables()
}