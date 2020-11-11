#' Get pathway annotations from MSigDB
#'
#' @param species `character` Which species to retrieve the pathway annotations
#'   for.
#' @param all `logical` Should the data for all species be retrieved?
#'
#' @return A `data.table` mapping gene symbol to pathway ids
#'
#' @md
#' @import data.table
#' @import msigdbr
#' @export
getMSigDBPathways <- function(species=c('Homo sapiens', 'Rattus norvegicus'),
    all=TRUE)
{
    .getMSigAsDT <- function(...) as.data.table(msigdbr(...))

    if (all) {
        dataL <- lapply(species, .getMSigAsDT)
        mSigDT <- rbindlist(dataL)
        return(mSigDT)
    } else if (length(species) > 1) {
        species <- match.arg(species)
    }
    mSigDT <- .getMSigAsDT(species=species)
    return(mSigDT)
}