#' Download pathway data from Comparative Toxicology Database
#'
#' Available data can be found here http://ctdbase.org/downloads/
#'
#' @param url `character`
#'
#' @md
#' @import data.table
#' @export
getCTDpathways <- function(url='http://ctdbase.org/reports/CTD_genes_pathways.csv.gz')
{
    CTD <- fread(url)
    setnames(CTD, c('GeneSymbol', ''))
    return(CTD)
}