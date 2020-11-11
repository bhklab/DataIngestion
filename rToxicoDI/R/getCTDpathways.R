#' Download pathway data from Comparative Toxicology Database
#'
#' Available data can be found here http://ctdbase.org/downloads/
#'
#' @param url `character` Url of the compressed .csv file
#'
#' @md
#' @import data.table
#' @export
getCTDpathways <- function(url='http://ctdbase.org/reports/CTD_genes_pathways.csv.gz')
{
    suppressWarnings({ CTD <- fread(url) })
    setnames(CTD, c('gene_symbol', 'ncbi_gene', 'pathway_name', 'pathway_id'))
    return(CTD)
}