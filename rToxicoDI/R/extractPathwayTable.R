#'
#'
#'
#' @import data.table
#' @import biomaRt
#' @import msigdbr
#' @importFrom BiocParallel bplapply
#' @export
extractPathwayTable <- function(path='metadata/pathways_raw.txt',
    outDir='procdata')
{

    # -- ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # -- download the raw data
    mSigPathways <- getMSigDBPathways()
    CTDpathways <- getCTDpathways()

    # -- extract the pathway data
    mSigDT <- mSigPathways[grepl('^GO_.*|^REACTOME.*', gs_name),
        .(gene_symbol, gs_name)]
    setnames(mSigDT, c('gs_name'), c('pathway_id'))
    pathway <- rbind(mSigDT, CTDpathways[, .(gene_symbol, pathway_id)])

    # -- connect to ensembl to map gene symbols
    require('biomaRt')  ## FIXME:: Remove this when using as package
    ensemblRat <- useEnsembl('genes', 'rnorvegicus_gene_ensembl')
    ensemblHuman <- useEnsembl('genes', 'hsapiens_gene_ensembl')
    ratAttrs <- as.data.table(listAttributes(ensemblRat))
    humanAttrs <- as.data.table(listAttributes(ensemblHuman))

    # get symbols available for each species
    humanSymbols <- grep('symbol', humanAttrs$name, value=TRUE)
    ratSymbols <- grep('symbol', ratAttrs$name, value=TRUE)

    # function to query biomaRt
    .getBM <- function(filter, attributes, values, mart)
        getBM(filter=filter, attributes=attributes,values=values, mart=mart)

    # parallelize queries on each symbol
    humanMapping <- rbindlist(bplapply(humanSymbols, .getBM,
        attributes=c(humanSymbols, 'ensembl_gene_id'),
        values=unique(pathway$gene_symbol), mart=ensemblHuman))
    ratMapping <- rbindlist(bplapply(ratSymbols, .getBM,
        attributes=c(ratSymbols, 'ensembl_gene_id'),
        values=unique(pathway$gene_symbol), mart=ensemblRat))

    # melt so the all the symbols are in one column
    humanMapping <- melt(humanMapping, id.vars='ensembl_gene_id',
        measure.vars=humanSymbols, variable.name='symbol_type',
        value.name='symbol')
    ratMapping <- melt(ratMapping, id.vars='ensembl_gene_id',
        measure.vars=ratSymbols, variable.name='symbol_type',
        value.name='symbol')

    # merge all the symbols into a single table mapping symbol to ensembl id
    geneSymbolDT <- unique(rbind(humanMapping, ratMapping))

    # -- join gene symbol mappings to pathway table
    setkeyv(pathway, 'gene_symbol')
    setkeyv(geneSymbolDT, 'symbol')
    pathway[geneSymbolDT, gene_id := i.ensembl_gene_id]

    for (table in c('pathway'))
        fwrite(pathway, file.path(outDir, paste0(table, '.csv')))
}

if (sys.nframe() == 0) {
    library(data.table)
    library(biomaRt)
    library(BiocParallel)

    outDir='procdata'

    extractPathwayTable()
}