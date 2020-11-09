#'
#'
#'
#' @import data.table
#' @import biomaRt
#' @importFrom BiocParallel bplapply
#' @export
extractPathwayTable <- function(path='metadata/pathways_raw.txt',
    outDir='procdata')
{

    # ensure the save directory exits
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    ## TODO:: Refactor function length!
    ## TODO:: Refactor old code for style

    # -- read in file lines a stings
    lines <- readLines(path)

    # -- split the lines into lists on spaces
    line_list <- strsplit(lines, split = ' ')

    # -- split the lines into vectors on tabs
    .strSplitTab <- function(x) strsplit(x, split = "\t")
    line_list <- sapply(unlist(line_list), .strSplitTab)
    # Name the list items by pathway
    names(line_list) <- vapply(seq_along(line_list),
                               function(x) {line_list[[x]][1] },
                               FUN.VALUE = character(1))
    # Assign the vectors to a list of pathways, dropping the pathway name and url
    pathway_vectors <- lapply(names(line_list),
        function(name) { line_list[[name]][c(-1, -2)] })
    # Name the pathways
    names(pathway_vectors) <- names(line_list)
    rm(line_list)

    ## TODO:: End of old code

    # -- add additional pathway data
    pathwayL <- lapply(pathway_vectors, as.data.table)
       .addColumn <- function(x, colName, value) {  # reference symantics
        x[,  newCol := value]
        setnames(x, 'newCol', colName)
        x
    }
    pathwayL <- mapply(.addColumn, x=pathwayL,
        colName='pathway', value=names(pathwayL), SIMPLIFY=FALSE)
    pathway <- rbindlist(pathwayL)
    setnames(pathway, 'V1', 'symbol')

    # -- add id column to pathway
    uniquePathway <- unique(pathway[, 'pathway'])
    uniquePathway[, id := seq_len(.N)]

    setkeyv(uniquePathway, 'pathway')
    setkeyv(pathway, 'pathway')
    pathway[uniquePathway, id := i.id]

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
        attributes=c(humanSymbols, 'ensembl_gene_id'), values=pathway$symbol,
        mart=ensemblHuman))
    ratMapping <- rbindlist(bplapply(ratSymbols, .getBM,
        attributes=c(ratSymbols, 'ensembl_gene_id'), values=pathway$symbol,
        mart=ensemblRat))

    # melt so the all the symbols are in one column
    humanMapping <- melt(humanMapping, id.vars='ensembl_gene_id',
        measure.vars=humanSymbols, variable.name='symbol_type',
        value.name='symbol')
    ratMapping <- melt(ratMapping, id.vars='ensembl_gene_id',
        measure.vars=ratSymbols, variable.name='symbol_type',
        value.name='symbol')

    # merge all the symbols into a single table mapping symbol to ensembl id
    geneSymbolDT <- rbind(humanMapping, ratMapping)

    # join with pathways to get the mapping
    setkeyv(pathway, 'symbol')
    setkeyv(geneSymbolDT, 'symbol')
    pathway[geneSymbolDT, ensembl_gene_id := i.ensembl_gene_id]

    ## TODO:: Parse drug bank pathways for use in rToxicoDI
    #drugBankPathways <- fread(file.path('metadata', 'drug_bank_pathways.csv'))

    for (table in c('pathway'))
        fwrite(pathway, file.path(outDir, paste0(table, '.csv')))
}

if (sys.nframe() == 0) {
    library(data.table)
    library(biomaRt)
    library(BiocParallel)
    path='metadata/pathways_raw.txt'
    outDir='procdata'

    extractPathwayTable()
}