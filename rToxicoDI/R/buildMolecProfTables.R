#' Build all tables related to molecular profiles for ToxicoDB and write to disk
#'
#' **WARNING**: This function uses a lot of RAM. It is not recommended to run
#'   this on systems with less the 32 GB of memory, otherwise it may crash
#'   your system.
#'
#' @param path `character` path to the directory where the ToxicoSet data was
#'  extracted.
#' @param outDir `character` path to save table csvs to.
#' @param ... `pairlist` force subsequent arguments to be named.
#' @param molecProf `character` Name of molecular profiles directory within
#'   the `path` directory. Default is parameter name (i.e., 'molecProf').
#' @param gene `character` Name of gene directory within the `path` directory.
#'   Default is parameter name.
#' @param analysis `character` Name of gene directory within the `path`
#'   directory. Default is parameter name.
#' @param geneSynonymPattern A `character` string with a valid regex pattern,
#'   passed to `pattern` of `list.files` to load in all gene synonym files.
#'
#' @return None. Writes to disk.
#'
#' @import data.table
#' @md
#' @export
buildMolecProfTables <- function(path='procdata', outDir='latest',
    geneSynonymPattern='GeneSynonyms.*csv', ...,  molecProf='molecProf',
    gene='gene', analysis='analysis')
{
    ## TODO:: Refactor this into at least two functions

    # -- ensure the output directory exists
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # -- load the required table for each tSet

    # molecular profiles
    molecProfTables <- loadTSetTables(path, molecProf)

    # genes
    geneTables <- loadTSetTables(path, gene)

    # differential expression analysis
    analysisTables <- loadTSetTables(path, analysis)

    # -- load table dependencies
    ## TODO:: handle errors, maybe refactor into a helper function?

    compound <- fread(file.path(outDir, 'compound.csv'))
    dataset <- fread(file.path(outDir, 'dataset.csv'))
    sample <- fread(file.path(outDir, 'sample.csv'))
    dataset_sample <- fread(file.path(outDir, 'dataset_sample.csv'))
    cell <- fread(file.path(outDir, 'cell.csv'))

    # -- build gene tables
    ## TODO:: Change the join assigns to merge.data.table to prevent slow coercions to 
    #  the correct column type.
    gene <- rbindlist(geneTables)
    gene[, gene_id := gsub('_at', '', gene_id)]
    gene$dataset_id <- unlist(mapply(rep, x=names(geneTables),
        times=vapply(geneTables, nrow, numeric(1)), SIMPLIFY=FALSE))
    #rm(geneTables)
    gene_dataset <- unique(gene[, .(gene_id, dataset_id)])
    gene[, dataset_id := NULL]
    gene <- unique(gene)
    # Dropping duplicate gene ids caused by 1:many relationship with transcripts
    gene <- gene[!duplicated(gene_id)]
    gene[, id := seq_len(.N)]

    # gene annotation table
    gene_annotation <- gene[, .(id, symbol, entrez_gid, transcript_name,
        transcript_id)]
    gene_annotation[, full_name := NA]
    setnames(gene_annotation, c('id', 'transcript_id'),
        c('gene_id', 'ensembl_tid'))
    setkeyv(gene_annotation, 'gene_id')

    # map gene to gene_id
    setkeyv(gene_dataset, 'gene_id')
    setkeyv(gene, 'gene_id')
    gene_dataset[gene, gene_id := i.id]

    # map dataset to dataset_id
    setkeyv(gene_dataset, 'dataset_id')
    setkeyv(dataset, 'name')
    gene_dataset[dataset, dataset_id := i.id]

    # convert columns to proper type
    gene_dataset <- gene_dataset[, lapply(.SD, as.numeric)]
    setorderv(gene_dataset, c('gene_id', 'dataset_id'))

    # gene
    gene <- unique(gene[, .(id, gene_id)])
    setnames(gene, 'gene_id', 'name')

    # sanity check for gene tables
    gene_idInGene <- gene_dataset$gene_id %in% gene$id
    if (!all(gene_idInGene)) stop(.errorMsg(.context(),
        "Some gene ids aren't in gene? Something has gone terribly wrong!"))

    # -- build compound gene response

    molecProfTables <- lapply(molecProfTables, `[`,  # fix column types
        j=sample_id := as.character(sample_id))
    compound_gene_response <- unique(rbindlist(molecProfTables))
    compound_gene_response[, gene_id := gsub('_at', '', gene_id)]
    #rm(molecProfTables)

    # map ids from other tables
    setkeyv(sample, 'name')
    setkeyv(compound_gene_response, 'sample_id')
    compound_gene_response[sample, sample_id := i.id]

    setkeyv(gene, 'name')
    setkeyv(compound_gene_response, 'gene_id')
    compound_gene_response[gene, gene_id := i.id]

    compound_gene_response[, `:=`(sample_id=as.integer(sample_id),
                                  gene_id=as.integer(gene_id))]

    # -- build analysis table
    .addColumn <- function(x, colName, value) {  # reference symantics
        x[,  newCol := value]
        setnames(x, 'newCol', colName)
        x
    }
    analysisTables <- mapply(.addColumn, x=analysisTables,
        colName='dataset_id', value=names(analysisTables), SIMPLIFY=FALSE)
    analysis <- rbindlist(analysisTables)
    analysis[, gene_id := gsub('_at', '', gene_id)]
    #rm(analysisTables)

    # map ids from other tables
    setkeyv(analysis, 'gene_id')
    analysis[gene, gene_id := i.id]

    setkeyv(analysis, 'compound_id')
    setkeyv(compound, 'name')
    analysis[compound, compound_id := i.id]

    setkeyv(analysis, 'dataset_id')
    analysis[dataset, dataset_id := i.id]

    # map cell_id to analysis
    setkeyv(analysis, 'cell_id')
    setkeyv(cell, 'name')
    analysis[cell, cell_id := i.id]

    # map sample_id to analysis
    sample <- merge.data.table(sample, dataset_sample, by.x='id',
        by.y='sample_id')
    ## FIXME:: Stop coercing column types it's slow, just add new column and delete old one
    analysis[, `:=`(id=seq_len(.N), gene_id=as.integer(gene_id),  # fix column types for join
        compound_id=as.integer(compound_id), dataset_id=as.integer(dataset_id),
        cell_id=as.integer(cell_id))]

    setkeyv(analysis, c('compound_id', 'dose', 'time', 'cell_id', 'dataset_id'))
    setkeyv(sample, c('compound_id', 'dose', 'time', 'cell_id', 'dataset_id'))
    analysis <- merge.data.table(analysis, 
        sample[, c('id', 'compound_id', 'dose', 'time', 'cell_id', 'dataset_id')], 
        allow.cartesian=TRUE)
    setnames(analysis, c('id.x', 'id.y'), c('id', 'sample_id'))

    analysis[, `:=`(sample_id=as.integer(sample_id),
        dataset_id=NULL, dose=NULL, time=NULL, cell_id=NULL, compound_id=NULL)]

    # map analysis_id to compoundGeneResponse
    setkeyv(compound_gene_response, c('gene_id', 'sample_id'))
    setkeyv(analysis,  c('gene_id', 'sample_id'))
    compound_gene_response[analysis, analysis_id := i.id]
    compound_gene_response[is.na(analysis_id), analysis_id := 0]

    nullRow <- data.table(id=0)
    analysis[, `:=`(gene_id=NULL, sample_id=NULL)]
    analysis <- rbindlist(list(analysis, nullRow), fill=TRUE)

    # remove sample id and drop duplicate rows
    analysis <- unique(analysis)

    # sanity checks
    if (any(duplicated(analysis$id)))
        stop("Duplicate primary key in analysis table after mapping to samples!")
    if (!all(analysis$id %in% compound_gene_response$analysis_id))
        stop("Dropped some analysis ids when joining with compound_gene_response. It is appropriate to PANIC!")
    if (any(is.na(compound_gene_response$analysis_id)))
        stop("NAs found in analysis_id column of compound_gene_response! EGATS!")
    if (!all(sample$id %in% unique(compound_gene_response$sample_id)))
        stop("Some samples are missing from compound_gene_response! Things have gone terribly wrong!")


    # -- sort tables before writing to disk
    setorderv(analysis, 'id')
    setorderv(gene, 'id')
    setorderv(compound_gene_response, 'gene_id')

    # make the synonym table
    synonymFiles <- list.files('metadata', pattern=geneSynonymPattern,
        full.names=TRUE)
    synonyms <- lapply(synonymFiles, fread)
    synonym <- rbindlist(synonyms)[, id := NULL]
    synonym[Synonyms == '-', Synonyms := NA]
    rm(synonyms)
    setkeyv(synonym, 'EnsemblID')
    setkeyv(gene, 'name')
    gene_synonym <- merge.data.table(synonym, gene, by.x='EnsemblID',
        by.y='name', allow.cartesian=TRUE)[, .(id, Synonyms)]
    gene_synonym <- na.omit(gene_synonym)
    setnames(gene_synonym, c('id', 'Synonyms'), c('gene_id', 'synonym'))

    # -- sanity checks
    if (!all(gene_synonym$id %in% gene$id))
        stop(.errorMsg(.context(), 'Gene id in gene_synonym not in gene,
            please proceed to panic!'))

    # sanity check
    if (any(duplicated(gene$name)))
        stop(.errorMsg(.context(), 'Duplicated gene names!'))

    # -- write to disk
    for (table in c('gene', 'gene_dataset', 'compound_gene_response',
        'analysis', 'gene_annotation', 'gene_synonym'))
    {
        fwrite(get(table), file.path(outDir, paste0(table, '.csv')))
    }
}



if (sys.nframe() == 0) {
    library(data.table)
    source('R/utilities.R')
    path='procdata'
    outDir='latest'
    geneSynonymPattern='GeneSynonyms.*csv'
    molecProf='molecProf'
    gene='gene'
    analysis='analysis'

}