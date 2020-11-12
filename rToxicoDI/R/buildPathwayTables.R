#' Build all tables related to samples for ToxicoDB and write to disk
#'
#' @param path `character` path to the directory where the ToxicoSet data was
#'  extracted
#' @param outDir `character` path to save table csvs to
#' @param ... `pairlist` force further arguments to be named
#'
#' @import data.table
#' @export
buildPathwayTables <- function(path='procdata', outDir='latest', ...)
{
    # -- ensure the output directory exists
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # -- load the compound table for each tSet
    files <- list.files(file.path(path, 'pathway'), pattern='csv',
        full.names=TRUE)

    # -- get the pathway table and remove it from the file path vectorzx
    whichPathwayPath <- grepl('pathway.csv', files)
    pathwayDT <- fread(files[whichPathwayPath])
    files <- files[!whichPathwayPath]

    # -- read in the stats data
    pathwayStats <- lapply(files, fread)
    names(pathwayStats) <- gsub( '.*/|.csv', '', files)

    # -- load dependent tables
    compound <- fread(file.path(outDir, 'compound.csv'))
    dataset <- fread(file.path(outDir, 'dataset.csv'))
    cell <- fread(file.path(outDir, 'cell.csv'))
    gene <- fread(file.path(outDir, 'gene.csv'))
    ## FIXME:: Remove this
    gene$name <- gsub('_at', '', gene$name)

    # drop non-existent genes
    pathwayDT <- pathwayDT[gene_id %in% gene$name]

    # -- build pathway table
    pathway <- unique(pathwayDT[, .(pathway_id)])
    pathway[, id := seq_len(.N)]
    setnames(pathway, 'pathway_id', 'name')

    # -- build the pathway gene table
    pathway_gene <- pathwayDT[, .(pathway_id, gene_id)]
    # get pathway_id
    setkeyv(pathway, 'name')
    setkeyv(pathway_gene, 'pathway_id')
    pathway_gene[pathway, pathway_id := i.id]
    # get gene_id
    setkeyv(gene, 'name')
    setkeyv(pathway_gene, 'gene_id')
    pathway_gene[gene, gene_id := i.id]
    pathway_gene[, `:=`(pathway_id=as.integer(pathway_id),
        gene_id=as.integer(gene_id))]

    # -- build pathway_stats table
    ## FIXME:: Use standardized file names to prevent having manually specify which TSet
    ## TODO:: Refactor this mess
    tghDT <- rbindlist(pathwayStats[grepl('TGH', names(pathwayStats))], fill=TRUE)
    tghDT[, `:=`(dataset_id='TGGATES_humanldh', cell_id='Hepatocyte')]
    tgrDT <- rbindlist(pathwayStats[grepl('TGR', names(pathwayStats))], fill=TRUE)
    tgrDT[, `:=`(dataset_id='TGGATES_ratldh', cell_id='Hepatocyte')]
    dmDT <- rbindlist(pathwayStats[grepl('DM', names(pathwayStats))], fill=TRUE)
    dmDT[, `:=`(dataset_id='drugMatrix_rat', cell_id='Hepatocyte')]
    hepargDT <- rbindlist(pathwayStats[grepl('RG', names(pathwayStats))], fill=TRUE)
    hepargDT[, `:=`(dataset_id='EMEXP2458', cell_id='HepaRG')]
    hepag2DT <- rbindlist(pathwayStats[grepl('G2', names(pathwayStats))], fill=TRUE)
    hepag2DT[, `:=`(dataset_id='EMEXP2458', cell_id='Hep-G2')]

    pathway_stats <- rbindlist(list(tghDT, tgrDT, dmDT, hepargDT, hepag2DT), fill=TRUE)
    newNames <- c('pathway_id', 'compound_id', 'genes_total', 'stat_dis', 'genes_up',
            'genes_down', 'p_value')
    setnames(pathway_stats,
        old=c('pathway', 'drugname', 'Genes..tot.', 'Stat..dist.dir.',
            'Genes..up.', 'Genes..down.', 'p'),
        new=newNames)
    pathway_stats <- pathway_stats[, .SD,
        .SDcols=c(newNames, 'dataset_id', 'cell_id', 'fdr')]
    # fix drugid
    pathway_stats[, compound_id := gsub('.fold_change', '',
        pathway_stats$compound_id)]

    # -- fix pathway ids
    pathway_stats[grepl('^KEGG', pathway_id),
        pathway_id := stringr::str_extract(pathway_id, '^KEGG:[^-]*')]
    pathway_stats[grepl('^REACT', pathway_id),
        pathway_id := stringr::str_extract(pathway_id, '^REACT:[^-]*-[^-]*-[^-]*')]

    # remove ids we couldn't map
    pathway_stats <- pathway_stats[pathway_id %in% pathway$name]


    # -- map ids to pathway_stats
    ## TODO:: Refactor into a helper function

    # -- rename old drugs
    drugRemappings <- fread(file.path('metadata', 'old_newDrugmapping.csv'))
    setkey(drugRemappings, 'Old name')
    setkey(pathway_stats, 'compound_id')
    pathway_stats[drugRemappings, compound_id := `New name`]
    if (any(drugRemappings$`Old name` %in% unique(pathway_stats$compound_id)))
        stop(.errorMsg(.context(), 'An old drug name is present in the
            pathway_stats table, something has gone wrong with remapping
            old drug names to new ones.'))

    # pathway_id
    setkeyv(pathway, 'name')
    setkeyv(pathway_stats, 'pathway_id')
    pathway_stats <- merge.data.table(pathway_stats, pathway,
        all.x=TRUE, by.x='pathway_id', by.y='name')
    pathway_stats[, pathway_id := NULL]
    setnames(pathway_stats, 'id', 'pathway_id')

    # compound_id
    setkeyv(compound, 'name')
    setkeyv(pathway_stats, 'compound_id')
    pathway_stats <- merge.data.table(pathway_stats, compound,
        by.x='compound_id', by.y='name', all.x=TRUE)
    pathway_stats[, compound_id := NULL]
    setnames(pathway_stats, 'id', 'compound_id')

    # cell_id
    setkeyv(cell, 'name')
    setkeyv(pathway_stats, 'cell_id')
    pathway_stats <- merge.data.table(pathway_stats, cell,
        by.x='cell_id', by.y='name', all.x=TRUE)
    pathway_stats[, `:=`(cell_id=NULL, tissue_id=NULL)]
    setnames(pathway_stats, 'id', 'cell_id')

    # dataset_id
    setkeyv(dataset, 'name')
    setkeyv(pathway_stats, 'dataset_id')
    pathway_stats[dataset, dataset_id := i.id][, dataset_id :=
        as.integer(dataset_id)]


    # -- build pathway_compound
    pathway_stats <- unique(pathway_stats)  # ensure there are no duplicate rows
    pathway_stats[, id := seq_len(.N)]
    pathway_compound <- unique(pathway_stats[, .(pathway_id, compound_id,
        cell_id, dataset_id)])

    # sanity check
    if (nrow(pathway_compound) != nrow(pathway_stats))
        stop(.errorMsg(.context(), 'The unique combination of ',
            'pathway_id, compound_id, dataset_id and cell_id has failed to ',
            'uniquely identify rows in pathway_stats!'))

    # -- build pathway_genes
    pathway_gene <- unique(pathwayDT[, .(pathway_id, gene_id)])

    # map gene id to table
    setkeyv(gene, 'name')
    setkeyv(pathway_gene, 'gene_id')
    pathway_gene[gene, gene_id := i.id]
    pathway_gene[, gene_id := as.integer(gene_id)]

    # map pathway id to table
    setkeyv(pathway_gene, 'pathway_id')
    setkeyv(pathway, 'name')
    pathway_gene[pathway, pathway_id := i.id]
    pathway_gene[, pathway_id := as.integer(pathway_id)]

    # map dataset id to table
    pathway_dataset <- unique(pathway_stats[, .(id, dataset_id)])
    setnames(pathway_dataset, 'id', 'dataset_id')

    # subset pathway_stats
    pathway_stats[, ontology := NA]
    pathway_stats <- pathway_stats[, .(pathway_id, ontology, genes_total,
        stat_dis, genes_up, genes_down, p_value, fdr)]

    for (table in c('pathway', 'pathway_stats', 'pathway_compound',
        'pathway_dataset'))
    {
        fwrite(get(table), file=file.path(outDir, paste0(table, '.csv')))
    }
}

if (sys.nframe() == 0) {
    library(data.table)
    outDir <- 'latest'
    path <- 'procdata'
}