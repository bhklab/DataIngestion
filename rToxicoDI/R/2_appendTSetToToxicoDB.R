##############################################
#### ADD NEW TSETS TO EXISTING DATABASE TABLES 
##############################################

## TODO::
# - The current implementation is fragile and dependent on the order in which
# TSets are appended together; fix this
# - Rewrite this method to use disk.frames to allow larger than memory joins


library(dplyr)
library(data.table)
library(ToxicoGx)

# Select which tables to regenerate 
tables <- c("compounds", "compounds_datasets", "datasets", "genes_datasets", "genes",
            "gene_annotations", "compound_annotations", "datasets_samples", "samples",
            "species", "tissues", "cells", "viabilities",
            "analysis", "pathways", "pathways_genes", "pathways_datasets",
            "compound_gene_response")

appendTSetToToxicoDB <- function(tables, lab_tables, lab_tSet, lab_out) {

    #### SET MULTITHREADING ####
    # #DTthreadsOld <- getDTthreads()
    # #setDTthreads(14)
    # #on.exit(setDTthreads(DTthreadsOld)) # Revert to system defaults after execution
    #
    #### READ IN TSET TABLES ####
    for (dt in tables) {
        assign(paste0(dt, 1), fread(paste0('results/', dt, lab_tSet, ".csv")))
    }

    #### READ IN CURRENT DB TABLES ####
    for (dt in tables) {
        assign(paste(dt), fread(paste0('results/', dt, lab_tables,".csv")))
    }

    #### MERGING TABLES AND SAVING TO DISK ####

    # Analysis table; no overlap so easy
    analysis1 <- na.omit(analysis1)
    analysis2 <- rbindlist(list(na.omit(analysis), analysis1[, id := (id + nrow(analysis))]))
    if (!(0 %in% analysis2$id)) { # Add zero row to analysis (to allow unmatched compound_gene_response rows)
        analysis2 <- rbindlist(list(data.table(id = 0, fold_change = NA, log_odds = NA, p_value = NA, fdr = NA, avg_expr = NA), analysis2))
    }

    if (
        any(analysis2[id %in% na.omit(analysis1)$id , ] != na.omit(analysis1)) ||
            any(analysis2[id %in% na.omit(analysis)$id, ] != na.omit(analysis))
    ) { stop('Error in analysis table mapping!') } else {
        fwrite(analysis2, file = paste0('results/', 'analysis', lab_out, '.csv'), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
        #rm(analysis2)
        gc()
    }

    ## Species table
    species3 <- species1[(which(!(species1$name %in% species$name))),]
    if (nrow(species3) > 0) {
        species2 <- rbindlist(list(species, species3[, dataset_id := (dataset_id + nrow(species))]))
    } else {
        species2 <- species
    }
    fwrite(species2, file = paste0('results/', "species", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")

    ## Tissues
    tissues3 <- tissues1[!(tissues1$name %in% tissues$name), 'id' := (id + nrow(species))]
    if (nrow(tissues3) > 0) {
        tissues2 <- unique(rbindlist(list(tissues, tissues3)))
    } else {
        tissues2 <- tissues
    }
    fwrite(tissues2, file = paste0('results/', "tissues", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
    #rm(tissues, tissues1, tissues2, tissues3)
    gc()

    ## Cells
    cells2 <- rbindlist(list(
        cells,
        cells1
    ))
    fwrite(unique(cells2), file = paste0('results/', "cells", lab_out, ".csv"),
           row.names = FALSE, na = "", eol = "\r\n", sep = ",")
    #rm(cells, cells1, cells2)
    gc()

    ## compounds
    compounds3 <- copy(compounds1)
    compounds1 <- compounds1[!(name %in% compounds$name), ][, id := (.I + nrow(compounds))]
    compounds2 <- rbindlist(list(
        compounds,
        compounds1
    ))



    if (
        any(compounds2[id %in% compounds$id, ] != compounds) ||
            any(compounds2[id %in% compounds1$id, ] != compounds1)
    ) { stop('Issue with compounds table mappings!') } else {
        fwrite(compounds2, file = paste0('results/', "compounds", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
        #rm(compounds2, compounds1, compounds)
        gc()
    }


    ## compound Annotations
    compound_annotations1 <- compound_annotations1[which(!(compounds3$name %in% compounds$name))][, compound_id := (.I + nrow(compound_annotations))]
    compound_annotations2 <- rbindlist(list(
        compound_annotations,
        compound_annotations1
    ))

    annots <- fread(paste0('metadata/', 'Drug_annotations_V2.1.csv'))[!duplicated(unique.drugid),]
    da_joined <- compound_annotations2[inchikey != "" & !duplicated(inchikey)][annots[inchikey != "" & !duplicated(inchikey),], on = 'inchikey', nomatch=NULL]
    if (any(compounds2[da_joined$compound_id,]$name != da_joined$unique.drugid)) {
        warning("Issue with compound annotations inchikey mapping!")
    }

    if (
        any(compound_annotations2[compound_annotations$compound_id, ] != compound_annotations, na.rm = TRUE) ||
            suppressWarnings(any(compound_annotations2[compound_annotations1$compound_id, -'symbol'] != compound_annotations1[, -'symbol'], na.rm = TRUE))
    ) { stop('Issue with compound_annotations table mappings!') } else {
        fwrite(compound_annotations2, file = paste0('results/', "compound_annotations", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
        #rm(compound_annotations, compound_annotations1, compound_annotations2)
        gc()
    }

    ## Datasets
    datasets1 <- datasets1[which(!(datasets1$name %in% datasets$name)),][, id := (.I + nrow(datasets))]
    datasets2 <- rbindlist(list(
        datasets,
        datasets1
    ))
    fwrite(datasets2, file = paste0('results/', "datasets", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
    #rm(datasets, datasets1, datasets2)
    gc()

    ## Genes
    genes3 <- copy(genes1)
    genes1 <- genes1[which(!(genes1$name %in% genes$name))][, id := (.I + nrow(genes))]
    genes2 <- rbindlist(list(
        genes,
        genes1
    ))

    if (
        any(genes2[id %in% genes$id, ] != genes) ||
            ifelse(nrow(genes1) != 0, any(genes1[id %in% genes1$id, ] != genes1), FALSE)
    ) { stop('Issue with genes table mappings!') } else {
        fwrite(genes2, file = paste0('results/', "genes", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
        #rm(genes, genes1, genes2)
        gc()
    }

    ## Gene Annotations
    gene_annotations1 <- gene_annotations1[which(!(ensembl_tid %in% gene_annotations$ensembl_tid))][, gene_id := (.I + nrow(gene_annotations))]
    gene_annotations2 <- rbindlist(list(
        gene_annotations,
        gene_annotations1
    ))

    if (
        any(gene_annotations2[gene_id %in% gene_annotations$gene_id, ] != gene_annotations, na.rm=TRUE) ||
            ifelse(nrow(gene_annotations1) != 0, any(gene_annotations2[gene_id %in% gene_annotations1$gene_id, ] != gene_annotations1, na.rm=TRUE), FALSE)
    ) { stop('Issue with gene annotation mapping') } else {
        fwrite(gene_annotations2, file = paste0('results/', "gene_annotations", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
        #rm(gene_annotations, gene_annotations1, gene_annotations2)
        gc()
    }

    ## Samples

    # Update compound_ids of samples table
    samples1_length <- nrow(samples1)
    compounds4 <- copy(compounds3)[, compound_id := id][, id := NULL]
    samples1 <- samples1[compounds4, on='compound_id',][, compound_id := NULL]
    compounds4 <- copy(compounds2)[, compound_id := id][, id := NULL]
    samples1 <- na.omit(samples1[compounds4, on='i.name==name'][order(id), .SD, .SDcols = colnames(samples)])
    rm(compounds4); gc()

    if (nrow(samples1) != samples1_length) {
        stop("Length of samples1 changed after joins to update compound_id!")
    }

    # Increment sample ids based on previous tSet
    samples1 <- samples1[, id := id + nrow(samples)]
    samples2 <- rbindlist(list(
        samples,
        samples1
    ))

    ##FIXME:: Why are the tables not ordered by ID?
    if (
        any(samples2[samples$id, ][order(id)] != samples[order(id)]) ||
            any(samples2[samples1$id, ][order(id)] != samples1[order(id)])
    ) { stop('Issue with samples table mappings!') } else {
        fwrite(samples2, file = paste0('results/', "samples", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
    }

    ## Species
    ## TODO:: Finish implementing cell id replacement based on tSet2 values
    if (!('name' %in% species1)) species1[, name := 'Human']
    species1 <- species1[which(!(name %in% species$name))][, dataset_id := (dataset_id + nrow(species))]
    species2 <- rbindlist(list(
        species,
        species1
    ))

    ## Viabilities
    if (lab_tSet %in% c('_drugMatrix_rat', '_EMEXP2458')) {
        viabilities1 <- viabilities
        viabilities2 <- viabilities
    } else {
        viabilities1 <- viabilities1[, sample_id := sample_id + nrow(viabilities)]
        viabilities2 <- rbindlist(list(
            viabilities,
            viabilities1
        ))
    }

    if (!(lab_tSet %in% c('_drugMatrix_rat', '_EMEXP2458'))) {
        if (
            any(viabilities2[sample_id %in% viabilities$sample_id, ] != viabilities, na.rm = TRUE) ||
            any(viabilities2[sample_id %in% viabilities1$sample_id, ] != viabilities1, na.rm = TRUE)
        ) {
            stop('Issue with viabilities table mappings!')
        }
    } else {
        fwrite(viabilities2, file = paste0('results/', "viabilities", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
        rm(viabilities, viabilities1, viabilities2)
        gc()
    }


    ## compound Gene Response



    # Update analyis_id, gene_id and sample_id for all values in compound_gene_response1
    genes4 <- copy(genes2)[, gene_id := id][, -'id']
    compound_gene_response1 <- compound_gene_response1[genes3, on=c(gene_id='id')][order(id), -'gene_id']
    compound_gene_response1 <- compound_gene_response1[genes4, on='name'][order(id), .SD, .SDcols=colnames(compound_gene_response)]
    compound_gene_response1 <- compound_gene_response1[analysis_id != 0, analysis_id := (analysis_id + nrow(analysis))]
    compound_gene_response1 <- na.omit(compound_gene_response1[, sample_id := (sample_id + nrow(samples))])

    if (
        any(analysis2[id %in% unique(compound_gene_response1[analysis_id > 0,]$analysis_id)] != analysis1)
    ) { stop('Issue with remapping analysis_ids in compound_gene_response')
    } else if (
        any(!(as.numeric(samples2[id %in% unique(compound_gene_response1$sample_id), ]$name) %in% as.numeric(samples1$name)))
    ) { stop('Issue with remapping sample_ids in compound_gene_response')
    } else if (
        any(!(ifelse(length(genes1$name > 0),
                     genes2[id %in% unique(compound_gene_response1$gene_id), ]$name %in% genes1$name,
                     TRUE)))
    ) { stop('Issue with remapping gene_ids in compound_gene_response') }

    # Set the new compound gene response
    compound_gene_response2 <- rbindlist(list(
        compound_gene_response,
        compound_gene_response1[, id := .I + nrow(compound_gene_response)]
    ))

    if (
        any(compound_gene_response2[id %in% compound_gene_response$id, ] != compound_gene_response) ||
            any(compound_gene_response2[id %in% compound_gene_response1$id, ] != compound_gene_response1)
    ) { stop('Issue with mapping between new and old compound_gene_response tables!')
    } else if (!all(analysis1$id %in% compound_gene_response2[analysis_id > 0]$analysis_id)) {
        stop('Issue with mapping between anaylsis1 and compound_gene_response2!')
    } else {
        fwrite(compound_gene_response2,
               file = paste0('results/', "compound_gene_response", lab_out, ".csv"),
               row.names = FALSE, na = "", eol = "\r\n", sep = ",")
        rm(compound_gene_response2, compound_gene_response1, compound_gene_response)
        gc()
    }

    ## Pathways
    ## TODO:: Genearlize this for including more pathways
    pathways1 <- pathways1[which(!(pathways1$gene_id %in% pathways$gene_id)),][, id := .I + nrow(pathways)]
    pathways2 <- rbindlist(list(
        pathways,
        pathways1
    ))
    fwrite(pathways2, file = paste0('results/', "pathways", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
    rm(pathways, pathways1)

    #### JOIN TABLES ####

    ## PATHWAYS GENES

    # Replace modified gene_ids
    new_genes <- which(!(genes1$name %in% genes$name))
    pathways_genes1 <- pathways_genes1[gene_id %in% new_genes,][, gene_id := (gene_id + nrow(genes))]
    pathways_genes2 <- rbindlist(list(
        pathways_genes,
        pathways_genes1
    ))
    fwrite(pathways_genes2, file = paste0('results/', "pathways_genes", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
    #rm(pathways_genes, pathways_genes1, pathways_genes2)
    gc()

    ## PATHWAYS DATASETS
    pathways_datasets1 <- pathways_datasets1[, dataset_id := which(datasets2$name %in% datasets1$name)]
    pathways_datasets2 <- rbindlist(list(
        pathways_datasets,
        pathways_datasets1
    ))
    fwrite(pathways_datasets2, file = paste0('results/', "pathways_datasets", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
    #rm(pathways_datasets, pathways_datasets1, pathways_datasets2)
    gc()

    ## DATASETS SAMPLES
    datasets_samples1 <- datasets_samples1[, dataset_id := which(datasets2$name %in% datasets1$name)][, sample_id := (sample_id + nrow(samples))]
    datasets_samples2 <- rbindlist(list(
        datasets_samples,
        datasets_samples1
    ))

    if (
        any(samples2[
                id %in% datasets_samples2[dataset_id == which(datasets2$name == datasets1$name)]$sample_id
            ] != samples1) ||
            any(samples2[
                    id %in% datasets_samples2[dataset_id %in% which(datasets2$name %in% datasets$name)]$sample_id
                ] != samples)
    ) {
        stop('Mismatch between data_sets samples and samples tables!')
    } else {
        fwrite(datasets_samples2, file = paste0('results/', "datasets_samples", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
        gc()
    }

    ## GENES DATASETS
    genes_datasets1 <- genes_datasets1[, dataset_id := which(datasets2$name %in% datasets1$name)][, gene_id := which(genes2$name %in% genes3$name)]
    genes_datasets2 <- rbindlist(list(
        genes_datasets,
        genes_datasets1
    ))

    if (
        !all(genes2[
                 id %in% genes_datasets2[dataset_id == which(datasets2$name == datasets1$name)]$gene_id
             ]$name %in% genes3$name) ||
            !all(genes2[
                     id %in% genes_datasets2[dataset_id %in% which(datasets2$name %in% datasets$name)]$gene_id
                 ]$name %in% genes$name)
    ) {
        stop('Mismatch between genes_datasets and genes tables!')
    } else {
        fwrite(genes_datasets2, file = paste0('results/', "genes_datasets", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
        gc()
    }

    ## TODO:: Determine if this test works?
    if (!all(genes1$name %in% genes2[genes_datasets2[dataset_id == which(datasets2$name %in% datasets1$name)]$gene_id]$name)) {
        stop(paste0('Mapping issue with genes_datasets - mismatches for: ',
                    genes2[which(!(genes2[genes_datasets2[dataset_id == which(datasets2$name %in% datasets1$name)]$gene_id]$name %in% genes1$name))]))
    }

    ## compoundS DATASETS
    compounds_datasets1 <- compounds_datasets1[, dataset_id := which(datasets2$name %in% datasets1$name)][, compound_id := merge(compounds2, compounds3, by="name")[order(id.y),]$id.x]

    ## TODO:: Remove this column remaning and include in tSetToDatabaseCSV
    colnames(compounds_datasets) <- c('compound_id', 'dataset_id', 'compound_uid')
    compounds_datasets2 <- rbindlist(list(
        compounds_datasets,
        compounds_datasets1
    ))

    if (
        !all(compounds2[
                 id %in% compounds_datasets2[dataset_id == which(datasets2$name == datasets1$name)]$compound_id
             ]$name %in% compounds3$name) ||
            !all(compounds2[
                     id %in% compounds_datasets2[dataset_id %in% which(datasets2$name %in% datasets$name)]$compound_id
                 ]$name %in% compounds$name)
    ) {
        stop('Mismatch between compounds_datasets table and compounds tables!')
    } else {
        fwrite(datasets_samples2, file = paste0('results/', "datasets_samples", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
        gc()
    }

    ## TODO:: Determine if this test works?
    if (!all(compounds1$name %in% compounds2[compounds_datasets2[dataset_id == which(datasets2$name %in% datasets1$name)]$compound_id]$name)) {
        stop(paste0('Mapping issue with compounds_datasets - mismatches for: ',
                    compounds2[which(!(compounds2[gene_datasets2[dataset_id == which(dataset2$name %in% datasets1$name)]$compound_id]$name %in% compounds1$name))]))
    }

    fwrite(compounds_datasets2, file = paste0('results/', "compounds_datasets", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
    rm(compounds_datasets, compounds_datasets1, compounds_datasets2)
    gc()
}

clean_env <- function() {
    rm(list = ls())
    gc()
}
