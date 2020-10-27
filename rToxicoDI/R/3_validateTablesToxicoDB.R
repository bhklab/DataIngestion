######################################################
#### VALIDATE MAPPING BETWEEN TABLES FOR TOXICODB ####
######################################################

# 2020/02/03
# by Christopher Eeles

# Dependencies ------------------------------------------------------------
library(data.table)
library(dplyr)
library(ToxicoGx)

# Load Data ---------------------------------------------------------------

## List of tables to load into memory
tables <- c("compounds", "compounds_datasets", "datasets", "genes_datasets", "genes", 
            "gene_annotations", "compound_annotations", "datasets_samples", "samples", 
            "species", "tissues", "cells", "viabilities",
            "analysis", "pathways", "pathways_genes", "pathways_datasets",
            "compound_gene_response")

validateTablesToxicoDB <- function(label) {
  ## Read in tables as data.table objects
  for (dt in tables) {
    assign(dt, fread(paste0('../results/', dt, label, ".csv")))
  }
  
  ## Read in tSets from .rds 
  TGGATES_humanldh <- readRDS('../tSets/TGGATES_humanldh.rds')
  TGGATES_ratldh <- readRDS('../tSets/TGGATES_ratldh.rds')
  drugMatrix_rat <- readRDS('../tSets/drugMatrix_rat.rds')
  
  # Test Mapping ------------------------------------------------------------
  
  ## compounds
  
  ### TGGATEShuman
  dn <- drugNames(TGGATES_humanldh)
  compound_ids <- compounds_datasets[dataset_id == 1, ]$compound_id
  if (!all(compounds[id %in% compound_ids,]$name %in% dn)) {
    stop("TGGATEShuman compounds mismatch")
  }
  
  ### TGGATESrat
  dn <- drugNames(TGGATES_ratldh)
  compound_ids <- compounds_datasets[dataset_id == 2, ]$compound_id
  if (!all(compounds[id %in% compound_ids,]$name %in% dn)) {
    stop("TGGATESrat compounds mismatch")
  }
  
  ### drugMatrix_rat
  dn <- drugNames(drugMatrix_rat)
  compound_ids <- compounds_datasets[dataset_id == 3, ]$compound_id
  if (!all(compounds[id %in% compound_ids,]$name %in% dn)) {
    stop("drugMatrix_rat compounds mismatch")
  }
  
  ## compounds Datasets
  
  uid <- unique(phenoInfo(TGGATES_humanldh, 'rna')$dataset_compoundid)
  unique_ids <- compounds_datasets[dataset_id == 1, ]$unique_id
  if (!all(uid == unique_ids)) {
    stop("TGGATEShuman unique)id mismatch")
  }
  
  ### TGGATESrat
  uid <- unique(phenoInfo(TGGATES_ratldh, 'rna')$dataset_compoundid)
  unique_ids <- compounds_datasets[dataset_id == 2, ]$unique_id
  if (!all(uid == unique_ids)) {
    stop("TGGATESrat unique_id mismatch")
  }
  
  ### drugMatrix_rat
  uid <- unique(phenoInfo(drugMatrix_rat, 'rna')$dataset_compoundid)
  unique_ids <- compounds_datasets[dataset_id == 3, ]$unique_id
  if (!all(uid == unique_ids)) {
    stop("drugMatrix_rat unique_id mismatch")
  }
    
  ## Samples
  
  ### TGGATEShuman
  sn <- phenoInfo(TGGATES_humanldh, 'rna')$samplename
  sample_ids <- datasets_samples[dataset_id == 1, ]$sample_id
  if (!all(as.character(samples[id %in% sample_ids]$name) %in% sn)) {
    stop('TGGATEShuman samples mismatch')
  }
  if (!all(samples[id %in% sample_ids]$compound_id) %in% compounds_datasets[dataset_id == 1, ]$compound_id) {
    stop("TGGATEShuman sample compound_ids mismatch")
  }
  if (!all(unique(compound_gene_response[sample_id %in% sample_ids]$gene_id) %in% genes_datasets[dataset_id == 1, ]$gene_id)) {
    stop("TGGATEShuamn compound_gene_response gene_ids mismatch")
  }
  
  ### TGGATESrat
  sn <- phenoInfo(TGGATES_ratldh, 'rna')$samplename
  sample_ids <- datasets_samples[dataset_id == 2, ]$sample_id
  if (!all(as.character(samples[id %in% sample_ids]$name) %in% sn)) {
    stop('TGGATESrat samples mismatch')
  }
  if (!all(samples[id %in% sample_ids]$compound_id %in% compounds_datasets[dataset_id == 2, ]$compound_id)) {
    stop("TGGATESrat sample compound_ids mismatch")
  }
  if (!all(unique(compound_gene_response[sample_id %in% sample_ids]$gene_id) %in% genes_datasets[dataset_id == 2, ]$gene_id)) {
    stop("TGGATESrat compound_gene_response gene_ids mismatch")
  }
  
  ### drugMatrix_rat
  sn <- phenoInfo(drugMatrix_rat, 'rna')$samplename
  sample_ids <- datasets_samples[dataset_id == 3, ]$sample_id
  if (!all(as.character(samples[id %in% sample_ids]$name) %in% sn)) {
    stop('drugMatrix_rat samples mismatch')
  }
  if (!all(samples[id %in% sample_ids]$compound_id) %in% compounds_datasets[dataset_id == 3, ]$compound_id) {
    stop("drugMatrix_rat sample compound_ids mismatch")
  }
  if (!all(unique(compound_gene_response[sample_id %in% sample_ids]$gene_id) %in% genes_datasets[dataset_id == 3, ]$gene_id)) {
    stop("drugMatrix_rat compound_gene_response gene_ids mismatch")
  }
  
  ## Genes
  gn <- fNames(TGGATES_humanldh, 'rna')
  gene_ids <- genes_datasets[dataset_id == 1, ]$name
  if (!all(genes[id %in% gene_ids,]$name %in% gn)) {
    stop("TGGATEShuman genes mismatch")
  }
  
  ### TGGATESrat
  gn <- fNames(TGGATES_ratldh, 'rna')
  gene_ids <- genes_datasets[dataset_id == 2, ]$name
  if (!all(genes[id %in% gene_ids,]$name %in% gn)) {
    stop("TGGATESrat genes mismatch")
  }
  
  ### drugMatrix_rat
  gn <- fNames(drugMatrix_rat, 'rna')
  gene_ids <- genes_datasets[dataset_id == 3, ]$name
  if (!all(genes[id %in% gene_ids,]$name %in% gn)) {
    stop("drugMatrix_rat genes mismatch")
  }
  
  ## Gene Annotations
  annots <- fread('../metadata/Drug_annotations_V2.1.csv')[!duplicated(unique.drugid),]
  da_joined <- compound_annotations[inchikey != "" & !duplicated(inchikey)][annots[inchikey != "" & !duplicated(inchikey),], on = 'inchikey', nomatch=NULL]
  if (any(compounds[da_joined$compound_id,]$name != da_joined$unique.compoundid)) {
    stop("Issue with compound annotations inchikey mapping!")
  }
  
  datasets_tissues <- data.table(expand.grid(datasets$id, tissues$id))
  colnames(datasets_tissues) <- c('dataset_id', 'tissue_id')
  fwrite(datasets_tissues, file = paste0('../results/', "datasets_tissues", label, ".csv"))
}

##FIXME:: What did de used to be? Replace with new results to check the validity of this test
## de is analysiscompoundsGenes; when was this test supposed to be run?
.analysisMappingTest <- function(tSet_label, final_tables_label) {
  
  # Load tSet data to compare against final tables
  analysis <- fread(file.path('..', 'results', paste0('analysiscompoundsGenes', final_tables_label)))
  compounds2 <- 
  genes2 <- 
  
  # Load final table data
  
  # Analysis Mapping Test
  # Take random sample of analysis ids, and evaluate that the mappings are correct
  for (i in sample(seq_len(nrow(analysis)), 150)) {
    compound_ids <- analysis[i,]$compound_id 
    gene_ids <- analysis[i,]$gene_id
    gene_ids2 <- compound_gene_response[analysis_id == i]$gene_id
    compound_ids2 <- samples[compound_gene_response[analysis_id == i]$sample_id,]$compound_id
    if (!all(compound_ids %in% compound_ids2)) {
      stop(paste0('Mapping issue for analysis ', i, ' mismatches at compounds', paste(which(!(compound_ids %in% compound_ids2)), collapse = ' ')))
    } else if (!all(gene_ids %in% gene_ids2)) {
      stop(paste0('Mapping issue for analysis ', i, ' mismatches at genes', paste(which(!(gene_ids %in% genes_ids2)), collapse = ' ')))
    }
    print(i)
  }
  

}

