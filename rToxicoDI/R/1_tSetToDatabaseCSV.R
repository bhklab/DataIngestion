####
# Data for ToxicoDB
####

library(ToxicoGx)
library(Biobase)
library(SummarizedExperiment)
library(S4Vectors)
library(tidyr)
library(dplyr)
library(data.table)
library(parallel)
library(BiocParallel)

#### PICK A tSET TO USE FOR THE DATABASE ####
tSetToDBtables <- function(tSet, lab_out) {

  #### POPULATE DATAFRAMES ####

  #### BASIC TABLES ####

  ## compound
  compounds <- data.table("id" = seq_along(drugNames(tSet)),
                      "name" = drugNames(tSet)
                      )
  compounds1 <- copy(compounds)
  
  annot <- fread("../metadata/Drug_annotations_V2.1.csv")
  colnames(annot)[2] <- "name"

  annot <- annot[compounds, on = "name"]
  annot <- annot[!duplicated(name), ]
  compounds1[an, on="name", ctd := i.ctd]
  
  
  whichAn <- which(an$name %in% compounds$name)
  
  ## compound ANNOTATION
  compound_annotations <- data.table(
    "compound_id" = compounds$id,
    "pubchem" = annot$cid,
    "chembl" = rep(NA, length(drugNames(tSet))),
    "compoundbank" = rep(NA, length(drugNames(tSet))),
    "targets" = rep(NA, length(drugNames(tSet))),
    'carcinogenicity' = annot$Carcinogenicity,
    "class_in_vivo" = annot$`Classif. in vivo`,
    "class_in_vitro" = annot$`Classif. in vitro`,
    "class_name" = rep(NA, length(drugNames(tSet))),
    "smiles" = annot$smiles,
    "inchikey" = annot$inchikey,
    "ctd" = compounds1$ctd
  )
  compound_annotations[carcinogenicity == "", carcinogenicity := NA,
                   by = carcinogenicity][class_in_vivo == "",
                                         class_in_vivo := NA,
                                         by = carcinogenicity]

  rm(compounds1)
  gc()
  
  ## DATASET
  datasets <- data.table("id" = c(1), "name" = name(tSet))

  ## GENE
  genes <- data.table("id" = seq_len(nrow(featureInfo(tSet, "rna"))),
                      "name" = featureInfo(tSet, "rna")$gene_id)

  ## GENE_ANNOTATION
  gene_annotations <- data.table(
    cbind(seq_len(nrow(featureInfo(tSet, "rna"))), # id
    as.data.frame(featureInfo(tSet, "rna")[, c("Symbol", "EntrezGene.ID", "transcript_name", "transcript_id")])) # ensembl, entrez
  )
  colnames(gene_annotations) <- c("gene_id", "Symbol", "entrez_gid", "transcript_name", "ensembl_tid")

  ## SPECIES
  species <- data.table("dataset_id" = datasets$id,
                        "name" = unique(phenoInfo(tSet, "rna")$species))
  ## TISSUES
  tissues <- data.table(
    "id" = seq_along(unique(phenoInfo(tSet, "rna")$organ_id)),
    "name" = unique(phenoInfo(tSet, "rna")$organ_id),
    "code" = rep(NA, length(unique(phenoInfo(tSet, "rna")$organ_id)))
  )

  ## TODO:: Consider adding Taxonomy table as three way join of cells, species and
  ## CELLS
  cells <- data.table(distinct(as.data.frame(phenoInfo(tSet, "rna")[, c("cellid", "organ_id")])))
  cells <- data.table(
    'id' = cells[,.I],
    cells,
    'tissue_id' = vapply(cells$organ_id, function(x) { which(tissues$name %in% x) }, FUN.VALUE = numeric(1))
    )
  cells <- cells[, organ_id := NULL]
  colnames(cells) <- c("id", "name", "tissue_id")

  ## SAMPLES
  samples <- data.table("id" = seq_along(phenoInfo(tSet, "rna")$samplename),
                        "compound_id" = vapply(phenoInfo(tSet, "rna")$drugid, function(x) which(compounds$name == x), FUN.VALUE = numeric(1)),
                        "cell_id" = vapply(phenoInfo(tSet, "rna")$cellid, function(x) which(cells$name == x), FUN.VALUE = numeric(1)),
                        "name" = phenoInfo(tSet, "rna")$samplename,
                        "dose" = phenoInfo(tSet, "rna")$dose_level,
                        "concentration" = phenoInfo(tSet, "rna")$concentration,
                        "time" = phenoInfo(tSet, "rna")$duration,
                        "replicate" = phenoInfo(tSet, "rna")$individual_id
  )

  # Subset samples to only samples with molecular profiless
  if (length(sensitivityInfo(tSet) > 0)) {
    s <- unlist(sensitivityInfo(tSet)[, c("Control", "Low", "Middle", "High")]) %>% na.omit() # Samples with viability measures
    samp <- s[which(s %in% samples$name)] # The samples with viability values which also have gene expression

    ## VIABILITY
    # Repeat each time value for control,
    times <- vapply(unlist(sensitivityInfo(tSet)$duration_h), function(x) { rep(x, 4) }, FUN.VALUE = character(4))
    times <- times[which(s %in% samp)]

    # Subset viabiltity measures to only the samples with molecular profiles
    viab <- unlist(as.data.frame(tSet@sensitivity$raw[,,2])) %>% na.omit()
    viab <- viab[which(samp %in% s)]

    ## TODO:: Write a script to make sure NAs line up between sensitivityInfo and viabilities
    viabilities <- data.table(
      "sample_id" = vapply(samp, function(x) which(samples$name %in% x), FUN.VALUE = numeric(1)),
      "viability" = viab
    )
  }
  ## PATHWAYS

  ## GENERATE PATHWAYS TABLE FROM FILE

  # Read in file lines a stings
  lines <- readLines("../metadata/pathways_raw.txt")

  # Split the lines in
  line_list <- strsplit(lines, split = ' ')

  # Split the lines into vectors on tabs
  line_list <- sapply(unlist(line_list), function(x) { strsplit(x, split = "\t")})

  # Name the list items by pathway
  names(line_list) <- vapply(seq_along(line_list), 
                             function(x) {line_list[[x]][1] }, 
                             FUN.VALUE = character(1))

  # Assign the vectors to a list of pathways, dropping the pathway name and url
  pathway_vectors <- lapply(names(line_list), function(name) { line_list[[name]][c(-1, -2)] })

  # Name the pathways
  names(pathway_vectors) <- names(line_list)
  rm(line_list)

  # Create the pathways table
  pathways <- data.table(
    "id" = seq_along(names(pathway_vectors)),
    "name" = names(pathway_vectors)
  )

  ## PATHWAYS_GENES
  pathway_gene_list <- lapply(pathways$name, function(name) {
    data.table(
      "pathway_id" = rep(which(pathways$name %in% name), length(which(gene_annotations$Symbol %in% pathway_vectors[[name]]))),
      "gene_id" = which(gene_annotations$Symbol %in% pathway_vectors[[name]])
    )
  })
  pathways_genes <- rbindlist(pathway_gene_list)

  ## PATHWAYS_DATASETS
  pathways_datasets <- as.data.table(expand_grid(unique(pathways_genes$pathway_id), unique(datasets$id)))
  colnames(pathways_datasets) <- c('pathway_id', 'dataset_id')

  ## TODO:: Fix this to not require all pathways map to a gene
  # ## VALIDATE RESULTS OF PATHWAYS
  # for (i in seq_along(pathway_vectors)) {
  #   if (all(!(genes[ pathways_genes %>% filter(id == i) %>% pull(gene_id), ]$name %in% pathway_vectors[[i]]))) {
  #     stop(paste(names(pathway_vectors)[i], paste0(pathway_vectors[[i]], collapse = ", ")))
  #   }
  # }

  ## compoundS_DATASETS
  compounds_datasets <- expand.grid(compounds$id, datasets$id)
  compounds_datasets$unique_id <- unique(phenoInfo(tSet, 'rna')$dataset_drugid)
  colnames(compounds_datasets) <- c("compound_id", "dataset_id", "compound_uid")

  ## GENES_DATASETS
  genes_datasets <- expand.grid(genes$id, datasets$id)
  colnames(genes_datasets) = c("gene_id", "dataset_id")

  ## DATASETS_SAMPLES
  datasets_samples <- expand.grid(datasets$id, samples$id)
  colnames(datasets_samples) <- c("dataset_id", "sample_id")

  ### LARGE TABLES ####
  ##TODO:: Determin what why there is a subset call here with no parameters?
  phenoInfo <- subset(phenoInfo(tSet, "rna"))
  molecularProf <- molecularProfiles(tSet, "rna")[, as.character(phenoInfo$samplename)]

  system.time({
    compound_gene_list  <- lapply(samples$id, function(i) {
      data.table(samples$id[rep(i, length(molecularProf[, i]))],
                 which(names(molecularProf[ , i]) %in% genes$name),
                 rep(0, length(molecularProf[, i])),
                 molecularProf[ , i]
      )
    })
    compound_gene_response <- rbindlist(compound_gene_list)
  })

  compound_gene_response <- data.table(seq_len(nrow(compound_gene_response)), compound_gene_response)
  colnames(compound_gene_response) <- c('id', 'sample_id', 'gene_id', 
                                    'analysis_id', 'expression')
  rownames(compound_gene_response) <- NULL

  #### DIFFENTIAL GENE EXPRESION ANALYSIS ####
  library(limma)

  #### ALL compoundS AT THE SAME TIME

  # Average replicate values for each compound gene pair
  SE <- tSet@molecularProfiles$rna
  eset <- as(SE, 'ExpressionSet')
  eset$drugid <- make.names(eset$drugid)
  
  # Samples, compounds, dose
  targets <- as.data.frame(pData(eset)[, c("samplename", "drugid", 
                                           "dose_level", "duration")])
  
  # Create combined design condition
  compound <- factor(targets$drugid)
  dose <- factor(targets$dose_level)
  duration <- factor(targets$duration)
  
  design <- model.matrix(~0 + compound:dose:duration)
  colnames(design) <- gsub(':', '.', colnames(design))
  
  # Fit a linear model based on design matrix
  fit <- lmFit(eset, design)

  # Set up parrallelization cluster to construct the contrast coefficients 
  #   we want to extract from the model
  ##TODO:: Implement parallelization with biocParallel
  #snowParam <- snowParam(workers=(parralel::detectCores() - 1))
  
  # Generate a vector of all valid dose-duration combinations for each compound
  coefs <- c()
  if (name(tSet) == 'drugMatrix_rat') {
    for (drg in levels(compound)) {
      doseDurations <- grep('High|Middle|Low',
                            paste0('compound', pData(eset)$drugid, 
                                   '.dose', pData(eset)$dose_level, 
                                   '.duration', pData(eset)$duration),
                            value=TRUE)
      coefs <- c(coefs,
                 paste0(doseDurations, '-', 
                        gsub('compound.*dose',
                             'compoundDMSO.dose',
                          gsub('doseHigh|doseMiddle|doseLow', 
                              'doseControl', doseDurations))
                 )
      )
    }
  } else {
    for (drg in levels(compound)) {
      doseDurations <- grep('High|Middle|Low',
                            paste0('compound', pData(eset)$drugid, 
                                   '.dose', pData(eset)$dose_level, 
                                   '.duration', pData(eset)$duration),
                            value=TRUE)
      coefs <- c(coefs,
                 paste0(doseDurations, '-', gsub('doseHigh|doseMiddle|doseLow', 
                                                 'doseControl', doseDurations))
      )
    }
  }
  
  # Remove duplicated contrasts
  coefs <- unique(coefs)
  
  # Make contrasts matrix for fitting LM
  contrasts <- makeContrasts(contrasts=coefs, levels=design)  
  contrFit <- contrasts.fit(fit, contrasts)
 
  # Predict coefficients using emperical Bayes moderation of SE
  # Generate a t-stat, moderated F-stat and log-odds of differential expression
  stats <- eBayes(contrFit)

  saveRDS(stats, file = paste0('../results/', 'limma_stats', lab_out, '.rds'))
  
  compoundNames <- make.names(compounds$name)
  resultList <- lapply(coefs, function(coef) {
    # Disassmble contrasts into annotations for this statistical test
    annotations <- unlist(strsplit(
      unlist(lapply(strsplit(coef, '-'), function(str) str[1])), 
      'compound|.dose|.duration'))
    # Get the compound_id from 
    compound_id <- which(compoundNames %in% annotations[2])
    data.table(
      "gene_id" = seq_len(nrow(stats)),
      "compound_id" = rep(compound_id, nrow(stats)),
      "dose" = rep(annotations[3], nrow(stats)), 
      "time" = rep(annotations[4], nrow(stats)),
      topTreat(stats, 
      coef=coef,
      sort.by="none",
      number='all',
      adjust.method="BH")[, c("logFC", "B", "P.Value", "adj.P.Val", "AveExpr")]
    )
  })
  analysis <- rbindlist(resultList)

  # Update column names
  colnames(analysis)[5:length(colnames(analysis))] <- c("fold_change", 
                                                        "log_odds", 
                                                        "p_value", 
                                                        "fdr", 
                                                        "avg_expr")
  analysis <- data.table("id" = seq_len(nrow(analysis)), analysis)
  
  ##TODO:: Move this to table creation
  # Fix data type for time in samples
  samples[, time:= as.character(time)]
  
  # Join analysis to samples to get sample_id for each compound/dose/time combination
  setkey(samples, compound_id, dose, time)
  setkey(analysis, gene_id, compound_id, dose, time)
  analysis[samples[replicate==1 & dose %in% unique(analysis)$dose, ], 
           on=c(compound_id='compound_id', dose='dose', time='time'), 
           sample_id := i.id]
  
  # Label compound_gene_response rows by matching analysis id
  # This is a join operation (due to the on parameter)
  setkey(compound_gene_response, sample_id, gene_id)
  compound_gene_response[analysis,
                     on=c(sample_id='sample_id', gene_id='gene_id'),
                     analysis_id := i.id]
  
  ## FIXME:: Write some checks for mapping of the analysis table!
  # Check that the mappings line up
  # if (any(compound_gene_response[analysis_id > 0, unique(analysis_id)] != analysis[order(id)]$id)) {
  #   stop("Mapping issue between compound_gene_response and analysis tables!")
  # }

  de <- merge(analysis, genes, by.x = 'gene_id', by.y = 'id')
  de <- merge(de, compounds, by.x = 'compound_id', by.y = 'id')
  colnames(de) <- c('compound_id', 'gene_id', 'analysis_id', 'dose', 'time', 
                    'fold_change', 'log_odds', 'p_value', 'fdr', 'avg_expr',
                    'sample_id', 'gene', 'compound')
  
  fwrite(de, file = paste0('../results/analysiscompoundsGenes', lab_out, '.csv'))
  
  # Drop extra columns used for matching
  analysis <- analysis[, .SD, .SDcols = !c('gene_id', 'compound_id', 'dose', 
                                           'time', 'sample_id')]

  if (name(tSet) == 'drugMatrix_rat') {
    viabilities <- data.table()
  }

  #### EXPORT AS CSV FILES ####
  for (df in c("compounds", "compounds_datasets", "datasets", "genes_datasets", "genes",
               "gene_annotations", "compound_annotations", "datasets_samples", "samples",
               "species", "tissues", "cells", "viabilities",
               "analysis", "pathways", "pathways_genes", "pathways_datasets")
  ) {
    fwrite(get(df), file = paste0('../results/',df, lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")
  }

  # Faster for writing to disk with large (> 1 million rows) data frames/tables
  fwrite(compound_gene_response, file = paste0("../results/", "compound_gene_response", lab_out, ".csv"), row.names = FALSE, na = "", eol = "\r\n", sep = ",")

  rm(list = c("compounds", "compounds_datasets", "datasets", "genes_datasets", "genes",
       "gene_annotations", "compound_annotations", "datasets_samples", "samples",
       "species", "tissues", "cells", "viabilities",
       "analysis", "pathways", "pathways_genes", "pathways_datasets", 'compound_gene_response'))
  gc()
}


#' Define an S4 Generic for the methods::as function
#'
#' This will allow creation of new definitons for object conversions
#'
#' @export
setGeneric('as', function(object, value) methods::as(object, value, ...))

#' Convert factor vector to numeric vector
#' 
#' @param f A \code{factor} vector with numeric levels
#' 
#' @return A \code{numeric} vector of the factor levels in the same order
#' 
#' @export
setMethod('as',
          signature('factor'),
          function(object, value) {
            switch(value,
                   'numeric'={as.numeric(levels(object))[object]},
                   'character'={as.character(object)},
                   'list'={as.list(object) },
                   stop(paste0('Cannot coerce to class ', value))
            )
          })

#' Coerce a SummarizedExperiment object to an ExpressionSet
#'
#' @warning This method assumes that all slots not present a SummarizedExperiment
#'   were moved into the metadata list from the original ExpressionSet
#'
#'
#' @importFrom SummarizedExperiment colData rowData assays assay 
#' @importFrom S4Vectors metadata
#' @export
##FIXME:: Subset meta data colums before making the eSet
setMethod('as',
          signature('SummarizedExperiment'),
          function(object, value) {
            if (value != 'ExpressionSet') 
              as(object, value)
            else
              Biobase::ExpressionSet(
                as.environment(as.list(assays(object))),
                phenoData=as(colData(object), 'AnnotatedDataFrame'),
                featureData=as(rowData(object), 'AnnotatedDataFrame'),
                experimentData=metadata(object)$experimentData,
                annotation=metadata(object)$annotation,
                protocolData=metadata(object)$protocolData
              )
          })

#' Coerce a DFrame object to an AnnotatedDataFrame
#'
#' @importFrom S4Vectors AnnotatedDataFrame
#' @export
setMethod('as',
          signature('DFrame'),
          function(object, value) {
            if (value != "AnnotatedDataFrame")
              as(object, value)
            else
              AnnotatedDataFrame(as.data.frame(object))
          })




