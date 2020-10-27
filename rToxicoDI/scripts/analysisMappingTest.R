# Analysis Mapping Test
# Take random sample of analysis ids, and evaluate that the mappings are correct
for (i in sample(seq_len(nrow(analysis)), 150)) {
  drug_ids <- analysis[i,]$drug_id 
  gene_ids <- analysis[i,]$gene_id
  gene_ids2 <- drug_gene_response[analysis_id == i]$gene_id
  drug_ids2 <- samples[drug_gene_response[analysis_id == i]$sample_id,]$drug_id
  if (!all(drug_ids %in% drug_ids2)) {
    stop(paste0('Mapping issue for analysis ', i, ' mismatches at drugs', paste(which(!(drug_ids %in% drug_ids2)), collapse = ' ')))
  } else if (!all(gene_ids %in% gene_ids2)) {
    stop(paste0('Mapping issue for analysis ', i, ' mismatches at genes', paste(which(!(gene_ids %in% genes_ids2)), collapse = ' ')))
  }
  print(i)
}

## Analysis Mapping Test
for (i in sample(nrow(de):nrow(analysis), 150)) {
  gene_names <- genes2[drug_gene_response[analysis_id == i]$gene_id,]$name
  drug_names <- drugs2[samples[drug_gene_response[analysis_id == 1]$sample_id]$drug_id, ]$name
  if (!all(gene_names %in% genes$name)) {
      stop(paste0('Mapping issue for analysis ', i, ' mismatches at genes ', paste(which(!(gene_names %in% genes$name)), collapse = ' ')))
  } else if (!all(drug_names %in% drugs$name)) {
      stop(paste0('Mapping issue for analysis ', i, ' mismatches at drugs', paste(which(!(drug_names %in% drugs$name)), collapse = ' ')))
  }
  print(i)
}