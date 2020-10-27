#############################################
######### REGENERATE TOXICODB TABLES ########
#############################################

# 2020/03/17
# by Christopher Eeles


# Load Dependencies -------------------------------------------------------

library(data.table)
library(ToxicoGx)

# Generate Individual tSet Tables -----------------------------------------

## Load functions
source('./1_tSetToDatabaseCSV.R')

## Get List of Availabel tSets
availableTSets <- ToxicoGx::availableTSets(saveDir="../tSets")
tSets <- availableTSets[[1]]

## Download any missing tSets and load all tSets
for (tSet in tSets) {
  tSetFiles <- list.files('../tSets', '*rds')
  if (!(paste0(tSet, '.rds') %in% tSetFiles)) {
   # stop('downloadTSet is currently not working, please download the tSets from Zenodo 
    #     at https://doi.org/10.5281/zenodo.3712423 and place them in the tSets folder!')
    ToxicoGx::downloadTSet(tSet, saveDir='../tSets')
  }
  # Load in each
  assign(tSet, readRDS(list.files('../tSets', paste0(tSet, '*'), full.names=TRUE)), 
         envir=globalenv())
}

## Generate database tables from individual tSets
for (tSet in tSets) {
  tSetToDBtables(get(tSet), paste0('_', tSet))
  rm(tSet); gc() # To save memory
}

# Append Together the tSet Tables -----------------------------------------

source('./2_appendTSetToToxicoDB.R')

##TODO:: Generalize this to work with N tSets of unkown name
##FIXME:: Find which tSet has '

appendTSetToToxicoDB(tables, paste0('_', tSets[1]), paste0('_', tSets[2]), paste0('_', tSets[1], '_', tSets[2]))
appendTSetToToxicoDB(tables, paste0('_', tSets[1], '_', tSets[2]), paste0('_', tSets[3]), "_latest")

# Validate the Mapping in the ToxicoDB Tables -----------------------------

source('./3_validateTablesToxicoDB.R')

## Check mappings in the final tables
validateTablesToxicoDB("_latest")