#############################################
######### REGENERATE TOXICODB TABLES ########
#############################################

# 2020/03/17
# by Christopher Eeles


# Load Dependencies -------------------------------------------------------

library(data.table)
library(ToxicoGx)

# Generate Individual tSet Tables -----------------------------------------



### Get List of Availabel tSets
#availableTSets <- ToxicoGx::availableTSets(saveDir="../tSets")
#tSets <- availableTSets[[1]]

### Download any missing tSets and load all tSets
#for (tSet in tSets) {
#  tSetFiles <- list.files('../tSets', '*rds')
#  if (!(paste0(tSet, '.rds') %in% tSetFiles)) {
#    # stop('downloadTSet is currently not working, please download the tSets from Zenodo
#    #     at https://doi.org/10.5281/zenodo.3712423 and place them in the tSets folder!')
#    ToxicoGx::downloadTSet(tSet, saveDir='../tSets')
#  }
#  # Load in each
#  assign(tSet, readRDS(list.files('../tSets', paste0(tSet, '*'), full.names=TRUE)),
#         envir=globalenv())
#}

## Load functions
source('R/1_tSetToDatabaseCSV.R')

tSetPaths <- list.files('tSets', pattern='EM.*rds|.*ldh.*rds|.*drugMatrix.*rds', full.names=TRUE)

## Generate database tables from individual tSets
#for (path in tSetPaths) {
#    tSet <- readRDS(path)
#    tSetToDBtables(tSet, paste0('_', name(tSet)))
#}

# Append Together the tSet Tables -----------------------------------------

source('R/2_appendTSetToToxicoDB.R')

##TODO:: Generalize this to work with N tSets of unkown name
##FIXME:: Find which tSet has '

tSets <- gsub('^.*/|\\.rds', '', tSetPaths)
tSetLabs <- paste0('_', tSets)

tSetCombos <- Reduce(c, tSets, accumulate=TRUE)
tableLabs <- unlist(lapply(tSetCombos, paste, collapse='-'))
tableLabs <- paste0('_', tableLabs)

for (i in seq_len(length(tableLabs) - 1)) {
    appendTSetToToxicoDB(tables, tableLabs[i], tSetLabs[i + 1], tableLabs[i + 1])
}


# Validate the Mapping in the ToxicoDB Tables -----------------------------

source('./3_validateTablesToxicoDB.R')

## Check mappings in the final tables
validateTablesToxicoDB("_latest")