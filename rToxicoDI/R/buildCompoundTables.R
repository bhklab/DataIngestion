#' Build all tables related to compound for ToxicoDB and write to disk
#'
#' @param path `character` path to the directory where the ToxicoSet data was
#'  extracted
#' @param outDir `character` path to save table csvs to
#' @param ... `pairlist` force further arguments to be named
#' @param annotColMap A `vector` named vector mapping the column names in
#'   the annotation .csv to the correct database table names. If we don't
#'   have data for some columns set them to NA and they will be added to the
#'   final .csv as NA columns.
#'
#'
#' @import data.table
#' @export
buildCompoundTables <- function(path='procdata',
    annotPath='metadata/drug_annotations.csv', outDir='latest', ...,
    updatedCompounds='metadata/old_newDrugmapping.csv',
    synonymPath='metadata/Compound_synonyms.csv',
    annotColMap=c(compound_id='compound_id', pubchem='pubchem', ctd='ctd',
        chembl='chembl', drugbank='drugbank', targets='targets',
        carcinogenicity='carcinogenicity', class_in_vivo='classif_in_vivo',
        class_in_vitro='classif_in_vitro', class_name='class_name',
        smiles='smiles', inchikey='inchikey', name='name',
        NTP='ntp', IARC='iarc', DILI_status='dili_status')
    )
{
    # ensure the output directory exists
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # load the compound table for each tSet
    files <- list.files(file.path(path, 'compound'), pattern='csv', full.names=TRUE)
    compoundTables <- lapply(files, fread)  # disable quoting input
    names(compoundTables) <-
        trimws(gsub('^.*/|.csv$', '', files))

    # load annotations file
    annotations <- fread(annotPath)[, drug_id := NULL]
    annotations <- unique(annotations)
    renameMap <- na.omit(annotColMap)
    setnames(annotations, renameMap, names(renameMap), skip_absent=TRUE)

    # build compound annotation table
    compoundTable <- unique(rbindlist(compoundTables))[, 'name']
    compoundTable[, name := gsub('\\\"', '', name)]
    compoundTable

    setkeyv(compoundTable, 'name')
    setkeyv(annotations, 'name')
    compound <- unique(merge.data.table(compoundTable, annotations, all.x=TRUE,
        sort=FALSE))

    # ensure nothing weird happened in the join
    if (!setequal(compound$name, compoundTable$name))
        stop(.context(), 'the compound table has more or less compound after
            joining with the annotations table. Something has gone wrong!')

    # ensure no duplciated
    if (any(duplicated(compound$name)))
        stop(.erroMsg(.context(), 'There are duplicated compound names in ',
            'the compound table!'))

    # build compound tables
    compound[, id := seq_len(.N)]
    missingAnnotCols <- annotColMap[is.na(annotColMap)]
    assignNAstring <- deparse(dput(missingAnnotCols), width.cutoff=500L)
    assignNAstring <- gsub('^c\\(', '`:=`(', assignNAstring)
    compound[, eval(str2lang(assignNAstring))]

    annotCols <- c(setdiff(names(annotColMap), c('name', 'compound_id')), 'id')
    compound_annotations <- compound[, ..annotCols]
    compound <- compound[, .(id, name)]

    # subset and sort columns
    setnames(compound_annotations, 'id', 'compound_id')
    inColNames <- names(annotColMap) %in% colnames(compound_annotations)
    setcolorder(compound_annotations, names(annotColMap)[inColNames])

    # build compound datasets tables
    compound_dataset <- mapply(cbind, compoundTables, names(compoundTables),
        SIMPLIFY=FALSE)
    compound_dataset <- unique(rbindlist(compound_dataset))
    setnames(compound_dataset, 'V2', 'dataset')
    setnames(compound_dataset, 'dataset_drugid', 'compound_uid')

    # build datasets table
    dataset <- data.table(
        id=seq_along(names(compoundTables)),
        name=names(compoundTables))

    setkeyv(dataset, 'name')
    setkeyv(compound_dataset, 'dataset')
    compound_dataset[dataset, `:=`(dataset_id=i.id)]
    setkeyv(compound, 'name')
    setkeyv(compound_dataset, 'name')
    compound_dataset <- merge(compound_dataset, compound)[, .(id, dataset_id, compound_uid)]
    setnames(compound_dataset, 'id', 'compound_id')

    synonym <- fread(synonymPath)
    synonym[, id := NULL]
    setkeyv(synonym, 'Drug')
    setkeyv(compound, 'name')
    compound_synonyms <- merge.data.table(compound, synonym, by.x='name',
        by.y='Drug', allow.cartesian=TRUE)
    compound_synonyms[, name := NULL]
    setnames(compound_synonyms, c('id', 'Synonym'), c('compound_id', 'synonym'))

    # filter for bad synonyms
    compound_synonyms <- compound_synonyms[nchar(synonym) < 30 & synonym != 'NO_PUBCHEM_ENTRY', ]

    for (table in c('compound', 'compound_annotations', 'compound_dataset',
        'dataset', 'compound_synonyms'))
    {
        fwrite(get(table), file.path(outDir, paste0(table, '.csv')))
    }
}

if (sys.nframe() == 0) {
    library(data.table)
    path='procdata'
    annotPath='metadata/drug_annotations.csv'
    outDir='latest'
    updatedCompounds='metadata/old_newDrugmapping.csv'
    synonymPath='metadata/Compound_synonyms.csv'
    annotColMap=c(compound_id='compound_id', pubchem='pubchem', ctd='ctd',
        chembl='chembl', drugbank='drugbank', targets='targets',
        carcinogenicity='carcinogenicity', class_in_vivo='classif_in_vivo',
        class_in_vitro='classif_in_vitro', class_name='class_name',
        smiles='smiles', inchikey='inchikey', name='name',
        NTP='ntp', IARC='iarc', DILI_status='dili_status')


    buildCompoundTables()
}