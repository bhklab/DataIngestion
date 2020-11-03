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
    annotPath='metadata/Drug_annotations_V2.1.csv', outDir='latest', ...,
    annotColMap=c(compound_id='id', pubchem=NA, cid='cid', chembl=NA,
    drugbank = NA, targets=NA, carcinogenicity='Carcinogenicity',
    class_in_vivo='Classif. in vivo', class_in_vitro='Classif. in vitro',
    class_name=NA, smiles='smiles', inchikey='inchikey', name='unique.drugid'))
{
    # ensure the output directory exists
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # load the compound table for each tSet
    files <- list.files(file.path(path, 'compound'), pattern='csv', full.names=TRUE)
    compoundTables <- lapply(files, fread, sep='\n')
    names(compoundTables) <-
        trimws(gsub('^.*/|.csv$', '', files))

    # load annotations file
    annotations <- fread(annotPath)
    renameMap <- na.omit(annotColMap)
    setnames(annotations, renameMap, names(renameMap), skip_absent=TRUE)

    # build compound annotation table
    compoundTable <- unique(rbindlist(compoundTables))
    setkeyv(compoundTable, 'name')
    setkeyv(annotations, 'name')
    compounds <- compoundTable[annotations]

    # ensure nothing weird happened in the join
    if (setequal(compounds$name, compoundTable$name))
        stop(.context(), 'the compounds table has more or less compounds after
            joining with the annotations table. Something has gone wrong!')

    # build compound tables
    compounds[, id := seq_len(.N)]
    missingAnnotCols <- annotColMap[is.na(annotColMap)]
    assignNAstring <- deparse(dput(missingAnnotCols), width.cutoff=500L)
    assignNAstring <- gsub('^c\\(', '`:=`(', assignNAstring)
    compounds[, eval(str2lang(assignNAstring))]

    annotCols <- c(setdiff(names(annotColMap), 'name'), 'id')
    compound_annotations <- compounds[, ..annotCols]
    compounds <- compounds[, .(id, name)]

    # build compounds datasets tables
    compound_dataset <- mapply(cbind, compoundTables, names(compoundTables),
        SIMPLIFY=FALSE)
    compound_dataset <- rbindlist(compound_dataset)
    setnames(compound_dataset, 'V2', 'dataset')

    # build datasets table if it doesn't already exist
    if (!file.exists(file.path(outDir, 'dataset.csv'))) {
        dataset <- data.table(
            id=seq_along(names(compoundTables)),
            name=names(compoundTables))
    } else {
        dataset <- fread(outDir, 'dataset.csv')
    }

    setkeyv(dataset, 'name')
    setkeyv(compound_dataset, 'dataset')
    compound_dataset[dataset, `:=`(dataset_id=i.id, compound_uid=NA)]
    setkeyv(compounds, 'name')
    setkeyv(compound_dataset, 'name')
    compound_dataset <- merge(compound_dataset, compounds)[, .(id, dataset_id, compound_uid)]
    setnames(compound_dataset, 'id', 'compound_id')

    for (table in c('compounds', 'compound_annotations'))
}