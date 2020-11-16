#' Build all tables related to samples for ToxicoDB and write to disk
#'
#' @param path `character` path to the directory where the ToxicoSet data was
#'  extracted
#' @param outDir `character` path to save table csvs to
#' @param ... `pairlist` force further arguments to be named
#'
#' @export
buildSampleTables <- function(path='procdata', outDir='latest', ...) {
    # -- ensure the output directory exists
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # -- load the compound table for each tSet
    files <- list.files(file.path(path, 'sample'), pattern='csv', full.names=TRUE)
    sampleTables <- lapply(files, fread)
    names(sampleTables) <-
        trimws(gsub('^.*/|.csv$', '', files))

    # -- check for required tables have been generated
    # build datasets table if it doesn't already exist
    ## TODO:: Refactor into helpers?
    if (!file.exists(file.path(outDir, 'dataset.csv'))) {
        dataset <- data.table(
            id=seq_along(names(sampleTables)),
            name=names(sampleTables))
    } else {
        dataset <- fread(file.path(outDir, 'dataset.csv'))
    }
    cellPath <- file.path(outDir, 'cell.csv')
    if (!file.exists(cellPath))
        stop(.errorMsg(.context(), 'Cell table now found in ', outDir, ', please ',
            'run buildCellTables() before this function!'))
    compoundPath <- file.path(outDir, 'compound.csv')
    if (!file.exists(compoundPath))
        stop(.errorMsg(.context(), 'Compund table now found in ', outDir, ', ',
            'please run buildCompunoundTables() before this function!'))

    # -- load table dependencies
    cell <- fread(cellPath)
    compound <- fread(compoundPath)

    # -- map samples to cell, compound and datset
    # deal with numeric only sample names
    .checkNameIsChar <- function(DT) typeof(DT[['name']]) == 'character'
    nameIsChar <- unlist(lapply(sampleTables, .checkNameIsChar))
    .nameColAsChar <-
        function(DT) DT[, name := as.character(name)]  # reference symantics
    sampleTables[!nameIsChar] <- lapply(sampleTables[!nameIsChar], .nameColAsChar)

    # build sample table
    sample <- rbindlist(sampleTables)
    sample$dataset <- unlist(mapply(rep, x=names(sampleTables),
        times=vapply(sampleTables, nrow, numeric(1)), SIMPLIFY=FALSE))

    compoundIDinCompound <- sample$compound_id %in% compound$name
    if (!all(compoundIDinCompound)) {
        stop(.errorMsg(.context(), 'The compounds ',
            paste0(sample$compound_id[!compoundIDinCompound], collapse=', '),
            'are not present in the compound table! Please check that',
            ' buildCompoundTables is working correctly.'))
    }

    # map cell
    setkeyv(sample, 'cell_id')
    setkeyv(cell, 'name')
    sample[cell, cell_id := i.id]
    #rm(cell)

    # map compund
    if (!(all(sample$compound_id %in% compound$name)))
        stop(.errorMsg(.collapse(), 'Not all sample compound_ids are in the ',
            'compound table!'))
    setkeyv(sample, 'compound_id')
    setkeyv(compound, 'name')
    sample[compound, compound_id := i.id]
    #rm(compound)

    # map dataset
    setkeyv(sample, 'dataset')
    setkeyv(dataset, 'name')
    sample[dataset, dataset_id := i.id]
    rm(dataset)

    # -- build sample tables
    sample[, id := seq_len(.N)]
    dataset_sample <- unique(sample[, .(id, dataset_id)])
    setnames(dataset_sample, 'id', 'sample_id')

    sample <- sample[, .SD, .SDcols=!grepl('dataset', colnames(sample))]
    setcolorder(sample, c('id', 'compound_id', 'cell_id', 'name', 'dose',
        'time', 'replicate', 'concentration'))

    for (table in c('sample', 'dataset_sample')) {
        fwrite(get(table), file.path(outDir, paste0(table, '.csv')))
    }
}


if (sys.nframe() == 0) {
    library(data.table)
    path='procdata'
    outDir='latest'
}