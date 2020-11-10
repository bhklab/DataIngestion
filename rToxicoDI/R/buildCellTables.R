#' Build all tables related to cells for ToxicoDB and write to disk
#'
#' @param path `character` path to the directory where the ToxicoSet data was
#'  extracted
#' @param outDir `character` path to save table csvs to
#' @param ... `pairlist` force further arguments to be named
#'
#' @import data.table
#' @export
buildCellTables <- function(path='procdata', outDir='latest', ...) {
    # -- ensure the output directory exists
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # -- load the compound table for each tSet
    files <- list.files(file.path(path, 'cell'), pattern='csv', full.names=TRUE)
    cellTables <- lapply(files, fread)
    names(cellTables) <-
        trimws(gsub('^.*/|.csv$', '', files))

    # -- join tSet cell tables
    cell <- rbindlist(cellTables)
    cell$dataset <- unlist(mapply(rep, x=names(cellTables),
        times=vapply(cellTables, nrow, numeric(1)), SIMPLIFY=FALSE))

    # -- tissues
    tissue <- unique(cell[, .(tissue_id)])
    tissue[, `:=`(id=seq_len(.N), code=NA)]
    setnames(tissue, 'tissue_id', 'name')
    setcolorder(tissue, c('id', 'name', 'code'))

    # map tissue_id to cell table
    setkeyv(cell, 'tissue_id')
    setkeyv(tissue, 'name')
    cell[tissue, tissue_id := i.id]

    # build datasets table if it doesn't already exist
    if (!file.exists(file.path(outDir, 'dataset.csv'))) {
        dataset <- data.table(
            id=seq_along(names(cellTables)),
            name=names(cellTables))
    } else {
        dataset <- fread(file.path(outDir, 'dataset.csv'))
    }

    # map dataset id to dataset
    setkeyv(cell, 'dataset')
    setkeyv(dataset, 'name')
    cell[dataset, dataset := i.id]

    # map tissue_id to tissue
    setkeyv(cell, 'tissue_id')
    setkeyv(tissue, 'name')
    cell[dataset, tissue_id := i.id]

    # -- species
    species <- unique(cell[, .(dataset, species)])
    setnames(species, c('species', 'dataset'), c('name', 'dataset_id'))
    setkeyv(species, 'dataset_id')

    # -- cell
    cell[, id := .GRP, by=name]

    # make join table
    dataset_cell <- unique(cell[, .(id, dataset)])
    setnames(dataset_cell, c('id', 'dataset'), c('cell_id', 'dataset_id'))

    cell <- unique(cell[, .(id, name, tissue_id)])
    setcolorder(cell, c('id', 'tissue_id', 'name'))
    setkeyv(cell, 'id')

    # -- write to disk
    for (table in c('tissue', 'species', 'cell', 'dataset_cell')) {
        fwrite(get(table), file.path(outDir, paste0(table, '.csv')))
    }
}


if (sys.nframe() == 0) {
    library(data.table)
    path='procdata'
    outDir='latest'
}