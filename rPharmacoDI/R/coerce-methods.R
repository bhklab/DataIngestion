#' Coerce a SummarizedExperiment object to a long data.table, retaining the data in all assays, rowData, colData
#'
#' @param from [`SummarizedExperiment`]
#' @param to [`character`]
#'
#' @return ['data.table]
#'
#' @importFrom SummarizedExperiment colData rowData
#' @import data.table
#'
#' @export
setMethod("coerce",
    signature(from="SummarizedExperiment", to="data.table"),
    function(from, to) {

    # Get the data to join on
    assaysDT <- .assaysToLongDT(assays(from), names(assays(from)))

    rowDataDT <- data.table(as(rowData(from), 'data.frame'), keep.rownames='features')
    colnames(rowDataDT) <- paste0('row_', colnames(rowDataDT))

    colDataDT <- data.table(as(colData(from), 'data.frame'), keep.rownames='samples')


    ## TODO:: Determine if there are any other items we need from `metadata`
    if('protocalData' %in% names(metadata(from))) {
                protocolDataDT <- data.table(as(metadata(from)$protocolData, 'data.frame'), keep.rownames='samples')
        colDataDT <- merge.data.table(colDataDT, protocolDataDT, by='samples')
    }

    colnames(colDataDT) <- paste0('col_', colnames(colDataDT))

    # Join colData and rowData to the long assaysDT
    longSummarizedExperiment <- merge.data.table(assaysDT, rowDataDT, by='features',
                                                 allow.cartesion=TRUE, fill.missing=TRUE)
    longSummarizedExperiment <- merge.data.table(longSummarizedExperiment, colDataDT, by='samples',
                                                 allow.cartesian=TRUE, fill.missing=TRUE)

    return(longSummarizedExperiment)
})


#' Converts each assay in a list of assays to a data.table, then iteratively merged the data.tables by the shared
#'    feature (rownames) and samples (colnames) columns.
#'
#' @param assays [`list`] A list of `matrix` objects, as returned by the `SummarizedExperiment::assays` function.
#' @param assayNames [`character`] Names for each assay, e.g. names of the list returned by
#'     `SummarizedExperiment::assays`
#'
#' @import data.table
#' @importFrom SummarizedExperiment assays
#'
#' @noRd
#' @keywords internal
#' @export
.assaysToLongDT <- function(assays, assayNames) {
        # Convert each assay to a data.table, return a list
        assaysDtL <- mapply(FUN=.assayToLongDT,
                            assays, assayNames,
                            SIMPLIFY=FALSE)

        # Metaprogram some R code to join an unknown number of data.tables by key
        codeDataTables <- paste0('assaysDtL[[', seq_along(assaysDtL), ']]')
        if (length(codeDataTables) - 2 > 0) {
            codeOperators <- c('[', rep('][', length(codeDataTables) - 2), ']')
        } else {
            codeOperators <- c('[', ']')
        }

        zippedStrings <- unlist(mapply(c, codeDataTables, codeOperators, SIMPLIFY=FALSE))
        chainedJoinExpression <- parse(text=paste0(zippedStrings, collapse=''))

        # Evaluate the R code, this allows merging n data.tables without using inefficient Reduce function
        assayDT <- eval(chainedJoinExpression)

        return(assayDT)
}

#' Coerce an assay matrix from a SummarizedExperiment to a long data.table
#'
#' @param assay [`matrix`] A numeric matrix with features as rows and samples as columns, as returend by
#'     `SummarizedExperiment::assay`
#' @param assayName [`character`] The name of the assay, this becomes the column name for the assay values
#'
#' @return [`data.table`] where the columns of the assay values are in column named `assayName` and the
#'     column names are in the `samples` column
#'
#' @noRd
#' @keywords internal
#' @export
.assayToLongDT <- function(assay, assayName) {
    DT <- data.table(assay, keep.rownames='features')
    DTmolten = melt.data.table(DT,  id.vars='features', value.vars=rownames(assay),
                                variable.name='samples', value.name=paste0('assay_', assayName)) # So we can identify assays vs annotations
    data.table::setkey(DTmolten, features, samples) # Set keys for joins
    return(DTmolten)
}