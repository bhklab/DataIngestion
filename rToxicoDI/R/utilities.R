#' Processes a tSet name into the appropriate associated file name
#'
#' @param fileName `character` The name of a tSet to parse into its file name.
#'
#' @return `character` The file name formatted approriately.
#'
#' @md
#' @keywords internal
#' @noRd
.preprocessFileName <- function(fileName) {

    # remove file format
    fileName <- gsub('.csv', '', fileName)

    # split on spaces
    fileName <- unlist(strsplit(fileName, ' '))

    # paste back together with all but the first item as lowercase
    if (length(fileName) > 1) {
        fileName <- c(fileName[1], tolower(fileName[2:length(fileName)]))
        fileName <- paste(fileName, collapse='_')
        fileName <- gsub('_ldh', 'ldh', fileName)
    }

    return(paste(fileName, '.csv'))
}