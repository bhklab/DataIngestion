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
    # remove any leading or trailed whitespace
    fileName <- trimws(fileName)

    return(paste(fileName, '.csv'))
}

## FIXME:: Move this into CoreGx
#' Return the name of the function and the name of the package that function
#'   is in when called within an R function.
#'
#' For providing context in user messages, warnings and errors
#'
#' @param n `integer` How far up the call stack to look for context. Defaults to
#'   2 since it is assumed this function will be used inside of `message`,
#'   `warning` or `stop`.
#'
#' @return `list`:
#' - fun: `character` The name of the function where `.getExectutionContext()`
#' was called
#' - pkg: `character` The name of the package `fun` is from, if applicable.
#'
#' @md
#' @keywords internal
#' @importFrom rlang trace_back
#' @noRd
#' @aliases .context
.getExecutionContext <- function(n=2) {

    # name of function which called this function
    callStack <- rlang::trace_back()$calls
    context <- deparse(callStack[[length(callStack) - n]])

    # remove function arguments
    context <- gsub('\\(.*\\)', '', context)

    return(paste0('\n[', context, '] ', collapse='::'))
}
#' @noRd
.context <- .getExecutionContext

#' @import CoreGx
.errorMsg <- CoreGx::.errorMsg

#' @import CoreGx
.warnMsg <- CoreGx::.warnMsg