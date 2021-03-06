% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{loadTSetTables}
\alias{loadTSetTables}
\title{Load in the files associated with a set of ToxicoDB tables for each TSet}
\usage{
loadTSetTables(path = "procdata", table)
}
\arguments{
\item{path}{\code{character} Where the table directories are located. Default
is 'procdata'.}

\item{table}{\code{character} The name of the table directory, as generated by
the associated \verb{extract*} method.}
}
\value{
A \code{list} of \code{data.table}s, with each item named for the TSet it was
generated from.
}
\description{
Load in the files associated with a set of ToxicoDB tables for each TSet
}
