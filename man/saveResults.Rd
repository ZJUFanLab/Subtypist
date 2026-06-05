% Generated manually from R/result_utils.R
\name{saveResults}
\alias{saveResults}
\title{Save Subtypist annotation results to disk}
\usage{
saveResults(result.table, path, name)
}
\arguments{
\item{result.table}{A data frame or tibble containing Subtypist results.}

\item{path}{Directory path to save the file.}

\item{name}{File name ending in ".csv", ".tsv", or ".xlsx".}
}
\value{
The processed results (invisibly).
}
\description{
Writes Subtypist result tables to CSV, TSV, or XLSX files.
}
\examples{
saveResults(res[[2]], path = "output", name = "subtype.csv")
saveResults(res[[2]], path = "output", name = "subtype.xlsx")
}
