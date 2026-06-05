% Generated manually from R/result_utils.R
\name{sortScore}
\alias{sortScore}
\title{Sort scoring results by resolution}
\usage{
sortScore(result.table, .f = mean)
}
\arguments{
\item{result.table}{Subtypist result table.}

\item{.f}{A function to summarize \code{Score} by \code{resolution} (default: \code{mean}).}
}
\value{
A tibble with \code{resolution} and summarized \code{value}.
}
\description{
Summarizes absolute cluster-level \code{Score} values by resolution.
}
\examples{
sortScore(res[[2]],mean)
}
