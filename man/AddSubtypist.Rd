% Generated manually from R/result_utils.R
\name{AddSubtypist}
\alias{AddSubtypist}
\title{Add Subtypist annotations to a Seurat object}
\usage{
AddSubtypist(
  object = NULL,
  resolution = NULL,
  result.table = NULL,
  result.list = NULL,
  prefix = "Subtypist",
  suffix = NULL,
  select_index = NULL,
  meta.prefix = "phenotypic_molecules_",
  value.suffix = NULL
)
}
\arguments{
\item{object}{A Seurat object.}

\item{resolution}{A vector of resolution values at which to assign annotations.}

\item{result.table}{A data frame containing Subtypist results with
\code{resolution}, \code{merged_cluster}, and \code{phenotypic_molecules}.}

\item{result.list}{Result list from \code{Subtypist_merge()}.}

\item{prefix}{Prefix for Subtypist cluster columns.}

\item{suffix}{A character suffix appended to Subtypist cluster columns.}

\item{select_index}{A named integer vector specifying which phenotype to
select for each \code{merged_cluster}.}

\item{meta.prefix}{Prefix used for new metadata annotation columns.}

\item{value.suffix}{Suffix appended to annotation values.}
}
\value{
A Seurat object with new metadata columns for Subtypist annotations.
}
\description{
Maps \code{phenotypic_molecules} from a Subtypist result table back to cells
through the selected Subtypist cluster column.
}
\examples{
res.table <- tibble::tibble(
  resolution = c(0.1, 0.1, 0.1),
  merged_cluster = c(0, 1, 2),
  initial_cluster = c(0, 1, 2),
  phenotypic_molecules = c("G1", "G2", "G3"),
  Score = c(1, 2, 4)
)
object <- AddSubtypist(
  object = res[["Object"]],
  result.table = res.table,
  resolution = 0.1
)
}

