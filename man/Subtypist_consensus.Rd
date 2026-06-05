% Generated manually from R/consensus_utils.R
\name{Subtypist_consensus}
\alias{Subtypist_consensus}
\title{Inter-resolution stability scoring for Subtypist results}
\usage{
Subtypist_consensus(
  object = NULL,
  result.table = NULL,
  result.list = NULL,
  selected.resolution = NULL,
  evaluate.all = TRUE,
  window.size = 1,
  prefix = "Subtypist",
  suffix = NULL,
  cluster.columns = NULL,
  score.summary = mean,
  verbose = TRUE
)
}
\arguments{
\item{object}{A Seurat object containing Subtypist cluster columns in
\code{object@meta.data}.}

\item{result.table}{A Subtypist result table. The canonical column names are
\code{resolution}, \code{merged_cluster}, \code{initial_cluster},
\code{phenotypic_molecules}, and \code{Score}.}

\item{result.list}{A result list returned by \code{Subtypist_merge()}. If
provided, \code{object} and \code{result.table} are taken from this list.}

\item{selected.resolution}{Optional resolution to assess when
\code{evaluate.all = FALSE}. If \code{NULL}, the function uses the resolution
with the highest summarized absolute \code{Score}.}

\item{evaluate.all}{Whether to evaluate every resolution in
\code{result.table}. Default is \code{TRUE}.}

\item{window.size}{Number of neighboring resolution ranks on each side of
each evaluated resolution used for stability scoring. Default is \code{1}.}

\item{prefix}{Prefix used by Subtypist metadata columns. Default is
\code{"Subtypist"}, matching columns such as \code{Subtypist_snn_res.0.5}.}

\item{suffix}{Optional suffix appended to Subtypist metadata columns.}

\item{cluster.columns}{Optional named character vector mapping resolution
values to metadata column names. Use this when cluster columns do not follow
the standard Subtypist naming pattern.}

\item{score.summary}{Function used to summarize cluster-level \code{Score}
values when \code{selected.resolution = NULL}. Default is \code{mean}.}

\item{verbose}{Whether to print the selected and evaluated resolutions.}
}
\value{
A list with augmented \code{result.table}, \code{consensus.match.table},
\code{selected.resolution}, and \code{evaluated.resolutions}.
}
\description{
\code{Subtypist_consensus()} adds inter-resolution stability scores to a
Subtypist result table.
}
