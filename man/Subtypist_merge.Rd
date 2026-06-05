% Generated manually from R/Subtypist.R
\name{Subtypist_merge}
\alias{Subtypist_merge}
\title{Reference-free identification of cell subtypes}
\usage{
Subtypist_merge(
  object,
  min.resolution = 0.3,
  max.resolution = 1.5,
  by = 0.1,
  max.steps = 100,
  marker_assay = "RNA",
  cluster_assay = "integrated",
  n_candidate_markers = 300,
  min.pct.1 = 0.1,
  min.diff = 0.4,
  min.avg_log2FC = 0.5,
  logfc.threshold = 0.1,
  prefix = "Subtypist",
  accelerated = FALSE,
  algorithm = 1,
  tua = 0,
  verbose = FALSE,
  regulation = c("up", "down", "both"),
  top_k = 3,
  termination.mode = c("default", "elbow"),
  elbow.metric = c("mean_score", "total_score"),
  elbow.min.points = 4,
  elbow.min.steps = 2,
  elbow.confirm.steps = 1,
  elbow.window = 1,
  elbow.stability.enabled = FALSE,
  elbow.stability.metric = c("ARI", "NMI"),
  elbow.stability.threshold = 0.85,
  elbow.subsample.fraction = 0.8,
  elbow.subsample.reps = 10,
  elbow.subsample.seed = 1L,
  return.diagnostics = FALSE
)
}
\arguments{
\item{object}{A Seurat object.}

\item{min.resolution}{Minimum clustering resolution.}

\item{max.resolution}{Maximum clustering resolution.}

\item{by}{Resolution increment.}

\item{max.steps}{Maximum number of merge steps per resolution.}

\item{marker_assay}{Assay used for marker detection.}

\item{cluster_assay}{Assay used for Seurat graph clustering in Seurat v4.}

\item{n_candidate_markers}{Number of candidate marker genes retained per cluster.}

\item{min.pct.1}{Minimum percent of cells expressing a marker in the cluster.}

\item{min.diff}{Minimum expression-percent difference for specificity scoring.}

\item{min.avg_log2FC}{Minimum average log2 fold change used by specificity scoring.}

\item{logfc.threshold}{Log fold-change threshold for marker testing.}

\item{prefix}{Prefix for Subtypist metadata columns.}

\item{accelerated}{Whether to use accelerated marker detection in Seurat v5.}

\item{algorithm}{Seurat clustering algorithm.}

\item{tua}{Threshold used by Subtypist to classify a cluster as specific.}

\item{verbose}{Whether to print progress.}

\item{regulation}{Marker regulation direction: up, down, or both.}

\item{top_k}{Number of marker genes used to summarize cluster scores.}

\item{termination.mode}{Termination rule for each resolution. \code{"default"}
keeps the standard Subtypist rule; \code{"elbow"} uses a segmented-regression elbow rule with optional
subsampling stability checks.}

\item{elbow.metric}{Score metric used when \code{termination.mode = "elbow"}.}

\item{elbow.min.points}{Minimum number of score points required to fit an elbow.}

\item{elbow.min.steps}{Minimum merge step before elbow calling.}

\item{elbow.confirm.steps}{Number of additional merge steps observed after the
estimated elbow before accepting it as a stopping point.}

\item{elbow.window}{Local step window around the estimated elbow used for
optional subsampling stability checks.}

\item{elbow.stability.enabled}{Whether to run optional subsampling stability
checks for the elbow mode.}

\item{elbow.stability.metric}{Similarity metric used for subsampling stability.}

\item{elbow.stability.threshold}{Minimum similarity required to accept the elbow.}

\item{elbow.subsample.fraction}{Fraction of cells sampled for each stability replicate.}

\item{elbow.subsample.reps}{Number of subsampling replicates.}

\item{elbow.subsample.seed}{Seed used for subsampling stability checks.}

\item{return.diagnostics}{Whether to return \code{step.history},
\code{elbow.table}, and termination settings.}
}
\value{
A list with \code{Object} and \code{result.table}. If
\code{return.diagnostics = TRUE}, the list also contains \code{step.history},
\code{elbow.table}, and termination settings.
}
\description{
Subtypist evaluates clustering results across multiple resolutions, merges
clusters with insufficient marker specificity, and reports subtype-associated
phenotypic molecules.
}
\examples{
Subtypist_merge(
  object = Seu,
  min.resolution = 0.1,
  max.resolution = 0.2,
  by = 0.1,
  marker_assay = "RNA",
  cluster_assay = "RNA"
)
}
