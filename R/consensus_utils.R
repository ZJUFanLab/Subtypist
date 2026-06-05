#' Inter-resolution stability scoring for Subtypist results
#'
#' `Subtypist_consensus()` adds inter-resolution stability scores to a
#' Subtypist result table.
#'
#' @param object A Seurat object containing Subtypist cluster columns in
#'   `object@meta.data`.
#' @param result.table A Subtypist result table. The canonical column names are
#'   `resolution`, `merged_cluster`, `initial_cluster`, `phenotypic_molecules`,
#'   and `Score`.
#' @param result.list A result list returned by `Subtypist_merge()`. If provided,
#'   `object` and `result.table` are taken from this list.
#' @param selected.resolution Optional resolution to assess when
#'   `evaluate.all = FALSE`. If `NULL`, the function uses the resolution with
#'   the highest summarized absolute `Score`.
#' @param evaluate.all Whether to evaluate every resolution in `result.table`.
#'   Default is `TRUE`.
#' @param window.size Number of neighboring resolution ranks on each side of
#'   each evaluated resolution used for stability scoring. Default is `1`.
#' @param prefix Prefix used by Subtypist metadata columns. Default is
#'   `"Subtypist"`, matching columns such as `Subtypist_snn_res.0.5`.
#' @param suffix Optional suffix appended to Subtypist metadata columns.
#' @param cluster.columns Optional named character vector mapping resolution
#'   values to metadata column names. Use this when cluster columns do not follow
#'   the standard Subtypist naming pattern.
#' @param score.summary Function used to summarize cluster-level `Score` values
#'   when `selected.resolution = NULL`. Default is `mean`.
#' @param verbose Whether to print the selected and evaluated resolutions.
#'
#' @return A list with:
#' \describe{
#'   \item{result.table}{The original result table with added consensus columns.}
#'   \item{consensus.match.table}{A detailed table containing one best-matched
#'   neighboring cluster per evaluated-resolution cluster and neighboring
#'   resolution.}
#'   \item{selected.resolution}{The resolution recommended by summarized
#'   `Score`.}
#'   \item{evaluated.resolutions}{The resolutions assessed by the function.}
#' }
#' @export
Subtypist_consensus <- function(object = NULL,
                                result.table = NULL,
                                result.list = NULL,
                                selected.resolution = NULL,
                                evaluate.all = TRUE,
                                window.size = 1,
                                prefix = "Subtypist",
                                suffix = NULL,
                                cluster.columns = NULL,
                                score.summary = mean,
                                verbose = TRUE) {
  if (!is.null(result.list)) {
    object <- result.list[["Object"]]
    result.table <- result.list[["result.table"]]
  }

  if (is.null(object)) {
    stop("Please provide a Seurat object or a Subtypist result.list.")
  }
  if (is.null(result.table)) {
    stop("Please provide a Subtypist result.table or a Subtypist result.list.")
  }
  if (!is.data.frame(result.table)) {
    stop("'result.table' must be a data.frame or tibble.")
  }
  if (is.null(object@meta.data) || nrow(object@meta.data) == 0) {
    stop("'object@meta.data' is empty.")
  }
  if (window.size < 1) {
    stop("'window.size' must be at least 1.")
  }

  result.table <- as.data.frame(result.table)
  required.columns <- c("resolution", "merged_cluster", "Score")
  missing.columns <- setdiff(required.columns, colnames(result.table))
  if (length(missing.columns) > 0) {
    stop(
      "The result.table is missing required columns: ",
      paste(missing.columns, collapse = ", ")
    )
  }

  result.table$.subtypist_row_id <- seq_len(nrow(result.table))
  result.table$.resolution_key <- .resolution_key(result.table$resolution)
  result.table$.cluster_key <- .cluster_key(result.table$merged_cluster)

  resolutions <- .ordered_resolutions(result.table$resolution)
  if (length(resolutions) < 2) {
    stop("At least two resolutions are required for inter-resolution stability scoring.")
  }

  recommended.resolution <- .select_resolution_by_score(
    result.table = result.table,
    score.summary = score.summary
  )

  if (evaluate.all) {
    target.resolutions <- resolutions
  } else {
    if (is.null(selected.resolution)) {
      selected.resolution <- recommended.resolution
    }
    target.resolutions <- .match_resolution(selected.resolution, resolutions)
  }

  score.results <- lapply(
    target.resolutions,
    FUN = function(current.resolution) {
      .score_one_resolution(
        object = object,
        result.table = result.table,
        resolutions = resolutions,
        current.resolution = current.resolution,
        window.size = window.size,
        prefix = prefix,
        suffix = suffix,
        cluster.columns = cluster.columns
      )
    }
  )

  match.table <- as.data.frame(dplyr::bind_rows(
    lapply(score.results, `[[`, "match.table")
  ))
  consensus.table <- as.data.frame(dplyr::bind_rows(
    lapply(score.results, `[[`, "consensus.table")
  ))

  output <- .append_consensus_columns(
    result.table = result.table,
    consensus.table = consensus.table,
    selected.resolution = recommended.resolution
  )

  output$.subtypist_row_id <- NULL
  output$.resolution_key <- NULL
  output$.cluster_key <- NULL

  if (verbose) {
    message("Recommended resolution: ", .resolution_key(recommended.resolution))
    if (evaluate.all) {
      message("Evaluated resolutions: ", paste(.resolution_key(target.resolutions), collapse = ", "))
    } else {
      message("Evaluated resolution: ", paste(.resolution_key(target.resolutions), collapse = ", "))
    }
  }

  consensus.result <- list(
    result.table = as.data.frame(output),
    consensus.match.table = as.data.frame(match.table),
    selected.resolution = recommended.resolution,
    evaluated.resolutions = target.resolutions
  )

  return(consensus.result)
}

.ordered_resolutions <- function(resolution) {
  resolution.key <- .resolution_key(resolution)
  resolution.unique <- resolution[!duplicated(resolution.key)]
  resolution.numeric <- suppressWarnings(as.numeric(as.character(resolution.unique)))

  if (all(!is.na(resolution.numeric))) {
    resolution.unique <- resolution.unique[order(resolution.numeric)]
  }

  resolution.unique
}

.resolution_key <- function(resolution) {
  as.character(resolution)
}

.match_resolution <- function(resolution, available.resolutions) {
  resolution.key <- .resolution_key(resolution)
  available.key <- .resolution_key(available.resolutions)

  if (resolution.key %in% available.key) {
    return(available.resolutions[match(resolution.key, available.key)])
  }

  resolution.numeric <- suppressWarnings(as.numeric(resolution.key))
  available.numeric <- suppressWarnings(as.numeric(available.key))
  if (!is.na(resolution.numeric) && all(!is.na(available.numeric))) {
    distance <- abs(available.numeric - resolution.numeric)
    if (min(distance) < 1e-8) {
      return(available.resolutions[which.min(distance)])
    }
  }

  stop(
    "selected.resolution = ", resolution.key,
    " was not found in result.table$resolution."
  )
}

.select_resolution_by_score <- function(result.table, score.summary) {
  split.score <- split(result.table$Score, result.table$.resolution_key)
  summary.score <- vapply(
    split.score,
    FUN.VALUE = numeric(1),
    FUN = function(x) {
      x <- abs(x)
      x <- x[!is.na(x)]
      if (length(x) == 0) {
        return(NA_real_)
      }
      value <- tryCatch(
        score.summary(x, na.rm = TRUE),
        error = function(e) score.summary(x)
      )
      as.numeric(value)[1]
    }
  )

  if (all(is.na(summary.score))) {
    stop("Cannot select a resolution because all summarized Score values are NA.")
  }

  selected.key <- names(summary.score)[which.max(summary.score)]
  result.table$resolution[match(selected.key, result.table$.resolution_key)]
}

.score_one_resolution <- function(object,
                                  result.table,
                                  resolutions,
                                  current.resolution,
                                  window.size,
                                  prefix,
                                  suffix = NULL,
                                  cluster.columns = NULL) {
  current.resolution <- .match_resolution(current.resolution, resolutions)
  current.key <- .resolution_key(current.resolution)
  current.index <- which(.resolution_key(resolutions) == current.key)
  window.index <- seq(
    from = max(1, current.index - window.size),
    to = min(length(resolutions), current.index + window.size)
  )
  resolution.window <- resolutions[window.index]
  neighbor.resolutions <- resolution.window[
    .resolution_key(resolution.window) != current.key
  ]

  if (length(neighbor.resolutions) == 0) {
    return(list(match.table = data.frame(), consensus.table = data.frame()))
  }

  current.column <- .find_cluster_column(
    object = object,
    resolution = current.resolution,
    prefix = prefix,
    suffix = suffix,
    cluster.columns = cluster.columns
  )
  neighbor.columns <- stats::setNames(
    vapply(
      neighbor.resolutions,
      FUN.VALUE = character(1),
      FUN = function(resolution) {
        .find_cluster_column(
          object = object,
          resolution = resolution,
          prefix = prefix,
          suffix = suffix,
          cluster.columns = cluster.columns
        )
      }
    ),
    .resolution_key(neighbor.resolutions)
  )

  current.labels <- .metadata_labels(object, current.column)
  current.rows <- result.table[result.table$.resolution_key == current.key, ]
  n.cells <- sum(!is.na(current.labels))

  match.table <- .calculate_consensus_matches(
    selected.rows = current.rows,
    selected.labels = current.labels,
    selected.column = current.column,
    object = object,
    neighbor.resolutions = neighbor.resolutions,
    neighbor.columns = neighbor.columns,
    n.cells = n.cells
  )
  consensus.table <- .summarize_consensus_matches(match.table)

  list(match.table = match.table, consensus.table = consensus.table)
}

.find_cluster_column <- function(object,
                                 resolution,
                                 prefix,
                                 suffix = NULL,
                                 cluster.columns = NULL) {
  resolution.key <- .resolution_key(resolution)

  if (!is.null(cluster.columns)) {
    if (is.null(names(cluster.columns))) {
      stop("'cluster.columns' must be a named character vector.")
    }
    if (resolution.key %in% names(cluster.columns)) {
      column <- cluster.columns[[resolution.key]]
      if (!column %in% colnames(object@meta.data)) {
        stop("Metadata column '", column, "' was not found in object@meta.data.")
      }
      return(column)
    }
  }

  candidates <- .subtypist_cluster_column(prefix, resolution.key, suffix)
  candidates <- candidates[!is.na(candidates)]
  matched <- candidates[candidates %in% colnames(object@meta.data)]

  if (length(matched) > 0) {
    return(matched[[1]])
  }

  stop(
    "No cluster column was found for resolution ", resolution.key,
    ". Tried: ", paste(candidates, collapse = ", "),
    ". You can provide a named 'cluster.columns' vector."
  )
}

.metadata_labels <- function(object, column) {
  labels <- object@meta.data[[column]]
  names(labels) <- rownames(object@meta.data)
  .cluster_key(labels)
}

.cluster_key <- function(cluster) {
  key <- as.character(cluster)
  names(key) <- names(cluster)
  numeric.key <- suppressWarnings(as.numeric(key))
  is.whole <- !is.na(numeric.key) & abs(numeric.key - round(numeric.key)) < 1e-8
  key[is.whole] <- as.character(as.integer(round(numeric.key[is.whole])))
  key
}

.calculate_consensus_matches <- function(selected.rows,
                                         selected.labels,
                                         selected.column,
                                         object,
                                         neighbor.resolutions,
                                         neighbor.columns,
                                         n.cells) {
  match.rows <- list()
  row.index <- 1

  for (i in seq_len(nrow(selected.rows))) {
    selected.cluster <- selected.rows$.cluster_key[[i]]
    selected.cells <- names(selected.labels)[
      !is.na(selected.labels) & selected.labels == selected.cluster
    ]
    selected.size <- length(selected.cells)

    for (neighbor.resolution in neighbor.resolutions) {
      neighbor.key <- .resolution_key(neighbor.resolution)
      neighbor.column <- neighbor.columns[[neighbor.key]]
      neighbor.labels <- .metadata_labels(object, neighbor.column)
      neighbor.clusters <- sort(unique(neighbor.labels[!is.na(neighbor.labels)]))

      best <- .best_neighbor_match(
        selected.cells = selected.cells,
        neighbor.labels = neighbor.labels,
        neighbor.clusters = neighbor.clusters
      )

      match.rows[[row.index]] <- tibble::tibble(
        resolution = selected.rows$resolution[[i]],
        merged_cluster = selected.rows$merged_cluster[[i]],
        current_column = selected.column,
        neighbor_resolution = neighbor.resolution,
        neighbor_column = neighbor.column,
        matched_cluster = best$cluster,
        selected_cluster_size = selected.size,
        matched_cluster_size = best$matched.size,
        intersection_size = best$intersection,
        union_size = best$union,
        matched_jaccard = best$jaccard,
        matched_preservation = best$preservation,
        cluster_fraction = if (n.cells == 0) NA_real_ else selected.size / n.cells
      )
      row.index <- row.index + 1
    }
  }

  if (length(match.rows) == 0) {
    return(data.frame())
  }

  as.data.frame(dplyr::bind_rows(match.rows))
}

.best_neighbor_match <- function(selected.cells,
                                 neighbor.labels,
                                 neighbor.clusters) {
  if (length(selected.cells) == 0 || length(neighbor.clusters) == 0) {
    return(list(
      cluster = NA_character_,
      matched.size = NA_integer_,
      intersection = NA_integer_,
      union = NA_integer_,
      jaccard = NA_real_,
      preservation = NA_real_
    ))
  }

  candidate.rows <- lapply(neighbor.clusters, function(neighbor.cluster) {
    neighbor.cells <- names(neighbor.labels)[
      !is.na(neighbor.labels) & neighbor.labels == neighbor.cluster
    ]
    intersection.size <- length(intersect(selected.cells, neighbor.cells))
    union.size <- length(union(selected.cells, neighbor.cells))
    jaccard <- if (union.size == 0) NA_real_ else intersection.size / union.size
    preservation <- intersection.size / length(selected.cells)

    tibble::tibble(
      cluster = neighbor.cluster,
      matched.size = length(neighbor.cells),
      intersection = intersection.size,
      union = union.size,
      jaccard = jaccard,
      preservation = preservation
    )
  })

  candidates <- dplyr::bind_rows(candidate.rows)
  candidates <- candidates[order(
    -candidates$jaccard,
    -candidates$preservation,
    candidates$cluster
  ), ]
  best <- candidates[1, ]

  list(
    cluster = as.character(best$cluster[[1]]),
    matched.size = as.integer(best$matched.size[[1]]),
    intersection = as.integer(best$intersection[[1]]),
    union = as.integer(best$union[[1]]),
    jaccard = as.numeric(best$jaccard[[1]]),
    preservation = as.numeric(best$preservation[[1]])
  )
}

.summarize_consensus_matches <- function(match.table) {
  if (nrow(match.table) == 0) {
    return(data.frame())
  }

  match.table <- match.table
  match.table$.resolution_key <- .resolution_key(match.table$resolution)
  match.table$.cluster_key <- .cluster_key(match.table$merged_cluster)

  split.matches <- split(
    match.table,
    paste(match.table$.resolution_key, match.table$.cluster_key, sep = "||")
  )

  summary.rows <- lapply(split.matches, function(x) {
    jaccard <- x$matched_jaccard[!is.na(x$matched_jaccard)]
    preservation <- x$matched_preservation[!is.na(x$matched_preservation)]
    consensus.score <- if (length(jaccard) == 0) NA_real_ else mean(jaccard)
    preservation.score <- if (length(preservation) == 0) NA_real_ else mean(preservation)

    tibble::tibble(
      .resolution_key = x$.resolution_key[[1]],
      .cluster_key = x$.cluster_key[[1]],
      consensus_score = consensus.score,
      preservation_score = preservation.score,
      cluster_fraction = x$cluster_fraction[[1]]
    )
  })

  as.data.frame(dplyr::bind_rows(summary.rows))
}

.append_consensus_columns <- function(result.table,
                                      consensus.table,
                                      selected.resolution) {
  selected.key <- .resolution_key(selected.resolution)
  output <- result.table

  output$selected_resolution <- selected.resolution
  output$is_selected_resolution <- output$.resolution_key == selected.key

  metric.columns <- c(
    "consensus_score",
    "preservation_score",
    "cluster_fraction"
  )
  for (column in metric.columns) {
    output[[column]] <- NA
  }

  if (nrow(consensus.table) > 0) {
    for (i in seq_len(nrow(consensus.table))) {
      row.id <- which(
        output$.resolution_key == consensus.table$.resolution_key[[i]] &
          output$.cluster_key == consensus.table$.cluster_key[[i]]
      )
      if (length(row.id) == 1) {
        for (column in metric.columns) {
          output[[column]][row.id] <- consensus.table[[column]][[i]]
        }
      }
    }
  }

  output
}
