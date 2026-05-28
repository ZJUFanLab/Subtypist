#' Subtypist merge with post-hoc plateau diagnostics
#'
#' This diagnostic wrapper records merge-step score trajectories for each
#' resolution. By default it keeps the original Subtypist termination rule.
#' Users can optionally set \code{termination.mode = "plateau"} to stop a
#' resolution when the selected score metric enters a plateau, or
#' \code{termination.mode = "elbow"} to use a two-stage elbow plus optional
#' subsampling-stability rule.
#'
#' @param object A Seurat object.
#' @param min.resolution Minimum clustering resolution.
#' @param max.resolution Maximum clustering resolution.
#' @param by Resolution increment.
#' @param max.steps Maximum number of merge steps per resolution.
#' @param use.assay Assay used for marker detection.
#' @param cluster_assay Assay used for Seurat graph clustering in Seurat v4.
#' @param n.top Number of top marker genes retained per cluster.
#' @param min.pct.1 Minimum percent of cells expressing a marker in the cluster.
#' @param min.diff Minimum expression-percent difference for specificity scoring.
#' @param min.avg_log2FC Minimum average log2 fold change used by specificity scoring.
#' @param logfc.threshold Log fold-change threshold for marker testing.
#' @param prefix Prefix for Subtypist metadata columns.
#' @param accelerated Whether to use accelerated marker detection in Seurat v5.
#' @param algorithm Seurat clustering algorithm.
#' @param tua Threshold used by Subtypist to classify a cluster as specific.
#' @param verbose Whether to print progress.
#' @param regulation Marker regulation direction: up, down, or both.
#' @param top_k Number of marker genes used to summarize cluster scores.
#' @param termination.mode Termination rule for each resolution. \code{"original"}
#'   keeps the original Subtypist rule; \code{"plateau"} uses the plateau rule as
#'   an alternative stopping criterion; \code{"elbow"} uses a segmented-regression
#'   elbow rule with optional subsampling stability checks.
#' @param plateau.metric Score metric used when \code{termination.mode =
#'   "plateau"}.
#' @param plateau.change Whether plateau is called from absolute score changes
#'   or relative score changes.
#' @param plateau.eps Minimum score change treated as meaningful improvement.
#' @param plateau.patience Number of consecutive low-change steps required for a plateau.
#' @param plateau.min.steps Minimum number of merge steps before plateau calling.
#' @param elbow.metric Score metric used when \code{termination.mode = "elbow"}.
#' @param elbow.min.points Minimum number of score points required to fit an elbow.
#' @param elbow.min.steps Minimum merge step before elbow calling.
#' @param elbow.confirm.steps Number of additional merge steps observed after
#'   the estimated elbow before accepting it as a stopping point.
#' @param elbow.window Local step window around the estimated elbow used for
#'   optional subsampling stability checks.
#' @param elbow.stability.enabled Whether to run optional subsampling stability
#'   checks for the elbow mode.
#' @param elbow.stability.metric Similarity metric used for subsampling stability.
#' @param elbow.stability.threshold Minimum similarity required to accept the elbow.
#' @param elbow.subsample.fraction Fraction of cells sampled for each stability replicate.
#' @param elbow.subsample.reps Number of subsampling replicates.
#' @param elbow.subsample.seed Seed used for subsampling stability checks.
#'
#' @return A list with \code{Object}, \code{result.table}, \code{step.history},
#'   \code{plateau.table}, and \code{elbow.table}.
#' @export
Subtypist_merge_diagnostics <- function(object,
                                        min.resolution = 0.3,
                                        max.resolution = 1.5,
                                        by = 0.1,
                                        max.steps = 100,
                                        use.assay = "RNA",
                                        cluster_assay = "integrated",
                                        n.top = 300,
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
                                        termination.mode = c("original", "plateau", "elbow"),
                                        plateau.metric = c("total_score", "mean_score"),
                                        plateau.change = c("absolute", "relative"),
                                        plateau.eps = 0.01,
                                        plateau.patience = 2,
                                        plateau.min.steps = 2,
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
                                        elbow.subsample.seed = 1L) {
  regulation <- match.arg(regulation)
  termination.mode <- match.arg(termination.mode, choices = c("original", "plateau", "elbow"))
  plateau.metric <- match.arg(plateau.metric)
  plateau.change <- match.arg(plateau.change)
  elbow.metric <- match.arg(elbow.metric)
  elbow.stability.metric <- match.arg(elbow.stability.metric)
  elbow.min.points <- max(4L, as.integer(elbow.min.points))
  elbow.min.steps <- max(1L, as.integer(elbow.min.steps))
  elbow.confirm.steps <- max(0L, as.integer(elbow.confirm.steps))
  elbow.window <- max(0L, as.integer(elbow.window))
  elbow.stability.threshold <- max(0, min(1, as.numeric(elbow.stability.threshold)))
  elbow.subsample.fraction <- max(0.05, min(1, as.numeric(elbow.subsample.fraction)))
  elbow.subsample.reps <- max(1L, as.integer(elbow.subsample.reps))

  if (!(is(object, "Seurat") || is(object, "SeuratObject"))) {
    stop("Error: Please input a Seurat or SeuratObject!")
  }
  if (is.null(object@graphs)) {
    stop("Error: Please run FindNeighbors before Subtypist_merge_diagnostics().")
  }
  if (!use.assay %in% names(object@assays)) {
    stop("Error: The assay selected by use.assay is not in the object.")
  }

  seurat.major <- substr(packageVersion("Seurat"), 1, 1)
  if (seurat.major == "4") {
    if (!length(object@assays[[object@active.assay]]@scale.data)) {
      stop("Error: Please provide a Seurat object with scaled data.")
    }
  } else if (seurat.major == "5") {
    if (!length(object@assays[[object@active.assay]]$scale.data)) {
      stop("Error: Please provide a Seurat object with scaled data.")
    }
  } else {
    stop("Unsupported Seurat major version: ", seurat.major)
  }

  obj2 <- object
  if (seurat.major == "4" && use.assay == "SCT") {
    obj2 <- Seurat::PrepSCTFindMarkers(object = obj2, verbose = FALSE)
  }

  clustertmp <- rep(FALSE, 100)
  results <- data.frame()
  step.history <- data.frame()
  elbow.table <- data.frame()

  for (i.resolution in seq(min.resolution, max.resolution, by = by)) {
    if (seurat.major == "4") {
      Seurat::DefaultAssay(obj2) <- cluster_assay
    }

    obj2 <- Seurat::FindClusters(
      object = obj2,
      resolution = i.resolution,
      verbose = FALSE,
      algorithm = algorithm
    )

    column <- paste(obj2@active.assay, "_snn_res.", as.character(i.resolution), sep = "")
    if (algorithm == 4) {
      obj2@meta.data[[column]] <- as.numeric(obj2@meta.data[[column]]) - 1
      Seurat::Idents(obj2) <- obj2@meta.data[[column]]
      obj2@meta.data[[column]] <- factor(
        obj2@meta.data[[column]],
        levels = 0:(length(unique(obj2@meta.data[[column]])) - 1)
      )
    }

    Newcolumn <- paste(prefix, "snn_res.", as.character(i.resolution), sep = "")
    clusterNum <- length(unique(obj2@meta.data[[column]]))
    initial.cluster.num <- clusterNum

    if (clusterNum == 1) next
    cluster_size <- table(obj2@meta.data[[column]])
    if (any(cluster_size < 3)) {
      message(
        "Skipping resolution ", i.resolution,
        ": cluster(s) with fewer than 3 cells: ",
        paste(names(cluster_size)[cluster_size < 3], collapse = ", ")
      )
      next
    }
    if (clustertmp[clusterNum]) {
      last.resolution <- c(last.resolution, i.resolution)
    }
    last.resolution <- list()
    last.resolution <- c(last.resolution, i.resolution)

    if (verbose) {
      cat(
        "Now calculate resolution:", i.resolution,
        "\nThere are", clusterNum, "subclusters at this resolution.",
        "\nThe results will be saved in ->", Newcolumn, "\n"
      )
    }

    clustertmp[clusterNum] <- TRUE
    tmp <- rep(FALSE, clusterNum)
    steps <- 0
    stop.reason <- NA_character_
    elbow.result <- NULL
    state.history <- list()

    while (TRUE) {
      if (steps == 0) {
        obj2@meta.data[[Newcolumn]] <- as.numeric(levels(obj2@meta.data[[column]]))[obj2@meta.data[[column]]]
        Seurat::Idents(obj2) <- Newcolumn
        Seurat::DefaultAssay(obj2) <- use.assay

        marker.res <- .diagnostic_find_initial_markers(
          obj = obj2,
          column = column,
          use.assay = use.assay,
          n.top = n.top,
          min.pct.1 = min.pct.1,
          min.diff = min.diff,
          min.avg_log2FC = min.avg_log2FC,
          logfc.threshold = logfc.threshold,
          regulation = regulation,
          accelerated = accelerated,
          seurat.major = seurat.major
        )
        Allmarkers_top <- marker.res$Allmarkers_top
        resMarker <- marker.res$resMarker

        M <- .getInitWeightedJaccardMatrix(clusterNum, Allmarkers_top, regulation = regulation)
        for (cluster in 0:(clusterNum - 1)) {
          tmp[cluster + 1] <- .check_standard(resMarker[resMarker$cluster == cluster, ], tua = tua)
        }
        mergedNodes <- vector("list", clusterNum)
        for (i in seq_len(clusterNum)) mergedNodes[[i]] <- list(i)
      }

      if (steps == 0) {
        step.history <- rbind(
          step.history,
          .diagnostic_step_row(
            resolution = i.resolution,
            step = steps,
            clusterNum = clusterNum,
            initial.cluster.num = initial.cluster.num,
            tmp = tmp,
            resMarker = resMarker,
            top_k = top_k,
            regulation = regulation,
            merged.cluster.min = NA_integer_,
            merged.cluster.max = NA_integer_,
            merge.jaccard = NA_real_,
            stop.reason = NA_character_
          )
        )
        state.history[[as.character(steps)]] <- .diagnostic_save_state(
          object = obj2,
          Newcolumn = Newcolumn,
          resMarker = resMarker,
          mergedNodes = mergedNodes,
          clusterNum = clusterNum,
          tmp = tmp
        )
      }

      if (sum(tmp) == clusterNum) {
        stop.reason <- "all_clusters_specific"
        if (verbose) cat("GOOD! get the perfect results!\n")
        break
      } else if (steps == max.steps) {
        stop.reason <- "max_steps"
        if (verbose) cat("Reach the max of steps!\n")
        break
      }

      max.col <- apply(M, 2, max)
      if (termination.mode == "original" &&
          sum(tmp) != clusterNum && sum(tmp) != 0 && min(max.col[tmp]) >= max(max.col[!tmp])) {
        stop.reason <- "original_nonspecific_boundary"
        break
      }

      IndexRes <- .Find_max_below_threshold(M, 1e9)
      firstMax <- IndexRes[[1]]
      Index.min <- IndexRes[[2]]
      Index.max <- IndexRes[[3]]
      if (Index.min == 0 || Index.max == 0) {
        stop.reason <- "no_merge_candidate"
        break
      }

      merged.max <- max(unlist(purrr::map(mergedNodes, .f = length)))
      if (tmp[Index.min] == FALSE || tmp[Index.max] == FALSE) {
        if ((purrr::map(mergedNodes, .f = length)[[Index.min]] > merged.max && tmp[[Index.min]] == FALSE) ||
            (purrr::map(mergedNodes, .f = length)[[Index.max]] > merged.max && tmp[Index.min] == FALSE) ||
            merged.max > length(unique(obj2@meta.data[[column]])) / 2) {
          stop.reason <- "original_merge_limit"
          break
        }
      }

      cluster.min <- Index.min - 1
      cluster.max <- Index.max - 1

      obj2 <- .updateSeuratObj(
        obj = obj2,
        column = Newcolumn,
        combined.min = cluster.min,
        combined.max = cluster.max,
        clusterNum = clusterNum
      )
      Seurat::Idents(obj2) <- Newcolumn

      if (clusterNum == 2) {
        stop.reason <- "single_cluster_after_merge"
        if (verbose) {
          cat(paste0("Resolution ", i.resolution, " produced a single cluster; stopping this resolution."))
        }
        break
      }

      newMarkers_top <- .diagnostic_find_merged_markers(
        obj = obj2,
        column = column,
        Newcolumn = Newcolumn,
        cluster.min = cluster.min,
        use.assay = use.assay,
        n.top = n.top,
        min.pct.1 = min.pct.1,
        min.diff = min.diff,
        min.avg_log2FC = min.avg_log2FC,
        logfc.threshold = logfc.threshold,
        regulation = regulation,
        accelerated = accelerated,
        seurat.major = seurat.major
      )

      if (accelerated && seurat.major == "5") {
        Allmarkers_top <- newMarkers_top
        resMarker <- data.frame()
        for (cluster in 0:(clusterNum - 1)) {
          cluster_top_with_score <- getSpecificity_score(
            Allmarkers_top[Allmarkers_top$cluster == cluster, ],
            min.pct.1 = min.pct.1,
            min.diff = min.diff,
            min.avg_log2FC = min.avg_log2FC,
            regulation = regulation
          )
          resMarker <- rbind(resMarker, cluster_top_with_score)
        }
      } else {
        resMarker <- .updateMarkerlist(
          resMarker,
          combined.min = cluster.min,
          combined.max = cluster.max,
          newMarkers_top,
          clusterNum
        )
      }

      tmp <- .updateState(tmp, Index.min, Index.max, resMarker[resMarker$cluster == cluster.min, ])
      M <- .updateDistanceMatrix(
        M = M,
        Index.min = Index.min,
        Index.max = Index.max,
        resMarker = resMarker,
        clusterNum = clusterNum,
        operation = .getWeightedJaccard,
        regulation = regulation
      )
      mergedNodes <- .mergeSteps(mergedNodes, combined.min = Index.min, combined.max = Index.max)

      steps <- steps + 1
      clusterNum <- clusterNum - 1

      step.history <- rbind(
        step.history,
        .diagnostic_step_row(
          resolution = i.resolution,
          step = steps,
          clusterNum = clusterNum,
          initial.cluster.num = initial.cluster.num,
          tmp = tmp,
          resMarker = resMarker,
          top_k = top_k,
          regulation = regulation,
          merged.cluster.min = cluster.min,
          merged.cluster.max = cluster.max,
          merge.jaccard = firstMax,
          stop.reason = NA_character_
        )
      )
      state.history[[as.character(steps)]] <- .diagnostic_save_state(
        object = obj2,
        Newcolumn = Newcolumn,
        resMarker = resMarker,
        mergedNodes = mergedNodes,
        clusterNum = clusterNum,
        tmp = tmp
      )

      if (termination.mode == "plateau") {
        current.history <- step.history[step.history$resolution == i.resolution, ]
        if (.diagnostic_plateau_reached(
          step.history = current.history,
          metric = plateau.metric,
          change = plateau.change,
          eps = plateau.eps,
          patience = plateau.patience,
          min.steps = plateau.min.steps
        )) {
          stop.reason <- paste0("plateau_", plateau.metric)
          if (verbose) {
            cat(
              "Plateau termination reached at resolution", i.resolution,
              "step", steps,
              "using", plateau.metric, "\n"
            )
          }
          break
        }
      }

      if (termination.mode == "elbow") {
        current.history <- step.history[step.history$resolution == i.resolution, ]
        elbow.result <- .diagnostic_elbow_reached(
          step.history = current.history,
          metric = elbow.metric,
          min.points = elbow.min.points,
          min.steps = elbow.min.steps,
          confirm.steps = elbow.confirm.steps,
          window = elbow.window,
          stability.enabled = elbow.stability.enabled,
          stability.metric = elbow.stability.metric,
          stability.threshold = elbow.stability.threshold,
          subsample.fraction = elbow.subsample.fraction,
          subsample.reps = elbow.subsample.reps,
          subsample.seed = elbow.subsample.seed,
          object = obj2,
          state.history = state.history,
          resolution = i.resolution,
          column = column,
          Newcolumn = Newcolumn,
          use.assay = use.assay,
          cluster_assay = cluster_assay,
          n.top = n.top,
          min.pct.1 = min.pct.1,
          min.diff = min.diff,
          min.avg_log2FC = min.avg_log2FC,
          logfc.threshold = logfc.threshold,
          regulation = regulation,
          accelerated = accelerated,
          seurat.major = seurat.major,
          max.steps = max.steps,
          algorithm = algorithm,
          tua = tua,
          top_k = top_k
        )
        if (isTRUE(elbow.result$should.stop)) {
          stop.reason <- elbow.result$stop.reason
          if (verbose) {
            cat(
              "Elbow termination reached at resolution", i.resolution,
              "step", steps,
              "using", elbow.metric, "\n"
            )
          }
          break
        }
      }
    }

    if (termination.mode == "elbow" &&
        !is.null(elbow.result) &&
        isTRUE(elbow.result$should.stop) &&
        !is.null(elbow.result$selected.step)) {
      selected.state <- state.history[[as.character(elbow.result$selected.step)]]
      if (!is.null(selected.state)) {
        obj2@meta.data[[Newcolumn]] <- selected.state$labels[rownames(obj2@meta.data)]
        Seurat::Idents(obj2) <- Newcolumn
        resMarker <- selected.state$resMarker
        mergedNodes <- selected.state$mergedNodes
        clusterNum <- selected.state$clusterNum
        tmp <- selected.state$tmp
        step.history <- step.history[
          !(step.history$resolution == i.resolution & step.history$step > elbow.result$selected.step),
          ,
          drop = FALSE
        ]
      }
    }

    if (nrow(step.history) > 0) {
      idx <- which(step.history$resolution == i.resolution & step.history$step == max(step.history$step[step.history$resolution == i.resolution]))
      if (length(idx) > 0) step.history$stop_reason[idx[length(idx)]] <- stop.reason
    }

    # Keep the final result table aligned with the original Subtypist_merge()
    # implementation. Plateau-specific finite handling is applied only to
    # step.history diagnostics, not to the canonical result.table.
    clu <- .setClulterInf(resMarker, mergedNodes, i.resolution, regulation = regulation, top_k = top_k)
    results <- rbind(results, as.data.frame(clu))

    if (!is.null(elbow.result) && !is.null(elbow.result$elbow.step)) {
      elbow.table <- rbind(
        elbow.table,
        data.frame(
          resolution = i.resolution,
          elbow_metric = elbow.metric,
          elbow_step = elbow.result$elbow.step,
          elbow_score = elbow.result$elbow.score,
          elbow_evaluated_step = elbow.result$evaluated.step,
          elbow_stability = elbow.result$elbow.stability,
          elbow_selected_step = elbow.result$selected.step,
          elbow_selected_score = elbow.result$selected.score,
          elbow_selected_stability = elbow.result$selected.stability,
          elbow_stability_metric = elbow.stability.metric,
          elbow_stability_threshold = elbow.stability.threshold,
          elbow_subsample_fraction = elbow.subsample.fraction,
          elbow_subsample_reps = elbow.subsample.reps,
          elbow_subsample_enabled = elbow.stability.enabled,
          elbow_stop_reason = stop.reason,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  plateau.table <- Subtypist_call_plateau(
    step.history = step.history,
    eps = plateau.eps,
    patience = plateau.patience,
    min.steps = plateau.min.steps,
    change = plateau.change
  )

  reslist <- list(
    Object = obj2,
    result.table = results,
    step.history = step.history,
    plateau.table = plateau.table,
    elbow.table = elbow.table,
    termination.mode = termination.mode,
    plateau.metric = plateau.metric,
    plateau.change = plateau.change,
    elbow.metric = elbow.metric,
    elbow.confirm.steps = elbow.confirm.steps,
    elbow.window = elbow.window,
    elbow.stability.enabled = elbow.stability.enabled,
    elbow.stability.metric = elbow.stability.metric,
    elbow.stability.threshold = elbow.stability.threshold
  )
  return(reslist)
}

#' Call post-hoc plateau windows from merge-step trajectories
#'
#' @param step.history The \code{step.history} data frame returned by
#'   \code{Subtypist_merge_diagnostics()}.
#' @param eps Minimum absolute score change treated as meaningful.
#' @param patience Number of consecutive low-change steps required.
#' @param min.steps Minimum step number at which plateau can be called.
#' @param change Whether plateau is called from absolute score changes or
#'   relative score changes.
#'
#' @return A data frame with one row per resolution and metric.
#' @export
Subtypist_call_plateau <- function(step.history,
                                   eps = 0.01,
                                   patience = 2,
                                   min.steps = 2,
                                   change = c("absolute", "relative")) {
  change <- match.arg(change)
  if (is.null(step.history) || nrow(step.history) == 0) {
    return(data.frame())
  }

  out <- list()
  i <- 1
  for (res in sort(unique(step.history$resolution))) {
    one <- step.history[step.history$resolution == res, ]
    one <- one[order(one$step), ]
    for (metric in c("total_score", "mean_score")) {
      score <- one[[metric]]
      delta <- .diagnostic_score_delta(score, change = change)
      low.change <- abs(delta) < eps
      low.change[is.na(low.change)] <- FALSE

      plateau.step <- NA_integer_
      if (nrow(one) >= max(min.steps, patience + 1)) {
        for (row.index in seq_len(nrow(one))) {
          current.step <- one$step[row.index]
          if (current.step < min.steps || row.index < patience + 1) next
          recent <- low.change[(row.index - patience + 1):row.index]
          if (all(recent)) {
            plateau.step <- current.step
            break
          }
        }
      }

      final.row <- one[nrow(one), ]
      out[[i]] <- data.frame(
        resolution = res,
        metric = metric,
        plateau_step = plateau.step,
        original_stop_step = .diagnostic_get_scalar(final.row, "step", NA_real_),
        original_stop_reason = .diagnostic_get_scalar(final.row, "stop_reason", NA_character_),
        final_cluster_num = .diagnostic_get_scalar(final.row, "cluster_num", NA_real_),
        final_total_score = .diagnostic_get_scalar(final.row, "total_score", NA_real_),
        final_mean_score = .diagnostic_get_scalar(final.row, "mean_score", NA_real_),
        plateau_change = change,
        plateau_called = !is.na(plateau.step),
        stringsAsFactors = FALSE
      )
      i <- i + 1
    }
  }
  if (length(out) == 0) return(data.frame())
  do.call(rbind, out)
}

.diagnostic_plateau_reached <- function(step.history,
                                        metric = c("total_score", "mean_score"),
                                        change = c("absolute", "relative"),
                                        eps = 0.01,
                                        patience = 2,
                                        min.steps = 2) {
  metric <- match.arg(metric)
  change <- match.arg(change)
  if (is.null(step.history) || nrow(step.history) == 0) return(FALSE)
  if (!metric %in% colnames(step.history)) return(FALSE)

  one <- step.history[order(step.history$step), ]
  if (nrow(one) < max(min.steps, patience + 1)) return(FALSE)

  current.step <- one$step[nrow(one)]
  if (current.step < min.steps) return(FALSE)

  score <- one[[metric]]
  delta <- .diagnostic_score_delta(score, change = change)
  low.change <- abs(delta) < eps
  low.change[is.na(low.change)] <- FALSE

  recent.index <- (length(low.change) - patience + 1):length(low.change)
  all(low.change[recent.index])
}

.diagnostic_elbow_reached <- function(step.history,
                                      metric = c("mean_score", "total_score"),
                                      min.points = 4,
                                      min.steps = 2,
                                      confirm.steps = 1,
                                      window = 1,
                                      stability.enabled = FALSE,
                                      stability.metric = c("ARI", "NMI"),
                                      stability.threshold = 0.85,
                                      subsample.fraction = 0.8,
                                      subsample.reps = 10,
                                      subsample.seed = 1L,
                                      object,
                                      state.history,
                                      resolution,
                                      column,
                                      Newcolumn,
                                      use.assay,
                                      cluster_assay,
                                      n.top,
                                      min.pct.1,
                                      min.diff,
                                      min.avg_log2FC,
                                      logfc.threshold,
                                      regulation,
                                      accelerated,
                                      seurat.major,
                                      max.steps,
                                      algorithm,
                                      tua,
                                      top_k) {
  metric <- match.arg(metric)
  stability.metric <- match.arg(stability.metric)
  if (is.null(step.history) || nrow(step.history) == 0) {
    return(list(should.stop = FALSE, stop.reason = NA_character_))
  }
  one <- step.history[order(step.history$step), ]
  if (!metric %in% colnames(one)) {
    return(list(should.stop = FALSE, stop.reason = NA_character_))
  }
  if (nrow(one) < min.points) {
    return(list(should.stop = FALSE, stop.reason = NA_character_))
  }

  elbow.step <- .diagnostic_estimate_elbow(one[[metric]], one$step, min.steps = min.steps)
  if (is.na(elbow.step)) {
    return(list(should.stop = FALSE, stop.reason = NA_character_))
  }
  current.step <- max(one$step, na.rm = TRUE)
  if (current.step < elbow.step + confirm.steps) {
    return(list(
      should.stop = FALSE,
      elbow.step = elbow.step,
      elbow.score = .diagnostic_get_scalar(one[one$step == elbow.step, , drop = FALSE], metric, NA_real_),
      evaluated.step = current.step,
      elbow.stability = NA_real_,
      selected.step = elbow.step,
      selected.score = NA_real_,
      selected.stability = NA_real_,
      stop.reason = "elbow_waiting_for_confirmation"
    ))
  }
  elbow.row <- one[one$step == elbow.step, , drop = FALSE]
  if (nrow(elbow.row) == 0) {
    elbow.row <- one[nrow(one), , drop = FALSE]
  }

  elbow.stability <- NA_real_
  elbow.score <- .diagnostic_get_scalar(elbow.row, metric, NA_real_)
  selected.step <- elbow.step
  selected.score <- elbow.score
  selected.stability <- NA_real_
  stop.reason <- paste0("elbow_", metric)

  if (stability.enabled) {
    candidate.steps <- sort(unique(one$step[abs(one$step - elbow.step) <= window]))
    candidate.steps <- candidate.steps[candidate.steps <= current.step]
    candidate.steps <- candidate.steps[vapply(
      candidate.steps,
      function(step) !is.null(state.history[[as.character(step)]]),
      logical(1)
    )]
    if (length(candidate.steps) == 0) {
      candidate.steps <- elbow.step
    }
    candidate.steps <- candidate.steps[order(abs(candidate.steps - elbow.step), candidate.steps)]

    stability.rows <- list()
    for (candidate.step in candidate.steps) {
      candidate.state <- state.history[[as.character(candidate.step)]]
      baseline.labels <- if (!is.null(candidate.state)) candidate.state$labels else NULL
      stability.result <- .diagnostic_elbow_stability(
        object = object,
        resolution = resolution,
        cluster_assay = cluster_assay,
        use.assay = use.assay,
        column = column,
        Newcolumn = Newcolumn,
        selected.step = candidate.step,
        baseline.labels = baseline.labels,
        stability.metric = stability.metric,
        stability.threshold = stability.threshold,
        subsample.fraction = subsample.fraction,
        subsample.reps = subsample.reps,
        subsample.seed = subsample.seed + as.integer(candidate.step),
        n.top = n.top,
        min.pct.1 = min.pct.1,
        min.diff = min.diff,
        min.avg_log2FC = min.avg_log2FC,
        logfc.threshold = logfc.threshold,
        regulation = regulation,
        accelerated = accelerated,
        seurat.major = seurat.major,
        max.steps = candidate.step,
        algorithm = algorithm,
        tua = tua,
        top_k = top_k
      )
      candidate.row <- one[one$step == candidate.step, , drop = FALSE]
      stability.rows[[length(stability.rows) + 1]] <- data.frame(
        step = candidate.step,
        score = .diagnostic_get_scalar(candidate.row, metric, NA_real_),
        stability = stability.result$selected.stability,
        passed = isTRUE(stability.result$passed),
        stop.reason = stability.result$stop.reason,
        stringsAsFactors = FALSE
      )
    }

    stability.table <- do.call(rbind, stability.rows)
    elbow.match <- stability.table[stability.table$step == elbow.step, , drop = FALSE]
    if (nrow(elbow.match) > 0) {
      elbow.stability <- elbow.match$stability[1]
    }
    passed.table <- stability.table[stability.table$passed, , drop = FALSE]
    if (nrow(passed.table) > 0) {
      passed.table <- passed.table[order(abs(passed.table$step - elbow.step), passed.table$step), ]
      selected.step <- passed.table$step[1]
      selected.score <- passed.table$score[1]
      selected.stability <- passed.table$stability[1]
      stop.reason <- paste0("elbow_", metric, "_stable_", stability.metric)
    } else {
      selected.step <- elbow.step
      selected.score <- .diagnostic_get_scalar(elbow.row, metric, NA_real_)
      selected.stability <- elbow.stability
      stop.reason <- paste0("elbow_", metric, "_stability_limited")
    }
  }

  list(
    should.stop = TRUE,
    elbow.step = elbow.step,
    elbow.score = elbow.score,
    evaluated.step = current.step,
    elbow.stability = elbow.stability,
    selected.step = selected.step,
    selected.score = selected.score,
    selected.stability = selected.stability,
    stop.reason = stop.reason
  )
}

.diagnostic_estimate_elbow <- function(score, steps, min.steps = 2) {
  score <- as.numeric(score)
  steps <- as.numeric(steps)
  ok <- is.finite(score) & is.finite(steps)
  score <- score[ok]
  steps <- steps[ok]
  if (length(score) < max(min.steps + 2, 4)) {
    return(NA_real_)
  }

  best.step <- NA_real_
  best.sse <- Inf
  candidate.idx <- seq.int(2, length(score) - 1)
  candidate.idx <- candidate.idx[steps[candidate.idx] >= min.steps]
  if (length(candidate.idx) == 0) {
    return(NA_real_)
  }
  for (k in candidate.idx) {
    x <- steps
    hinge <- pmax(0, steps - steps[k])
    fit <- stats::lm(score ~ x + hinge)
    sse <- sum(stats::resid(fit)^2, na.rm = TRUE)
    if (is.finite(sse) && sse < best.sse) {
      best.sse <- sse
      best.step <- steps[k]
    }
  }
  best.step
}

.diagnostic_elbow_stability <- function(object,
                                        resolution,
                                        cluster_assay,
                                        use.assay,
                                        column,
                                        Newcolumn,
                                        selected.step,
                                        baseline.labels = NULL,
                                        stability.metric = c("ARI", "NMI"),
                                        stability.threshold = 0.85,
                                        subsample.fraction = 0.8,
                                        subsample.reps = 10,
                                        subsample.seed = 1L,
                                        n.top,
                                        min.pct.1,
                                        min.diff,
                                        min.avg_log2FC,
                                        logfc.threshold,
                                        regulation,
                                        accelerated,
                                        seurat.major,
                                        max.steps,
                                        algorithm,
                                        tua,
                                        top_k) {
  stability.metric <- match.arg(stability.metric)
  target.col <- if (stability.metric == "ARI") "ARI" else "NMI"
  all.cells <- rownames(object@meta.data)
  if (is.null(all.cells) || length(all.cells) == 0) {
    return(list(
      elbow.stability = NA_real_,
      selected.step = selected.step,
      selected.score = NA_real_,
      selected.stability = NA_real_,
      stop.reason = "elbow_no_cells",
      passed = FALSE
    ))
  }

  set.seed(as.integer(subsample.seed))
  reps <- max(1, as.integer(subsample.reps))
  frac <- max(min(as.numeric(subsample.fraction), 1), 0.1)
  subset.size <- max(2L, floor(length(all.cells) * frac))
  scores <- numeric(reps)
  scores[] <- NA_real_

  for (b in seq_len(reps)) {
    chosen <- sample(all.cells, size = subset.size, replace = FALSE)
    sub.obj <- subset(object, cells = chosen)
    sub.cells <- rownames(sub.obj@meta.data)
    baseline <- NULL
    if (!is.null(baseline.labels)) {
      baseline <- baseline.labels[sub.cells]
    }
    if (is.null(baseline)) baseline <- sub.obj@meta.data[[Newcolumn]]
    if (is.null(baseline)) baseline <- sub.obj@meta.data[[column]]
    if (is.null(baseline)) next
    names(baseline) <- sub.cells

    if (seurat.major == "4" && use.assay == "SCT") {
      sub.obj <- Seurat::PrepSCTFindMarkers(object = sub.obj, verbose = FALSE)
    }
    stability.prefix <- paste0("elbow_stability_", b, "_")
    stability.col <- paste(stability.prefix, "snn_res.", as.character(resolution), sep = "")
    result <- tryCatch(Subtypist_merge(
      object = sub.obj,
      min.resolution = resolution,
      max.resolution = resolution,
      by = 1,
      max.steps = max(1, as.integer(max.steps)),
      use.assay = use.assay,
      cluster_assay = cluster_assay,
      n.top = n.top,
      min.pct.1 = min.pct.1,
      min.diff = min.diff,
      min.avg_log2FC = min.avg_log2FC,
      logfc.threshold = logfc.threshold,
      prefix = stability.prefix,
      accelerated = accelerated,
      algorithm = algorithm,
      tua = tua,
      verbose = FALSE,
      regulation = regulation,
      top_k = top_k
    ), error = function(e) e)
    if (inherits(result, "error")) next
    pred <- result$Object@meta.data[[stability.col]]
    if (is.null(pred)) next
    pred.cells <- rownames(result$Object@meta.data)
    truth <- baseline[pred.cells]
    ok <- !is.na(truth) & !is.na(pred)
    if (sum(ok) == 0) next
    if (target.col == "ARI") {
      scores[b] <- .diagnostic_adjusted_rand_index(truth[ok], pred[ok])
    } else {
      scores[b] <- .diagnostic_normalized_mutual_info(truth[ok], pred[ok])
    }
  }

  mean.score <- mean(scores, na.rm = TRUE)
  if (is.na(mean.score)) {
    return(list(
      elbow.stability = NA_real_,
      selected.step = selected.step,
      selected.score = NA_real_,
      selected.stability = NA_real_,
      stop.reason = "elbow_stability_unavailable",
      passed = FALSE
    ))
  }

  if (mean.score < stability.threshold) {
    return(list(
      elbow.stability = mean.score,
      selected.step = selected.step,
      selected.score = NA_real_,
      selected.stability = mean.score,
      stop.reason = paste0("elbow_stability_below_", stability.metric),
      passed = FALSE
    ))
  }

  list(
    elbow.stability = mean.score,
    selected.step = selected.step,
    selected.score = NA_real_,
    selected.stability = mean.score,
    stop.reason = paste0("elbow_", stability.metric),
    passed = TRUE
  )
}

.diagnostic_save_state <- function(object,
                                   Newcolumn,
                                   resMarker,
                                   mergedNodes,
                                   clusterNum,
                                   tmp) {
  labels <- object@meta.data[[Newcolumn]]
  names(labels) <- rownames(object@meta.data)
  list(
    labels = labels,
    resMarker = resMarker,
    mergedNodes = mergedNodes,
    clusterNum = clusterNum,
    tmp = tmp
  )
}

.diagnostic_adjusted_rand_index <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) stop("Arguments must be vectors of the same length.")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  expected <- (a + b) * (a + c) / (a + b + c + d)
  max.index <- ((a + b) + (a + c)) / 2
  if (max.index == expected) return(NA_real_)
  (a - expected) / (max.index - expected)
}

.diagnostic_normalized_mutual_info <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  px <- rowSums(tab) / n
  py <- colSums(tab) / n
  pxy <- tab / n
  nz <- pxy > 0
  mi <- sum(pxy[nz] * log(pxy[nz] / (px[row(pxy)[nz]] * py[col(pxy)[nz]])))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 && hy == 0) return(1)
  if (hx == 0 || hy == 0) return(0)
  mi / sqrt(hx * hy)
}

.diagnostic_score_delta <- function(score, change = c("absolute", "relative")) {
  change <- match.arg(change)
  score <- as.numeric(score)
  if (change == "absolute") {
    return(c(NA_real_, diff(score)))
  }
  previous <- score[-length(score)]
  current <- score[-1]
  c(NA_real_, (current - previous) / (abs(previous) + .Machine$double.eps))
}

.diagnostic_get_scalar <- function(x, name, default) {
  if (!name %in% colnames(x)) return(default)
  value <- x[[name]][1]
  if (length(value) == 0) return(default)
  value
}

.diagnostic_step_row <- function(resolution,
                                 step,
                                 clusterNum,
                                 initial.cluster.num,
                                 tmp,
                                 resMarker,
                                 top_k,
                                 regulation,
                                 merged.cluster.min,
                                 merged.cluster.max,
                                 merge.jaccard,
                                 stop.reason) {
  score.summary <- .diagnostic_score_summary(resMarker, top_k = top_k, regulation = regulation)
  data.frame(
    resolution = resolution,
    step = step,
    initial_cluster_num = initial.cluster.num,
    cluster_num = clusterNum,
    specific_cluster_num = sum(tmp),
    nonspecific_cluster_num = clusterNum - sum(tmp),
    total_score = score.summary$total_score,
    mean_score = score.summary$mean_score,
    median_score = score.summary$median_score,
    min_score = score.summary$min_score,
    max_score = score.summary$max_score,
    merged_cluster_min = merged.cluster.min,
    merged_cluster_max = merged.cluster.max,
    merge_jaccard = merge.jaccard,
    stop_reason = stop.reason,
    stringsAsFactors = FALSE
  )
}

.diagnostic_score_summary <- function(resMarker, top_k, regulation = c("up", "down", "both")) {
  regulation <- match.arg(regulation)
  if (is.null(resMarker) || nrow(resMarker) == 0) {
    return(list(
      total_score = NA_real_,
      mean_score = NA_real_,
      median_score = NA_real_,
      min_score = NA_real_,
      max_score = NA_real_
    ))
  }

  resMarker$specificity_score_for_summary <- .diagnostic_numeric_score(
    resMarker$specificity_score
  )

  marker.ranked <- resMarker %>%
    dplyr::arrange(
      cluster,
      dplyr::desc(specificity_score_for_summary),
      dplyr::desc(abs(avg_log2FC))
    )

  top.score <- marker.ranked %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = specificity_score_for_summary, n = top_k, with_ties = FALSE) %>%
    dplyr::summarise(cluster_score = sum(abs(specificity_score_for_summary), na.rm = TRUE), .groups = "drop")

  scores <- top.score$cluster_score
  if (length(scores) == 0) {
    return(list(
      total_score = NA_real_,
      mean_score = NA_real_,
      median_score = NA_real_,
      min_score = NA_real_,
      max_score = NA_real_
    ))
  }
  list(
    total_score = sum(scores, na.rm = TRUE),
    mean_score = mean(scores, na.rm = TRUE),
    median_score = stats::median(scores, na.rm = TRUE),
    min_score = min(scores, na.rm = TRUE),
    max_score = max(scores, na.rm = TRUE)
  )
}

.diagnostic_numeric_score <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  x
}

.diagnostic_find_initial_markers <- function(obj,
                                             column,
                                             use.assay,
                                             n.top,
                                             min.pct.1,
                                             min.diff,
                                             min.avg_log2FC,
                                             logfc.threshold,
                                             regulation,
                                             accelerated,
                                             seurat.major) {
  only.pos <- (regulation == "up")
  idents.all <- sort(unique(Seurat::Idents(object = obj)))

  if (seurat.major == "5" && accelerated) {
    obj_join_layers <- SeuratObject::JoinLayers(obj)
    all.markers <- presto::wilcoxauc(X = obj_join_layers, group.by = column, verbose = FALSE)
    all.markers <- all.markers %>% dplyr::filter(pct_in >= 0.1)
    all.markers$pct_in <- all.markers$pct_in / 100
    all.markers$pct_out <- all.markers$pct_out / 100
    colnames(all.markers) <- c(
      "gene", "cluster", "avgExpr", "avg_log2FC", "statistic", "auc",
      "p_val", "p_val_adj", "pct.1", "pct.2"
    )
    all.markers <- all.markers[, c("cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    if (only.pos) all.markers <- subset(all.markers, avg_log2FC > 0)
  } else {
    obj_for_markers <- obj
    if (seurat.major == "5") {
      obj_for_markers <- SeuratObject::JoinLayers(obj)
    }
    all.markers <- data.frame()
    for (ident in idents.all) {
      if (seurat.major == "4") {
        i.markers <- Seurat::FindMarkers(
          obj_for_markers,
          ident.1 = ident,
          only.pos = only.pos,
          min.pct = 0.1,
          verbose = FALSE,
          logfc.threshold = 0.1
        )
      } else {
        i.markers <- Seurat::FindMarkers(
          obj_for_markers,
          ident.1 = ident,
          only.pos = only.pos,
          assay = use.assay,
          min.pct = 0.1,
          verbose = FALSE,
          logfc.threshold = logfc.threshold
        )
      }
      i.markers$cluster <- ident
      i.markers$gene <- rownames(i.markers)
      all.markers <- rbind(all.markers, i.markers)
    }
  }

  Allmarkers_top <- all.markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = abs(avg_log2FC), n = n.top)

  if (seurat.major == "4") {
    Allmarkers_top$cluster <- suppressWarnings(as.numeric(as.character(Allmarkers_top$cluster)))
    if (any(is.na(Allmarkers_top$cluster))) {
      Allmarkers_top$cluster <- as.numeric(levels(Allmarkers_top$cluster))[Allmarkers_top$cluster]
    }
  } else {
    Allmarkers_top$cluster <- as.numeric(as.character(Allmarkers_top$cluster))
  }

  resMarker <- data.frame()
  for (cluster in sort(unique(Allmarkers_top$cluster))) {
    cluster_top_with_score <- getSpecificity_score(
      Allmarkers_top[Allmarkers_top$cluster == cluster, ],
      min.pct.1 = min.pct.1,
      min.diff = min.diff,
      min.avg_log2FC = min.avg_log2FC,
      regulation = regulation
    )
    resMarker <- rbind(resMarker, cluster_top_with_score)
  }

  list(Allmarkers_top = Allmarkers_top, resMarker = resMarker)
}

.diagnostic_find_merged_markers <- function(obj,
                                            column,
                                            Newcolumn,
                                            cluster.min,
                                            use.assay,
                                            n.top,
                                            min.pct.1,
                                            min.diff,
                                            min.avg_log2FC,
                                            logfc.threshold,
                                            regulation,
                                            accelerated,
                                            seurat.major) {
  only.pos <- (regulation == "up")

  if (seurat.major == "5" && accelerated) {
    obj_join_layers <- SeuratObject::JoinLayers(obj)
    all.markers <- presto::wilcoxauc(X = obj_join_layers, group.by = Newcolumn, verbose = FALSE)
    all.markers <- all.markers %>% dplyr::filter(pct_in >= 0.1)
    all.markers$pct_in <- all.markers$pct_in / 100
    all.markers$pct_out <- all.markers$pct_out / 100
    colnames(all.markers) <- c(
      "gene", "cluster", "avgExpr", "avg_log2FC", "statistic", "auc",
      "p_val", "p_val_adj", "pct.1", "pct.2"
    )
    all.markers <- all.markers[, c("cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    if (only.pos) all.markers <- subset(all.markers, avg_log2FC > 0)
    all.markers <- all.markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(order_by = abs(avg_log2FC), n = n.top)
    all.markers$cluster <- as.numeric(all.markers$cluster)
    return(all.markers)
  }

  obj_for_markers <- obj
  if (seurat.major == "5") {
    obj_for_markers <- SeuratObject::JoinLayers(obj)
  }
  Seurat::Idents(obj_for_markers) <- Newcolumn
  if (seurat.major == "4") {
    newCluster_marker <- Seurat::FindMarkers(
      obj_for_markers,
      ident.1 = cluster.min,
      ident.2 = NULL,
      only.pos = only.pos,
      min.pct = 0,
      verbose = FALSE
    )
  } else {
    newCluster_marker <- Seurat::FindMarkers(
      obj_for_markers,
      ident.1 = cluster.min,
      ident.2 = NULL,
      only.pos = only.pos,
      assay = use.assay,
      min.pct = 0,
      verbose = FALSE,
      logfc.threshold = logfc.threshold
    )
  }
  newMarkers_top <- newCluster_marker %>%
    dplyr::slice_max(order_by = abs(avg_log2FC), n = n.top)
  newMarkers_top <- getSpecificity_score(
    newMarkers_top,
    min.pct.1 = min.pct.1,
    min.diff = min.diff,
    min.avg_log2FC = min.avg_log2FC,
    regulation = regulation
  )
  newMarkers_top$gene <- rownames(newMarkers_top)
  newMarkers_top
}
