#' @importFrom dplyr %>%
NULL

.subtypist_cluster_column <- function(prefix, resolution, suffix = NULL) {
  prefix <- sub("_+$", "", prefix)
  paste0(prefix, "_snn_res.", resolution, suffix)
}

#' Sort scoring results by resolution
#'
#' @param result.table Subtypist result table.
#' @param .f A function to summarize `Score` by `resolution`. Default is `mean`.
#'
#' @return A tibble with `resolution` and summarized `value`.
#' @export
sortScore <- function(result.table, .f = mean) {
  rank <- result.table %>%
    dplyr::group_by(resolution) %>%
    dplyr::summarise(value = .f(abs(Score)))
  return(rank)
}

#' Add Subtypist annotations to a Seurat object
#'
#' @param object A Seurat object.
#' @param resolution A vector of resolution values at which to assign annotations.
#' @param result.table A data frame containing Subtypist results with
#'   `resolution`, `merged_cluster`, and `phenotypic_molecules`.
#' @param result.list Result list from `Subtypist_merge()`.
#' @param prefix Prefix for Subtypist cluster columns.
#' @param suffix A character suffix appended to Subtypist cluster columns.
#' @param select_index A named integer vector specifying which phenotype to
#'   select for each `merged_cluster`.
#' @param meta.prefix Prefix used for new metadata annotation columns.
#' @param value.suffix Suffix appended to annotation values.
#'
#' @return A Seurat object with new metadata columns for Subtypist annotations.
#' @export
AddSubtypist <- function(object = NULL,
                         resolution = NULL,
                         result.table = NULL,
                         result.list = NULL,
                         prefix = "Subtypist",
                         suffix = NULL,
                         select_index = NULL,
                         meta.prefix = "phenotypic_molecules_",
                         value.suffix = NULL) {
  if (!is.null(result.list)) {
    object <- result.list[["Object"]]
    result.table <- result.list[["result.table"]]
  } else {
    if (is.null(object)) {
      stop("There is no Seurat object provided")
    }
    if (is.null(result.table)) {
      stop("Please provide the results")
    }
  }
  if (is.null(resolution)) {
    stop("Please provide the resolution at which annotations need to be added!")
  }

  result <- result.table[c("resolution", "merged_cluster", "phenotypic_molecules")]
  result$resolution <- as.character(result$resolution)
  result <- result[result$resolution %in% as.character(resolution), ]

  if (!is.null(select_index)) {
    if (length(resolution) != 1) {
      stop("'select_index' can only be used when a single resolution is specified.")
    }
    if (length(select_index) != nrow(result[result$resolution == as.character(resolution), ])) {
      stop("Length of 'select_index' must match the number of clusters at resolution ", resolution, ".")
    }
    result <- result %>%
      dplyr::group_by(merged_cluster) %>%
      dplyr::mutate(
        selected_index = select_index[as.character(merged_cluster)],
        phenotypic_molecules = purrr::map2_chr(
          phenotypic_molecules,
          selected_index,
          function(x, idx) {
            if (length(x) >= idx) x[[idx]] else NA_character_
          }
        )
      ) %>%
      dplyr::ungroup()
  } else if (is.list(result$phenotypic_molecules)) {
    result$phenotypic_molecules <- purrr::map_chr(
      result[["phenotypic_molecules"]],
      ~ paste(.x, collapse = " / ")
    )
  }

  if (!is.null(value.suffix)) {
    result$phenotypic_molecules <- paste0(result$phenotypic_molecules, value.suffix)
  }

  Addmeta <- lapply(
    X = resolution,
    FUN = function(x) {
      combined.column.name <- .subtypist_cluster_column(prefix, x, suffix)
      if (!combined.column.name %in% colnames(object@meta.data)) {
        stop(paste0("Column ", combined.column.name, " not found in meta.data"))
      }
      resmeta <- object@meta.data[combined.column.name]
      names(resmeta) <- "Selected_resolution_Column"
      resMarkersTable <- result[result$resolution == as.character(x), ]
      if (nrow(resMarkersTable) == 0) {
        warning(paste0("No annotations found for resolution ", x, " in result.table"))
      }
      resmeta <- dplyr::left_join(
        resmeta,
        resMarkersTable,
        by = c(Selected_resolution_Column = "merged_cluster")
      )
      resmeta <- resmeta[c("phenotypic_molecules")]
      colnames(resmeta) <- paste0(meta.prefix, x)
      return(resmeta)
    }
  )
  Addmeta_df <- do.call(cbind, Addmeta)
  object@meta.data <- cbind(object@meta.data, Addmeta_df)
  return(object)
}

#' Save Subtypist annotation results to disk
#'
#' @param result.table A data frame or tibble containing Subtypist results.
#' @param path Directory path to save the file.
#' @param name File name ending in `.csv`, `.tsv`, or `.xlsx`.
#'
#' @return The processed results invisibly.
#' @export
saveResults <- function(result.table, path, name) {
  if (!is.data.frame(result.table)) {
    stop("'results' must be a data.frame or tibble.")
  }
  if (missing(path) || !nzchar(path)) {
    stop("'path' must be a non-empty string.")
  }
  if (!dir.exists(path)) {
    message("Path does not exist. Creating directory: ", path)
    dir.create(path, recursive = TRUE)
  }
  if (missing(name) || !nzchar(name)) {
    stop("'name' must be a non-empty string.")
  }

  if ("phenotypic_molecules" %in% colnames(result.table) && is.list(result.table$phenotypic_molecules)) {
    result.table$phenotypic_molecules <- purrr::map_chr(
      result.table$phenotypic_molecules,
      ~ paste(.x, collapse = " / ")
    )
  }

  if ("initial_cluster" %in% colnames(result.table) && is.list(result.table$initial_cluster)) {
    result.table$initial_cluster <- purrr::map_chr(
      result.table$initial_cluster,
      ~ paste(.x, collapse = ", ")
    )
  }

  file_path <- file.path(path, name)
  if (grepl("\\.csv$", name, ignore.case = TRUE)) {
    write.csv(result.table, file = file_path, row.names = FALSE)
  } else if (grepl("\\.tsv$", name, ignore.case = TRUE)) {
    readr::write_tsv(result.table, file = file_path)
  } else if (grepl("\\.xlsx$", name, ignore.case = TRUE)) {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Package 'writexl' is required to write .xlsx files. Please install it.")
    }
    writexl::write_xlsx(result.table, path = file_path)
  } else {
    stop("Unsupported file extension. Please use '.csv', '.tsv', or '.xlsx'.")
  }

  message("Results saved to: ", file_path)
  return(invisible(result.table))
}
