#' @importFrom dplyr %>%
NULL
#' Title
#'
#' @param object a Seurat object
#' @param min.resolution the minimum value of resolution
#' @param max.resolution the maximum value of resolution
#' @param by increment of the sequence.
#' @param max.steps
#' @param use.assay Name of assay to use
#' @param cluster.assay Name of the assay in the Seurat object to use for clustering
#' @param n.top the number of top genes
#' @param min.pct.1 the minimum value of pct.1
#' @param min.diff the minimum difference between pct.1 and pct.2
#' @param logfc.threshold the threshold of the log2 Flod Change
#' @param prefix String prefix for naming metadata columns. Each column name will be constructed as "<prefix>.snn_res.<resolution>"
#'
#' @return list of objects and merged results
#' @export
#'
#' @examples Subtypist_merge(object=Seu,min.resolution=0.1,max.resolution=0.2,by=0.1,use.assay="RNA",cluster_assay = "RNA")
Subtypist_merge <- function(object,
                             min.resolution=0.3,
                             max.resolution=1.5,
                             by=0.1,
                             max.steps = 100,
                             use.assay="RNA",
                             cluster_assay = "integrated",
                             n.top=500,
                             min.pct.1=0.1,
                             min.diff=0.1,
                             min.avg_log2FC=0.5,
                             logfc.threshold = 0.1,
                             prefix = 'Subtypist')
{
  if(!(is(object, "Seurat") || is(object, "SeuratObject"))){
    stop("Error: Please input a Seurat or SeuratObject!")
  }
  if(is.null(object@graphs)) {
    stop("Error: Before running FindClusters, you need to compute the neighborhood graph using the FindNeighbors function. FindNeighbors builds the neighborhood graph, such as the k-nearest neighbors (KNN) graph, which is essential for cluster identification. Please run FindNeighbors with your data first to create the required neighborhood graph.")
  }
  if(!length(object@assays[[object@active.assay]]@scale.data)){
    stop("Error: Please provide a Seurat object with scaled data. Make sure you have already performed data scaling before using this function.")
  }
  if(!use.assay %in% names(object@assays)){
    stop("Error: The assay you selected is not in the Object! ")
  }
  INF = 1e9
  obj2 <- object
  obj2@active.assay <- use.assay

  if(use.assay == "SCT"){
    obj2 <- Seurat::PrepSCTFindMarkers(object=obj2,verbose=FALSE)
  }
  # Global parameters for each function call
  clustertmp <- rep(FALSE,100)
  results <- tibble::tibble()
  mergedNodeslist <- list()
  # Traverse different resolutions
  resolution.list <- c()
  for(i.resolution in seq(min.resolution,max.resolution,by=by)){
    Seurat::DefaultAssay(obj2) <- cluster_assay
    obj2 <- Seurat::FindClusters(object = obj2, resolution = i.resolution,verbose=FALSE)
    column <- paste(obj2@active.assay,"_snn_res.",as.character(i.resolution),sep="")
    Newcolumn <- paste(prefix,"snn_res.",as.character(i.resolution),sep="")
    clusterNum <- length(unique(obj2@meta.data[[column]]))
    if(clustertmp[clusterNum]){
      # results[unlist(map(results$resolution, ~identical(.,last.resolution))),]$resolution <- map(results[unlist(map(results$resolution, ~identical(.,last.resolution))),]$resolution,~c(.,i.resolution))
      # clu$resolution <- i.resolution
      # results <- rbind(results,clu)
      last.resolution = c(last.resolution,i.resolution)
      next # The number of clusters corresponding to this resolution has already appeared, omitting this merging process
    }
    last.resolution = list()
    last.resolution = c(last.resolution,i.resolution)
    resolution.list <- c(resolution.list,i.resolution)

    cat("Now caculate the resolution: ",i.resolution,"\nThere are",clusterNum,"subclusters at this resolution.","\nThe results will be saved in ->",Newcolumn,"\n")
    clustertmp[clusterNum] <- TRUE
    tmp <- rep(FALSE,clusterNum)
    steps = 0
    init.false <- length(tmp) - sum(tmp)
    while(TRUE){

      if(steps == 0){
        # obj2[[Newcolumn]] <- obj2[[column]]
        obj2@meta.data[[Newcolumn]] <- as.numeric(levels(obj2@meta.data[[column]]))[obj2@meta.data[[column]]]
        Seurat::Idents(obj2) <- Newcolumn
        Seurat::DefaultAssay(obj2) <- use.assay
        idents.all <- sort(x = unique(x = Seurat::Idents(object = obj2)))
        markers.list <- list()
        all.markers.RNA <- tibble::tibble()
        for(i.ident in 1:length(idents.all)){
          i.markers <- Seurat::FindMarkers(obj2,ident.1=idents.all[i.ident],only.pos=T,min.pct=min.pct.1,assay=use.assay,verbose = FALSE,logfc.threshold=logfc.threshold) # other parameter logfc.threshold
          i.markers$cluster <- idents.all[i.ident]
          i.markers$gene <- rownames(i.markers)
          all.markers.RNA <- rbind(all.markers.RNA,i.markers)
        }
        # all.markers.RNA$gene <- rownames(all.markers.RNA)
        Allmarkers_top <- all.markers.RNA %>% dplyr::group_by(cluster) %>%  dplyr::top_n(n=n.top,wt=avg_log2FC)
        Allmarkers_top$cluster <- as.numeric(levels(Allmarkers_top$cluster))[Allmarkers_top$cluster]
        # Distance Matrix
        M <- getInitWeightedJaccardMatrix(clusterNum,Allmarkers_top)
        # markers with specificial score and tmp
        resMarker <- tibble::tibble()
        for(cluster in 0:(clusterNum-1)){
          cluster_top_with_score <- getSpecificity_score(Allmarkers_top[Allmarkers_top$cluster==cluster,],min.pct.1= min.pct.1,min.diff = 0.4)# ,min.gap =
          resMarker <- rbind(resMarker,cluster_top_with_score)
          tmp[cluster + 1] <- check_standard(cluster_top_with_score)
        }
        selected.merge.times <- 0
        mergedNodes <- list()
        for(i in 1:clusterNum){
          mergedNodes[[i]] <- list(i)
        }
      }


      if(sum(tmp)==clusterNum){
        cat("GOOD! get the prefect results!\n")
        break
      }else if(steps==max.steps){
        cat("Reach the max of steps!\n")
        break}
      max.col <- apply(M,2,max)
      # Find the two clusters that need to be merged
      # Control that the maximum number of merges for clusters containing False is always less than the maximum number of merges for the entire data


      if(sum(tmp) != clusterNum & sum(tmp)!=0 & min(max.col[tmp]) >= max(max.col[!tmp])){
        break
      }

      # Find the two clusters that need to be merged
      # Control that the maximum number of merges for clusters containing False is always less than the maximum number of merges for the entire data
      IndexRes <- Find_max_below_threshold(M,1e9)
      firstMax <- IndexRes[[1]]
      Index.min <- IndexRes[[2]]
      Index.max <- IndexRes[[3]]
      # find the Max of the merging times
      merged.max = max(unlist(purrr::map(mergedNodes,.f=length)))
      tmpMax <- firstMax

      if(tmp[Index.min]==FALSE | tmp[Index.max]==FALSE){
        if((purrr::map(mergedNodes,.f=length)[[Index.min]] > merged.max & tmp[[Index.min]]==FALSE)
           |(purrr::map(mergedNodes,.f=length)[[Index.max]] > merged.max & tmp[Index.min]==FALSE)
           | merged.max > length(unique(obj2@meta.data[[column]]))/2){break}
      }

      cluster.min <- Index.min - 1
      cluster.max <- Index.max - 1
      # It's time to merge
      # Update the SeuratObject column
      obj2 <- updateSeuratObj(obj=obj2,column = Newcolumn,combined.min=cluster.min,combined.max=cluster.max,clusterNum=clusterNum)
      # Update thr MarkersList
      Seurat::Idents(obj2) <- Newcolumn
      newCluster_marker <- Seurat::FindMarkers(obj2, ident.1 = cluster.min, ident.2 = NULL,only.pos=T,min.pct=0,verbose = FALSE,logfc.threshold=0)
      newMarkers_top <- newCluster_marker %>% dplyr::top_n(n=n.top,wt=avg_log2FC)
      newMarkers_top <- getSpecificity_score(newMarkers_top,min.pct.1= min.pct.1,min.diff=min.diff)

      newMarkers_top$gene <- rownames(newMarkers_top)
      resMarker <- updateMarkerlist(resMarker,combined.min=cluster.min,combined.max=cluster.max,newMarkers_top,clusterNum) # newMarkers_top with cluster
      # Update the condition list
      tmp <- updateState(tmp,Index.min,Index.max,resMarker[resMarker$cluster==cluster.min,])
      # Update DistanceMatrix
      M <- updateDistanceMatrix(M=M,Index.min=Index.min,Index.max=Index.max,resMarker=resMarker,clusterNum=clusterNum,operation=getWeightedJaccard)
      # Update the merged list
      mergedNodes <- mergeSteps(mergedNodes,combined.min=Index.min,combined.max=Index.max)
      # do it for everysteps
      steps = steps + 1
      clusterNum = clusterNum - 1
    }

    # printSteps(mergedNodes,column)
    clu <- setClulterInf(resMarker,mergedNodes,i.resolution)
    results <- rbind(results,clu)
  }
  reslist <- list(obj2,results)
  names(reslist) <- c('Object','result.table')
  return(reslist)
}

#' Title Sort scoring results by resolution
#'
#' @param result.table results  data frame.
#' @param .f A function to summarize `Score` by `resolution` (default: `mean`).
#'
#' @return a sorted tibble with `resolution` and summarized `value`.
#' @export
#'
#' @examples sortScore(res[[2]],mean)
sortScore <- function(result.table,.f=mean){
  rank <- result.table %>% dplyr::group_by(resolution) %>%
    dplyr::summarise(value = .f(Score))
  return(rank)
}

#' Title Add Subtypist Annotations to a Seurat Object
#'
#' @param object A Seurat object.
#' @param resolution A vector of resolution values at which to assign annotations.
#' @param result.table A data frame containing Subtypist results with `resolution`, `merge_cluster`, and `molecular_phenotype`.
#' @param result.list result list from Subtypist_merge
#' @param prefix Prefix for new metadata columns. Default is `"Subtypist"`.
#' @param suffix A character suffix to append to the new annotation column names
#' @param select_index A named integer vector specifying which phenotype to select for each class (merge_cluster). Names should match cluster IDs in the `merge_cluster` column of the result table.
#'
#' @return A Seurat object with new metadata columns for Subtypist annotations.
#' @export
#'
#' @examples
#' library(Seurat)
#' res.table = tibble::tibble(resolution =c(0.1,0.1,0.1),merge_cluster = c(0,1,2),initial_cluster = c(0,1,2),molecular_phenotype = c('G1','G2','G3'),Score = c(1,2,4))
#' object <- AddSubtypist(object =res[[1]],result.table=res.table,resolution = c(0.1))
AddSubtypist <- function(object=NULL,resolution=NULL,result.table=NULL,result.list=NULL,prefix='Subtypist',suffix = NULL,select_index = NULL,...){
  if(!is.null(result.list)){
    object <- result.list[['Object']]
    result.table <- result.list[['result.table']]

  }else{
    if(is.null(object)){
      stop("There is no Seurat object provided")
    }
    if(is.null(result.table)){
      stop("Please provide the results")
    }
  }
  if(is.null(resolution)){
    stop("Please provide the resolution at which annotations need to be added!")
  }
  result <- result.table[c('resolution','merge_cluster','molecular_phenotype')]
  result$resolution = as.character(result$resolution)
  # If select_index is provided, select the specified molecular phenotype.
  if (!is.null(select_index)) { # Exclude quantity mismatches
    if (length(resolution) != 1) {
      stop("'select_index' can only be used when a single resolution is specified.")
    }
    if (length(select_index) != nrow(result[result$resolution == as.character(resolution), ])) {
      stop("Length of 'select_index' must match the number of clusters at resolution ", x, ".")
    }
    # select_index should be a named vector where names correspond to merge_cluster identifiers.
    result <- result %>%
      dplyr::group_by(merge_cluster) %>%
      dplyr::mutate(
        selected_index = select_index[as.character(merge_cluster)],
        molecular_phenotype = purrr::map2_chr(molecular_phenotype, selected_index, function(x, idx) {
          if (length(x) >= idx) x[[idx]] else NA_character_
        })
      ) %>% dplyr::ungroup()
  } else {
    # By default, concatenate all molecular phenotypes.
    if (is.list(result$molecular_phenotype)) {
      result$molecular_phenotype <- purrr::map_chr(result$molecular_phenotype, ~paste(.x, collapse = " / "))
    }
  }

  Addmeta <- lapply(
    X = resolution,
    FUN = function(x){
      combined_column_name <- paste0(prefix,"snn_res.",x,suffix)
      if (!combined_column_name %in% colnames(object@meta.data)) {
        stop(paste0("Column ", combined_column_name, " not found in meta.data"))
      }
      resmeta <- object@meta.data[combined_column_name]
      names(resmeta) <- 'Selected_resolution_Column'
      resMarkersTable <- result[result$resolution == as.character(x), ]
      if (nrow(resMarkersTable) == 0) {
        warning(paste0("No annotations found for resolution ", x, " in result.table"))
      }
      resmeta <- dplyr::left_join(resmeta,resMarkersTable,by=c(Selected_resolution_Column='merge_cluster'))
      resmeta <- resmeta[c("molecular_phenotype")]
      colnames(resmeta) <- paste0("Molecular_phenotype_", x)
      return(resmeta)
    }
  )
  Addmeta_df <- do.call(cbind, Addmeta)
  object@meta.data <- cbind(object@meta.data, Addmeta_df)
  return(object)
}


#' Title
#'
#' @param result.table results A data.frame or tibble containing scSubAnn results.
#' @param path Directory path to save the file.
#' @param name File name, should end with ".csv", ".tsv", or ".xlsx".
#'
#' @return The processed results (invisibly).
#' @export
#'
#' @examples saveResults(res[[2]], path = "output", name = "subtype.csv")
#' saveResults(res[[2]], path = "output", name = "subtype.xlsx")
saveResults <- function(result.table,path,name)
{
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
    # Format columns if present and are lists
  if ("molecular_phenotype" %in% colnames(result.table) && is.list(result.table$molecular_phenotype)) {
    result.table$molecular_phenotype <- purrr::map_chr(
      result.table$molecular_phenotype, ~ paste(.x, collapse = " / ")
    )
  }

  if ("initial_cluster" %in% colnames(result.table) && is.list(result.table$initial_cluster)) {
    result.table$initial_cluster <- purrr::map_chr(
      result.table$initial_cluster, ~ paste(.x, collapse = ", ")
    )
  }
  file_path <- file.path(path, name)
  # Determine format and write
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
