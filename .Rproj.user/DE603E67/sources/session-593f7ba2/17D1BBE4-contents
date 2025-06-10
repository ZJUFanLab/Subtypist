#' Title
#'
#' @param object A Seurat object.
#' @param result.table A data frame containing annotation results.
#' @param results.list Subtypist result list
#' @param reduction Name of the dimensionality reduction to use.
#' @param resolution A vector of resolution values for which to retrieve annotation columns.
#' @param show Type of annotation to visualize. One of: 'Result_clusterIndex', 'initial.clusterIndex', 'molecular_phenotype'.
#' @param label Whether to label clusters on the plot.
#' @param label.size Size of the cluster labels.
#' @param label.color Color of the cluster labels.
#' @param label.box Whether to draw a box around cluster labels.
#' @param repel Whether to repel text labels.
#' @param cells.highlight Cells to highlight on the plot.
#' @param cols.highlight Color(s) for highlighted cells.
#' @param sizes.highlight Size(s) for highlighted cells.
#' @param na.value Color for NA values.
#' @param ncol Number of columns if combining multiple plots.
#' @param combine Whether to combine plots.
#' @param raster Whether to use rasterized plotting.
#' @param raster.dpi DPI for rasterized output.
#' @param prefix Prefix for new metadata columns. Default is `"Subtypist"`.
#' @param suffix A character suffix to append to the new annotation column names
#' @param select_index A named integer vector specifying which phenotype to select for each class (merge_cluster). Names should match cluster IDs in the `merge_cluster` column of the result table.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples Subtypist_Dimplot(object = res[[1]], result.table = res[[2]], resolution = c(0.1), show = "molecular_phenotype")
Subtypist_Dimplot <- function(object=NULL,result.table=NULL,results.list=NULL,
                                reduction=NULL,resolution=NULL,show=NULL,
                                label = FALSE,label.size = 4,label.color = 'black',
                                label.box = FALSE,repel = FALSE,cells.highlight = NULL,
                                cols.highlight = '#DE2D26',sizes.highlight = 1,na.value = 'grey50',
                                ncol = NULL,combine = TRUE,raster = NULL,raster.dpi = c(512, 512),prefix='Subtypist',suffix=NULL,select_index= NULL){
  if(!is.null(results.list)){
    object <- results.list[['Object']]
    result.table <- results.list[['result.table']]
  }else{
    if(is.null(object)){
      stop("There is no Seurat object provided")
    }
  }
  if(is.null(resolution)){
    stop("Please provide the resolution at which annotations need to be added!")
  }
  if (!is.null(result.table)) {
    keep_res <- c()
    for (index in seq_along(resolution)) {
      res_x <- resolution[index]
      if (nrow(result.table[result.table$resolution == res_x, ]) == 0) {
        warning(paste("The resolution", res_x, "has no result and will be skipped."))
      } else {
        keep_res <- c(keep_res, res_x)
      }
    }
    resolution <- keep_res
  }
  if(is.null(result.table)){
    if(show == "Result_clusterIndex"){column.name <- unlist(lapply(X=resolution,FUN=paste0("Combined_snn_res.",x)))}
    if(show =="initial.clusterIndex"){column.name <- unlist(lapply(X=resolution,FUN=paste0("initial_cluster_",x)))}
    if(show =="molecular_phenotype"){column.name <- unlist(lapply(X=resolution,FUN=paste0(prefix,'snn_res.',x)))}
    for(index in 1:length(column.name)){
      if(!column.name[index] %in% colnames(object@meta.data)){
        cat("The resoluton of",x,"has no result.")
        resolution <- resolution[-index]
      }
    }
  }

  Combined_Comlumn.name <- unlist(lapply(X=resolution,FUN=function(x){paste0(prefix,'snn_res.',x)}))
  for(index in 1:length(Combined_Comlumn.name)){
    if (!Combined_Comlumn.name[index] %in% colnames(object@meta.data)){
      cat("There is no result of the resolution of",Combined_Comlumn.name[index],"in this object!")
      Combined_Comlumn.name <- Combined_Comlumn.name[-index]
    }
  }
  plot <- NULL
  if(show == "Result_clusterIndex"){
    plot <- DimPlot(object=object,group.by=Combined_Comlumn.name,ncol=3,label = TRUE,label.size = 4,label.color = 'black',label.box = FALSE,repel = FALSE,cells.highlight = NULL,cols.highlight = '#DE2D26',sizes.highlight = 1,na.value = 'grey50',ncol = NULL,combine = TRUE,raster = NULL,raster.dpi = c(512, 512))
  }

  if(show =="initial.clusterIndex"){
    if(is.list(result.table$initial_cluster)){
      result.table$initial_cluster <- purrr::map_chr(result.table$initial_cluster,.f=~paste(.x,collapse=", "))
    }
    result.table <- result.table[c('resolution','merge_cluster','initial_cluster')]
    Addmeta <- apply(
      X = resolution,
      FUN = function(x){
        Combined_Comlumn.name <- paste0("Combined_snn_res.",x)
        resmeta <- object@meta.data[Combined_Comlumn.name]
        names(resmeta) <- 'Selected_resolution_Column'
        resmeta$Selected_resolution_Column <- as.numeric(resmeta$Selected_resolution_Column)
        resMarkersTable <- result.table[result.table$resolution==x,]
        resmeta <- left_join(resmeta,resMarkersTable,by=c(Selected_resolution_Column='merge_cluster'))
        names(resmeta[3]) <- paste0("initial_cluster_",x)
        resmeta <- resmeta[3]
        return(resmeta)
      }
    )
    object@meta.data <- cbind(object@meta.data,Addmeta)
    plot <- DimPlot(object=object,group.by=colnames(Addmeta),ncol=ncol,label = TRUE,label.size = 4,label.color = 'black',label.box = FALSE,repel = FALSE,cells.highlight = NULL,cols.highlight = '#DE2D26',sizes.highlight = 1,na.value = 'grey50',combine = TRUE,raster = NULL,raster.dpi = c(512, 512))
  }
  if(show =="molecular_phenotype"){
    if(is.list(result.table$molecular_phenotype)){
      result.table$molecular_phenotype <- purrr::map_chr(result.table$molecular_phenotype,.f=~paste(.x, collapse = "/"))
    }
    molecular_phenotype.name <- lapply(X=resolution,FUN=function(x){return(paste0("Molecular_phenotype_",x))})
    Addmeta <- lapply(X=resolution,FUN=function(x){
      resmeta <- object@meta.data
      if(!paste0("Molecular_phenotype_",x) %in% colnames(object) & !is.null(result.table) & x %in% result.table$resolution){
        object <- AddSubtypist(object = object,resolution = x,result.table = result.table,prefix=prefix,suffix = suffix,select_index = select_index)
        resmeta <- cbind(resmeta,object[[paste0("Molecular_phenotype_",x)]])
      }
      return(resmeta)
    })
    object@meta.data <- cbind(object@meta.data,Addmeta)
    plot <- Seurat::DimPlot(object = object,group.by=unlist(molecular_phenotype.name),ncol=ncol,label = TRUE,label.size = 5)
  }
  return(plot)
}

