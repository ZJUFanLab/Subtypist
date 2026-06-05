
.check_standard <- function(resMarker,tua=0){
  temp <- FALSE
  if(sum(abs(resMarker$specificity_score)) > tua) {
    temp = TRUE
  }
  return(temp)
}

.cap_nonfinite_logfc <- function(x) {
  x <- as.numeric(x)
  finite.idx <- is.finite(x)
  if (all(finite.idx)) return(x)

  finite.max <- 0
  if (any(finite.idx)) {
    finite.max <- max(abs(x[finite.idx]), na.rm = TRUE)
    if (!is.finite(finite.max)) finite.max <- 0
  }

  inf.idx <- is.infinite(x)
  x[inf.idx] <- sign(x[inf.idx]) * finite.max
  x[!is.finite(x)] <- 0
  x
}

.getInitState <- function(resMarker,clusterNum,tmp,min.score=0){
  for(i in 1:clusterNum){
    tmp[i] <- .check_standard(resMarkers[resMarkers$cluster==i,])
  }
  return(tmp)
}


.updateState <- function(tmp,combined.min,combined.max,resMarker){
  tmp[combined.min] <- .check_standard(resMarker[resMarker$cluster==(combined.min-1),])
  tmp <- tmp[-combined.max]
  return(tmp)
}

.getJaccard <- function(cluster_top,a,b){
  intersections <- length(intersect(cluster_top[cluster_top$cluster==a,]$gene,cluster_top[cluster_top$cluster==b,]$gene))
  unions <- length(union(cluster_top[cluster_top$cluster==a,]$gene,cluster_top[cluster_top$cluster==b,]$gene))
  return(intersections / unions)
}

.getWeightedJaccard <- function(cluster_top, a, b) {
  genes_a <- cluster_top$gene[cluster_top$cluster == a]
  genes_b <- cluster_top$gene[cluster_top$cluster == b]

  unions <- length(union(genes_a, genes_b))
  intersectGene <- intersect(genes_a, genes_b)

  if (length(intersectGene) == 0) return(0)

  intersectInf <- cluster_top[
    cluster_top$cluster %in% c(a, b) & cluster_top$gene %in% intersectGene,]

  return(sum(abs(.cap_nonfinite_logfc(intersectInf$avg_log2FC))) / unions)
}

.getWeightedJaccardByRegulation <- function(cluster_top, a, b, regulation = c("both", "up", "down")) {
  regulation <- match.arg(regulation)
  if (regulation == "both") {
    up <- cluster_top[cluster_top$avg_log2FC > 0, ]
    down <- cluster_top[cluster_top$avg_log2FC < 0, ]
    return(mean(c(
      .getWeightedJaccard(up, a, b),
      .getWeightedJaccard(down, a, b)
    ), na.rm = TRUE))
  }
  .getWeightedJaccard(cluster_top, a, b)
}

.getInitWeightedJaccardMatrix <- function(clusterNum, resMarker, regulation = c("both", "up", "down")) {
  regulation <- match.arg(regulation)
  M <- matrix(0, nrow = clusterNum, ncol = clusterNum)
  for(i in 1:(clusterNum-1)){
    for(j in (i+1):clusterNum){
      JaccardScore <- .getWeightedJaccardByRegulation(resMarker, i - 1, j - 1, regulation = regulation)
      M[i,j] <- M[j,i] <- JaccardScore
    }
  }
  diag(M) <- 0
  return(M)
}


# .getInitWeightedJaccardMatrix <- function(clusterNum,resMarker){
#   M <- matrix(nrow=clusterNum,ncol=clusterNum)
#   for(i in 1:(clusterNum-1)){
#     M[i,i] <- 0
#     for(j in (i+1):clusterNum){
#       JaccardScore <- .getWeightedJaccard(resMarker,i-1,j-1)
#       #break
#       M[i,j] <- M[j,i] <- JaccardScore
#     }
#     M[j,j] <- 0
#   }
#   return(M)
# }

.getInitJaccardMatrix <- function(clusterNum,resMarker){
  M <- matrix(nrow=clusterNum,ncol=clusterNum)
  for(i in 1:(clusterNum-1)){
    M[i,i] <- 0
    for(j in (i+1):clusterNum){
      genes.1 <- resMarker[resMarker$cluster==i-1,]$gene
      genes.2 <- resMarker[resMarker$cluster==j-1,]$gene
      JaccardScore <- getJaccard(resMarker,i-1,j-1)
      #break
      M[i,j] <- M[j,i] <- JaccardScore
    }
    M[j,j] <- 0
  }
  return(M)
}



.updateDistanceMatrix <- function(M,Index.min,Index.max,resMarker,clusterNum,operation,regulation = NULL){ # operation传入函数
  # delete combined.max row
  M <- M[-Index.max,-Index.max]
  gene.new <- resMarker[resMarker$cluster==(Index.min-1),]$gene
  # update the information of combined.min row
  for(index in 1:(clusterNum-1)){
    if(index == Index.min){
      M[Index.min,index] <- 0
      next
    }
    if (is.null(regulation)) {
      JaccardScore <- operation(resMarker,Index.min-1,index-1)
    } else {
      JaccardScore <- .getWeightedJaccardByRegulation(resMarker, Index.min - 1, index - 1, regulation = regulation)
    }
    M[Index.min,index] <- JaccardScore
    M[index,Index.min] <- JaccardScore
  }
  return(M)
}

.initialClusterLabels <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  y <- suppressWarnings(as.numeric(x))
  if (all(!is.na(y))) return(y)
  as.character(x)
}

.updateSeuratObj <- function(obj,column,combined.min,combined.max,clusterNum){
  labels <- obj@meta.data[[column]]
  if (is.factor(labels)) labels <- as.character(labels)
  labels <- suppressWarnings(as.numeric(labels))

  labels[labels == combined.max] <- combined.min
  if(combined.max <= clusterNum-2){
    for(i in (combined.max+1):(clusterNum-1)){
      labels[labels == i] <- i - 1
    }
  }
  obj@meta.data[[column]] <- labels
  return(obj)
}

.updateMarkerlist <- function(cluster_top,combined.min,combined.max,newMarkers_top,clusterNum){ # newMarkers_top with cluster

  cluster_top <- subset(cluster_top,subset=cluster!=combined.min & cluster!=combined.max)

  newMarkers_top$cluster <- combined.min
  cols <- colnames(cluster_top)

  if(combined.max <= clusterNum-2){
    for(index in (combined.max+1):(clusterNum-1)){
      cluster_top$cluster[cluster_top$cluster==index] <- index - 1
    }
  }
  cluster_top <- rbind(cluster_top,newMarkers_top[,cols])
  return(cluster_top)
}

.mergeSteps <- function(mergedNodes,combined.min,combined.max){
  mergedNodes[[combined.min]] <- c(mergedNodes[[combined.min]],mergedNodes[[combined.max]])
  mergedNodes <- mergedNodes[-combined.max]
  return(mergedNodes)
}

.printSteps <- function(mergedNodes,column){
  for(i in 1:length(mergedNodes)){
    cat("The New group",(i-1),"is composed of :\n")
    for(y in 1:length(mergedNodes[[i]])){
      cat("inital cluster in",column," --> ",(mergedNodes[[i]][[y]]-1),"\n")
    }
    cat("\n")
  }
}

.Find_max_below_threshold <- function(mat,threshold){
  tmpmax <- -1
  first <- -1
  sceond <- -1
  for(i in 1:(ncol(mat)-1)){
    for(j in (i+1):ncol(mat)){
      if(i !=j){
        if(mat[i,j] > tmpmax && mat[i,j] < threshold){
          tmpmax <- mat[i,j]
          first <- i
          sceond <- j
        }
      }
    }
  }
  if(tmpmax == -1){
    return(list(threshold,0,0))
  }
  return(list(tmpmax,first,sceond))
}

.which.max.False <- function(M,tmp,max.col){
  first <- which(max.col==max(max.col[!tmp]))
  second <- which(M[,first] == max(max.col[!tmp]))
  return(list(first,second))
}

.setClulterInf <- function(resMarker, mergedNodes, resolution, regulation = c("up", "down", "both"),top_k = 3) {
  regulation <- match.arg(regulation)

  resMarker <- resMarker %>%
    dplyr::arrange(cluster, dplyr::desc(specificity_score), dplyr::desc(abs(avg_log2FC)))

  # Select top n marker genes per cluster
  top_genes <- resMarker %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = specificity_score, n = top_k, with_ties = FALSE)

  merge_cluster <- unique(top_genes$cluster)

  molecular_phenotype <- aggregate(top_genes$gene, by = list(type = top_genes$cluster), list)[-1]
  molecular_phenotype <- tibble::tibble(purrr::map(molecular_phenotype$x, ~ head(.x, top_k)))

  genes_score <- resMarker %>% dplyr::group_by(cluster)
  genes_score.top <- genes_score %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = specificity_score, n = top_k, with_ties = FALSE)

  score <- aggregate(genes_score.top$specificity_score, by = list(type = genes_score.top$cluster), sum)[,-1]

  clu <- cbind(
    resolution = rep(resolution, length(merge_cluster)),
    tibble::tibble(merge_cluster = sort(merge_cluster, decreasing = FALSE)),
    initial_cluster = tibble::tibble(purrr::map(mergedNodes, .f = function(x) purrr::map(x, .f = function(y) y - 1))),
    molecular_phenotype = molecular_phenotype,
    Score = score
  )

  colnames(clu) <- c("resolution", "merged_cluster", "initial_cluster", "phenotypic_molecules", "Score")
  return(clu)
}
#' Calculate marker specificity scores
#'
#' @param markers_top marker genes per cluster.
#' @param min.pct.1 Minimum percent of cells expressing the gene in the cluster.
#' @param min.diff Minimum absolute difference in expression percent between
#'   the target cluster and other clusters.
#' @param min.avg_log2FC Minimum average log2 fold change retained for backward
#'   compatibility.
#' @param regulation Marker regulation direction: up, down, or both.
#'
#' @return return a specificity score
#' @export
#'
getSpecificity_score <- function(markers_top,
                                 min.pct.1 = 0.1,
                                 min.diff = 0.5,
                                 min.avg_log2FC = 0.5,
                                 regulation = c("up", "down", "both")) {
  regulation <- match.arg(regulation)

  x1 <- markers_top$pct.1 - markers_top$pct.2
  avg_log2FC <- .cap_nonfinite_logfc(markers_top$avg_log2FC)
  markers_top$specificity_score <- 0
  delta <- max(c(as.numeric(min.diff), 0.25), na.rm = TRUE)

  if (regulation == "up") {
    valid <- avg_log2FC > 0 & x1 > 0
  } else if (regulation == "down") {
    valid <- avg_log2FC < 0 & x1 < 0
  } else {
    valid <- sign(avg_log2FC) == sign(x1) & sign(avg_log2FC) != 0
  }

  idx <- which(valid & abs(x1) >= delta)
  if (length(idx) > 0) {
    markers_top$specificity_score[idx] <- x1[idx] * avg_log2FC[idx]
  }

  markers_top$avg_log2FC <- avg_log2FC

  markers_top <- dplyr::arrange(markers_top, dplyr::desc(specificity_score), dplyr::desc(abs(avg_log2FC)))
  return(markers_top)
}
