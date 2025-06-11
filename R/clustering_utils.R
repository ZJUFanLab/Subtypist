
.check_standard <- function(resMarker){
  temp <- FALSE
  if(sum(resMarker$specificity_score) > 0) {
    temp = TRUE
  }
  return(temp)
}

.getInitState <- function(resMarker,clusterNum,tmp,min.score=0){
  for(i in 1:clusterNum){
    tmp[i] <- check_standard(resMarkers[resMarkers$cluster==i,])
  }
  return(tmp)
}


.updateState <- function(tmp,combined.min,combined.max,resMarker){
  tmp[combined.min] <- check_standard(resMarker[resMarker$cluster==(combined.min-1),])
  tmp <- tmp[-combined.max]
  return(tmp)
}

.getJaccard <- function(cluster_top,a,b){
  intersections <- length(intersect(cluster_top[cluster_top$cluster==a,]$gene,cluster_top[cluster_top$cluster==b,]$gene))
  unions <- length(union(cluster_top[cluster_top$cluster==a,]$gene,cluster_top[cluster_top$cluster==b,]$gene))
  return(intersections / unions)
}

.getWeightedJaccard <- function(cluster_top,a,b){
  unions <- length(union(cluster_top[cluster_top$cluster==a,]$gene,cluster_top[cluster_top$cluster==b,]$gene))

  intersectGene <- intersect(cluster_top[cluster_top$cluster==a,]$gene,cluster_top[cluster_top$cluster==b,]$gene)
  # get avg_log2FC

  if(length(intersectGene) == 0){return(0)}
  intersectInf <- cluster_top[(cluster_top$cluster==b | cluster_top$cluster==a ) & cluster_top$gene %in% intersectGene,]
  return(sum(intersectInf$avg_log2FC) / unions)
}

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

.getInitWeightedJaccardMatrix <- function(clusterNum,resMarker){
  M <- matrix(nrow=clusterNum,ncol=clusterNum)
  for(i in 1:(clusterNum-1)){
    M[i,i] <- 0
    for(j in (i+1):clusterNum){
      JaccardScore <- .getWeightedJaccard(resMarker,i-1,j-1)
      #break
      M[i,j] <- M[j,i] <- JaccardScore
    }
    M[j,j] <- 0
  }
  return(M)
}

.updateDistanceMatrix <- function(M,Index.min,Index.max,resMarker,clusterNum,operation){ # operation传入函数
  # delete combined.max row
  M <- M[-Index.max,-Index.max]
  gene.new <- resMarker[resMarker$cluster==(Index.min-1),]$gene
  # update the information of combined.min row
  for(index in 1:(clusterNum-1)){
    if(index == Index.min){
      M[Index.min,index] <- 0
      next
    }
    JaccardScore <- operation(resMarker,Index.min-1,index-1)
    M[Index.min,index] <- JaccardScore
    M[index,Index.min] <- JaccardScore
  }
  return(M)
}

.updateSeuratObj <- function(obj,column,combined.min,combined.max,clusterNum){
  sc_meta <- obj@meta.data
  sc_meta[sc_meta[[column]]==combined.max,][[column]] <- combined.min
  if(combined.max <= clusterNum-2){
    for(i in (combined.max+1):(clusterNum-1)){
      new <- i-1
      sc_meta[sc_meta[[column]]==i,][[column]] <- new
    }
  }
  obj@meta.data <- sc_meta
  return(obj)
}

.updateMarkerlist <- function(cluster_top,combined.min,combined.max,newMarkers_top,clusterNum){ # newMarkers_top with cluster

  cluster_top <- subset(cluster_top,subset=cluster!=combined.min & cluster!=combined.max)

  newMarkers_top$cluster <- combined.min
  cols <- colnames(cluster_top)

  if(combined.max <= clusterNum-2){
    for(index in (combined.max+1):(clusterNum-1)){
      cluster_top[cluster_top$cluster==index,]$cluster <- index - 1
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

.setClulterInf <- function(resMarker,mergedNodes,resolution){
  top_genes <- resMarker %>% dplyr::arrange(desc(avg_log2FC)) %>% dplyr::arrange(desc(specificity_score))%>% dplyr::arrange(desc(cluster)) %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=3,wt=specificity_score)
  merge_cluster <- unique(top_genes$cluster)
  molecular_phenotype <- aggregate(top_genes$gene, by=list(type=top_genes$cluster),list)[-1]
  molecular_phenotype <- tibble::tibble(purrr::map(molecular_phenotype$x,.f=function(y){return(y[0:3])}))
  genes_score <- resMarker %>% dplyr::arrange(desc(avg_log2FC)) %>% dplyr::arrange(desc(specificity_score))%>% dplyr::arrange(desc(cluster)) %>% dplyr::group_by(cluster)
  # Allscore <- genes_score %>% group.by(cluster) %>% top_n(n=3,wt=Score)
  # score <- aggregate(AllScore,by=list(type=Allscore$cluster),sum)[,-1]
  genes_score.top <- genes_score %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=3,wt=avg_log2FC)
  score<-aggregate(genes_score.top$specificity_score, by=list(type=genes_score.top$cluster),sum)[,-1]
  clu <- cbind(resolution=rep(resolution,length(merge_cluster)),
               tibble::tibble(merge_cluster=sort(merge_cluster,decreasing = FALSE)),
               initial_cluster=tibble::tibble(purrr::map(mergedNodes,.f=function(x){return(purrr::map(x,.f=function(y){return(y-1)}))})),
               molecular_phenotype=molecular_phenotype,
               Score=score)
  colnames(clu)[1] <- "resolution"
  colnames(clu)[2] <- "merge_cluster"
  colnames(clu)[3] <- 'initial_cluster'
  colnames(clu)[4] <- "molecular_phenotype"
  colnames(clu)[5] <- "Score"
  #clu$resolution <- map(clu$resolution,list)
  return(clu)
}

getSpecificity_score <- function(markers_top,min.pct.1=0.1,min.diff=0,min.avg_log2FC=0.5){ # pct.1 findMarkers(min.pct),
  b <- max(min.diff,0.5)
  x1 <- markers_top$pct.1 - markers_top$pct.2
  markers_top$specificity_score <- 0
  if(nrow(markers_top[x1>=b,]) > 0){markers_top[x1>=b,]$specificity_score <- (markers_top[x1>=b,]$pct.1 - markers_top[x1>=b,]$pct.2) * markers_top[x1>=b,]$avg_log2FC}
  if(nrow(markers_top[x1<0.25,]) > 0){markers_top[x1<0.25,]$specificity_score <- 0}
  if(nrow(markers_top[x1>=0.25
                      & x1<b
                      & markers_top$pct.1<=0.5,]) > 0)
  {markers_top[x1>=0.25
               & x1<b
               & markers_top$pct.1<=0.5,]$specificity_score <- (markers_top[x1>=0.25 & x1<b & markers_top$pct.1<=0.5,]$pct.1-markers_top[x1>=0.25 & x1<b & markers_top$pct.1<=0.5,]$pct.2) * markers_top[x1>=0.25 & x1<b & markers_top$pct.1<=0.5,]$avg_log2FC
  }
  if(nrow(markers_top[x1>=0.25
                      & x1<b
                      & markers_top$pct.1>0.5
                      & markers_top$avg_log2FC>=min.avg_log2FC,]) > 0)
  {markers_top[x1>=0.25
               & x1<b
               & markers_top$pct.1>0.5
               & markers_top$avg_log2FC>=min.avg_log2FC,]$specificity_score <- (markers_top[x1>=0.25 & x1<b & markers_top$pct.1>0.5 & markers_top$avg_log2FC>=min.avg_log2FC,]$pct.1-markers_top[x1>=0.25 & x1<b & markers_top$pct.1>0.5 & markers_top$avg_log2FC>=min.avg_log2FC,]$pct.2) * markers_top[x1>=0.25 & x1<b & markers_top$pct.1>0.5 & markers_top$avg_log2FC>=min.avg_log2FC,]$avg_log2FC}
  if(nrow(markers_top[x1>=0.25
                      & x1<b
                      & markers_top$pct.1>0.5
                      & markers_top$avg_log2FC<min.avg_log2FC,]) > 0)
  {markers_top[x1>=0.25
               & x1<b
               & markers_top$pct.1>0.5
               & markers_top$avg_log2FC<min.avg_log2FC,]$specificity_score <- (markers_top[x1>=0.25 & x1<b & markers_top$pct.1>0.5 & markers_top$avg_log2FC < min.avg_log2FC,]$pct.1-markers_top[x1>=0.25 & x1<b & markers_top$pct.1>0.5 & markers_top$avg_log2FC <  min.avg_log2FC,]$pct.2) * markers_top[x1>=0.25 & x1<b & markers_top$pct.1>0.5 & markers_top$avg_log2FC <  min.avg_log2FC,]$avg_log2FC}
  markers_top <- dplyr::arrange(markers_top,desc(specificity_score))
  return(markers_top)
}
