# Subtypist


![R >= 4.0](https://img.shields.io/badge/R-%3E=4.0-blue) ![Seurat >= 4.0.0](https://img.shields.io/badge/Seurat-%3E=4.0-green)

#### **Reference-free identification of cell subtypes for single-cell transcriptomic data**

![curation](https://github.com/ZJUFanLab/Subtypist/blob/main/img/curation.png)

Identification of novel cell subtypes is critically crucial in revealing the pathogenesis and heterogeneity of disease, which provides unprecedented insights into the development of therapeutic strategies. Although numerous cell-type identification methods exist, these methods heavily rely on the reference with the fixed cell labels, which fail to uncover new cell subtypes marked with phenotypic molecules within a specific disease context. To fill this gap, we propose a pioneering reference-free annotation method, Subtypist, to identify disease-associated cell subtypes expressing phenotypic features using an ensemble-based strategy. 

# Install

- install dependent packages `devtools` 

```R
> install.packages(pkgs = 'devtools')
```

- then install Subtypist

```R
> devtools::install_github('ZJUFanLab/Subtypist')

# or download the repository as ZIP
> devtools::install_local("/path/to/Subtypist-main.zip")
```

# **usage**

1. **Loading or subsetting specific cell type**

``` R
# Load Seurat object from a full dataset and subset specific cell types
> FullObject <- Load("FullObject.RData")
> Seu <- subset(FullObject, idents = c("T cells")

# Alternatively, load a pre-processed Seurat object with annotated cell types, derived by subsetting the full dataset
> # Seu <- readRDS("Seu.rds")
```

2. **Cell type identification using an ensemble strategy without reference**

```R
# Object: a Loading or subsetting specific cell type Seurat object
# min.resolution: the minimum value of resolution
# max.resolution: the maxium value of resolution
# use.assay: Name of assay to use
# cluster.assay: Name of the assay in the Seurat object to use for clustering
> result <- Subtypist_merge(object=Seu,min.resolution=0.3,max.resolution=1.5,by=0.1,use.assay="RNA",cluster_assay = "RNA")

# Show results
> print(result)
$Object
An object of class Seurat 
1000 features across 1000 samples within 1 assay 
Active assay: RNA (1000 features, 985 variable features)
 3 dimensional reductions calculated: pca, umap, tsne
$result.table
   resolution merge_cluster initial_cluster       molecular_phenotype     Score
1         0.4             0               0   Gene7, Gene679, Gene990 2.0322883
2         0.4             1               1 Gene570, Gene807, Gene470 0.6984284
3         0.4             2               2 Gene871, Gene559, Gene247 3.3095874
4         0.4             3               3 Gene474, Gene746, Gene790 3.5678922
5         0.4             4               4 Gene776, Gene507, Gene323 2.9614237
6         1.1             0               0   Gene7, Gene679, Gene990 2.0341814
7         1.1             1               1 Gene871, Gene559, Gene247 3.3095874
8         1.1             2             2,5 Gene470, Gene807, Gene243 0.0000000
9         1.1             3               3 Gene474, Gene746, Gene577 3.4186799
10        1.1             4             4,7 Gene807, Gene801, Gene566 0.0000000
11        1.1             5               6 Gene776, Gene507, Gene323 2.9614237

```

3. **Evaluating clustering resolutions and annotating subtypes with specific phenotypic markers**

```R
# To evaluate and rank clustering resolutions based on their corresponding subtype identification results.
> sortScore(result$result.table) ## resolution = 0.4: highest 
# A tibble: 5 × 2
  resolution value
       <dbl> <dbl>
1        0.1  1.62
2        0.2  1.85
3        0.4  2.51
4        1.1  1.95
5        1.2  1.81

> # Add the result to the object 
> Seu <- AddSubtypist(result$Object,result.table=result$result.table,prefix='Subtypist')
> # To assign more specific phenotypic molecules to each subtype, 
> # the `select_index` parameter can be used to specify which gene to select 
> Seu <- Subtypist::AddSubtypist(result$Object,resolution=c(0.4),result.table=result$result.table,prefix = 'Subtypist',meta.prefix = 'phenotypic melocules_',value.suffix='+ B',select_index=c('0'=1,'1'=1,'2'=1,'3'=2,'4'=3))
> print(unique(Seu@meta.data['phenotypic melocules_0.4'])
> [1] "Gene570+ B" "Gene7+ B"   "Gene746+ B" "Gene871+ B" "Gene323+ B"

```

## Plot

- Visualize subtype-level distributions across dimensionality reduction space (e.g., UMAP) using Subtypist_Dimplot(). This function overlays the specified phenotypic molecules annotations—derived at selected clustering resolutions—onto the Seurat object. For example:

```R
> p <- Subtypist::Subtypist_Dimplot(Seu,result.table = result$result.table,resolution = c(0.4,1.1), show = "molecular_phenotype_",prefix = 'Subtypist')
```

# About

Subtypist was developed by Yue Yao. Should you have any questions, please contact Yue Yao at [yuey@zju.edu.cn](mailto:yueyo@zju.edu.cn)

