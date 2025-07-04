---
title: "Introduction to Single Cell RNA-Seq Part 8: Integration"
author: "UCD Bioinformatics Core"
date: "2025-07-04"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---
# Introduction to Single Cell RNA-Seq Part 8: Integration

More and more experiments involve a large number of samples/datasets, that may have been prepared in separate batches. Or in the case where one would like to include or integrate publicly available datasets. It is important to properly integrate these datasets, and we will see the effect the integration has at the end of this documentation.

Most of the methods that were developed to integrate single cell datasets fall into two categories. The first is the "anchor" based approach. In this approach, the first step is to select a batch as the "anchor" and convert other batches to the "anchor" batch. Among these approaches are [MNN](https://github.com/MarioniLab/MNN2017), [iMAP](https://github.com/Svvord/iMAP), and [SCALEX](https://github.com/jsxlei/SCALEX). The advantage of the anchor-based approach is that different batches of cells can be studied under the same experimental conditions, and the disadvantage is that it is not possible to fully combine the features of each batch because the cell types contained in each batch are unknown.

The second approach is to transform all batches of data to a low-dimensional space to correct batch effects, such as implemented in [Harmony](https://github.com/immunogenomics/harmony), [DESC](https://www.nature.com/articles/s41467-020-15851-3), [BBKNN](https://github.com/Teichlab/bbknn), [STACAS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8098019/) and [Seurat's integration](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8). This second approach has the advantage of extracting biologically relevant latent features and reducing the impact of noise, but it cannot be used for differential gene expression analysis. Many of these existing methods work well when the batches of datasets have the same cell types, however, they fail when there are different cell types involved in different datasets. [Scanorama](https://github.com/brianhie/scanorama) uses similar approach, but it allows integration of datasets that don't always share a common cell type among all and the batch-corrected gene expression data can be used for differential gene expression analysis.

[scVI](https://www.nature.com/articles/s41592-018-0229-2) is based on a hierarchical Baysian model with conditional distributions specified by deep neural networks. The expression of a gene in a cell is modeled using a zero-inflated negative binomial distribution, conditioned on batch annotaton (if available), as well as two unoberserved random variable. One is the variation due to differences in capture efficiency and sequencing depth and it serves as a cell-specific scaling factor. The other captures the biological differences. This frame work used by scVI allows for integration of datasets from different experiment, and permits differential expression analysis on estimated expression data.

Recently, [IMGG](https://www.mdpi.com/1422-0067/23/4/2082) has been developed that uses connected graphs and generative adversarial networks (GAN) to achieve the goal of eliminating nonbiological noise between batches of datasets. This new method has been demonstrated to work well both in the situation where datasets have the same cell types and in the situation where datasets may have different cell types.

In this workshop, we are going to look at Seurat's integration approach using reciprocal PCA, which is supurior to its first integration approach using canonical correlation analysis. The basic idea is to identify cross-dataset pairs of cells that are in a matched biological state ("anchors"), and use them to correct technical differences between datasets. The integration method we use has been implemented in Seurat and you can find the details of the method in [its publication](https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub).

## Set up workspace

``` r
library(Seurat)
library(ggplot2)
set.seed(12345)
experiment.aggregate <- readRDS("scRNA_workshop-02.rds") # filtered object
```

## Prepare data for integration
Prior to integration samples should be processed independently. First, we split the filtered object by sample to create a list of Seurat objects.

``` r
experiment.split <- SplitObject(experiment.aggregate, split.by = "orig.ident")
rm(experiment.aggregate)
```

Each object is then normalized and variable features detected. CellCycleScoring is run here for later visualization purpose, and is not required for integration purpose.

``` r
experiment.split <- lapply(experiment.split, function(sce){
  sce = NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce = CellCycleScoring(sce, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
  sce = FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
})
```

## Select integration features
Integration features are genes that are repeatedly variable across the objects to integrate. These are used to scale the data and run the PCA.

``` r
features <- SelectIntegrationFeatures(object.list = experiment.split)
```

## Scale data and run PCA
Once integration features have been identified, we can scale the data and run the PCA, which is required for finding the integration anchors.

``` r
experiment.split <- lapply(experiment.split,function(sce){
  sce = ScaleData(sce, features = features, vars.to.regress = c("S.Score", "G2M.Score", "percent_MT", "nFeature_RNA"))
  RunPCA(sce, features = features)
})
```

## Idenfity integration anchors
The integration anchors are pairs of cells that are mutual nearest neighbors on the shared low-dimensional representation. These may be calculated either for each Seurat object relative to a reference object, or pairwise between all objects if no reference is provided.

``` r
anchors <- FindIntegrationAnchors(object.list = experiment.split, anchor.features = features, reduction = "rpca")
```

## Integrate

``` r
experiment.integrated <- IntegrateData(anchorset = anchors)
experiment.integrated$group <- factor(experiment.integrated$group, levels=c("Normal", "Polyp", "Colorectal Cancer"))
```
The new experiment.integrated object has two assays: RNA and integrated. The RNA assay contains the normalized, scaled data from the individual experiment.split objects merged into a single table, while the data in the integrated assay has been scaled in such a way that it is no longer appropriate to use this assay for differential expression.

The authors recommend using the integrated assay for clustering and visualization (UMAP plots).

## Impact of integration
In the dimensionality reduction section we performed PCA on the un-integrated experiment.aggregate object, where we used the vars.to.regress argument of the ScaleData function to adjust for cell cycle, nucleus integrity, and sequencing depth. The PCA biplot looked like this:

![Previous PCA plot](04-dimensionality_reduction_files/figure-html/plot_pca-1.png)

After integration, the appearance of the PCA biplot has changed; cells no longer separate by group.

``` r
DefaultAssay(experiment.integrated) <- "integrated"
experiment.integrated <- ScaleData(experiment.integrated, assay="integrated")
experiment.integrated <- RunPCA(experiment.integrated, assay="integrated")
DimPlot(experiment.integrated,
        group.by = "group",
        reduction = "pca",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "mako")
```

![](08-integration_files/figure-html/PCA-1.png)<!-- -->

A similar effect can be seen in the UMAP. Previously, the un-integrated UMAP plot had this appearance:

![Previous UMAP plot](04-dimensionality_reduction_files/figure-html/unnamed-chunk-1-1.png)

After integration, the polyp and colorectal cancer cells are more co-localized on the biplot.

``` r
experiment.integrated <- RunUMAP(experiment.integrated, dims = 1:50)
DimPlot(experiment.integrated,
        reduction = "umap",
        group.by = "group",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "mako")
```

![](08-integration_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

In the next section, we will use the integrated data to perform clustering.

### Visualize metadata


``` r
lapply(c("nCount_RNA", "nFeature_RNA", "percent_MT"), function(feature){
  FeaturePlot(experiment.integrated,
              reduction = "umap",
              features = feature)
})
```

```
## [[1]]
```

![](08-integration_files/figure-html/meta-1.png)<!-- -->

```
## 
## [[2]]
```

![](08-integration_files/figure-html/meta-2.png)<!-- -->

```
## 
## [[3]]
```

![](08-integration_files/figure-html/meta-3.png)<!-- -->

### Visualize cell cycle phase


``` r
DimPlot(experiment.integrated,
        reduction = "umap",
        group.by = "Phase",
        shuffle = TRUE) +
  scale_color_viridis_d()
```

![](08-integration_files/figure-html/phase-1.png)<!-- -->

### Clusters using the integrated data


``` r
experiment.integrated <- FindNeighbors(experiment.integrated, reduction = "pca", dims = 1:50)
experiment.integrated <- FindClusters(experiment.integrated,
                                     resolution = seq(0.04, 0.07, 0.01))
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6315
## Number of edges: 271793
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9787
## Number of communities: 10
## Elapsed time: 0 seconds
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6315
## Number of edges: 271793
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9761
## Number of communities: 10
## Elapsed time: 0 seconds
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6315
## Number of edges: 271793
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9736
## Number of communities: 10
## Elapsed time: 0 seconds
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6315
## Number of edges: 271793
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9710
## Number of communities: 10
## Elapsed time: 0 seconds
```

``` r
cluster.resolutions <- grep("res", colnames(experiment.integrated@meta.data), value = TRUE)
sapply(cluster.resolutions, function(res){
  length(levels(experiment.integrated@meta.data[,res]))
})
```

```
## integrated_snn_res.0.04 integrated_snn_res.0.05 integrated_snn_res.0.06 
##                      10                      10                      10 
## integrated_snn_res.0.07 
##                      10
```

``` r
cluster.resolutions
```

```
## [1] "integrated_snn_res.0.04" "integrated_snn_res.0.05"
## [3] "integrated_snn_res.0.06" "integrated_snn_res.0.07"
```

``` r
lapply(cluster.resolutions, function(res){
  DimPlot(experiment.integrated,
          group.by = res,
          reduction = "umap",
          shuffle = TRUE) +
    scale_color_viridis_d(option = "turbo")
})
```

```
## [[1]]
```

![](08-integration_files/figure-html/clusters.integrated-1.png)<!-- -->

```
## 
## [[2]]
```

![](08-integration_files/figure-html/clusters.integrated-2.png)<!-- -->

```
## 
## [[3]]
```

![](08-integration_files/figure-html/clusters.integrated-3.png)<!-- -->

```
## 
## [[4]]
```

![](08-integration_files/figure-html/clusters.integrated-4.png)<!-- -->

``` r
lapply(cluster.resolutions, function(res){
         tmp = experiment.integrated@meta.data[,c(res, "group")]
         colnames(tmp) = c("cluster", "group")
         ggplot(tmp, aes(x = cluster, fill = group)) +
           geom_bar() +
           scale_fill_viridis_d(option = "mako") +
           theme_classic()
})
```

```
## [[1]]
```

![](08-integration_files/figure-html/clusters.integrated-5.png)<!-- -->

```
## 
## [[2]]
```

![](08-integration_files/figure-html/clusters.integrated-6.png)<!-- -->

```
## 
## [[3]]
```

![](08-integration_files/figure-html/clusters.integrated-7.png)<!-- -->

```
## 
## [[4]]
```

![](08-integration_files/figure-html/clusters.integrated-8.png)<!-- -->

``` r
FeaturePlot(experiment.integrated,
            reduction = "umap",
            features = "KCNMA1")
```

![](08-integration_files/figure-html/clusters.integrated-9.png)<!-- -->

``` r
VlnPlot(experiment.integrated,
        group.by = "integrated_snn_res.0.07",
        features = "KCNMA1",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](08-integration_files/figure-html/clusters.integrated-10.png)<!-- -->

``` r
Idents(experiment.integrated) <- "integrated_snn_res.0.07"
experiment.integrated <- BuildClusterTree(experiment.integrated, dims = 1:50)
PlotClusterTree(experiment.integrated)
```

![](08-integration_files/figure-html/clusters.integrated-11.png)<!-- -->



``` r
experiment.aggregate <- readRDS("scRNA_workshop-05.rds")
identical(rownames(experiment.aggregate@meta.data), rownames(experiment.integrated@meta.data))
```

```
## [1] TRUE
```

``` r
experiment.aggregate <- AddMetaData(experiment.aggregate,
                                    metadata = experiment.integrated$integrated_snn_res.0.07,
                                    col.name = "integrated_snn_res.0.07")
DimPlot(experiment.aggregate,
        reduction = "umap",
        group.by = "integrated_snn_res.0.07",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "turbo")
```

![](08-integration_files/figure-html/AddMetaData-1.png)<!-- -->

In this case, integration appears to have interfered with forming easily-interpreted clusters. There is very little relationship between location on the UMAP and cluster identity, which makes it harder to identify possible cell populations at a glance. We can add the integrated object's cluster identities to the un-integrated object if we choose, though projecting the integrated clustering onto the un-integrated UMAP is unlikely to be useful.

#### Save the Seurat object


``` r
saveRDS(experiment.integrated, file="scRNA_workshop-08.rds")
```


#### Session information

``` r
sessionInfo()
```

```
## R version 4.4.3 (2025-02-28)
## Platform: aarch64-apple-darwin20
## Running under: macOS Ventura 13.7.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/Los_Angeles
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_3.5.1      Seurat_5.2.1       SeuratObject_5.0.2 sp_2.1-4          
## 
## loaded via a namespace (and not attached):
##   [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3         
##   [4] rlang_1.1.3            magrittr_2.0.3         RcppAnnoy_0.0.22      
##   [7] spatstat.geom_3.2-9    matrixStats_1.3.0      ggridges_0.5.6        
##  [10] compiler_4.4.3         png_0.1-8              vctrs_0.6.5           
##  [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3       
##  [16] fastmap_1.2.0          labeling_0.4.3         utf8_1.2.4            
##  [19] promises_1.3.0         rmarkdown_2.27         ggbeeswarm_0.7.2      
##  [22] purrr_1.0.2            xfun_0.44              cachem_1.1.0          
##  [25] jsonlite_1.8.8         goftest_1.2-3          highr_0.11            
##  [28] later_1.3.2            spatstat.utils_3.1-2   irlba_2.3.5.1         
##  [31] parallel_4.4.3         cluster_2.1.8          R6_2.5.1              
##  [34] ica_1.0-3              bslib_0.7.0            stringi_1.8.4         
##  [37] RColorBrewer_1.1-3     spatstat.data_3.0-4    reticulate_1.39.0     
##  [40] parallelly_1.37.1      lmtest_0.9-40          jquerylib_0.1.4       
##  [43] scattermore_1.2        Rcpp_1.0.12            knitr_1.47            
##  [46] tensor_1.5             future.apply_1.11.2    zoo_1.8-12            
##  [49] sctransform_0.4.1      httpuv_1.6.15          Matrix_1.7-2          
##  [52] splines_4.4.3          igraph_2.0.3           tidyselect_1.2.1      
##  [55] abind_1.4-5            rstudioapi_0.16.0      yaml_2.3.8            
##  [58] spatstat.random_3.2-3  codetools_0.2-20       miniUI_0.1.1.1        
##  [61] spatstat.explore_3.2-7 listenv_0.9.1          lattice_0.22-6        
##  [64] tibble_3.2.1           plyr_1.8.9             withr_3.0.0           
##  [67] shiny_1.8.1.1          ROCR_1.0-11            ggrastr_1.0.2         
##  [70] evaluate_0.23          Rtsne_0.17             future_1.33.2         
##  [73] fastDummies_1.7.3      survival_3.8-3         polyclip_1.10-6       
##  [76] fitdistrplus_1.1-11    pillar_1.9.0           KernSmooth_2.23-26    
##  [79] plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0        
##  [82] munsell_0.5.1          scales_1.3.0           globals_0.16.3        
##  [85] xtable_1.8-4           glue_1.7.0             lazyeval_0.2.2        
##  [88] tools_4.4.3            data.table_1.15.4      RSpectra_0.16-1       
##  [91] RANN_2.6.1             dotCall64_1.1-1        cowplot_1.1.3         
##  [94] grid_4.4.3             ape_5.8                tidyr_1.3.1           
##  [97] colorspace_2.1-0       nlme_3.1-167           patchwork_1.2.0       
## [100] beeswarm_0.4.0         vipor_0.4.7            cli_3.6.2             
## [103] spatstat.sparse_3.0-3  spam_2.10-0            fansi_1.0.6           
## [106] viridisLite_0.4.2      dplyr_1.1.4            uwot_0.2.2            
## [109] gtable_0.3.5           sass_0.4.9             digest_0.6.35         
## [112] progressr_0.14.0       ggrepel_0.9.5          farver_2.1.2          
## [115] htmlwidgets_1.6.4      htmltools_0.5.8.1      lifecycle_1.0.4       
## [118] httr_1.4.7             mime_0.12              MASS_7.3-64
```
