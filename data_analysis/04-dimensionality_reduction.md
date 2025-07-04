---
title: "Introduction to Single Cell RNA-Seq Part 4: Dimensionality reduction"
author: "UCD Bioinformatics Core"
date: "2025-07-03"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Introduction to Single Cell RNA-Seq Part 4: Dimensionality reduction
Single cell (or nucleus) data are extremely high-dimensional. In order to reduce the complexity of analysis and remove sources of noise, dimensionality reduction is an important step in the analysis workflow. In this section, we will be using two dimension reduction methods: PCA and UMAP.


## Set up workspace

``` r
library(Seurat)
library(ggplot2)
experiment.aggregate <- readRDS(file="scRNA_workshop-03.rds")
experiment.aggregate
```

```
## An object of class Seurat 
## 11474 features across 6315 samples within 1 assay 
## Active assay: RNA (11474 features, 7110 variable features)
##  3 layers present: counts, data, scale.data
```

``` r
set.seed(12345)
```

## Perform dimensionality reduction with PCA
Principal Components Analysis (PCA) is a widely-used dimension reduction method. Each PC is a vector in the reduced-dimensional space that is orthogonal to all preceding PCs. The first of these explains the largest amount of variation and each subsequent PC explains slightly less than the preceding component. PCA is performed on the scaled data, and therefore uses only the variable features. 

``` r
?RunPCA
```


``` r
experiment.aggregate <- RunPCA(experiment.aggregate, npcs = 100)
```

While it is theoretically possible to calculate as many PCs as there are features in the data, typically 100 PCs is more than sufficient. In fact, many of these PCs may explain negligible amounts of variation. Seurat provides a number of ways to visualize the PCA results.

### Principal components plot
The PCA biplot is a scatter plot showing the placement of each cell on two selected components. By default, the first and second PC are used, but any two calculted PCs may be specified.

At this point in the analysis, since we are no longer performing QA and filtering, we can move to examining relationships between cells on a per-group rather than per-sample basis.

``` r
DimPlot(experiment.aggregate,
        group.by = "group",
        reduction = "pca",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "mako")
```

![](04-dimensionality_reduction_files/figure-html/plot_pca-1.png)<!-- -->

The axes are unit-less; points (cells or nuclei) that are farther apart are more dissimilar on the displayed PC than points that are closer together.

### PCA loadings
Each PC can be imagined as a sort of meta-gene (a linear combination of the genes included in creating the PCs: the variable features) for which every cell has an expression value. The top genes associated with a PC component can be plotted for a selected dimension(s) using the `VizDimLoadings` function.

``` r
VizDimLoadings(experiment.aggregate,
               dims = 1,
               nfeatures = 25,
               reduction = "pca",
               ncol = 1) +
  theme_minimal(base_size = 8)
```

![](04-dimensionality_reduction_files/figure-html/viz_pca-1.png)<!-- -->

### Heat map
Heat maps can be used to visualize the expressions of the top associated gene with a PC component. On the x-axis, cells are ordered by their embeddings ("expression" of the PC), while on the y-axis, genes are ordered by PC loading. When fewer than the total number of cells is selected, this results in selection of the cells with the largest absolute value embeddings, which emphasizes variation on the PC.

``` r
DimHeatmap(experiment.aggregate,
           dims = 1,
           nfeatures = 25,
           cells = 500,
           reduction = "pca",
           balanced = TRUE,
           slot = "scale.data")
```

![](04-dimensionality_reduction_files/figure-html/heatmap_pca-1.png)<!-- -->

#### Explore
Re-import the original data and try modifying the ScaleData vars.to.regress argument. You could remove some variables, or add others. What happens? See how choices effect the plots.

``` r
experiment.explore <- readRDS("scRNA_workshop-03.rds")
experiment.explore <- ScaleData(experiment.explore) # make changes here to explore the data
experiment.explore <- RunPCA(experiment.explore) # what happens if you adjust npcs?
VizDimLoadings(experiment.explore, dims = 1:2)
DimPlot(experiment.explore, reduction = "pca")
DimHeatmap(experiment.explore, dims = 1:6, cells = 500, balanced = TRUE) # adjust parameters
```

## Selecting PCs to use
To overcome the extensive technical noise in any single gene, Seurat clusters cells based on their PCA scores, with each PC essentially representing a meta-gene that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step.

### Elbow plot
An elbow plot displays the standard deviations (or approximate singular values if running PCAFast) of the principle components for easy identification of an elbow in the graph. This elbow often corresponds well with the significant PCs and is much faster to run.  This is the traditional approach to selecting principal components.

The appearance of elbow plots tends to be highly consistent across single cell / single nucleus experiments. Generally, the line approaches zero at around PC 50. This is a reasonable number of PCs to use for the downstream steps.

``` r
ElbowPlot(experiment.aggregate, ndims = 100)
```

![](04-dimensionality_reduction_files/figure-html/elbow-1.png)<!-- -->

### JackStraw
The JackStraw function randomly permutes a subset of data, and calculates projected PCA scores for these genes. The PCA scores for these randomly permuted genes are then compared with the observed PCA scores to determine statistical significance. The end result is a p-value for each gene's association with each principal component.

PCs with a strong enrichment of low p-value genes are identified as significant components.

**The JackStraw permutation is computationally intensive and can be quite slow. Consider skipping this step and exploring the function when you have some extra time.**

``` r
experiment.aggregate <- JackStraw(experiment.aggregate, dims = 100)
experiment.aggregate <- ScoreJackStraw(experiment.aggregate, dims = 1:100)
JackStrawPlot(object = experiment.aggregate, dims = 1:100) +
  scale_color_viridis_d() +
  theme(legend.position="bottom")
```

![](04-dimensionality_reduction_files/figure-html/jackstraw-1.png)<!-- -->

Let's use the first 50 PCs.


``` r
use.pcs <- 1:50
```

## UMAP
[Uniform Manifold Approximation and Projection](https://arxiv.org/pdf/1802.03426v3.pdf) (UMAP) is a dimensionality reduction method that is commonly used in single cell RNA-Seq analysis. Single cell data is extremely high-dimensional; UMAP calculates a nearest neighbor network describing the relationships between cells as measured by the PC loadings of variable genes and creates a low-dimensional space that preserves these relationships.

``` r
# calculate UMAP
experiment.aggregate <- RunUMAP(experiment.aggregate,
                                dims = use.pcs)
```

While UMAP can be a general non-linear dimensionality reduction approach, it's most frequently used as a visualization technique. A UMAP biplot offers a very useful graphical representation of the relationships captured by the nearest neighbor graph.


``` r
# UMAP colored by sample identity
DimPlot(experiment.aggregate,
        group.by = "group",
        reduction = "umap",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "mako")
```

![](04-dimensionality_reduction_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

## Prepare for the next section

#### Save object

``` r
saveRDS(experiment.aggregate, file="scRNA_workshop-04.rds")
```

#### Download Rmd

``` r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2025-July-Single-Cell-RNA-Seq-Analysis/main/data_analysis/05-clustering_celltype.Rmd", "05-clustering_celltype.Rmd")
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
##  [19] promises_1.3.0         rmarkdown_2.27         purrr_1.0.2           
##  [22] xfun_0.44              cachem_1.1.0           jsonlite_1.8.8        
##  [25] goftest_1.2-3          highr_0.11             later_1.3.2           
##  [28] spatstat.utils_3.1-2   irlba_2.3.5.1          parallel_4.4.3        
##  [31] cluster_2.1.8          R6_2.5.1               ica_1.0-3             
##  [34] bslib_0.7.0            stringi_1.8.4          RColorBrewer_1.1-3    
##  [37] spatstat.data_3.0-4    reticulate_1.39.0      parallelly_1.37.1     
##  [40] lmtest_0.9-40          jquerylib_0.1.4        scattermore_1.2       
##  [43] Rcpp_1.0.12            knitr_1.47             tensor_1.5            
##  [46] future.apply_1.11.2    zoo_1.8-12             sctransform_0.4.1     
##  [49] httpuv_1.6.15          Matrix_1.7-2           splines_4.4.3         
##  [52] igraph_2.0.3           tidyselect_1.2.1       abind_1.4-5           
##  [55] rstudioapi_0.16.0      yaml_2.3.8             spatstat.random_3.2-3 
##  [58] codetools_0.2-20       miniUI_0.1.1.1         spatstat.explore_3.2-7
##  [61] listenv_0.9.1          lattice_0.22-6         tibble_3.2.1          
##  [64] plyr_1.8.9             withr_3.0.0            shiny_1.8.1.1         
##  [67] ROCR_1.0-11            evaluate_0.23          Rtsne_0.17            
##  [70] future_1.33.2          fastDummies_1.7.3      survival_3.8-3        
##  [73] polyclip_1.10-6        fitdistrplus_1.1-11    pillar_1.9.0          
##  [76] KernSmooth_2.23-26     plotly_4.10.4          generics_0.1.3        
##  [79] RcppHNSW_0.6.0         munsell_0.5.1          scales_1.3.0          
##  [82] globals_0.16.3         xtable_1.8-4           glue_1.7.0            
##  [85] lazyeval_0.2.2         tools_4.4.3            data.table_1.15.4     
##  [88] RSpectra_0.16-1        RANN_2.6.1             dotCall64_1.1-1       
##  [91] cowplot_1.1.3          grid_4.4.3             tidyr_1.3.1           
##  [94] colorspace_2.1-0       nlme_3.1-167           patchwork_1.2.0       
##  [97] cli_3.6.2              spatstat.sparse_3.0-3  spam_2.10-0           
## [100] fansi_1.0.6            viridisLite_0.4.2      dplyr_1.1.4           
## [103] uwot_0.2.2             gtable_0.3.5           sass_0.4.9            
## [106] digest_0.6.35          progressr_0.14.0       ggrepel_0.9.5         
## [109] farver_2.1.2           htmlwidgets_1.6.4      htmltools_0.5.8.1     
## [112] lifecycle_1.0.4        httr_1.4.7             mime_0.12             
## [115] MASS_7.3-64
```
