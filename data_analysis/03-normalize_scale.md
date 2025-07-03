---
title: "Introduction to Single Cell RNA-Seq Part 3: Normalize and scale"
author: "UCD Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

# Introduction to Single Cell RNA-Seq Part 3: Normalize and scale


## Set up workspace
First, we need to load the required libraries.

``` r
library(Seurat)
library(kableExtra)
experiment.aggregate <- readRDS("scRNA_workshop-02.rds")
experiment.aggregate
```

```
## An object of class Seurat 
## 11474 features across 6315 samples within 1 assay 
## Active assay: RNA (11474 features, 0 variable features)
##  1 layer present: counts
```

``` r
set.seed(12345)
```

## Normalize the data
At the end of filtering step, the data should be at good quality for analysis. The next step is to normalize the data. We employ a global-scaling normalization method, LogNormalize, that normalizes the gene expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and then log-transforms the data.

``` r
?NormalizeData
```


``` r
experiment.aggregate <- NormalizeData(
  object = experiment.aggregate,
  normalization.method = "LogNormalize",
  scale.factor = 10000)
```

## Cell cycle assignment
Cell cycle phase can be a significant source of variation in single cell and single nucleus experiments. There are a number of automated cell cycle stage detection methods available for single cell data. For this workshop, we will be using the built-in Seurat cell cycle function, `CellCycleScoring`. This tool compares gene expression in each cell to a list of cell cycle marker genes and scores each barcode based on marker expression. The phase with the highest score is selected for each barcode. Seurat includes a list of cell cycle genes in human single cell data.

``` r
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
```

For other species, a user-provided gene list may be substituted, or the orthologs of the human gene list used instead.

**Do not run the code below for human experiments!**

``` r
# mouse code DO NOT RUN for human data
library(biomaRt)
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useEnsembl("ensembl",
                     dataset = "hsapiens_gene_ensembl",
                     mirror = "uswest")
  mouse = useEnsembl("ensembl",
                     dataset = "mmusculus_gene_ensembl",
                     mirror = "uswest")
  genes = getLDS(attributes = c("hgnc_symbol"),
                 filters = "hgnc_symbol",
                 values = x ,
                 mart = human,
                 attributesL = c("mgi_symbol"),
                 martL = mouse,
                 uniqueRows=T)
  humanx = unique(genes[, 2])
  print(head(humanx)) # print first 6 genes found to the screen
  return(humanx)
}
# convert lists to mouse orthologs
s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
```

Once an appropriate gene list has been identified, the `CellCycleScoring` function can be run.

``` r
experiment.aggregate <- CellCycleScoring(experiment.aggregate,
                                         s.features = s.genes,
                                         g2m.features = g2m.genes,
                                         set.ident = TRUE)
table(experiment.aggregate@meta.data$Phase) %>%
  kable(caption = "Number of Cells in each Cell Cycle Stage",
        col.names = c("Stage", "Count"),
        align = "c") %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Number of Cells in each Cell Cycle Stage</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> Stage </th>
   <th style="text-align:center;"> Count </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> G1 </td>
   <td style="text-align:center;"> 3803 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> G2M </td>
   <td style="text-align:center;"> 1169 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> S </td>
   <td style="text-align:center;"> 1343 </td>
  </tr>
</tbody>
</table>

Because the "set.ident" argument was set to TRUE (this is also the default behavior), the active identity of the Seurat object was changed to the phase. To return the active identity to the sample identity, use the `Idents` function.

``` r
table(Idents(experiment.aggregate))
```

```
## 
##    S  G2M   G1 
## 1343 1169 3803
```

``` r
Idents(experiment.aggregate) <- "orig.ident"
table(Idents(experiment.aggregate))
```

```
## 
## A001-C-007 A001-C-104 B001-A-301 
##       1037       1848       3430
```

## Identify variable genes
The function FindVariableFeatures identifies the most highly variable genes (default 2000 genes) by fitting a line to the relationship of log(variance) and log(mean) using loess smoothing, uses this information to standardize the data, then calculates the variance of the standardized data.  This helps avoid selecting genes that only appear variable due to their expression level.

``` r
?FindVariableFeatures
```


``` r
experiment.aggregate <- FindVariableFeatures(
  object = experiment.aggregate,
  selection.method = "vst")
length(VariableFeatures(experiment.aggregate))
```

```
## [1] 2000
```

``` r
top10 <- head(VariableFeatures(experiment.aggregate), 10)
top10
```

```
##  [1] "BEST4"           "CLCA4"           "TPH1"            "SMOC2"          
##  [5] "NRG1"            "TRPM3"           "IRAG2"           "PTPRR"          
##  [9] "ENSG00000260573" "CACNA1A"
```

``` r
var.feat.plot <- VariableFeaturePlot(experiment.aggregate)
var.feat.plot <- LabelPoints(plot = var.feat.plot, points = top10, repel = TRUE)
var.feat.plot
```

![](03-normalize_scale_files/figure-html/find_variable_genes-1.png)<!-- -->

**How do the results change if you use selection.method = "dispersion" or selection.method = "mean.var.plot"?**

FindVariableFeatures isn't the only way to set the "variable features" of a Seurat object. Another reasonable approach is to select a set of "minimally expressed" genes.

``` r
min.value <- 2
min.cells <- 10

num.cells <- Matrix::rowSums(GetAssayData(experiment.aggregate, slot = "count") > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])
length(genes.use)
```

```
## [1] 7110
```

``` r
VariableFeatures(experiment.aggregate) <- genes.use
```

## Scale the data
The `ScaleData` function scales and centers genes in the dataset. If variables are provided with the "vars.to.regress" argument, they are individually regressed against each gene, and the resulting residuals are then scaled and centered unless otherwise specified. We regress out cell cycle results S.Score and G2M.Score, mitochondrial RNA level (percent_MT), and the number of features (nFeature_RNA) as a proxy for sequencing depth.

``` r
experiment.aggregate <- ScaleData(experiment.aggregate,
                                  vars.to.regress = c("S.Score", "G2M.Score", "percent_MT", "nFeature_RNA"))
```

## Prepare for the next section

#### Save object

``` r
saveRDS(experiment.aggregate, file = "scRNA_workshop-03.rds")
```

#### Download Rmd

``` r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2025-July-Single-Cell-RNA-Seq-Analysis/main/data_analysis/04-dimensionality_reduction.Rmd", "04-dimensionality_reduction.Rmd")
```

#### Session Information

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
## [1] kableExtra_1.4.0   Seurat_5.2.1       SeuratObject_5.0.2 sp_2.1-4          
## 
## loaded via a namespace (and not attached):
##   [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3         
##   [4] rlang_1.1.3            magrittr_2.0.3         RcppAnnoy_0.0.22      
##   [7] spatstat.geom_3.2-9    matrixStats_1.3.0      ggridges_0.5.6        
##  [10] compiler_4.4.3         systemfonts_1.1.0      png_0.1-8             
##  [13] vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1         
##  [16] pkgconfig_2.0.3        fastmap_1.2.0          labeling_0.4.3        
##  [19] utf8_1.2.4             promises_1.3.0         rmarkdown_2.27        
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
##  [67] shiny_1.8.1.1          ROCR_1.0-11            evaluate_0.23         
##  [70] Rtsne_0.17             future_1.33.2          fastDummies_1.7.3     
##  [73] survival_3.8-3         polyclip_1.10-6        xml2_1.3.6            
##  [76] fitdistrplus_1.1-11    pillar_1.9.0           KernSmooth_2.23-26    
##  [79] plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0        
##  [82] ggplot2_3.5.1          munsell_0.5.1          scales_1.3.0          
##  [85] globals_0.16.3         xtable_1.8-4           glue_1.7.0            
##  [88] lazyeval_0.2.2         tools_4.4.3            data.table_1.15.4     
##  [91] RSpectra_0.16-1        RANN_2.6.1             dotCall64_1.1-1       
##  [94] cowplot_1.1.3          grid_4.4.3             tidyr_1.3.1           
##  [97] colorspace_2.1-0       nlme_3.1-167           patchwork_1.2.0       
## [100] cli_3.6.2              spatstat.sparse_3.0-3  spam_2.10-0           
## [103] fansi_1.0.6            viridisLite_0.4.2      svglite_2.1.3         
## [106] dplyr_1.1.4            uwot_0.2.2             gtable_0.3.5          
## [109] sass_0.4.9             digest_0.6.35          progressr_0.14.0      
## [112] ggrepel_0.9.5          farver_2.1.2           htmlwidgets_1.6.4     
## [115] htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7            
## [118] mime_0.12              MASS_7.3-64
```
