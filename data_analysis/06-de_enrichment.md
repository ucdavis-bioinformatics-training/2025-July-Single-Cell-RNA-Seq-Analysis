---
title: "Introduction to Single Cell RNA-Seq Part 6: Enrichment and model-based differential expression"
author: "Bioinformatics Core"
date: "2025-07-10"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Introduction to Single Cell RNA-Seq Part 6: Enrichment and model-based differential expression


## Set up workspace

``` r
library(Seurat)
library(limma)
library(topGO)
library(dplyr)
library(kableExtra)
set.seed(12345)
experiment.aggregate <- readRDS("scRNA_workshop-05.rds")
Idents(experiment.aggregate) <- "finalcluster"
```

## 1. Gene Ontology (GO) Enrichment of Genes Expressed in a Cluster
[Gene Ontology](http://geneontology.org/docs/ontology-documentation/) provides a controlled vocabulary for describing gene products.  Here we use enrichment analysis to identify GO terms that are over-represented among the gene expressed in cells in a given cluster. 

``` r
cluster10 <- subset(experiment.aggregate, idents = '10')
expr <- as.matrix(GetAssayData(cluster10))

# Select genes that are expressed > 0 in at least half of cells
n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.5)]
all.genes <- rownames(expr)

# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
	GOdata <- new("topGOdata",
		ontology = "BP", # use biological process ontology
		allGenes = geneList,
		geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
# Test for enrichment using Fisher's Exact Test
	resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
	GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> GO.ID </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Term </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Annotated </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Significant </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Expected </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Fisher </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GO:0000381 </td>
   <td style="text-align:left;"> regulation of alternative mRNA splicing, via spliceosome </td>
   <td style="text-align:right;"> 41 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 1.85 </td>
   <td style="text-align:left;"> 6.9e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006338 </td>
   <td style="text-align:left;"> chromatin remodeling </td>
   <td style="text-align:right;"> 918 </td>
   <td style="text-align:right;"> 66 </td>
   <td style="text-align:right;"> 41.41 </td>
   <td style="text-align:left;"> 0.00028 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0035176 </td>
   <td style="text-align:left;"> social behavior </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.99 </td>
   <td style="text-align:left;"> 0.00033 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:2000650 </td>
   <td style="text-align:left;"> negative regulation of sodium ion transmembrane transporter ... </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.18 </td>
   <td style="text-align:left;"> 0.00035 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0051660 </td>
   <td style="text-align:left;"> establishment of centrosome localization </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.41 </td>
   <td style="text-align:left;"> 0.00043 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006895 </td>
   <td style="text-align:left;"> Golgi to endosome transport </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.77 </td>
   <td style="text-align:left;"> 0.00072 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0045053 </td>
   <td style="text-align:left;"> protein retention in Golgi apparatus </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.23 </td>
   <td style="text-align:left;"> 0.00085 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0032534 </td>
   <td style="text-align:left;"> regulation of microvillus assembly </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.23 </td>
   <td style="text-align:left;"> 0.00085 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0030036 </td>
   <td style="text-align:left;"> actin cytoskeleton organization </td>
   <td style="text-align:right;"> 482 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 21.74 </td>
   <td style="text-align:left;"> 0.00093 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0044331 </td>
   <td style="text-align:left;"> cell-cell adhesion mediated by cadherin </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1.22 </td>
   <td style="text-align:left;"> 0.00107 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0009653 </td>
   <td style="text-align:left;"> anatomical structure morphogenesis </td>
   <td style="text-align:right;"> 1541 </td>
   <td style="text-align:right;"> 107 </td>
   <td style="text-align:right;"> 69.52 </td>
   <td style="text-align:left;"> 0.00111 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0045944 </td>
   <td style="text-align:left;"> positive regulation of transcription by RNA polymerase II </td>
   <td style="text-align:right;"> 825 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 37.22 </td>
   <td style="text-align:left;"> 0.00118 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:1902275 </td>
   <td style="text-align:left;"> regulation of chromatin organization </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.90 </td>
   <td style="text-align:left;"> 0.00161 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0042060 </td>
   <td style="text-align:left;"> wound healing </td>
   <td style="text-align:right;"> 261 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 11.77 </td>
   <td style="text-align:left;"> 0.00162 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0007043 </td>
   <td style="text-align:left;"> cell-cell junction assembly </td>
   <td style="text-align:right;"> 91 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 4.11 </td>
   <td style="text-align:left;"> 0.00172 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0090257 </td>
   <td style="text-align:left;"> regulation of muscle system process </td>
   <td style="text-align:right;"> 128 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 5.77 </td>
   <td style="text-align:left;"> 0.00184 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0031175 </td>
   <td style="text-align:left;"> neuron projection development </td>
   <td style="text-align:right;"> 591 </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 26.66 </td>
   <td style="text-align:left;"> 0.00185 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0070830 </td>
   <td style="text-align:left;"> bicellular tight junction assembly </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1.80 </td>
   <td style="text-align:left;"> 0.00185 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0048669 </td>
   <td style="text-align:left;"> collateral sprouting in absence of injury </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:left;"> 0.00203 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0051674 </td>
   <td style="text-align:left;"> localization of cell </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:left;"> 0.00203 </td>
  </tr>
</tbody>
</table>
* Annotated: number of genes (out of all.genes) that are annotated with that GO term
* Significant: number of genes that are annotated with that GO term and meet our criteria for "expressed"
* Expected: Under random chance, number of genes that would be expected to be annotated with that GO term and meeting our criteria for "expressed"
* Fisher: (Raw) p-value from Fisher's Exact Test

### Organisms other than human
Many organisms have a database equivalent to org.Hs.eg.db on Bioconductor: https://www.bioconductor.org/packages/release/data/annotation/

If not, and your organism is on NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) it's straightforward to create your own using the AnnotationForge R package and the function makeOrgPackageFromNCBI (https://bioconductor.org/packages/devel/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html)

## 2. Model-based DE analysis in limma
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html) is an R package for differential expression analysis of bulk RNASeq and microarray data.  We apply it here to single cell data.

Limma can be used to fit any linear model to expression data and is useful for analyses that go beyond two-group comparisons.  A detailed tutorial of model specification in limma is available [here](https://ucdavis-bioinformatics-training.github.io/2025-June-RNA-Seq-Analysis/data_analysis/DE_Analysis_mm_with_quizzes) and in the [limma User's Guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf).


``` r
# filter genes to those expressed in at least 10% of cells
keep <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.1)]
expr2 <- expr[keep,]

# Set up "design matrix" with statistical model
cluster10$proper.group <- make.names(cluster10$group)
mm <- model.matrix(~0 + proper.group + S.Score + G2M.Score + percent_MT + nFeature_RNA, data = cluster10[[]])
head(mm) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupColorectal.Cancer </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupPolyp </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent_MT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> AAACCCAAGTTATGGA_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0554074 </td>
   <td style="text-align:right;"> -0.1461764 </td>
   <td style="text-align:right;"> 0.5717008 </td>
   <td style="text-align:right;"> 1513 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAGCCATCAAGACCTT_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0660183 </td>
   <td style="text-align:right;"> -0.0729290 </td>
   <td style="text-align:right;"> 0.4841313 </td>
   <td style="text-align:right;"> 1265 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAGTTCGGTACCTATG_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0202042 </td>
   <td style="text-align:right;"> -0.0451182 </td>
   <td style="text-align:right;"> 0.5148005 </td>
   <td style="text-align:right;"> 1148 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AATGAAGTCAGCGTCG_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.1225441 </td>
   <td style="text-align:right;"> 0.0612193 </td>
   <td style="text-align:right;"> 1.0398614 </td>
   <td style="text-align:right;"> 1272 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACAGAAATCCAGCTCT_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0769269 </td>
   <td style="text-align:right;"> -0.0707690 </td>
   <td style="text-align:right;"> 0.6334125 </td>
   <td style="text-align:right;"> 1637 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACGGTTACAAATCGGG_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0550882 </td>
   <td style="text-align:right;"> 0.1286419 </td>
   <td style="text-align:right;"> 0.7523716 </td>
   <td style="text-align:right;"> 1934 </td>
  </tr>
</tbody>
</table>

``` r
tail(mm) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupColorectal.Cancer </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupPolyp </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent_MT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> TTCCACGAGTATTGCC_A001-C-104 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.4520556 </td>
   <td style="text-align:right;"> 0.0266469 </td>
   <td style="text-align:right;"> 0.3015580 </td>
   <td style="text-align:right;"> 4562 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTCTAGTAGATAACAC_A001-C-104 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0218701 </td>
   <td style="text-align:right;"> -0.0793631 </td>
   <td style="text-align:right;"> 1.3125000 </td>
   <td style="text-align:right;"> 1140 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTGCCTGTCCGCGATG_A001-C-104 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0612575 </td>
   <td style="text-align:right;"> -0.0292726 </td>
   <td style="text-align:right;"> 0.3594352 </td>
   <td style="text-align:right;"> 2241 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTGGGATTCGTTCCCA_A001-C-104 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0709637 </td>
   <td style="text-align:right;"> -0.0806137 </td>
   <td style="text-align:right;"> 0.4050223 </td>
   <td style="text-align:right;"> 1562 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTTCAGTTCGCACTCT_A001-C-104 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.1047877 </td>
   <td style="text-align:right;"> -0.0540112 </td>
   <td style="text-align:right;"> 1.8099548 </td>
   <td style="text-align:right;"> 1309 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GACTCAACACACACGC_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0481511 </td>
   <td style="text-align:right;"> -0.0900060 </td>
   <td style="text-align:right;"> 0.1574803 </td>
   <td style="text-align:right;"> 990 </td>
  </tr>
</tbody>
</table>

``` r
# Fit model in limma
fit <- lmFit(expr2, mm)
head(coef(fit)) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupColorectal.Cancer </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupPolyp </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent_MT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000291215 </td>
   <td style="text-align:right;"> -0.2357714 </td>
   <td style="text-align:right;"> -0.2248263 </td>
   <td style="text-align:right;"> -0.0954233 </td>
   <td style="text-align:right;"> 0.5209505 </td>
   <td style="text-align:right;"> -0.3817891 </td>
   <td style="text-align:right;"> 0.0510813 </td>
   <td style="text-align:right;"> 0.0001589 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C1orf159 </td>
   <td style="text-align:right;"> -0.0141535 </td>
   <td style="text-align:right;"> -0.1325994 </td>
   <td style="text-align:right;"> -0.0259810 </td>
   <td style="text-align:right;"> -0.1442860 </td>
   <td style="text-align:right;"> 0.2517446 </td>
   <td style="text-align:right;"> 0.1333100 </td>
   <td style="text-align:right;"> 0.0001426 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SDF4 </td>
   <td style="text-align:right;"> 0.3616255 </td>
   <td style="text-align:right;"> 0.1839501 </td>
   <td style="text-align:right;"> 0.4912284 </td>
   <td style="text-align:right;"> 0.2506248 </td>
   <td style="text-align:right;"> 1.0914474 </td>
   <td style="text-align:right;"> -0.1494113 </td>
   <td style="text-align:right;"> -0.0000750 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CCNL2 </td>
   <td style="text-align:right;"> 0.4402786 </td>
   <td style="text-align:right;"> 1.6154478 </td>
   <td style="text-align:right;"> 0.0118319 </td>
   <td style="text-align:right;"> 0.4012170 </td>
   <td style="text-align:right;"> -1.8906665 </td>
   <td style="text-align:right;"> 0.2768585 </td>
   <td style="text-align:right;"> 0.0003573 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MIB2 </td>
   <td style="text-align:right;"> 0.0480968 </td>
   <td style="text-align:right;"> -0.0562502 </td>
   <td style="text-align:right;"> 0.2279206 </td>
   <td style="text-align:right;"> -0.1502399 </td>
   <td style="text-align:right;"> -0.4402191 </td>
   <td style="text-align:right;"> -0.1676190 </td>
   <td style="text-align:right;"> 0.0000508 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CDK11B </td>
   <td style="text-align:right;"> 0.6889952 </td>
   <td style="text-align:right;"> -0.1659776 </td>
   <td style="text-align:right;"> -0.0648067 </td>
   <td style="text-align:right;"> -0.7018856 </td>
   <td style="text-align:right;"> -0.2333263 </td>
   <td style="text-align:right;"> 0.1621756 </td>
   <td style="text-align:right;"> 0.0001548 </td>
  </tr>
</tbody>
</table>

``` r
# Test 'Normal' - 'Colorectal.Cancer'
contr <- makeContrasts(proper.groupNormal - proper.groupColorectal.Cancer, levels = colnames(coef(fit)))
contr %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal - proper.groupColorectal.Cancer </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> proper.groupColorectal.Cancer </td>
   <td style="text-align:right;"> -1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> proper.groupNormal </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> proper.groupPolyp </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S.Score </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G2M.Score </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> percent_MT </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nFeature_RNA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

``` r
fit2 <- contrasts.fit(fit, contrasts = contr)
fit2 <- eBayes(fit2)
out <- topTable(fit2, n = Inf, sort.by = "P")
head(out, 30) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> logFC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> AveExpr </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> t </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> P.Value </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> adj.P.Val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> B </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> UTP6 </td>
   <td style="text-align:right;"> 2.809558 </td>
   <td style="text-align:right;"> 0.1665493 </td>
   <td style="text-align:right;"> 5.729581 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0002398 </td>
   <td style="text-align:right;"> 7.9550347 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S100PBP </td>
   <td style="text-align:right;"> 2.764765 </td>
   <td style="text-align:right;"> 0.2078133 </td>
   <td style="text-align:right;"> 4.898211 </td>
   <td style="text-align:right;"> 0.0000024 </td>
   <td style="text-align:right;"> 0.0057367 </td>
   <td style="text-align:right;"> 4.4838247 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FDPS </td>
   <td style="text-align:right;"> 2.194920 </td>
   <td style="text-align:right;"> 0.1608721 </td>
   <td style="text-align:right;"> 4.481281 </td>
   <td style="text-align:right;"> 0.0000142 </td>
   <td style="text-align:right;"> 0.0127450 </td>
   <td style="text-align:right;"> 2.8889818 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CRKL </td>
   <td style="text-align:right;"> 2.062359 </td>
   <td style="text-align:right;"> 0.1502419 </td>
   <td style="text-align:right;"> 4.431309 </td>
   <td style="text-align:right;"> 0.0000175 </td>
   <td style="text-align:right;"> 0.0127450 </td>
   <td style="text-align:right;"> 2.7049563 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NMRAL1 </td>
   <td style="text-align:right;"> 2.284191 </td>
   <td style="text-align:right;"> 0.1670087 </td>
   <td style="text-align:right;"> 4.430493 </td>
   <td style="text-align:right;"> 0.0000175 </td>
   <td style="text-align:right;"> 0.0127450 </td>
   <td style="text-align:right;"> 2.7019633 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SOAT1 </td>
   <td style="text-align:right;"> 2.026978 </td>
   <td style="text-align:right;"> 0.1554562 </td>
   <td style="text-align:right;"> 4.415095 </td>
   <td style="text-align:right;"> 0.0000187 </td>
   <td style="text-align:right;"> 0.0127450 </td>
   <td style="text-align:right;"> 2.6455817 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CMSS1 </td>
   <td style="text-align:right;"> 2.179765 </td>
   <td style="text-align:right;"> 0.1843519 </td>
   <td style="text-align:right;"> 4.398089 </td>
   <td style="text-align:right;"> 0.0000200 </td>
   <td style="text-align:right;"> 0.0127450 </td>
   <td style="text-align:right;"> 2.5834884 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TXNL4A </td>
   <td style="text-align:right;"> 2.044588 </td>
   <td style="text-align:right;"> 0.1549221 </td>
   <td style="text-align:right;"> 4.379329 </td>
   <td style="text-align:right;"> 0.0000216 </td>
   <td style="text-align:right;"> 0.0127450 </td>
   <td style="text-align:right;"> 2.5152036 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PROX1 </td>
   <td style="text-align:right;"> 2.223377 </td>
   <td style="text-align:right;"> 0.1608351 </td>
   <td style="text-align:right;"> 4.344625 </td>
   <td style="text-align:right;"> 0.0000249 </td>
   <td style="text-align:right;"> 0.0127450 </td>
   <td style="text-align:right;"> 2.3894760 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ADAMTSL1 </td>
   <td style="text-align:right;"> 3.378881 </td>
   <td style="text-align:right;"> 0.4051661 </td>
   <td style="text-align:right;"> 4.315871 </td>
   <td style="text-align:right;"> 0.0000280 </td>
   <td style="text-align:right;"> 0.0127450 </td>
   <td style="text-align:right;"> 2.2858848 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LNX2 </td>
   <td style="text-align:right;"> 2.802389 </td>
   <td style="text-align:right;"> 0.3089559 </td>
   <td style="text-align:right;"> 4.306615 </td>
   <td style="text-align:right;"> 0.0000290 </td>
   <td style="text-align:right;"> 0.0127450 </td>
   <td style="text-align:right;"> 2.2526522 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SEC23B </td>
   <td style="text-align:right;"> 2.024936 </td>
   <td style="text-align:right;"> 0.1510704 </td>
   <td style="text-align:right;"> 4.256737 </td>
   <td style="text-align:right;"> 0.0000355 </td>
   <td style="text-align:right;"> 0.0142802 </td>
   <td style="text-align:right;"> 2.0745228 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCP11L1 </td>
   <td style="text-align:right;"> 2.194426 </td>
   <td style="text-align:right;"> 0.1880198 </td>
   <td style="text-align:right;"> 4.112117 </td>
   <td style="text-align:right;"> 0.0000629 </td>
   <td style="text-align:right;"> 0.0218800 </td>
   <td style="text-align:right;"> 1.5672384 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNF3 </td>
   <td style="text-align:right;"> 2.241657 </td>
   <td style="text-align:right;"> 0.1841910 </td>
   <td style="text-align:right;"> 4.097439 </td>
   <td style="text-align:right;"> 0.0000666 </td>
   <td style="text-align:right;"> 0.0218800 </td>
   <td style="text-align:right;"> 1.5165269 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBL1 </td>
   <td style="text-align:right;"> 2.107244 </td>
   <td style="text-align:right;"> 0.1733018 </td>
   <td style="text-align:right;"> 4.092354 </td>
   <td style="text-align:right;"> 0.0000680 </td>
   <td style="text-align:right;"> 0.0218800 </td>
   <td style="text-align:right;"> 1.4989902 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000287292 </td>
   <td style="text-align:right;"> 2.649530 </td>
   <td style="text-align:right;"> 0.2658202 </td>
   <td style="text-align:right;"> 4.044865 </td>
   <td style="text-align:right;"> 0.0000817 </td>
   <td style="text-align:right;"> 0.0231668 </td>
   <td style="text-align:right;"> 1.3360639 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WDTC1 </td>
   <td style="text-align:right;"> 2.184789 </td>
   <td style="text-align:right;"> 0.1787432 </td>
   <td style="text-align:right;"> 4.032124 </td>
   <td style="text-align:right;"> 0.0000858 </td>
   <td style="text-align:right;"> 0.0231668 </td>
   <td style="text-align:right;"> 1.2926107 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MED17 </td>
   <td style="text-align:right;"> 2.135271 </td>
   <td style="text-align:right;"> 0.1824321 </td>
   <td style="text-align:right;"> 4.009965 </td>
   <td style="text-align:right;"> 0.0000935 </td>
   <td style="text-align:right;"> 0.0231668 </td>
   <td style="text-align:right;"> 1.2172980 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RXRA </td>
   <td style="text-align:right;"> 2.060384 </td>
   <td style="text-align:right;"> 0.1753848 </td>
   <td style="text-align:right;"> 4.005170 </td>
   <td style="text-align:right;"> 0.0000952 </td>
   <td style="text-align:right;"> 0.0231668 </td>
   <td style="text-align:right;"> 1.2010447 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CDK18 </td>
   <td style="text-align:right;"> 2.168277 </td>
   <td style="text-align:right;"> 0.1831622 </td>
   <td style="text-align:right;"> 4.003186 </td>
   <td style="text-align:right;"> 0.0000959 </td>
   <td style="text-align:right;"> 0.0231668 </td>
   <td style="text-align:right;"> 1.1943236 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ADPGK </td>
   <td style="text-align:right;"> 2.058125 </td>
   <td style="text-align:right;"> 0.1664844 </td>
   <td style="text-align:right;"> 3.989304 </td>
   <td style="text-align:right;"> 0.0001012 </td>
   <td style="text-align:right;"> 0.0232690 </td>
   <td style="text-align:right;"> 1.1473755 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNF480 </td>
   <td style="text-align:right;"> 2.236793 </td>
   <td style="text-align:right;"> 0.2006418 </td>
   <td style="text-align:right;"> 3.965854 </td>
   <td style="text-align:right;"> 0.0001107 </td>
   <td style="text-align:right;"> 0.0242922 </td>
   <td style="text-align:right;"> 1.0683681 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPIL4 </td>
   <td style="text-align:right;"> 2.226664 </td>
   <td style="text-align:right;"> 0.2090412 </td>
   <td style="text-align:right;"> 3.892616 </td>
   <td style="text-align:right;"> 0.0001460 </td>
   <td style="text-align:right;"> 0.0306594 </td>
   <td style="text-align:right;"> 0.8240318 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IRF3 </td>
   <td style="text-align:right;"> 2.687759 </td>
   <td style="text-align:right;"> 0.3212877 </td>
   <td style="text-align:right;"> 3.873166 </td>
   <td style="text-align:right;"> 0.0001571 </td>
   <td style="text-align:right;"> 0.0316070 </td>
   <td style="text-align:right;"> 0.7597636 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SARS1 </td>
   <td style="text-align:right;"> 2.368320 </td>
   <td style="text-align:right;"> 0.2254349 </td>
   <td style="text-align:right;"> 3.835478 </td>
   <td style="text-align:right;"> 0.0001808 </td>
   <td style="text-align:right;"> 0.0346093 </td>
   <td style="text-align:right;"> 0.6359795 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TMEM51-AS1 </td>
   <td style="text-align:right;"> 1.994662 </td>
   <td style="text-align:right;"> 0.1726518 </td>
   <td style="text-align:right;"> 3.808954 </td>
   <td style="text-align:right;"> 0.0001995 </td>
   <td style="text-align:right;"> 0.0346093 </td>
   <td style="text-align:right;"> 0.5494559 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PRAG1 </td>
   <td style="text-align:right;"> 2.223802 </td>
   <td style="text-align:right;"> 0.2200442 </td>
   <td style="text-align:right;"> 3.793651 </td>
   <td style="text-align:right;"> 0.0002111 </td>
   <td style="text-align:right;"> 0.0346093 </td>
   <td style="text-align:right;"> 0.4997596 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNF13 </td>
   <td style="text-align:right;"> 2.212740 </td>
   <td style="text-align:right;"> 0.2147929 </td>
   <td style="text-align:right;"> 3.778006 </td>
   <td style="text-align:right;"> 0.0002237 </td>
   <td style="text-align:right;"> 0.0346093 </td>
   <td style="text-align:right;"> 0.4491215 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KALRN </td>
   <td style="text-align:right;"> 2.052896 </td>
   <td style="text-align:right;"> 0.1794657 </td>
   <td style="text-align:right;"> 3.775951 </td>
   <td style="text-align:right;"> 0.0002254 </td>
   <td style="text-align:right;"> 0.0346093 </td>
   <td style="text-align:right;"> 0.4424844 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST6GAL1 </td>
   <td style="text-align:right;"> 2.169592 </td>
   <td style="text-align:right;"> 0.1866848 </td>
   <td style="text-align:right;"> 3.770198 </td>
   <td style="text-align:right;"> 0.0002302 </td>
   <td style="text-align:right;"> 0.0346093 </td>
   <td style="text-align:right;"> 0.4239141 </td>
  </tr>
</tbody>
</table>

**Output columns:**

* logFC: log fold change (since we are working with Seurat's natural log transformed data, will be natural log fold change)
* AveExpr: Average expression across all cells in expr2
* t: logFC divided by its standard error
* P.Value: Raw p-value (based on t) from test that logFC differs from 0
* adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
* B: log-odds that gene is DE 

## Save files

``` r
write.csv(GenTable(GOdata, Fisher = resultFisher), file = "cluster10_GOdata.csv")
write.csv(out, file = "cluster10_Normal-Colorectal.Cancer_topTable.csv")
```

## A note on pseudobulk DE

Pseudobulk differential expression uses count data summed across all cells in each sample (typically within each cell type or cluster).  Unlike cell-level DE, pseudobulk DE *requires biological replicates* so we won't perform it on this dataset.

Once counts are summed, pseudobulk data are analyzed like bulk RNASeq data.

Pseudobulk DE may result in better false discovery rate control than cell-level DE, as shown [here](https://www.nature.com/articles/s41467-021-25960-2).

The Seurat function `AggregateExpression()` can be used to sum counts as described [here](https://satijalab.org/seurat/articles/de_vignette).

A tutorial on using limma for bulk RNASeq is available [here](https://ucdavis-bioinformatics-training.github.io/2025-June-RNA-Seq-Analysis/data_analysis/DE_Analysis_mm_with_quizzes).

## Prepare for the next section

#### Download Rmd document

``` r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2025-July-Single-Cell-RNA-Seq-Analysis/main/data_analysis/07-doublet_detection.Rmd", "07-doublet_detection.Rmd")
```

#### Session Information

``` r
sessionInfo()
```

```
## R version 4.5.0 (2025-04-11 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 26100)
## 
## Matrix products: default
##   LAPACK version 3.12.1
## 
## locale:
## [1] LC_COLLATE=English_United States.utf8 
## [2] LC_CTYPE=English_United States.utf8   
## [3] LC_MONETARY=English_United States.utf8
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.utf8    
## 
## time zone: America/Los_Angeles
## tzcode source: internal
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] org.Hs.eg.db_3.21.0  kableExtra_1.4.0     dplyr_1.1.4         
##  [4] topGO_2.60.1         SparseM_1.84-2       GO.db_3.21.0        
##  [7] AnnotationDbi_1.70.0 IRanges_2.42.0       S4Vectors_0.46.0    
## [10] Biobase_2.68.0       graph_1.86.0         BiocGenerics_0.54.0 
## [13] generics_0.1.4       limma_3.64.1         Seurat_5.3.0        
## [16] SeuratObject_5.1.0   sp_2.2-0            
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_2.0.0         
##   [4] magrittr_2.0.3          spatstat.utils_3.1-4    farver_2.1.2           
##   [7] rmarkdown_2.29          vctrs_0.6.5             ROCR_1.0-11            
##  [10] memoise_2.0.1           spatstat.explore_3.4-3  htmltools_0.5.8.1      
##  [13] sass_0.4.10             sctransform_0.4.2       parallelly_1.45.0      
##  [16] KernSmooth_2.23-26      bslib_0.9.0             htmlwidgets_1.6.4      
##  [19] ica_1.0-3               plyr_1.8.9              plotly_4.11.0          
##  [22] zoo_1.8-14              cachem_1.1.0            igraph_2.1.4           
##  [25] mime_0.13               lifecycle_1.0.4         pkgconfig_2.0.3        
##  [28] Matrix_1.7-3            R6_2.6.1                fastmap_1.2.0          
##  [31] GenomeInfoDbData_1.2.14 fitdistrplus_1.2-4      future_1.58.0          
##  [34] shiny_1.11.1            digest_0.6.37           patchwork_1.3.1        
##  [37] tensor_1.5.1            RSpectra_0.16-2         irlba_2.3.5.1          
##  [40] textshaping_1.0.1       RSQLite_2.4.1           progressr_0.15.1       
##  [43] spatstat.sparse_3.1-0   httr_1.4.7              polyclip_1.10-7        
##  [46] abind_1.4-8             compiler_4.5.0          bit64_4.6.0-1          
##  [49] DBI_1.2.3               fastDummies_1.7.5       MASS_7.3-65            
##  [52] tools_4.5.0             lmtest_0.9-40           httpuv_1.6.16          
##  [55] future.apply_1.20.0     goftest_1.2-3           glue_1.8.0             
##  [58] nlme_3.1-168            promises_1.3.3          grid_4.5.0             
##  [61] Rtsne_0.17              cluster_2.1.8.1         reshape2_1.4.4         
##  [64] gtable_0.3.6            spatstat.data_3.1-6     tidyr_1.3.1            
##  [67] data.table_1.17.6       xml2_1.3.8              XVector_0.48.0         
##  [70] spatstat.geom_3.4-1     RcppAnnoy_0.0.22        ggrepel_0.9.6          
##  [73] RANN_2.6.2              pillar_1.11.0           stringr_1.5.1          
##  [76] spam_2.11-1             RcppHNSW_0.6.0          later_1.4.2            
##  [79] splines_4.5.0           lattice_0.22-7          survival_3.8-3         
##  [82] bit_4.6.0               deldir_2.0-4            tidyselect_1.2.1       
##  [85] Biostrings_2.76.0       miniUI_0.1.2            pbapply_1.7-2          
##  [88] knitr_1.50              gridExtra_2.3           svglite_2.2.1          
##  [91] scattermore_1.2         xfun_0.52               statmod_1.5.0          
##  [94] matrixStats_1.5.0       UCSC.utils_1.4.0        stringi_1.8.7          
##  [97] lazyeval_0.2.2          yaml_2.3.10             evaluate_1.0.4         
## [100] codetools_0.2-20        tibble_3.3.0            cli_3.6.5              
## [103] uwot_0.2.3              systemfonts_1.2.3       xtable_1.8-4           
## [106] reticulate_1.42.0       jquerylib_0.1.4         GenomeInfoDb_1.44.0    
## [109] Rcpp_1.1.0              globals_0.18.0          spatstat.random_3.4-1  
## [112] png_0.1-8               spatstat.univar_3.1-3   parallel_4.5.0         
## [115] ggplot2_3.5.2           blob_1.2.4              dotCall64_1.2          
## [118] listenv_0.9.1           viridisLite_0.4.2       scales_1.4.0           
## [121] ggridges_0.5.6          crayon_1.5.3            purrr_1.0.4            
## [124] rlang_1.1.6             KEGGREST_1.48.1         cowplot_1.2.0
```
