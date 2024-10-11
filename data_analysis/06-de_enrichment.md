---
title: "Introduction to Single Cell RNA-Seq Part 6: Enrichment and model-based differential expression"
author: "Bioinformatics Core"
date: "2024-08-29"
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
```

We will be continuing the work from Part 1 and so need to load in the RDS file

``` r
experiment.aggregate <- readRDS("scRNA_workshop-05.rds")
experiment.aggregate
```

```
## An object of class Seurat 
## 38606 features across 14018 samples within 1 assay 
## Active assay: RNA (38606 features, 2000 variable features)
##  3 layers present: counts, data, scale.data
##  2 dimensional reductions calculated: pca, umap
```

Lets go ahead and set that common seed for everyone

``` r
set.seed(12345)
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
   <td style="text-align:left;"> GO:0002181 </td>
   <td style="text-align:left;"> cytoplasmic translation </td>
   <td style="text-align:right;"> 159 </td>
   <td style="text-align:right;"> 134 </td>
   <td style="text-align:right;"> 38.18 </td>
   <td style="text-align:left;"> &lt; 1e-30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0042776 </td>
   <td style="text-align:left;"> proton motive force-driven mitochondrial ATP synthesis </td>
   <td style="text-align:right;"> 57 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 13.69 </td>
   <td style="text-align:left;"> 1.6e-28 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0043161 </td>
   <td style="text-align:left;"> proteasome-mediated ubiquitin-dependent protein catabolic pr... </td>
   <td style="text-align:right;"> 388 </td>
   <td style="text-align:right;"> 190 </td>
   <td style="text-align:right;"> 93.18 </td>
   <td style="text-align:left;"> 2.0e-24 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0051301 </td>
   <td style="text-align:left;"> cell division </td>
   <td style="text-align:right;"> 640 </td>
   <td style="text-align:right;"> 284 </td>
   <td style="text-align:right;"> 153.70 </td>
   <td style="text-align:left;"> 1.7e-20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006120 </td>
   <td style="text-align:left;"> mitochondrial electron transport, NADH to ubiquinone </td>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 10.09 </td>
   <td style="text-align:left;"> 2.4e-18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0032543 </td>
   <td style="text-align:left;"> mitochondrial translation </td>
   <td style="text-align:right;"> 130 </td>
   <td style="text-align:right;"> 76 </td>
   <td style="text-align:right;"> 31.22 </td>
   <td style="text-align:left;"> 4.2e-17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0000398 </td>
   <td style="text-align:left;"> mRNA splicing, via spliceosome </td>
   <td style="text-align:right;"> 304 </td>
   <td style="text-align:right;"> 200 </td>
   <td style="text-align:right;"> 73.01 </td>
   <td style="text-align:left;"> 5.1e-15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006338 </td>
   <td style="text-align:left;"> chromatin remodeling </td>
   <td style="text-align:right;"> 655 </td>
   <td style="text-align:right;"> 285 </td>
   <td style="text-align:right;"> 157.30 </td>
   <td style="text-align:left;"> 9.6e-14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0032981 </td>
   <td style="text-align:left;"> mitochondrial respiratory chain complex I assembly </td>
   <td style="text-align:right;"> 55 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 13.21 </td>
   <td style="text-align:left;"> 2.6e-13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0045893 </td>
   <td style="text-align:left;"> positive regulation of DNA-templated transcription </td>
   <td style="text-align:right;"> 1677 </td>
   <td style="text-align:right;"> 556 </td>
   <td style="text-align:right;"> 402.74 </td>
   <td style="text-align:left;"> 9.4e-13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:1903241 </td>
   <td style="text-align:left;"> U2-type prespliceosome assembly </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 6.00 </td>
   <td style="text-align:left;"> 9.6e-13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0070936 </td>
   <td style="text-align:left;"> protein K48-linked ubiquitination </td>
   <td style="text-align:right;"> 82 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 19.69 </td>
   <td style="text-align:left;"> 2.3e-12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006886 </td>
   <td style="text-align:left;"> intracellular protein transport </td>
   <td style="text-align:right;"> 671 </td>
   <td style="text-align:right;"> 307 </td>
   <td style="text-align:right;"> 161.14 </td>
   <td style="text-align:left;"> 1.1e-11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006281 </td>
   <td style="text-align:left;"> DNA repair </td>
   <td style="text-align:right;"> 598 </td>
   <td style="text-align:right;"> 296 </td>
   <td style="text-align:right;"> 143.61 </td>
   <td style="text-align:left;"> 1.7e-11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006406 </td>
   <td style="text-align:left;"> mRNA export from nucleus </td>
   <td style="text-align:right;"> 71 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 17.05 </td>
   <td style="text-align:left;"> 3.0e-11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006325 </td>
   <td style="text-align:left;"> chromatin organization </td>
   <td style="text-align:right;"> 798 </td>
   <td style="text-align:right;"> 358 </td>
   <td style="text-align:right;"> 191.64 </td>
   <td style="text-align:left;"> 8.3e-11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0050821 </td>
   <td style="text-align:left;"> protein stabilization </td>
   <td style="text-align:right;"> 211 </td>
   <td style="text-align:right;"> 91 </td>
   <td style="text-align:right;"> 50.67 </td>
   <td style="text-align:left;"> 6.6e-10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0090148 </td>
   <td style="text-align:left;"> membrane fission </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 10.57 </td>
   <td style="text-align:left;"> 7.0e-10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006446 </td>
   <td style="text-align:left;"> regulation of translational initiation </td>
   <td style="text-align:right;"> 84 </td>
   <td style="text-align:right;"> 57 </td>
   <td style="text-align:right;"> 20.17 </td>
   <td style="text-align:left;"> 1.2e-09 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0008380 </td>
   <td style="text-align:left;"> RNA splicing </td>
   <td style="text-align:right;"> 451 </td>
   <td style="text-align:right;"> 280 </td>
   <td style="text-align:right;"> 108.31 </td>
   <td style="text-align:left;"> 1.5e-09 </td>
  </tr>
</tbody>
</table>

* Annotated: number of genes (out of all.genes) that are annotated with that GO term
* Significant: number of genes that are annotated with that GO term and meet our criteria for "expressed"
* Expected: Under random chance, number of genes that would be expected to be annotated with that GO term and meeting our criteria for "expressed"
* Fisher: (Raw) p-value from Fisher's Exact Test

## 2. Model-based DE analysis in limma
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html) is an R package for differential expression analysis of bulk RNASeq and microarray data.  We apply it here to single cell data.

Limma can be used to fit any linear model to expression data and is useful for analyses that go beyond two-group comparisons.  A detailed tutorial of model specification in limma is available [here](https://ucsf-cat-bioinformatics.github.io/2021-June-RNA-Seq-Analysis/data_analysis/DE_Analysis_mm_with_quizzes) and in the [limma User's Guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf).


``` r
# filter genes to those expressed in at least 10% of cells
keep <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.1)]
expr2 <- expr[keep,]

# Set up "design matrix" with statistical model
cluster10$proper.group <- make.names(cluster10$orig.ident)
mm <- model.matrix(~0 + proper.group + S.Score + G2M.Score + percent.mito + nFeature_RNA, data = cluster10[[]])
head(mm) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent.mito </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> AAACCAGCACATTAGC+LRTI_WRK1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.1586883 </td>
   <td style="text-align:right;"> 1.0371359 </td>
   <td style="text-align:right;"> 3.202255 </td>
   <td style="text-align:right;"> 5479 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAACCCATCAATGTAC+LRTI_WRK1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.7347704 </td>
   <td style="text-align:right;"> 0.6561995 </td>
   <td style="text-align:right;"> 3.305649 </td>
   <td style="text-align:right;"> 6910 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAACCGTGTCGTACCT+LRTI_WRK1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.8925801 </td>
   <td style="text-align:right;"> 0.4061249 </td>
   <td style="text-align:right;"> 3.362421 </td>
   <td style="text-align:right;"> 4621 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAACGGCCAGGTACAA+LRTI_WRK1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.5515437 </td>
   <td style="text-align:right;"> 0.8579274 </td>
   <td style="text-align:right;"> 4.211121 </td>
   <td style="text-align:right;"> 4638 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAACTTATCGGTTCGG+LRTI_WRK1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.5691890 </td>
   <td style="text-align:right;"> 0.9509525 </td>
   <td style="text-align:right;"> 5.031228 </td>
   <td style="text-align:right;"> 7667 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAAGGGTTCATCAGGC+LRTI_WRK1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.8256540 </td>
   <td style="text-align:right;"> 0.3402400 </td>
   <td style="text-align:right;"> 6.204631 </td>
   <td style="text-align:right;"> 6241 </td>
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
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent.mito </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GCCCATTGTCACGGGT+LRTI_WRK2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0361409 </td>
   <td style="text-align:right;"> -0.0412531 </td>
   <td style="text-align:right;"> 2.910205 </td>
   <td style="text-align:right;"> 2253 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GGCTTAATCGGCGATT+LRTI_WRK2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0639132 </td>
   <td style="text-align:right;"> -0.0778717 </td>
   <td style="text-align:right;"> 4.119639 </td>
   <td style="text-align:right;"> 1131 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GGTACTATCGTAGCGT+LRTI_WRK2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0276898 </td>
   <td style="text-align:right;"> 0.0288580 </td>
   <td style="text-align:right;"> 2.727273 </td>
   <td style="text-align:right;"> 2167 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTCGCCAGTACTATGA+LRTI_WRK3 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0426008 </td>
   <td style="text-align:right;"> -0.0503667 </td>
   <td style="text-align:right;"> 2.546757 </td>
   <td style="text-align:right;"> 6586 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AATCCGTCAAAGCGTA+LRTI_WRK4 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.4301881 </td>
   <td style="text-align:right;"> 0.2368817 </td>
   <td style="text-align:right;"> 2.721987 </td>
   <td style="text-align:right;"> 7958 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CACACTCGTACCACCT+LRTI_WRK4 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0814435 </td>
   <td style="text-align:right;"> -0.0556294 </td>
   <td style="text-align:right;"> 5.963636 </td>
   <td style="text-align:right;"> 836 </td>
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
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent.mito </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> LINC01128 </td>
   <td style="text-align:right;"> 0.0104547 </td>
   <td style="text-align:right;"> 0.1828878 </td>
   <td style="text-align:right;"> 0.0116083 </td>
   <td style="text-align:right;"> -0.0300109 </td>
   <td style="text-align:right;"> 0.0698463 </td>
   <td style="text-align:right;"> 0.0062522 </td>
   <td style="text-align:right;"> 0.0064823 </td>
   <td style="text-align:right;"> -3.80e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NOC2L </td>
   <td style="text-align:right;"> -0.0822753 </td>
   <td style="text-align:right;"> 0.1396044 </td>
   <td style="text-align:right;"> -0.0935686 </td>
   <td style="text-align:right;"> -0.0203252 </td>
   <td style="text-align:right;"> 0.1025678 </td>
   <td style="text-align:right;"> -0.0104968 </td>
   <td style="text-align:right;"> 0.0069699 </td>
   <td style="text-align:right;"> 5.26e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KLHL17 </td>
   <td style="text-align:right;"> 0.0713555 </td>
   <td style="text-align:right;"> 0.0281549 </td>
   <td style="text-align:right;"> 0.0366542 </td>
   <td style="text-align:right;"> 0.1126129 </td>
   <td style="text-align:right;"> 0.0452909 </td>
   <td style="text-align:right;"> -0.0207237 </td>
   <td style="text-align:right;"> -0.0065486 </td>
   <td style="text-align:right;"> -2.90e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ISG15 </td>
   <td style="text-align:right;"> 1.6378853 </td>
   <td style="text-align:right;"> 0.6244244 </td>
   <td style="text-align:right;"> 0.2260651 </td>
   <td style="text-align:right;"> 1.9633809 </td>
   <td style="text-align:right;"> -0.3275326 </td>
   <td style="text-align:right;"> -0.0656142 </td>
   <td style="text-align:right;"> -0.0763379 </td>
   <td style="text-align:right;"> 3.31e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C1orf159 </td>
   <td style="text-align:right;"> 0.0279527 </td>
   <td style="text-align:right;"> 0.0208602 </td>
   <td style="text-align:right;"> 0.8150750 </td>
   <td style="text-align:right;"> -0.0640302 </td>
   <td style="text-align:right;"> 0.1261940 </td>
   <td style="text-align:right;"> 0.0248038 </td>
   <td style="text-align:right;"> 0.0134276 </td>
   <td style="text-align:right;"> -6.60e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TNFRSF18 </td>
   <td style="text-align:right;"> 1.1823790 </td>
   <td style="text-align:right;"> 1.4462475 </td>
   <td style="text-align:right;"> 1.4931473 </td>
   <td style="text-align:right;"> 0.1115185 </td>
   <td style="text-align:right;"> -0.2780431 </td>
   <td style="text-align:right;"> -0.2435467 </td>
   <td style="text-align:right;"> -0.0291735 </td>
   <td style="text-align:right;"> 2.46e-05 </td>
  </tr>
</tbody>
</table>

``` r
# Test 'Normal' - 'Colorectal.Cancer'
contr <- makeContrasts(proper.groupLRTI_WRK1	 - proper.groupLRTI_WRK2	, levels = colnames(coef(fit)))
contr %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupLRTI_WRK1 - proper.groupLRTI_WRK2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> proper.groupLRTI_WRK1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> proper.groupLRTI_WRK2 </td>
   <td style="text-align:right;"> -1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> proper.groupLRTI_WRK3 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> proper.groupLRTI_WRK4 </td>
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
   <td style="text-align:left;"> percent.mito </td>
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
   <td style="text-align:left;"> VAV3 </td>
   <td style="text-align:right;"> -1.7508024 </td>
   <td style="text-align:right;"> 0.3110052 </td>
   <td style="text-align:right;"> -13.236840 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 62.96266 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TBC1D4 </td>
   <td style="text-align:right;"> -1.3900664 </td>
   <td style="text-align:right;"> 0.2833838 </td>
   <td style="text-align:right;"> -11.568702 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 48.94408 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HDAC9 </td>
   <td style="text-align:right;"> -1.2517491 </td>
   <td style="text-align:right;"> 0.1753274 </td>
   <td style="text-align:right;"> -10.333661 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 39.11020 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IL1R1 </td>
   <td style="text-align:right;"> -0.9995498 </td>
   <td style="text-align:right;"> 0.1300682 </td>
   <td style="text-align:right;"> -10.226861 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 38.28639 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CCR8 </td>
   <td style="text-align:right;"> -0.8353974 </td>
   <td style="text-align:right;"> 0.1163754 </td>
   <td style="text-align:right;"> -10.036890 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 36.83248 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EML4 </td>
   <td style="text-align:right;"> 1.5052329 </td>
   <td style="text-align:right;"> 1.7012424 </td>
   <td style="text-align:right;"> 9.561668 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 33.26267 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ARPC2 </td>
   <td style="text-align:right;"> 1.1548512 </td>
   <td style="text-align:right;"> 2.4380833 </td>
   <td style="text-align:right;"> 9.464631 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 32.54606 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MT-ND2 </td>
   <td style="text-align:right;"> -1.0490250 </td>
   <td style="text-align:right;"> 2.3930152 </td>
   <td style="text-align:right;"> -8.952967 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 28.84073 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LAYN </td>
   <td style="text-align:right;"> -0.7911447 </td>
   <td style="text-align:right;"> 0.1180913 </td>
   <td style="text-align:right;"> -8.823339 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 27.92232 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RPS4Y1 </td>
   <td style="text-align:right;"> 1.0730097 </td>
   <td style="text-align:right;"> 1.1304680 </td>
   <td style="text-align:right;"> 8.153916 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 23.31908 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000227240 </td>
   <td style="text-align:right;"> -1.2591894 </td>
   <td style="text-align:right;"> 0.2575334 </td>
   <td style="text-align:right;"> -7.784149 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 20.88238 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAPZB </td>
   <td style="text-align:right;"> 1.0079794 </td>
   <td style="text-align:right;"> 1.8103635 </td>
   <td style="text-align:right;"> 7.756462 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 20.70311 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TXNIP </td>
   <td style="text-align:right;"> -1.2963896 </td>
   <td style="text-align:right;"> 0.6787166 </td>
   <td style="text-align:right;"> -7.580633 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 19.57527 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAPZA1 </td>
   <td style="text-align:right;"> 0.9017427 </td>
   <td style="text-align:right;"> 1.2495491 </td>
   <td style="text-align:right;"> 7.516820 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 19.17053 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FTL </td>
   <td style="text-align:right;"> -1.2627828 </td>
   <td style="text-align:right;"> 2.1523660 </td>
   <td style="text-align:right;"> -7.179795 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 17.07449 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FRYL </td>
   <td style="text-align:right;"> 0.9189332 </td>
   <td style="text-align:right;"> 0.8948729 </td>
   <td style="text-align:right;"> 7.026423 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 16.14428 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UTY </td>
   <td style="text-align:right;"> 0.9005716 </td>
   <td style="text-align:right;"> 0.5991670 </td>
   <td style="text-align:right;"> 6.983129 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 15.88444 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SGPP1 </td>
   <td style="text-align:right;"> -0.6447435 </td>
   <td style="text-align:right;"> 0.2190291 </td>
   <td style="text-align:right;"> -6.961139 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 15.75292 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IL2RG </td>
   <td style="text-align:right;"> 1.0240127 </td>
   <td style="text-align:right;"> 1.8196712 </td>
   <td style="text-align:right;"> 6.898402 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 15.37944 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GIMAP7 </td>
   <td style="text-align:right;"> 1.4491604 </td>
   <td style="text-align:right;"> 1.4094892 </td>
   <td style="text-align:right;"> 6.831560 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 14.98436 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PTPRM </td>
   <td style="text-align:right;"> -0.7061117 </td>
   <td style="text-align:right;"> 0.1497160 </td>
   <td style="text-align:right;"> -6.752201 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 14.51912 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IRF1 </td>
   <td style="text-align:right;"> 0.9660324 </td>
   <td style="text-align:right;"> 1.1877894 </td>
   <td style="text-align:right;"> 6.728966 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 14.38370 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MCF2L2 </td>
   <td style="text-align:right;"> -0.7044998 </td>
   <td style="text-align:right;"> 0.2238412 </td>
   <td style="text-align:right;"> -6.688666 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 14.14968 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DEF6 </td>
   <td style="text-align:right;"> 0.8446568 </td>
   <td style="text-align:right;"> 0.9440441 </td>
   <td style="text-align:right;"> 6.655879 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 13.96009 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TUT7 </td>
   <td style="text-align:right;"> 0.8658212 </td>
   <td style="text-align:right;"> 0.7549669 </td>
   <td style="text-align:right;"> 6.647806 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 13.91351 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTN </td>
   <td style="text-align:right;"> -0.7275733 </td>
   <td style="text-align:right;"> 0.1723645 </td>
   <td style="text-align:right;"> -6.544020 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1e-07 </td>
   <td style="text-align:right;"> 13.31871 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UTRN </td>
   <td style="text-align:right;"> 1.0748347 </td>
   <td style="text-align:right;"> 1.2238449 </td>
   <td style="text-align:right;"> 6.437446 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1e-07 </td>
   <td style="text-align:right;"> 12.71557 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CELF2 </td>
   <td style="text-align:right;"> 1.1100942 </td>
   <td style="text-align:right;"> 1.7269740 </td>
   <td style="text-align:right;"> 6.408570 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2e-07 </td>
   <td style="text-align:right;"> 12.55350 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPP1R18 </td>
   <td style="text-align:right;"> 0.7727838 </td>
   <td style="text-align:right;"> 1.0353629 </td>
   <td style="text-align:right;"> 6.363937 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2e-07 </td>
   <td style="text-align:right;"> 12.30411 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SYNRG </td>
   <td style="text-align:right;"> 0.8090564 </td>
   <td style="text-align:right;"> 0.8666973 </td>
   <td style="text-align:right;"> 6.330475 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2e-07 </td>
   <td style="text-align:right;"> 12.11806 </td>
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

## Prepare for the next section

#### Session Information

``` r
sessionInfo()
```

```
## R version 4.4.1 (2024-06-14)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sonoma 14.6.1
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
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] org.Hs.eg.db_3.19.1  kableExtra_1.4.0     dplyr_1.1.4         
##  [4] topGO_2.56.0         SparseM_1.84-2       GO.db_3.19.1        
##  [7] AnnotationDbi_1.66.0 IRanges_2.38.1       S4Vectors_0.42.1    
## [10] Biobase_2.64.0       graph_1.82.0         BiocGenerics_0.50.0 
## [13] limma_3.60.4         Seurat_5.1.0         SeuratObject_5.0.2  
## [16] sp_2.1-4            
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3      rstudioapi_0.16.0       jsonlite_1.8.8         
##   [4] magrittr_2.0.3          spatstat.utils_3.1-0    rmarkdown_2.28         
##   [7] zlibbioc_1.50.0         vctrs_0.6.5             ROCR_1.0-11            
##  [10] memoise_2.0.1           spatstat.explore_3.3-2  htmltools_0.5.8.1      
##  [13] sass_0.4.9              sctransform_0.4.1       parallelly_1.38.0      
##  [16] KernSmooth_2.23-24      bslib_0.8.0             htmlwidgets_1.6.4      
##  [19] ica_1.0-3               plyr_1.8.9              plotly_4.10.4          
##  [22] zoo_1.8-12              cachem_1.1.0            igraph_2.0.3           
##  [25] mime_0.12               lifecycle_1.0.4         pkgconfig_2.0.3        
##  [28] Matrix_1.7-0            R6_2.5.1                fastmap_1.2.0          
##  [31] GenomeInfoDbData_1.2.12 fitdistrplus_1.2-1      future_1.34.0          
##  [34] shiny_1.9.1             digest_0.6.37           colorspace_2.1-1       
##  [37] patchwork_1.2.0         tensor_1.5              RSpectra_0.16-2        
##  [40] irlba_2.3.5.1           RSQLite_2.3.7           progressr_0.14.0       
##  [43] fansi_1.0.6             spatstat.sparse_3.1-0   httr_1.4.7             
##  [46] polyclip_1.10-7         abind_1.4-5             compiler_4.4.1         
##  [49] bit64_4.0.5             DBI_1.2.3               fastDummies_1.7.4      
##  [52] highr_0.11              MASS_7.3-61             tools_4.4.1            
##  [55] lmtest_0.9-40           httpuv_1.6.15           future.apply_1.11.2    
##  [58] goftest_1.2-3           glue_1.7.0              nlme_3.1-166           
##  [61] promises_1.3.0          grid_4.4.1              Rtsne_0.17             
##  [64] cluster_2.1.6           reshape2_1.4.4          generics_0.1.3         
##  [67] gtable_0.3.5            spatstat.data_3.1-2     tidyr_1.3.1            
##  [70] data.table_1.16.0       xml2_1.3.6              XVector_0.44.0         
##  [73] utf8_1.2.4              spatstat.geom_3.3-2     RcppAnnoy_0.0.22       
##  [76] ggrepel_0.9.5           RANN_2.6.2              pillar_1.9.0           
##  [79] stringr_1.5.1           spam_2.10-0             RcppHNSW_0.6.0         
##  [82] later_1.3.2             splines_4.4.1           lattice_0.22-6         
##  [85] bit_4.0.5               survival_3.7-0          deldir_2.0-4           
##  [88] tidyselect_1.2.1        Biostrings_2.72.1       miniUI_0.1.1.1         
##  [91] pbapply_1.7-2           knitr_1.48              gridExtra_2.3          
##  [94] svglite_2.1.3           scattermore_1.2         xfun_0.47              
##  [97] statmod_1.5.0           matrixStats_1.3.0       UCSC.utils_1.0.0       
## [100] stringi_1.8.4           lazyeval_0.2.2          yaml_2.3.10            
## [103] evaluate_0.24.0         codetools_0.2-20        tibble_3.2.1           
## [106] cli_3.6.3               uwot_0.2.2              systemfonts_1.1.0      
## [109] xtable_1.8-4            reticulate_1.38.0       munsell_0.5.1          
## [112] jquerylib_0.1.4         GenomeInfoDb_1.40.1     Rcpp_1.0.13            
## [115] globals_0.16.3          spatstat.random_3.3-1   png_0.1-8              
## [118] spatstat.univar_3.0-0   parallel_4.4.1          blob_1.2.4             
## [121] ggplot2_3.5.1           dotCall64_1.1-1         listenv_0.9.1          
## [124] viridisLite_0.4.2       scales_1.3.0            ggridges_0.5.6         
## [127] crayon_1.5.3            leiden_0.4.3.1          purrr_1.0.2            
## [130] rlang_1.1.4             KEGGREST_1.44.1         cowplot_1.1.3
```
