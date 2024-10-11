---
title: "Basics of Single Cell RNA-Seq Part 2: QA and filtering"
author: "UCSF CAT Bioinformatics"
date: "2024-08-27"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Introduction to Single Cell RNA-Seq Part 2: QA and filtering


## Set up workspace
First, we need to load the required libraries.

``` r
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(kableExtra)
<div class='r_output'>
If you are continuing directly from part 1, the experiment.aggregate object is likely already in your workspace. In case you cleared your workspace at the end of the previous section, or are working on this project at a later date after re-starting R, you can use the `readRDS` function to read your saved Seurat object from part 1.

``` r
experiment.aggregate <- readRDS("scRNA_workshop-01.rds")
experiment.aggregate
</div>
<div class='r_output'> An object of class Seurat 
 38606 features across 39196 samples within 1 assay 
 Active assay: RNA (38606 features, 0 variable features)
  1 layer present: counts
</div>
#### First lets look at the number of valid cells in each batch

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>Cell counts per sample</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Sample </th>
   <th style="text-align:right;"> Number of Cells </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> LRTI_WRK1 </td>
   <td style="text-align:right;"> 24478 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LRTI_WRK2 </td>
   <td style="text-align:right;"> 14718 </td>
  </tr>
</tbody>
</table>

A seed is used to initialize pseudo-random functions. Some of the functions we will be using have pseudo-random elements. Setting a common seed ensures that all of us will get the same results, and that the results will remain stable when re-run.

``` r
set.seed(12345)
<div class='r_output'>
 Display metadata QA/QC

Using a few nested functions, we can produce prettier, more detailed, versions of the simple exploratory summary statistics we generated for the available metadata in the last section. In the code below sections below,

Further, Seurat has a number of convenient built-in functions for visualizing metadata. These functions produce ggplot objects, which can easily be modified using ggplot2. Of course, all of these visualizations can be reproduced with custom code as well, and we will include some examples of both modifying Seurat plots and generating plots from scratch as the analysis continues.

1) 10% quantile tables are produced for each metadata value, separated by sample identity.
2) Ridge Plot
3) Violin Plot

Each is a different way of looking at the data as to get a better understanding of how each sample compares

 Feature counts (genes) per cell
<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>Feature count distribution by sample</caption>
 <thead>
  <tr>
   <th style="text-align:left;">  </th>
   <th style="text-align:right;"> LRTI_WRK1 </th>
   <th style="text-align:right;"> LRTI_WRK2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0% </td>
   <td style="text-align:right;"> 300.00 </td>
   <td style="text-align:right;"> 300.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5% </td>
   <td style="text-align:right;"> 307.00 </td>
   <td style="text-align:right;"> 387.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10% </td>
   <td style="text-align:right;"> 317.00 </td>
   <td style="text-align:right;"> 551.70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15% </td>
   <td style="text-align:right;"> 330.00 </td>
   <td style="text-align:right;"> 1135.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20% </td>
   <td style="text-align:right;"> 350.00 </td>
   <td style="text-align:right;"> 1643.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25% </td>
   <td style="text-align:right;"> 377.00 </td>
   <td style="text-align:right;"> 2090.25 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 30% </td>
   <td style="text-align:right;"> 409.00 </td>
   <td style="text-align:right;"> 2762.40 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 35% </td>
   <td style="text-align:right;"> 449.00 </td>
   <td style="text-align:right;"> 3359.90 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 40% </td>
   <td style="text-align:right;"> 510.00 </td>
   <td style="text-align:right;"> 3818.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 45% </td>
   <td style="text-align:right;"> 665.00 </td>
   <td style="text-align:right;"> 4176.65 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 50% </td>
   <td style="text-align:right;"> 983.00 </td>
   <td style="text-align:right;"> 4484.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 55% </td>
   <td style="text-align:right;"> 1380.00 </td>
   <td style="text-align:right;"> 4736.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 60% </td>
   <td style="text-align:right;"> 1795.00 </td>
   <td style="text-align:right;"> 4966.20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 65% </td>
   <td style="text-align:right;"> 2407.00 </td>
   <td style="text-align:right;"> 5182.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 70% </td>
   <td style="text-align:right;"> 3005.90 </td>
   <td style="text-align:right;"> 5360.90 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 75% </td>
   <td style="text-align:right;"> 3598.75 </td>
   <td style="text-align:right;"> 5548.75 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 80% </td>
   <td style="text-align:right;"> 4168.60 </td>
   <td style="text-align:right;"> 5738.60 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 85% </td>
   <td style="text-align:right;"> 4608.00 </td>
   <td style="text-align:right;"> 5934.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 90% </td>
   <td style="text-align:right;"> 5050.30 </td>
   <td style="text-align:right;"> 6185.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 95% </td>
   <td style="text-align:right;"> 5614.30 </td>
   <td style="text-align:right;"> 6588.15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 100% </td>
   <td style="text-align:right;"> 12390.00 </td>
   <td style="text-align:right;"> 10160.00 </td>
  </tr>
</tbody>
</table>

![](02-filtering_files/figure-html/nFeature-1.png)<!-- -->![](02-filtering_files/figure-html/nFeature-2.png)<!-- -->

 UMI counts per cell

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>UMI count distribution by sample</caption>
 <thead>
  <tr>
   <th style="text-align:left;">  </th>
   <th style="text-align:right;"> LRTI_WRK1 </th>
   <th style="text-align:right;"> LRTI_WRK2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0% </td>
   <td style="text-align:right;"> 329.00 </td>
   <td style="text-align:right;"> 369.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5% </td>
   <td style="text-align:right;"> 376.00 </td>
   <td style="text-align:right;"> 637.85 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10% </td>
   <td style="text-align:right;"> 394.00 </td>
   <td style="text-align:right;"> 1029.70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15% </td>
   <td style="text-align:right;"> 418.00 </td>
   <td style="text-align:right;"> 2169.65 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20% </td>
   <td style="text-align:right;"> 451.00 </td>
   <td style="text-align:right;"> 3280.20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25% </td>
   <td style="text-align:right;"> 494.00 </td>
   <td style="text-align:right;"> 4612.50 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 30% </td>
   <td style="text-align:right;"> 547.00 </td>
   <td style="text-align:right;"> 7444.30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 35% </td>
   <td style="text-align:right;"> 625.00 </td>
   <td style="text-align:right;"> 10452.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 40% </td>
   <td style="text-align:right;"> 782.00 </td>
   <td style="text-align:right;"> 13295.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 45% </td>
   <td style="text-align:right;"> 1127.65 </td>
   <td style="text-align:right;"> 15955.55 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 50% </td>
   <td style="text-align:right;"> 1698.50 </td>
   <td style="text-align:right;"> 18491.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 55% </td>
   <td style="text-align:right;"> 2396.00 </td>
   <td style="text-align:right;"> 20944.15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 60% </td>
   <td style="text-align:right;"> 3338.40 </td>
   <td style="text-align:right;"> 23336.40 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 65% </td>
   <td style="text-align:right;"> 5220.00 </td>
   <td style="text-align:right;"> 25543.05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 70% </td>
   <td style="text-align:right;"> 7494.90 </td>
   <td style="text-align:right;"> 27622.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 75% </td>
   <td style="text-align:right;"> 10529.75 </td>
   <td style="text-align:right;"> 29866.75 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 80% </td>
   <td style="text-align:right;"> 14225.80 </td>
   <td style="text-align:right;"> 32275.60 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 85% </td>
   <td style="text-align:right;"> 18083.00 </td>
   <td style="text-align:right;"> 34763.45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 90% </td>
   <td style="text-align:right;"> 21817.60 </td>
   <td style="text-align:right;"> 38190.50 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 95% </td>
   <td style="text-align:right;"> 26494.45 </td>
   <td style="text-align:right;"> 44213.30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 100% </td>
   <td style="text-align:right;"> 150158.00 </td>
   <td style="text-align:right;"> 100628.00 </td>
  </tr>
</tbody>
</table>

![](02-filtering_files/figure-html/nCount-1.png)<!-- -->![](02-filtering_files/figure-html/nCount-2.png)<!-- -->

 Percentage of Mitochondria per cell 

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>Mitochondrial percentage distribution by sample</caption>
 <thead>
  <tr>
   <th style="text-align:left;">  </th>
   <th style="text-align:right;"> LRTI_WRK1 </th>
   <th style="text-align:right;"> LRTI_WRK2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0% </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5% </td>
   <td style="text-align:right;"> 0.328 </td>
   <td style="text-align:right;"> 0.561 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10% </td>
   <td style="text-align:right;"> 0.694 </td>
   <td style="text-align:right;"> 1.517 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15% </td>
   <td style="text-align:right;"> 1.028 </td>
   <td style="text-align:right;"> 2.015 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20% </td>
   <td style="text-align:right;"> 1.340 </td>
   <td style="text-align:right;"> 2.311 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25% </td>
   <td style="text-align:right;"> 1.617 </td>
   <td style="text-align:right;"> 2.534 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 30% </td>
   <td style="text-align:right;"> 1.860 </td>
   <td style="text-align:right;"> 2.723 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 35% </td>
   <td style="text-align:right;"> 2.072 </td>
   <td style="text-align:right;"> 2.894 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 40% </td>
   <td style="text-align:right;"> 2.277 </td>
   <td style="text-align:right;"> 3.047 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 45% </td>
   <td style="text-align:right;"> 2.473 </td>
   <td style="text-align:right;"> 3.193 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 50% </td>
   <td style="text-align:right;"> 2.687 </td>
   <td style="text-align:right;"> 3.340 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 55% </td>
   <td style="text-align:right;"> 2.936 </td>
   <td style="text-align:right;"> 3.489 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 60% </td>
   <td style="text-align:right;"> 3.245 </td>
   <td style="text-align:right;"> 3.647 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 65% </td>
   <td style="text-align:right;"> 3.695 </td>
   <td style="text-align:right;"> 3.810 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 70% </td>
   <td style="text-align:right;"> 4.327 </td>
   <td style="text-align:right;"> 3.992 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 75% </td>
   <td style="text-align:right;"> 5.345 </td>
   <td style="text-align:right;"> 4.180 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 80% </td>
   <td style="text-align:right;"> 6.809 </td>
   <td style="text-align:right;"> 4.434 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 85% </td>
   <td style="text-align:right;"> 8.838 </td>
   <td style="text-align:right;"> 4.780 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 90% </td>
   <td style="text-align:right;"> 12.828 </td>
   <td style="text-align:right;"> 5.396 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 95% </td>
   <td style="text-align:right;"> 27.174 </td>
   <td style="text-align:right;"> 8.453 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 100% </td>
   <td style="text-align:right;"> 93.319 </td>
   <td style="text-align:right;"> 85.343 </td>
  </tr>
</tbody>
</table>

![](02-filtering_files/figure-html/pMito-1.png)<!-- -->![](02-filtering_files/figure-html/pMito-2.png)<!-- -->

 Percentage of Ribosomal (protein) per cell 

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>Ribosomal percentage distribution by sample</caption>
 <thead>
  <tr>
   <th style="text-align:left;">  </th>
   <th style="text-align:right;"> LRTI_WRK1 </th>
   <th style="text-align:right;"> LRTI_WRK2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0% </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5% </td>
   <td style="text-align:right;"> 0.873 </td>
   <td style="text-align:right;"> 0.646 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10% </td>
   <td style="text-align:right;"> 1.789 </td>
   <td style="text-align:right;"> 1.221 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15% </td>
   <td style="text-align:right;"> 2.857 </td>
   <td style="text-align:right;"> 4.471 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20% </td>
   <td style="text-align:right;"> 3.678 </td>
   <td style="text-align:right;"> 5.520 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25% </td>
   <td style="text-align:right;"> 4.263 </td>
   <td style="text-align:right;"> 5.966 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 30% </td>
   <td style="text-align:right;"> 4.737 </td>
   <td style="text-align:right;"> 6.323 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 35% </td>
   <td style="text-align:right;"> 5.166 </td>
   <td style="text-align:right;"> 6.623 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 40% </td>
   <td style="text-align:right;"> 5.562 </td>
   <td style="text-align:right;"> 6.904 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 45% </td>
   <td style="text-align:right;"> 5.964 </td>
   <td style="text-align:right;"> 7.170 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 50% </td>
   <td style="text-align:right;"> 6.413 </td>
   <td style="text-align:right;"> 7.425 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 55% </td>
   <td style="text-align:right;"> 6.908 </td>
   <td style="text-align:right;"> 7.701 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 60% </td>
   <td style="text-align:right;"> 7.465 </td>
   <td style="text-align:right;"> 7.978 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 65% </td>
   <td style="text-align:right;"> 8.184 </td>
   <td style="text-align:right;"> 8.269 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 70% </td>
   <td style="text-align:right;"> 9.135 </td>
   <td style="text-align:right;"> 8.622 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 75% </td>
   <td style="text-align:right;"> 10.345 </td>
   <td style="text-align:right;"> 9.060 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 80% </td>
   <td style="text-align:right;"> 11.635 </td>
   <td style="text-align:right;"> 9.665 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 85% </td>
   <td style="text-align:right;"> 13.014 </td>
   <td style="text-align:right;"> 10.719 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 90% </td>
   <td style="text-align:right;"> 14.680 </td>
   <td style="text-align:right;"> 12.943 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 95% </td>
   <td style="text-align:right;"> 16.986 </td>
   <td style="text-align:right;"> 17.130 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 100% </td>
   <td style="text-align:right;"> 45.979 </td>
   <td style="text-align:right;"> 40.878 </td>
  </tr>
</tbody>
</table>

![](02-filtering_files/figure-html/pRibo-1.png)<!-- -->![](02-filtering_files/figure-html/pRibo-2.png)<!-- -->

 Modifying Seurat plots

Modifying the ggplot objects produced by a Seurat plotting function works best on individual panels. Therefore, to recreate the function above with modifications, we can use `lapply` to create a list of plots. In some cases it may be more appropriate to create the plots individually so that different modifications can be applied to each plot.


``` r
VlnPlot(experiment.aggregate, features = "nCount_RNA", pt.size = 0.01) + 
  scale_y_continuous(trans = "log10") +
  scale_fill_viridis_d(option = "mako") +
  ggtitle("log10(nCount_RNA)")
</div>
![](02-filtering_files/figure-html/violins_list-1.png)<!-- -->

These can later be stitched together with another library, like patchwork, or cowplot.

### Custom plots
The Seurat built-in functions are useful and easy to interact with, but sometimes you may wish to visualize something for which a plotting function does not already exist. For example, we might want to see how many cells are expressing each gene over some UMI threshold.

The code below produces a ranked plot similar to the barcode inflection plots from the last section. On the x-axis are the genes arranged from least ubiquitously expressed to most. In a single cell dataset, many genes are expessed in a relatively small number of cells, or not at all. The y-axis displays the number of cells in which each gene is expressed.

**Note: This function is SLOW. You may want to skip this code block or run it while you take a break for a few minutes.**

``` r
plot(sort(Matrix::rowSums(GetAssayData(experiment.aggregate) >= 3)) , xlab="gene rank", ylab="number of cells", main="Cells per genes (reads/gene >= 3 )")
<div class='r_output'>
![](02-filtering_files/figure-html/gene_range-1.png)<!-- -->

# Scatter plots
Scatter plots allow us to visualize the relationships between the metadata variables.

Gene Plot, scatter plot of gene expression across cells, (colored by sample)
![](02-filtering_files/figure-html/gene_plot-1.png)<!-- -->![](02-filtering_files/figure-html/gene_plot-2.png)<!-- -->![](02-filtering_files/figure-html/gene_plot-3.png)<!-- -->

# Cell filtering

We use the information above to filter out cells. Here we choose those that have percent mitochondrial genes max of 10% and unique UMI counts under 50,000 or greater than 500. Further requiring the number of features persent per cell to be 1000 genes.

![](02-filtering_files/figure-html/filterviz-1.png)<!-- -->

 Cell Complexity

The standard way of calculating this is log10(genes)/log10(counts) however this gives absolute values which are difficult to judge. A possibly better approach is to fit a line through the cloud and then calculate the difference from the observed value to the expected.


</div>## 
## Call:
## lm(formula = log10(qc.metrics$nFeature_RNA) ~ log10(qc.metrics$nCount_RNA))
## 
## Coefficients:
##                  (Intercept)  log10(qc.metrics$nCount_RNA)  
##                       0.8027                        0.6678
<div class='r_output'>
![](02-filtering_files/figure-html/complex-1.png)<!-- -->
And we can add the information to the scatter plot

![](02-filtering_files/figure-html/complexplot-1.png)<!-- -->
  
  
# Histograms of metadata

Histograms can also be useful to look at the distribution over all the cells

![](02-filtering_files/figure-html/count_hist-1.png)<!-- -->

Displaying the number of genes per cell

![](02-filtering_files/figure-html/features_hist-1.png)<!-- -->


With Mitochondrial expression

![](02-filtering_files/figure-html/mito_hist-1.png)<!-- -->


 Cell filtering
The goal of cell filtering is to remove cells with anomolous expression profiles, typically low UMI cells, which may correspond to low-quality cells or background barcodes. It may also be appropriate to remove outlier cells with extremely high UMI counts.
In this case, the proposed cut-offs on the high end of the distributions are quite conservative, in part to reduce the size of the object and speed up analysis during the workshop.


These filters can be put in place with the `subset` function.

 Prefiltered data

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>Cell countse per sample</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Sample </th>
   <th style="text-align:right;"> Number of Cells </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> LRTI_WRK1 </td>
   <td style="text-align:right;"> 24478 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LRTI_WRK2 </td>
   <td style="text-align:right;"> 14718 </td>
  </tr>
</tbody>
</table>




 Post filtered data

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>Cell countse per sample</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Sample </th>
   <th style="text-align:right;"> Number of Cells </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> LRTI_WRK1 </td>
   <td style="text-align:right;"> 11727 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LRTI_WRK2 </td>
   <td style="text-align:right;"> 12405 </td>
  </tr>
</tbody>
</table>

**Play with the filtering parameters, and see how the results change. Is there a set of parameters you feel is more appropriate? Why?**

 Feature filtering
When creating the base Seurat object, we had the opportunity filter out some genes using the "min.cells" argument. At the time, we set the min feature to keep a cell to 300. Since we didn't filter out any features then (set to 0), we can apply a filter at this point. If we had filtered when the object was created, this would be an opportunity to be more aggressive. The custom code below provides a function that filters genes requiring a min.umi in at least min.cells, or takes a user-provided list of genes.



``` r
# define function
FilterGenes <- function(object, min.umi = NA, min.cells = NA, genes = NULL) {
  genes.use = NA
  if (!is.null(genes)) {
    genes.use = intersect(rownames(object), genes)
    } else if (min.cells & min.umi) {
      num.cells = Matrix::rowSums(GetAssayData(object) >= min.umi)
      genes.use = names(num.cells[which(num.cells >= min.cells)])
    }
  object = object[genes.use,]
  object = LogSeuratCommand(object = object)
  return(object)
}
# apply filter
experiment.filter <- FilterGenes(object = experiment.aggregate.filtered, min.umi = 2, min.cells = 10)
# filtering results
experiment.filter
</div>
<div class='r_output'> An object of class Seurat 
 16019 features across 24132 samples within 1 assay 
 Active assay: RNA (16019 features, 0 variable features)
  1 layer present: counts
</div>
``` r
experiment.aggregate.filtered <- experiment.filter
<div class='r_output'>

![](02-filtering_files/figure-html/filteredviz-1.png)<!-- -->

 Prepare for the next section

 Save object

``` r
saveRDS(experiment.aggregate.filtered, file="scRNA_workshop-02.rds")
</div>
#### Download Rmd

``` r
download.file("https://raw.githubusercontent.com/ucsf-cat-bioinformatics/2024-08-SCRNA-Seq-Analysis/main/data_analysis/03-normalize_scale.Rmd", "03-normalize_scale.Rmd")
<div class='r_output'>
 Session Information

``` r
sessionInfo()
</div>
<div class='r_output'> R version 4.4.1 (2024-06-14)
 Platform: aarch64-apple-darwin20
 Running under: macOS Sonoma 14.6.1
 
 Matrix products: default
 BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
 LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
 
 locale:
 [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
 
 time zone: America/Los_Angeles
 tzcode source: internal
 
 attached base packages:
 [1] stats     graphics  grDevices utils     datasets  methods   base     
 
 other attached packages:
 [1] kableExtra_1.4.0   dplyr_1.1.4        tidyr_1.3.1        ggplot2_3.5.1     
 [5] Seurat_5.1.0       SeuratObject_5.0.2 sp_2.1-4          
 
 loaded via a namespace (and not attached):
   [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3         
   [4] rlang_1.1.4            magrittr_2.0.3         RcppAnnoy_0.0.22      
   [7] spatstat.geom_3.3-2    matrixStats_1.3.0      ggridges_0.5.6        
  [10] compiler_4.4.1         systemfonts_1.1.0      png_0.1-8             
  [13] vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1         
  [16] pkgconfig_2.0.3        fastmap_1.2.0          labeling_0.4.3        
  [19] utf8_1.2.4             promises_1.3.0         rmarkdown_2.28        
  [22] purrr_1.0.2            xfun_0.47              cachem_1.1.0          
  [25] jsonlite_1.8.8         goftest_1.2-3          highr_0.11            
  [28] later_1.3.2            spatstat.utils_3.1-0   irlba_2.3.5.1         
  [31] parallel_4.4.1         cluster_2.1.6          R6_2.5.1              
  [34] ica_1.0-3              spatstat.data_3.1-2    bslib_0.8.0           
  [37] stringi_1.8.4          RColorBrewer_1.1-3     reticulate_1.38.0     
  [40] spatstat.univar_3.0-0  parallelly_1.38.0      lmtest_0.9-40         
  [43] jquerylib_0.1.4        scattermore_1.2        Rcpp_1.0.13           
  [46] knitr_1.48             tensor_1.5             future.apply_1.11.2   
  [49] zoo_1.8-12             sctransform_0.4.1      httpuv_1.6.15         
  [52] Matrix_1.7-0           splines_4.4.1          igraph_2.0.3          
  [55] tidyselect_1.2.1       abind_1.4-5            rstudioapi_0.16.0     
  [58] yaml_2.3.10            spatstat.random_3.3-1  codetools_0.2-20      
  [61] miniUI_0.1.1.1         spatstat.explore_3.3-2 listenv_0.9.1         
  [64] lattice_0.22-6         tibble_3.2.1           plyr_1.8.9            
  [67] withr_3.0.1            shiny_1.9.1            ROCR_1.0-11           
  [70] evaluate_0.24.0        Rtsne_0.17             future_1.34.0         
  [73] fastDummies_1.7.4      survival_3.7-0         polyclip_1.10-7       
  [76] xml2_1.3.6             fitdistrplus_1.2-1     pillar_1.9.0          
  [79] KernSmooth_2.23-24     plotly_4.10.4          generics_0.1.3        
  [82] RcppHNSW_0.6.0         munsell_0.5.1          scales_1.3.0          
  [85] globals_0.16.3         xtable_1.8-4           glue_1.7.0            
  [88] lazyeval_0.2.2         tools_4.4.1            data.table_1.15.4     
  [91] RSpectra_0.16-2        RANN_2.6.2             leiden_0.4.3.1        
  [94] dotCall64_1.1-1        cowplot_1.1.3          grid_4.4.1            
  [97] colorspace_2.1-1       nlme_3.1-166           patchwork_1.2.0       
 [100] cli_3.6.3              spatstat.sparse_3.1-0  spam_2.10-0           
 [103] fansi_1.0.6            viridisLite_0.4.2      svglite_2.1.3         
 [106] uwot_0.2.2             gtable_0.3.5           sass_0.4.9            
 [109] digest_0.6.37          progressr_0.14.0       ggrepel_0.9.5         
 [112] farver_2.1.2           htmlwidgets_1.6.4      htmltools_0.5.8.1     
 [115] lifecycle_1.0.4        httr_1.4.7             mime_0.12             
 [118] MASS_7.3-61
</div>