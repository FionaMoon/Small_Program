# Small_Program
 pieces of codes

## Run Azimuth locally
[Azimuth](https://satijalab.org/azimuth/) is an annotation tool which is recommended by Satiji Lab (Seurat team).
 It also offers R packages. However, running Azimuth in R typically requires an internet connection, which can be impractical when using a computer cluster.

I've modified the source code so that you can now run Azimuth locally. You only need to download the Azimuth reference data or create your own reference in the required format.

```r
rm(list=ls())
gc() ## release memory
## load packages
library(Seurat)
library(dplyr)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(progress)

source("./RunAzimuth_locally.R")
BM_refDir <- "./celltype_ref/Seurat_ref/bonemarrowref.SeuratData/inst/azimuth" ## Azimuth reference for bone marrow
pbmc_refDir <- "./celltype_ref/Seurat_ref/pbmcref.SeuratData/inst/azimuth" ## Azimuth reference for PBMC

## annotate by azimuth
## query/reference: annotate data
cat("\n ********* annotate data ************ \n")
SeuratObj <- RunAzimuth.Seurat(query = SeuratObj, reference = BM_refDir, assay = "RNA")
head(SeuratObj[[]])
table(SeuratObj$predicted.celltype.l2)
table(SeuratObj$predicted.celltype.l1)
