# Small_Program
 pieces of codes

## Run Azimuth locally
[Azimuth](https://satijalab.org/azimuth/) is an annotation tool which is recommended by Satiji Lab (Seurat team).
It also has R packages. However, running Azimuth in R requires internet connection，which is impossible when performing on computer cluster.

I've changed the source code which allowed you running Azimuth locally.
The only thing you have to do is download Azimuth reference data or make your own reference in required form.

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
