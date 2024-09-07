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
```
## Generate Alluvial Plots with a Single Line of Code
Many R packages can be used to make Alluvial plots, such as [easyalluvial](https://erblast.github.io/easyalluvial/), [ggforce](https://ggforce.data-imaginist.com/reference/geom_parallel_sets.html) and [ggalluvial](https://corybrunson.github.io/ggalluvial/).
With the sc_Alluvial_pl function, generating Alluvial plots becomes simpler. Here’s a quick example:
```R
source("./Alluvial_stack_barplot.R")

# example for you to learn
data(majors)
majors$curriculum <- as.factor(majors$curriculum)
p <- sc_Alluvial_pl(majors, x_key = "semester", 
    order_x = c("CURR1", "CURR3", "CURR5", "CURR7", "CURR9", "CURR11", "CURR13", "CURR15"),
stratum_group = "curriculum")
p
```

The function also supports frequency and percentage inputs and is not limited to categorical data, making it particularly useful for FACS (Fluorescence-Activated Cell Sorting) data. 
For example, if you have pre-calculated percentages, you can use the following code:
```R
# if we already calculate percentage
df <- read.table(text = "6  0.5 4   0.95    0.001   0.001   0.0012  0.01    0.0013  0.0025  91  0.75    0.95
8.1 2.1 2.9 0.7 5   2   4.5 1   0.5 0.7 90  0.9 0.2
8.1 2.1 2.9 0.7 5   2   4.5 1   0.5 0.7 90  0.9 0.2
5.3 2.9 2.2 0.3 3.2 3.8 2.1 25  16  7   56  0.03    0.001
")
rownames(df) <- c("7w", "8w", "11w", "17w")
colnames(df) <- c("CD33", "CD14", "Granulocytes", "Neutrophils", "DCs", "M1", "M2", "CD3", "CD4", "CD8", "B cell", "NK cell", "Plasma cell")
library(tidyr)
library(dplyr)
df <- df %>% mutate(time = rownames(df))
huNSGS <- df %>%
  pivot_longer(!time, names_to = "Celltypes", values_to = "percentage")

p <- sc_Alluvial_pl(huNSGS,x_key = "time", order_x = c("7w", "8w", "11w", "17w"),
stratum_group = "Celltypes", percent = "percentage")
p
```
