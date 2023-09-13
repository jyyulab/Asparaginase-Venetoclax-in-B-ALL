##Coded by Xin Huang(xhuang@stjude.org)
##R version 3.6.3 (2020-02-29)

library(NetBID2)
library(Seurat)
library("HCAData")
library(dplyr)
setwd('./B_HCA')

HCAData()
suppressPackageStartupMessages({
  library("ExperimentHub")
  library("SingleCellExperiment")
})

eh <- ExperimentHub()
query(eh, "HCAData")

# EH2047 | Human Cell Atlas - Census of Immune Cells, Bone marrow, 'dense matrix' format              
# EH2048 | Human Cell Atlas - Census of Immune Cells, Bone marrow, sample (column) annotation         
# EH2049 | Human Cell Atlas - Census of Immune Cells, Bone marrow, gene (row) annotation              
# EH2050 | Human Cell Atlas - Census of Immune Cells, Umbilical cord blood, 'dense matrix' format     
# EH2051 | Human Cell Atlas - Census of Immune Cells, Umbilical cord blood, sample (column) annotation
# EH2052 | Human Cell Atlas - Census of Immune Cells, Umbilical cord blood, gene (row) annotation  

# these three are the components to the bone marrow dataset
bonemarrow_h5densematrix <- eh[["EH2047"]]
bonemarrow_coldata <- eh[["EH2048"]]
bonemarrow_rowdata <- eh[["EH2049"]]
# and are put together when calling...
sce_bonemarrow <- HCAData("ica_bone_marrow")
save(sce_bonemarrow, file = 'sce_bonemarrow.rdata')


