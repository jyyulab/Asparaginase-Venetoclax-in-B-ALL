##Coded by Xin Huang(xhuang@stjude.org)
##R version 3.6.3 (2020-02-29)

scminer
suerrat
setwd('./HCA/scMINer')

# 1_seurat2eset #######
load('./B_sub.rdata')
#expression_data
expression_data<- B_sub@assays$RNA@data 
#cell_metadata: 
data<-B_sub@meta.data
data$cellName<-data$orig.ident
cell_metadata<-new('AnnotatedDataFrame', data = data)
#gene_annotation: 
gene_annotation <- data.frame(gene_short_name =as.character (row.names(expression_data)), geneSymbol=as.character (row.names(expression_data)))
rownames(gene_annotation)<-gene_annotation$gene_short_name
gene_annotation$geneSymbol<-gene_annotation$gene_short_name
gene_annotation$nCells<-Matrix::rowSums(expression_data!=0)

#new eset
eset <- CreateSparseEset(data=expression_data, meta.data = cell_metadata@data, feature.data = gene_annotation, add.meta = F)
# filter out genes expressed in less than 50cells
eset.filter.cell.log2 <- eset[fData(eset)$nCells>50,] 

# 2_generateSJARACNeInput #######

# ref data is tf_sigs_hg.RData
generateSJARACNeInput(
  input_eset = eset.filter.cell.log2, funcType = NULL, 
  ref = "hg",  #human
  wd.src = "SJARACNE",  #Output directory
  group_name = "celltype")







