

##Coded by Xin Huang(xhuang@stjude.org)
##R version 3.6.3 (2020-02-29)

library(NetBID2)

## 1.load.exp.RNASeq.demoSalmon ####
Gene.eset <-load.exp.RNASeq.demoSalmon(salmon_dir = './salmon', tx2gene = NULL,
                                       use_phenotype_info = pd, use_sample_col = "RNAseq_SJID",
                                       use_design_col = "subtype", return_type = "eset", merge_level = "gene")
save(Gene.eset, file = "Gene.eset")

## 2.ID conversion ######
db.preload() 
table<- get_IDtransfer(from_type = 'ensembl_gene_id_version', to_type = 'external_gene_name', use_genes = fData(eset)$gene, ignore_version = TRUE)
eset.update<- update_eset.feature(use_eset = eset, use_feature_info = fData(eset),from_feature = 'gene',to_feature = 'symbol',merge_method = 'median')

## 3.QC ######
fData(eset)$IQR <- apply(exprs(eset), 1, IQR) ## add IQR to fd for further filtration
eset.sel.0.75 <- eset[fData(eset)$IQR >= quantile(fData(eset)$IQR, 0.75),] ## remove genes with small IQR
draw.eset.QC (eset.sel.0.75, outdir = "QC/IQR0.75",
              intgroup = c("subtype"),choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), prefix = '2D')
