##Coded by Xin Huang(xhuang@stjude.org)
##R version 3.6.3 (2020-02-29)

library(NetBID2)

###################################### ZNF384 #########################################
project_main_dir <- './networks'

project_name <- 'ZNF384'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)
net_eset <- BALL.eset[,BALL.eset$primary.subtype=='ZNF384']
dim(net_eset)
49
# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
# 72 19595 
###################################### PAX5alt #########################################
project_name <- 'PAX5alt'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='PAX5alt']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
#   1 19666 

###################################### PAX5_P80R #########################################
project_name <- 'PAX5_P80R'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='PAX5 P80R']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
# 22 19645 
##################################### NUTM1 #########################################
project_name <- 'NUTM1'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='NUTM1']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

# phe <- pData(network.par$net.eset)
# use.samples <- rownames(phe) ## use all samples, or choose to use some samples
# # 
# use_gene_type <- 'external_gene_name' # user-defined
# use_genes <- rownames(fData(network.par$net.eset))
# use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
#
# SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
#                  TF_list=use_list$tf,SIG_list=use_list$sig,
#                  IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
#                  SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# # 
# # use_vec
# # FALSE  TRUE 
# # 31203 19166  

#################################### Near haploid  #########################################
project_name <- 'Near_haploid '
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='Near haploid']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
# 253 19414 

#################################### MEF2D #########################################
project_name <- 'MEF2D '
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='MEF2D']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
#  16 19651 

#################################### Low hypodiploid  #########################################
project_name <- 'Low_hypodiploid'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='Low hypodiploid']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
# 848 18819 

#################################### Low hyperdiploid  #########################################
project_name <- 'Low_hyperdiploid'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='Low hyperdiploid']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
#  66 19601 

#################################### allKMT2A #########################################
project_name <- 'allKMT2A'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype%in%c('KMT2A','KMT2A-like')]
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
# 20 19647 

#################################### IKZF1 N159Y #########################################
project_name <- 'IKZF1_N159Y'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='IKZF1 N159Y']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

# phe <- pData(network.par$net.eset)
# use.samples <- rownames(phe) ## use all samples, or choose to use some samples
# # 
# use_gene_type <- 'external_gene_name' # user-defined
# use_genes <- rownames(fData(network.par$net.eset))
# use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
#
# SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
#                  TF_list=use_list$tf,SIG_list=use_list$sig,
#                  IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
#                  SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# # 
# # use_vec
# # FALSE  TRUE 
# # 30800 19548 

#################################### iAMP21 #########################################
project_name <- 'iAMP21'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='iAMP21']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
#   131 19536 


#################################### HLF #########################################
project_name <- 'HLF'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='HLF']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

# phe <- pData(network.par$net.eset)
# use.samples <- rownames(phe) ## use all samples, or choose to use some samples
# # 
# use_gene_type <- 'external_gene_name' # user-defined
# use_genes <- rownames(fData(network.par$net.eset))
# use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
#
# SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
#                  TF_list=use_list$tf,SIG_list=use_list$sig,
#                  IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
#                  SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# # 
# # use_vec
# # FALSE  TRUE 
# # 29735 18985

#################################### High hyperdiploid  #########################################
project_name <- 'High_hyperdiploid'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='High hyperdiploid']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
# 135 19532 

#################################### allETV6-RUNX1  #########################################
project_name <- 'allETV6-RUNX1'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype%in%c('ETV6-RUNX1','ETV6-RUNX1-like')]
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
#    13 19654 


#################################### DUX4  #########################################
project_name <- 'DUX4'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='DUX4']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
#     1 19666 


#################################### CRLF2(non-Ph-like)  #########################################
project_name <- 'CRLF2-non-Ph-like'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='CRLF2(non-Ph-like)']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
#   123 19544 


#################################### BCL2/MYC  #########################################
project_name <- 'BCL2_MYC'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

net_eset <- BALL.eset[,BALL.eset$primary.subtype=='BCL2/MYC']
dim(net_eset)

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
# 
# use_vec
# FALSE  TRUE 
#   949 18718 

#################################### TCF3_PBX1 #########################################
project_name <- 'TCF3_PBX1'
## network.par (is an object) is very essential in the analysis!! main dir have several projects
if(exists('network.par')==TRUE) rm(network.par)

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

table(BALL.eset$primary.subtype)
BCL2/MYC CRLF2(non-Ph-like)               DUX4         ETV6-RUNX1 
18                 16                106                186 
ETV6-RUNX1-like  High hyperdiploid                HLF             iAMP21 
42                279                  9                 40 
IKZF1 N159Y              KMT2A         KMT2A-like   Low hyperdiploid 
8                135                  5                 51 
Low hypodiploid              MEF2D       Near haploid              NUTM1 
78                 43                 29                 11 
Other          PAX5 P80R            PAX5alt                 Ph 
125                 44                148                123 
Ph-like          TCF3-PBX1             ZNF384        ZNF384-like 
359                 78                 49                  4 
net_eset <- BALL.eset[,BALL.eset$primary.subtype=='TCF3-PBX1']
dim(net_eset)
19667       78 

# add the variable into network.par, very essential !
network.par$net.eset <- net_eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup='RNA-seq.library',do.logtransform=FALSE,prefix='beforeQC_',
             pca_plot_type='2D.interactive', choose_plot = c( "pca", "density"))

phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) ## use all samples, or choose to use some samples

use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0,IQR.loose_thre = 0,  ##user defined
                 SJAR.project_name='',SJAR.main_dir=network.par$out.dir.SJAR)
use_vec
FALSE  TRUE 
53 19614
