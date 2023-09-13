##Coded by Xin Huang(xhuang@stjude.org)
##R version 3.6.3 (2020-02-29)

library(NetBID2)


## CRLF2_non_Ph_like #############################################

network.dir <- './CRLF2-non-Ph-like'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='CRLF2(non-Ph-like)']
dim(analysis.par$cal.eset)
1
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

CRLF2_non_Ph_like.network <- analysis.par$merge.network
CRLF2_non_Ph_like.ac.eset <- analysis.par$merge.ac.eset 
save(CRLF2_non_Ph_like.network,file = 'CRLF2_non_Ph_like.network')
save(CRLF2_non_Ph_like.ac.eset,file = 'CRLF2_non_Ph_like.ac.eset')


## DUX4 #############################################

network.dir <- './DUX4'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='DUX4']
dim(analysis.par$cal.eset)
11
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

DUX4.network <- analysis.par$merge.network
DUX4.ac.eset <- analysis.par$merge.ac.eset 
save(DUX4.network,file = 'DUX4.network')
save(DUX4.ac.eset,file = 'DUX4.ac.eset')

## allETV6_RUNX1 #############################################

network.dir <- './allETV6-RUNX1'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='ETV6-RUNX1'|ASP.eset$subtypes=='ETV6-RUNX1-like(ETV6r)']
dim(analysis.par$cal.eset)
20
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

allETV6_RUNX1.network <- analysis.par$merge.network
allETV6_RUNX1.ac.eset <- analysis.par$merge.ac.eset 
save(allETV6_RUNX1.network,file = 'allETV6_RUNX1.network')
save(allETV6_RUNX1.ac.eset,file = 'allETV6_RUNX1.ac.eset')

## Hyperdiploid #############################################

network.dir <- './High_hyperdiploid'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='Hyperdiploid']
dim(analysis.par$cal.eset)
20
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

Hyperdiploid.network <- analysis.par$merge.network
Hyperdiploid.ac.eset <- analysis.par$merge.ac.eset 
save(Hyperdiploid.network,file = 'Hyperdiploid.network')
save(Hyperdiploid.ac.eset,file = 'Hyperdiploid.ac.eset')


## iAMP21 #############################################

network.dir <- './iAMP21'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='iAMP21']
dim(analysis.par$cal.eset)
3
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

iAMP21.network <- analysis.par$merge.network
iAMP21.ac.eset <- analysis.par$merge.ac.eset 
save(iAMP21.network,file = 'iAMP21.network')
save(iAMP21.ac.eset,file = 'iAMP21.ac.eset')

## allKMT2A #############################################

network.dir <- './allKMT2A'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='KMT2A']
dim(analysis.par$cal.eset)

32
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

KMT2A.network <- analysis.par$merge.network
KMT2A.ac.eset <- analysis.par$merge.ac.eset 
save(KMT2A.network,file = 'KMT2A.network')
save(KMT2A.ac.eset,file = 'KMT2A.ac.eset')


## Low_hypodiploid #############################################

network.dir <- './Low_hypodiploid'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='Low hypodiploid']
dim(analysis.par$cal.eset)
7
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

Low_hypodiploid.network <- analysis.par$merge.network
Low_hypodiploid.ac.eset <- analysis.par$merge.ac.eset 
save(Low_hypodiploid.network,file = 'Low_hypodiploid.network')
save(Low_hypodiploid.ac.eset,file = 'Low_hypodiploid.ac.eset')


## MEF2D #############################################

network.dir <- './MEF2D'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='MEF2D']
dim(analysis.par$cal.eset)
5
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean',memory_constrain = T)
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

MEF2D.network <- analysis.par$merge.network
MEF2D.ac.eset <- analysis.par$merge.ac.eset 
save(MEF2D.network,file = 'MEF2D.network')
save(MEF2D.ac.eset,file = 'MEF2D.ac.eset')


## Near_haploid #############################################

network.dir <- './Near_haploid'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='Near haploid']
dim(analysis.par$cal.eset)
3
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean',memory_constrain = T)
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

Near_haploid.network <- analysis.par$merge.network
Near_haploid.ac.eset <- analysis.par$merge.ac.eset 
save(Near_haploid.network,file = 'Near_haploid.network')
save(Near_haploid.ac.eset,file = 'Near_haploid.ac.eset')

## PAX5_P80R #############################################

network.dir <- './PAX5_P80R'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='PAX5 P80R']
dim(analysis.par$cal.eset)
2
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

PAX5_P80R.network <- analysis.par$merge.network
PAX5_P80R.ac.eset <- analysis.par$merge.ac.eset 
save(PAX5_P80R.network,file = 'PAX5_P80R.network')
save(PAX5_P80R.ac.eset,file = 'PAX5_P80R.ac.eset')

## PAX5alt #############################################

network.dir <- './PAX5alt'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='PAX5alt']
dim(analysis.par$cal.eset)
7
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

PAX5alt.network <- analysis.par$merge.network
PAX5alt.ac.eset <- analysis.par$merge.ac.eset 
save(PAX5alt.network,file = 'PAX5alt.network')
save(PAX5alt.ac.eset,file = 'PAX5alt.ac.eset')

## ZNF384 #############################################

network.dir <- './ZNF384'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='ZNF384']
dim(analysis.par$cal.eset)
4
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

ZNF384.network <- analysis.par$merge.network
ZNF384.ac.eset <- analysis.par$merge.ac.eset 
save(ZNF384.network,file = 'ZNF384.network')
save(ZNF384.ac.eset,file = 'ZNF384.ac.eset')


## TCF3_PBX1 #############################################

network.dir <- './TCF3_PBX1'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes=='TCF3-PBX1']
dim(analysis.par$cal.eset)
4
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

TCF3_PBX1.network <- analysis.par$merge.network
TCF3_PBX1.ac.eset <- analysis.par$merge.ac.eset 
save(TCF3_PBX1.network,file = 'TCF3_PBX1.network')
save(TCF3_PBX1.ac.eset,file = 'TCF3_PBX1.ac.eset')


## allPh #############################################

network.dir <- './allPh'
network.project.name <- ''
rm(analysis.par)
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

analysis.par$cal.eset <- ASP.eset[,ASP.eset$subtypes%in%c('Ph','Ph-like_CRLF2 ','Ph-like_non_CRLF2')]
dim(analysis.par$cal.eset)
67
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

allPh.network <- analysis.par$merge.network
allPh.ac.eset <- analysis.par$merge.ac.eset 
save(allPh.network,file = 'allPh.network')
save(allPh.ac.eset,file = 'allPh.ac.eset')




