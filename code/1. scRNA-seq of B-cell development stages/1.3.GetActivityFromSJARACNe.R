##Coded by Xin Huang(xhuang@stjude.org)
##R version 3.6.3 (2020-02-29)

library(NetBID2)
library(scMINER)

setwd('./B_HCA')
load('eset.filter.cell.log2')
acs_sc <- GetActivityFromSJARACNe( 
  SJARACNe_output_path ="./SJARACNE/",
  SJARACNe_input_eset = eset.filter.cell.log2,
  activity.method="unweighted", # we highly recommend using 'unweighted' as activity calculation method
  activity.norm=TRUE, 
  group_name = "celltype", # which group was used to partition expression profiles
  save_network_file=T, # whether or not save network for each group
  save_path="./networks") 