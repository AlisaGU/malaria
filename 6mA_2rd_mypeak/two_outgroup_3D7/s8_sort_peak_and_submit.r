#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
peak_xls_name <- Args[6]
output_dir<-Args[7]
# output_dir<-"/home/gushanshan/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
peak_xls<-fread(cmd=paste0("grep \"^Pf3D7\""," ",peak_xls_name,"|sort -k1,1 -k2,2n -k7,7nr"),stringsAsFactors = F,header=F)
colnames(peak_xls)<-c("chr","start","end","length","abs_submit","pileup","log_p","fold_enrichment","log_q","name")
one_peak_one_submit<-peak_xls %>%
    group_by(chr, start, end) %>%
    filter(log_p==max(log_p))
one_peak_one_submit<-as.data.frame(one_peak_one_submit)

sort_peak<-one_peak_one_submit[order(one_peak_one_submit$log_p,decreasing = T),]
write.table(sort_peak[,c(1,2,3,5)],paste0(output_dir,"/sorted_peak.txt"),sep="\t",row.names=F,col.names=F,quote=F)
