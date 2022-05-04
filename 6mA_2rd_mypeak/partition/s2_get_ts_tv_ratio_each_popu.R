#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------


# 2. functions ------------------------------------------------------------


# 3. input ----------------------------------------------------------------
Args <- commandArgs()
popu_symbol<-Args[6]

# 4. variable setting of test module---------------------------------------


# 5. process --------------------------------------------------------------
working_dir<-paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/2_2rd_partition/intergenic/",popu_symbol,"/mutation_density")
setwd(working_dir)
ts<-read.table("transition_mean_variant_count",as.is=T)
tv<-read.table("transversion_mean_variant_count",as.is=T)
result<-sapply(1:ncol(ts),function(i){
    if(i==2|i==7){
        return(rep(NA,14))
    }else{
        ts[,i]/tv[,i]
    }
})
write.table(result,"ts_tv_ratio.txt",row.names=F,col.names=F,quote=F,sep="\t")