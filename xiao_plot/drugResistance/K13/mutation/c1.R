#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library("VariantAnnotation")
library("GenomicFeatures")
library("dplyr")
library("data.table")
# 2. functions ------------------------------------------------------------ TODO:
my_readvcf<-function(vcf_filename=NULL,header_filename=NULL){
  vcf<-fread(vcf_filename,header=FALSE,stringsAsFactors = F)
  header<-unlist(fread(header_filename,header=FALSE,stringsAsFactors = F))
  header[1]<-"CHROM"
  colnames(vcf)<-header
  vcf<-vcf %>% select(!c(QUAL,FILTER,INFO,FORMAT))
  
  GT<-lapply(vcf[,6:ncol(vcf)],function(x){
    k<-sapply(x,function(y){unlist(strsplit(unlist(y),":"))[1]})
    k1<-t(sapply(k,function(y){unlist(strsplit(y,"/"))}))
    return(k1)
  })
  
  GT<-do.call(cbind,GT)
  colnames(GT)<-paste(rep(colnames(vcf)[6:ncol(vcf)],each=2),rep(1:2,times=ncol(vcf)-5),sep=".")
  result<-data.frame(vcf$POS,GT)
  return(result)
}

# 3. input ---------------------------------------------------------------- TODO:


# 4. variable setting of test module--------------------------------------- TODO:
vcf_filename<-"I:/projects/malaria/dataAndResult/xiao_plot/drugResistance/K13/variant_annotation/K13.vcf"
header_filename<-"I:/projects/malaria/dataAndResult/xiao_plot/drugResistance/K13/variant_annotation/header"
feizhou_samplelist<-"I:/projects/malaria/dataAndResult/1_2rd_initial_evaluation/WAF_CAF_EAF/WAF_CAF_EAF.sample.list"

WSEA_samplelist<-"I:/projects/malaria/dataAndResult/1_2rd_initial_evaluation/WSEA/WSEA.sample.list"
ESEA_samplelist<-"I:/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA/ESEA.sample.list"
SEA_samplelist<-"I:/projects/malaria/dataAndResult/1_2rd_initial_evaluation/SEA.sample.list"
# 5. process -------------------------------------------------------------- TODO:
## build pf BSgenome database

d<-makeTxDbFromGFF("I:/projects/malaria/dataAndResult/xiao_plot/drugResistance/K13/variant_annotation/K13.anno.gff")
b<-readVcf("I:/projects/malaria/dataAndResult/xiao_plot/drugResistance/K13/variant_annotation/K13.complete.vcf")

a<-FaFile("I:\\projects\\malaria\\dataAndResult\\ref_genome\\pf_3D7\\genome\\PlasmoDB-36_Pfalciparum3D7_Genome.fasta")
FaFileList(a)
open(a)
coding <- predictCoding(b,d, a)
close(a)

coding_df<-data.frame(seqname=seqnames(coding),start=start(coding),end=end(coding),
                      ref_dna_inGenomePosStrand=coding$REF,
                      alt_dna_inGenomePosStrand=sapply(coding$ALT,function(x){paste(x,collapse =", ")}),filter_ornot=coding$FILTER,
                      varAllele_inTranscriptDirec=sapply(coding$varAllele,function(x){paste(x,collapse =", ")}),
                      proteinLOC=sapply(coding$PROTEINLOC,function(x){paste(x,collapse =", ")}),
                      consequnece=coding$CONSEQUENCE,
                      refcodon=sapply(coding$REFCODON,function(x){paste(x,collapse =", ")}),
                      varcodon=sapply(coding$VARCODON,function(x){paste(x,collapse =", ")}),refAA=coding$REFAA,varAA=coding$VARAA)

coding_df_highQual<-coding_df %>% filter(filter_ornot=="PASS")%>% select(!c(end,filter_ornot)) 
coding_df_highQual$id<-sapply(1:nrow(coding_df_highQual),function(x){which(x==which(coding_df_highQual$start[x]==coding_df_highQual$start))})

GT<-my_readvcf(vcf_filename =vcf_filename,header_filename=header_filename )
## 非洲
sample<-unlist(read.table(feizhou_samplelist,header=F,as.is=T))
sample<-paste(gsub("-",".",sample),".",sep="")
o<-sapply(sample,function(x){grep(x,colnames(GT),fixed=T)})
GT_feizhou<-GT[,c(1,unlist(o))]
coding_df_highQual$Africa<-sapply(1:nrow(coding_df_highQual),function(x){
  id<-coding_df_highQual$id[x]
  start<-coding_df_highQual$start[x]
  fraction<-length(which(GT_feizhou[which(GT_feizhou$vcf.POS==start),-1]==id))/(ncol(GT_feizhou)-1)
  return(fraction)
})

## WSEA
sample<-unlist(read.table(WSEA_samplelist,header=F,as.is=T))
sample<-paste(gsub("-",".",sample),".",sep="")
o<-sapply(sample,function(x){grep(x,colnames(GT),fixed=T)})
GT_WSEA<-GT[,c(1,unlist(o))]
coding_df_highQual$WSEA<-sapply(1:nrow(coding_df_highQual),function(x){
  id<-coding_df_highQual$id[x]
  start<-coding_df_highQual$start[x]
  fraction<-length(which(GT_WSEA[which(GT_WSEA$vcf.POS==start),-1]==id))/(ncol(GT_WSEA)-1)
  return(fraction)
})


## ESEA
sample<-unlist(read.table(ESEA_samplelist,header=F,as.is=T))
sample<-paste(gsub("-",".",sample),".",sep="")
o<-sapply(sample,function(x){grep(x,colnames(GT),fixed=T)})
GT_ESEA<-GT[,c(1,unlist(o))]
coding_df_highQual$ESEA<-sapply(1:nrow(coding_df_highQual),function(x){
  id<-coding_df_highQual$id[x]
  start<-coding_df_highQual$start[x]
  fraction<-length(which(GT_ESEA[which(GT_ESEA$vcf.POS==start),-1]==id))/(ncol(GT_ESEA)-1)
  return(fraction)
})

## SEA
sample<-unlist(read.table(SEA_samplelist,header=F,as.is=T))
sample<-paste(gsub("-",".",sample),".",sep="")
o<-sapply(sample,function(x){grep(x,colnames(GT),fixed=T)})
GT_SEA<-GT[,c(1,unlist(o))]
coding_df_highQual$SEA<-sapply(1:nrow(coding_df_highQual),function(x){
  id<-coding_df_highQual$id[x]
  start<-coding_df_highQual$start[x]
  fraction<-length(which(GT_SEA[which(GT_SEA$vcf.POS==start),-1]==id))/(ncol(GT_SEA)-1)
  return(fraction)
})
## summary
write.table(coding_df_highQual,"I:/projects/malaria/dataAndResult/xiao_plot/drugResistance/K13/variant_annotation/variation_fraction_in_population.txt",sep="\t",col.names=T,row.names = F,quote=F)
