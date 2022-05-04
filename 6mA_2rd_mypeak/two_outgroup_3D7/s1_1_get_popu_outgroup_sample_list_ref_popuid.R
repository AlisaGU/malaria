#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
ori_popu_symbol <- Args[6]
outgroup_symbol <- Args[7]

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd")
sample_metadata <- fread("sample_provenance_and_sequencing_metadata.txt", header = T, stringsAsFactors = F)
# 4. variable setting of test module--------------------------------------- TODO:
# ori_popu_symbol<-"OCE"
# outgroup_symbol<-"ref"
# 5. process -------------------------------------------------------------- TODO:
QC_pass_sample_data <- sample_metadata[sample_metadata$QCpass, ]

popu_symbol <- unlist(strsplit(ori_popu_symbol, "_"))
popu_QC_pass <- lapply(popu_symbol, function(x) {
    QC_pass_sample_data[QC_pass_sample_data$Population == x, ]
})
popu_QC_pass <- do.call(rbind, popu_QC_pass)

popu_sample <- popu_QC_pass$Sample

outgroup_sample <- NULL
if (outgroup_symbol != "ref") {
    outgroup_sample <- QC_pass_sample_data$Sample[which(QC_pass_sample_data$Population == outgroup_symbol)[1]]
}
write.table(c(outgroup_sample, popu_sample),
    quote = F, row.names = F, col.names = F,
    file = paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", ori_popu_symbol, "/", ori_popu_symbol, ".sample.list")
)
