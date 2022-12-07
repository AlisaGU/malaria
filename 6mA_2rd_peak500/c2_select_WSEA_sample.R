#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd")
sample_metadata <- fread("sample_provenance_and_sequencing_metadata.txt", header = T, stringsAsFactors = F)

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
QC_pass_sample_data <- sample_metadata[sample_metadata$QCpass, ]
WESA_QC_pass <- QC_pass_sample_data[QC_pass_sample_data$Population == "WSEA", ]
write.table(WESA_QC_pass$Sample,
    quote = F, row.names = F, col.names = F,
    file = "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_motif/WSEA.sample.list"
)