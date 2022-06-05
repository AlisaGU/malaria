#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(data.table)

# 2. functions ------------------------------------------------------------


# 3. input ----------------------------------------------------------------
Args <- commandArgs()
intermediate_file <- Args[6]
motif_file <- Args[7]
central_motif_file <- Args[8]
# 4. variable setting of test module---------------------------------------
# intermediate_file <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/central_motif/Pf3D7_01_v3.intermediate"
# motif_file <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/motif_loc/Pf3D7_01_v3.motif.bed"
# central_motif_file <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/central_motif/Pf3D7_01_v3.centralMotif.bed"
# 5. process --------------------------------------------------------------
intermediate <- fread(intermediate_file, header = F, stringsAsFactors = F, select = 4)
motif <- fread(motif_file, header = F, stringsAsFactors = F)
row_index <- sapply(intermediate, function(x) {
    gsub("motif_", "", x)
})
result <- motif[as.numeric(row_index), ]
write.table(result, central_motif_file, row.names = F, col.names = F, sep = "\t", quote = F)
