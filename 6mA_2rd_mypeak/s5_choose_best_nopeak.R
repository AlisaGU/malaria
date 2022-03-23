#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
filename <- Args[6]
outfile_name <- Args[7]
# 4. variable setting of test module--------------------------------------- TODO:
# filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/a"

# 5. process -------------------------------------------------------------- TODO:
file <- fread(filename, stringsAsFactors = F, sep = "\t")
colnames(file)[c(4, ncol(file))] <- c("name", "width")

splited_peak <- split(file, f = file$name)
result <- lapply(splited_peak, function(data) {
    max_width_index <- which(data$width == max(data$width))
    return(data[max_width_index, ])
})
result <- do.call(rbind, result)
write.table(result, outfile_name, quote = F, sep = "\t", col.names = F, row.names = F)