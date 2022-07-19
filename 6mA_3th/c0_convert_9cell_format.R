#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- fread("3D7-6mA-9Cell.txt", header = T, stringsAsFactors = F)
attributes <- t(sapply(data$Attributes, function(x) {
    a <- unlist(strsplit(x, ";"))
    b <- unlist(lapply(a, function(y) {
        u <- unlist(strsplit(y, "="))
        result <- u[2]
        names(result) <- u[1]
        return(result)
    }))
    result <- c(b["coverage"], b["context"], b["IPDRatio"], b["frac"], b["fracLow"], b["fracUp"], b["identificationQv"])
    return(result)
}))
rownames(attributes) <- NULL
result <- cbind(data[, -9], attributes)
result1 <- result[, -c(2, 3)]
result1$Start <- result1$Start - 1
# write.table(result1, "3D7-6mA-9Cell.bed", sep = "\t", row.names = F, col.names = F, quote = F)
# write.table(names(result1), "column_names_of_3D7-6mA-9Cell.txt", sep = "\t", row.names = F, col.names = F, quote = F)

result2 <- result1[as.numeric(result1$identificationQv) > 20 & as.numeric(result1$coverage) > 10, ]
write.table(result2, "3D7-6mA-9Cell_QVgt20_COVgt10.bed", sep = "\t", row.names = F, col.names = F, quote = F)

result3 <- result1[as.numeric(result1$coverage) >= 25, ]
write.table(result3, "3D7-6mA-9Cell_COVgt25.bed", sep = "\t", row.names = F, col.names = F, quote = F)

result4 <- result1[as.numeric(result1$coverage) >= 50, ]
write.table(result4, "3D7-6mA-9Cell_COVgt50.bed", sep = "\t", row.names = F, col.names = F, quote = F)

result5 <- result1[as.numeric(result1$frac) >= 0.2 & as.numeric(result1$coverage) >= 25, ]
write.table(result5, "3D7-6mA-9Cell_FRACgt20_COVgt25.bed", sep = "\t", row.names = F, col.names = F, quote = F)
