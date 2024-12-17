#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
suppressWarnings(suppressMessages(library(preprocessCore)))

# 2. functions ------------------------------------------------------------ TODO:
detect_operating_system <- function() {
    switch(Sys.info()[["sysname"]],
        Windows = {
            return("Windows")
        },
        Linux = {
            return("Linux")
        }
    )
}

path_convert <- function(path = NULL) {
    result <- paste0("I:", gsub("\\picb\\evolgen\\users\\gushanshan", "", gsub("/", "\\", path, fixed = T), fixed = T))
    return(result)
}

data_download_parse <- function(GSEid = NULL, output = NULL) {
    suppressWarnings(suppressMessages(library("GEOquery")))
    suppressWarnings(suppressMessages(library("AnnoProbe")))

    eSet <- geoChina(GSEid)

    pdf(gsub("Rdata", "pdf", output), width = 7, height = 7)

    GSEdata <- lapply(eSet, function(sub_eSet) {
        exp <- exprs(sub_eSet)
        exp_range <- range(exp, na.rm = TRUE)
        boxplot(exp)
        pd <- pData(sub_eSet)
        metaInfo_exp_order_identical <- identical(rownames(pd), colnames(exp))
        array_platform <- sub_eSet@annotation
        gpl_probe_map <- read.delim(dir()[grepl(array_platform, dir())], comment.char = "#", header = TRUE)
        result <- list(
            exp = exp, exp_range = exp_range, metaInfo = pd,
            metaInfo_exp_order_identical = metaInfo_exp_order_identical,
            array_platform = array_platform, gpl_probe_map = gpl_probe_map
        )
        return(result)
    })
    save(GSEdata, file = output)

    dev.off()
}
# 3. input ---------------------------------------------------------------- TODO:
server_path <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/drugResistance/paper_SachelMok"
GSEid <- "GSE25883"
data_repository <- paste0(GSEid, ".exp_meta_probe.Rdata")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd(ifelse(detect_operating_system() == "Windows", path_convert(server_path), server_path))

# 下载和解析数据只运行一次，之后注释掉。这一步在windows上进行，因此将专属于这一块的R包载入封装到函数里
# exp_meta_probe <- data_download_parse(GSEid = GSEid, output = data_repository)

# exp数据的过滤和probe-gene的映射
load(data_repository)
exp <- GSEdata[[1]]$exp
normExp <- normalize.quantiles(exp)
rownames(normExp) <- rownames(exp)
colnames(normExp) <- colnames(exp)

pdf(gsub("Rdata", "normalized.pdf", data_repository), width = 7, height = 7)
boxplot(normExp)
dev.off()
