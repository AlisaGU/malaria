#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
suppressWarnings(suppressMessages(library(preprocessCore)))
library(ggplot2)
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
GSEid <- "GSE59099"
data_repository <- paste0(GSEid, ".exp_meta_probe.Rdata")
probes <- c("aMAL13P1.255_0", "aMAL13P1.255_1") # 在GPL18893平台上，PF3D7_1350700基因对应两个探针：aMAL13P1.255_0和aMAL13P1.255_1
actinII_probes <- c("aPF14_0124_0", "aPF14_0124_1")
actinI_probes <- c("aPFL2215w_0", "aPFL2215w_1")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd(ifelse(detect_operating_system() == "Windows", path_convert(server_path), server_path))

# 1. 下载和解析数据只运行一次，之后注释掉。这一步在windows上进行，因此将专属于这一块的R包载入封装到函数里
# exp_meta_probe <- data_download_parse(GSEid = GSEid, output = data_repository)

# 2. 分析
load(data_repository)

exp <- GSEdata[[1]]$exp
normExp <- exp
# normExp <- normalize.quantiles(exp)
# rownames(normExp) <- rownames(exp)
# colnames(normExp) <- colnames(exp)

metainfo <- GSEdata[[1]]$metaInfo
write.table(metainfo, paste0(GSEid, ".metainfo.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
## 3. 根据样本来源将clearance time分组：invivo和patient。这两个分组来自临床信息。暂时不考虑第三小节的分析
{
    parasite_clearance_halflife <- metainfo[, c(
        "geo_accession",
        "parasite clearance halflife upon artemisinin treatment (h):ch1",
        "parasite clearance halflife upon artemisinin treatment in patient:ch1"
    )]

    colnames(parasite_clearance_halflife) <- c("geo_accession", "halflife_inVivo", "halflife_patient")
    ## 3.1 in vitro
    inVivo_sample_index <- which((!is.na(parasite_clearance_halflife[, "halflife_inVivo"])) & (parasite_clearance_halflife[, "halflife_inVivo"] != "NA")) # 1025个个体
    inVivo_meta <- parasite_clearance_halflife[inVivo_sample_index, c("geo_accession", "halflife_inVivo")]
    inVivo_meta$"halflife_inVivo" <- as.numeric(gsub("h", "", inVivo_meta$"halflife_inVivo"))
    inVivo_meta$Art <- sapply(inVivo_meta$halflife_inVivo, function(x) {
        ifelse(x <= 5, "sensitive", "resistant")
    }) # resistant:298, sensitive:727
    inVivo_meta$geo_loc <- metainfo[match(inVivo_meta$geo_accession, metainfo$geo_accession), "characteristics_ch1.3"]

    normExp_inVivo_sensitive <- normExp[
        match(probes, rownames(normExp)),
        match(inVivo_meta$geo_accession[inVivo_meta$Art == "sensitive"], colnames(normExp))
    ]
    normExp_inVivo_sensitive_gene <- colMeans(normExp_inVivo_sensitive)

    normExp_inVivo_resistant <- normExp[
        match(probes, rownames(normExp)),
        match(inVivo_meta$geo_accession[inVivo_meta$Art == "resistant"], colnames(normExp))
    ]
    normExp_inVivo_resistant_gene <- colMeans(normExp_inVivo_resistant)
    t.test(normExp_inVivo_sensitive_gene, normExp_inVivo_resistant_gene, alternative = "less")
    ### 3.1.1
    geo_locs <- unique(inVivo_meta$geo_loc)
    geo_loc_info <- c()
    for (geo_loc in geo_locs) {
        if (geo_loc != "geographic origin: Ramu, Bangladesh") {
            inVivo_meta_geo <- inVivo_meta[which(geo_loc == inVivo_meta$geo_loc), ]
            normExp_inVivo_geo_sensitive <- normExp[
                match(probes, rownames(normExp)),
                match(inVivo_meta_geo$geo_accession[inVivo_meta_geo$Art == "sensitive"], colnames(normExp))
            ]
            normExp_inVivo_geo_sensitive_gene <- colMeans(normExp_inVivo_geo_sensitive)

            normExp_inVivo_geo_resistant <- normExp[
                match(probes, rownames(normExp)),
                match(inVivo_meta_geo$geo_accession[inVivo_meta_geo$Art == "resistant"], colnames(normExp))
            ]
            normExp_inVivo_geo_resistant_gene <- colMeans(normExp_inVivo_geo_resistant)

            p <- t.test(normExp_inVivo_geo_sensitive_gene, normExp_inVivo_geo_resistant_gene, alternative = "less")$p.value
            geo_loc_info <- rbind(geo_loc_info, c(gsub("geographic origin: ", "", geo_loc), ncol(normExp_inVivo_geo_sensitive), ncol(normExp_inVivo_geo_resistant), p))
        }
    }
    geo_loc_info <- as.data.frame(geo_loc_info)
    colnames(geo_loc_info) <- c("geographical_location", "count of sensitive sample", "count of resistant sample", "p value")
    ## 3.2 patient
    patient_sample_index <- which((!is.na(parasite_clearance_halflife[, "halflife_patient"])) & (parasite_clearance_halflife[, "halflife_patient"] != "NA")) # 110个个体
    patient_meta <- parasite_clearance_halflife[patient_sample_index, c("geo_accession", "halflife_patient")]
    patient_meta$"halflife_patient" <- as.numeric(gsub("h", "", patient_meta$"halflife_patient"))
    patient_meta$Art <- sapply(patient_meta$halflife_patient, function(x) {
        ifelse(x <= 5, "sensitive", "resistant")
    }) # resistant:82, sensitive:28

    normExp_patient_sensitive <- normExp[
        match(probes, rownames(normExp)),
        match(patient_meta$geo_accession[patient_meta$Art == "sensitive"], colnames(normExp))
    ]
    normExp_patient_sensitive_gene <- colMeans(normExp_patient_sensitive)

    normExp_patient_resistant <- normExp[
        match(probes, rownames(normExp)),
        match(patient_meta$geo_accession[patient_meta$Art == "resistant"], colnames(normExp))
    ]
    normExp_patient_resistant_gene <- colMeans(normExp_patient_resistant)
    t.test(normExp_patient_sensitive_gene, normExp_patient_resistant_gene, alternative = "less")
}

## 4.根据文章中给出的group信息分组，我自己的操作。暂时不考虑
{
    group_info <- read.table("GSE59099_group.txt", header = T, as.is = T, sep = "\t")
    splitted_info <- lapply(split(group_info, f = as.factor(group_info[, "Kmeans.Grp"])), function(x) {
        split(x, f = x[, "Asexual.stage..hpi."])
    })

    metainfo[, "characteristics_ch1.2"] <- gsub("patient sample id: ", "", metainfo[, "characteristics_ch1.2"])
    metainfo[, "parasite clearance halflife"] <- apply(metainfo, 1, function(x) {
        invivo <- x["parasite clearance halflife upon artemisinin treatment (h):ch1"]
        patient <- x["parasite clearance halflife upon artemisinin treatment in patient:ch1"]
        if (!is.na(invivo)) {
            return(invivo)
        } else if (!is.na(patient)) {
            return(patient)
        } else {
            return(NA)
        }
    })
    metainfo[, "parasite clearance halflife"] <- gsub("h", "", metainfo[, "parasite clearance halflife"])
    a <- sapply(splitted_info, function(y) {
        sapply(y, function(x) {
            Patient_code <- x[, "Patient.Code"]
            group_i_info <- metainfo[match(Patient_code, metainfo[, "characteristics_ch1.2"]), c("geo_accession", "parasite clearance halflife upon artemisinin treatment (h):ch1", "parasite clearance halflife upon artemisinin treatment in patient:ch1", "parasite clearance halflife")]
            group_i_info$Art <- sapply(group_i_info[, "parasite clearance halflife"], function(x) {
                ifelse(x <= 5, "sensitive", "resistant")
            })

            normExp_sensitive <- normExp[
                match(probes, rownames(normExp)),
                match(group_i_info$geo_accession[group_i_info$Art == "sensitive"], colnames(normExp))
            ]

            normExp_resistant <- normExp[
                match(probes, rownames(normExp)),
                match(group_i_info$geo_accession[group_i_info$Art == "resistant"], colnames(normExp))
            ]

            if (is.matrix(normExp_sensitive) & is.matrix(normExp_resistant) & NCOL(normExp_sensitive) > 0 & NCOL(normExp_resistant) > 0) {
                normExp_sensitive_gene <- colMeans(normExp_sensitive)
                normExp_resistant_gene <- colMeans(normExp_resistant)
                if (length(which(!is.na(normExp_sensitive_gene))) > 1 & length(which(!is.na(normExp_resistant_gene))) > 1) {
                    p <- t.test(normExp_sensitive_gene, normExp_resistant_gene, alternative = "less")$p.value
                    print(c(unique(x[, "Kmeans.Grp"]), unique(x[, "Asexual.stage..hpi."]), ncol(normExp_sensitive), ncol(normExp_resistant), p))
                }
            }
        })
    })
}

## 5. 按照XB的操作重新做一遍
sample_clearanceTime <- metainfo[, c("geo_accession", "characteristics_ch1.6")]
colnames(sample_clearanceTime) <- c("sample", "clearanceTime")

a <- gsub("parasite clearance halflife upon artemisinin treatment (h): ", "", sample_clearanceTime$clearanceTime, fixed = T)
b <- gsub("parasite clearance halflife upon artemisinin treatment in patient: ", "", a, fixed = T)
sample_clearanceTime$clearanceTime <- gsub("h", "", b, fixed = T)
na_index <- which(sample_clearanceTime$clearanceTime == "NA") # 18个样本的清除时间为NA，排除这18个样本
sample_clearanceTime <- sample_clearanceTime[-na_index, ]
sample_clearanceTime$clearanceTime <- as.numeric(sample_clearanceTime$clearanceTime)
sample_clearanceTime$reponse <- sapply(sample_clearanceTime$clearanceTime, function(x) {
    ifelse(x <= 5, "sensitive", "resistant")
})

namt_normExp <- normExp[match(probes, rownames(normExp)), match(sample_clearanceTime$sample, colnames(normExp))]
namt_normExp_means <- colMeans(namt_normExp)

actinI_normExp <- normExp[match(actinI_probes, rownames(normExp)), match(sample_clearanceTime$sample, colnames(normExp))]
actinI_normExp_means <- colMeans(actinI_normExp)

# log2ratio<-log2(namt_normExp_means/actinI_normExp_means)
data_for_plot <- data.frame(t(rbind(namt_normExp, namt_normExp_means, actinI_normExp, actinI_normExp_means)))
data_for_plot$response <- sample_clearanceTime$reponse[match(rownames(sample_clearanceTime), sample_clearanceTime$sample)]
data_for_plot$response <- factor(data_for_plot$response, levels = c("resistant", "sensitive"))

data_for_plot$namt_plus_actin <- data_for_plot$namt_normExp_means - data_for_plot$actinI_normExp_means
data_for_plot$log2ratio <- log2(data_for_plot$namt_normExp_means / data_for_plot$actinI_normExp_means)

ggplot(data_for_plot, aes(x = response, y = namt_normExp_means)) +
    geom_violin() +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 36, color = "black"),
        legend.position = "none"
    )
t.test(data_for_plot$namt_normExp_means[data_for_plot$response == "resistant"],
    data_for_plot$namt_normExp_means[data_for_plot$response == "sensitive"],
    alternative = "greater"
)

ggplot(data_for_plot, aes(x = response, y = actinI_normExp_means)) +
    geom_violin() +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 36, color = "black"),
        legend.position = "none"
    )
t.test(data_for_plot$actinI_normExp_means[data_for_plot$response == "resistant"],
    data_for_plot$actinI_normExp_means[data_for_plot$response == "sensitive"],
    alternative = "greater"
)

ggplot(data_for_plot, aes(x = response, y = namt_plus_actin)) +
    geom_violin() +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 36, color = "black"),
        legend.position = "none"
    )
t.test(data_for_plot$namt_plus_actin[data_for_plot$response == "resistant"],
    data_for_plot$namt_plus_actin[data_for_plot$response == "sensitive"],
    alternative = "greater"
)

ggplot(data_for_plot, aes(x = response, y = log2ratio)) +
    geom_violin() +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 36, color = "black"),
        legend.position = "none"
    )
t.test(data_for_plot$log2ratio[data_for_plot$response == "resistant"],
    data_for_plot$log2ratio[data_for_plot$response == "sensitive"],
    alternative = "greater"
)
