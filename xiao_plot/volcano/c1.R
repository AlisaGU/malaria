#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(openxlsx)
library(dplyr)
library(ggrepel)
library(data.table)
# 2. functions ------------------------------------------------------------ TODO:


# 3. variable setting of test module--------------------------------------- TODO:


# 4. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/volcano")

# 5. process -------------------------------------------------------------- TODO:
NAMT_OE23vsWT810_R <- fread("plotdata_NAMT-OE23vsWT810-R.txt", header = T, stringsAsFactors = F, sep = "\t")
KD_T <- read.xlsx("6mAKD RNA seq-T.xlsx")
T_NamtOEvs3D7_TPM <- fread("T_NamtOEvs3D7_dotplotdata.txt", header = T, stringsAsFactors = F, sep = "\t")
MMSL_RNA_2rep_T_230327_V2 <- read.xlsx("MMSL RNA 2rep T_230327 V2.xlsx")
NAMT_OE12vsWT78_S <- read.table("plotdata_NAMT-OE12vsWT78-S.txt", header = T, as.is = T, stringsAsFactors = F, sep = "\t")
R_6mAKDvs3D7_dotplotdata <- read.table("R_6mAKDvs3D7_dotplotdata.txt", header = T, sep = "\t", as.is = T)
plotdata_S_6mAKDvs3D7_TPM <- fread("plotdata_S_6mAKDvs3D7_TPM.txt", header = T, stringsAsFactors = F, sep = "\t")

interested_gene <- read.xlsx("2023 plot list for gss.xlsx")
DDRs <- c(interested_gene$HDR, interested_gene$MMEJ, interested_gene$BER, interested_gene$NER, interested_gene$MMR)
DDRs <- DDRs[!is.na(DDRs)]
Rhoptry <- interested_gene$Rhoptry
Rhoptry <- Rhoptry[!is.na(Rhoptry)]

### 1. MMSL
{
    colnames(MMSL_RNA_2rep_T_230327_V2)[1] <- "ID"
    textGene <- data.frame(ID = c("PF3D7_0107800", "PF3D7_0505500", "PF3D7_0803400", "PF3D7_1452000", "PF3D7_1350700"), genename = c("MRE11", "MSH6", "RAD54", "RON2", "PfNamt"))


    data <- MMSL_RNA_2rep_T_230327_V2
    data$colorType <- NA
    # data$colorType[data$pvalue < 0.05 & data$log2FoldChange > 1] <- "up"
    # data$colorType[data$pvalue < 0.05 & data$log2FoldChange < -1] <- "down"
    # data$colorType[data$ID %in% Rhoptry] <- "Rhoptry"
    data$colorType[data$ID %in% DDRs] <- "DDRs"
    data$colorType[data$ID == "PF3D7_1350700"] <- "PfNamt"

    data$text <- NA
    textGene1 <- textGene %>% filter(genename != "RON2")
    data$text[match(textGene1$ID, data$ID)] <- textGene1$genename

    ggplot(data %>% filter(significant != "LowTPM"), aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(color = "#b8ada7") +
        geom_point(data = data %>% filter(!is.na(colorType)), aes(color = colorType)) +
        geom_text_repel(aes(label = text), vjust = -1, size = 4) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        scale_x_continuous(breaks = c(-8, -4, 0, 4, 8)) +
        scale_color_manual(values = c(
            "DDRs" = "#214b9c", "PfNamt" = "#aa2423"
        )) +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            plot.title = element_text(
                colour = "black", face = "bold",
                size = 14, vjust = 1, hjust = 0.5
            )
        )
    ggsave("MMSL_RNA_2rep_T_230327_V2.pdf", height = 7, width = 7)
}


### 2. T KD
{
    colnames(KD_T)[1] <- "ID"
    textGene <- data.frame(ID = c("PF3D7_0107800", "PF3D7_0505500", "PF3D7_0803400", "PF3D7_1452000"), genename = c("MRE11", "MSH6", "RAD54", "RON2"))

    data <- KD_T
    data$colorType <- NA
    data$colorType[data$significant == "up"] <- "up"
    data$colorType[data$significant == "down"] <- "down"
    data$colorType[data$ID %in% Rhoptry] <- "Rhoptry"
    # data$colorType[data$ID %in% DDRs] <- "DDRs"
    # data$colorType[data$ID == "PF3D7_1350700"] <- "PfNamt"

    data$text <- NA
    textGene1 <- textGene %>% filter(genename == "RON2")
    data$text[match(textGene1$ID, data$ID)] <- textGene1$genename

    ggplot(data %>% filter(significant != "LowTPM"), aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(color = "#b8ada7") +
        geom_point(data = data %>% filter(!is.na(colorType)), aes(color = colorType)) +
        geom_text_repel(aes(label = text), vjust = -1, size = 4) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        scale_x_continuous(breaks = c(-8, -4, 0, 4, 8), limits = c(-8.5, 8.5)) +
        scale_color_manual(values = c(
            "up" = "#db8c6b", "down" = "#4494b0", "Rhoptry" = "purple", "DDRs" = "green"
        )) +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            plot.title = element_text(
                colour = "black", face = "bold",
                size = 14, vjust = 1, hjust = 0.5
            )
        )
    ggsave("6mAKD RNA seq-T.pdf", height = 7, width = 7)
}


### 3. T OE
{
    colnames(T_NamtOEvs3D7_TPM)[1] <- "ID"
    textGene <- data.frame(ID = c("PF3D7_0107800", "PF3D7_0505500", "PF3D7_0803400", "PF3D7_1452000"), genename = c("MRE11", "MSH6", "RAD54", "RON2"))

    data <- T_NamtOEvs3D7_TPM
    data$colorType <- NA
    data$colorType[data$significant == "up"] <- "up"
    data$colorType[data$significant == "down"] <- "down"
    data$colorType[data$ID %in% Rhoptry] <- "Rhoptry"
    # data$colorType[data$ID %in% DDRs] <- "DDRs"
    # data$colorType[data$ID == "PF3D7_1350700"] <- "PfNamt"

    data$text <- NA
    textGene1 <- textGene %>% filter(genename == "RON2")
    data$text[match(textGene1$ID, data$ID)] <- textGene1$genename

    ggplot(data %>% filter(significant != "LowTPM"), aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(color = "#b8ada7") +
        geom_point(data = data %>% filter(!is.na(colorType)), aes(color = colorType)) +
        geom_text_repel(aes(label = text), vjust = -1, size = 4) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        scale_x_continuous(breaks = c(-8, -4, 0, 4, 8), limits = c(-8.5, 8.5)) +
        scale_color_manual(values = c(
            "up" = "#db8c6b", "down" = "#4494b0", "Rhoptry" = "purple", "DDRs" = "green"
        )) +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            plot.title = element_text(
                colour = "black", face = "bold",
                size = 14, vjust = 1, hjust = 0.5
            )
        )
    ggsave("T_NamtOEvs3D7_TPM.pdf", height = 7, width = 7)
}


### 4. R OE
{
    textGene <- data.frame(ID = c("PF3D7_0107800", "PF3D7_0505500", "PF3D7_0803400", "PF3D7_1452000"), genename = c("MRE11", "MSH6", "RAD54", "RON2"))

    data <- NAMT_OE23vsWT810_R
    data$colorType <- NA
    data$colorType[data$significant == "up"] <- "up"
    data$colorType[data$significant == "down"] <- "down"
    # data$colorType[data$ID %in% Rhoptry] <- "Rhoptry"
    data$colorType[data$ID %in% DDRs] <- "DDRs"
    # data$colorType[data$ID == "PF3D7_1350700"] <- "PfNamt"

    data$text <- NA
    textGene1 <- textGene %>% filter(genename != "RON2")
    data$text[match(textGene1$ID, data$ID)] <- textGene1$genename

    ggplot(data %>% filter(significant != "LowTPM"), aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(color = "#b8ada7") +
        geom_point(data = data %>% filter(!is.na(colorType)), aes(color = colorType)) +
        geom_text_repel(aes(label = text), vjust = -1, size = 4) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        scale_x_continuous(breaks = c(-8, -4, 0, 4, 8), limits = c(-8.5, 8.5)) +
        scale_color_manual(values = c(
            "up" = "#db8c6b", "down" = "#4494b0", "Rhoptry" = "purple", "DDRs" = "green"
        )) +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            plot.title = element_text(
                colour = "black", face = "bold",
                size = 14, vjust = 1, hjust = 0.5
            )
        )
    ggsave("NAMT_OE23vsWT810_R.pdf", height = 7, width = 7)
}

### 5. S OE
{
    textGene <- data.frame(ID = c("PF3D7_0107800", "PF3D7_0505500", "PF3D7_0803400", "PF3D7_1452000"), genename = c("MRE11", "MSH6", "RAD54", "RON2"))

    data <- NAMT_OE12vsWT78_S
    # data$colorType <- NA
    # data$colorType[data$pvalue < 0.05 & data$log2FoldChange > 1] <- "up"
    # data$colorType[data$pvalue < 0.05 & data$log2FoldChange < -1] <- "down"
    # data$colorType[data$ID %in% Rhoptry] <- "Rhoptry"
    # data$colorType[data$ID %in% DDRs] <- "DDRs"
    # data$colorType[data$ID == "PF3D7_1350700"] <- "PfNamt"

    # data$text <- NA
    # textGene1 <- textGene %>% filter(genename != "RON2")
    # data$text[match(textGene$ID, data$ID)] <- textGene$genename

    ggplot(data %>% filter(significant != "LowTPM"), aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(color = "#b8ada7") +
        geom_point(data = data %>% filter(significant %in% c("up", "down")), aes(color = significant)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        scale_x_continuous(breaks = c(-8, -4, 0, 4, 8), limits = c(-8.5, 8.5)) +
        scale_color_manual(values = c(
            "up" = "#db8c6b", "down" = "#4494b0", "Rhoptry" = "purple", "DDRs" = "green"
        )) +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            plot.title = element_text(
                colour = "black", face = "bold",
                size = 14, vjust = 1, hjust = 0.5
            )
        )
    ggsave("NAMT_OE12vsWT78_S.pdf", height = 7, width = 7)
}


### 6. R KD
{
    colnames(R_6mAKDvs3D7_dotplotdata)[1] <- "ID"

    textGene <- data.frame(ID = c("PF3D7_0107800", "PF3D7_0505500", "PF3D7_0803400", "PF3D7_1452000"), genename = c("MRE11", "MSH6", "RAD54", "RON2"))

    data <- R_6mAKDvs3D7_dotplotdata
    # data$colorType <- NA
    # data$colorType[data$pvalue < 0.05 & data$log2FoldChange > 1] <- "up"
    # data$colorType[data$pvalue < 0.05 & data$log2FoldChange < -1] <- "down"
    # data$colorType[data$ID %in% Rhoptry] <- "Rhoptry"
    # data$colorType[data$ID %in% DDRs] <- "DDRs"
    # data$colorType[data$ID == "PF3D7_1350700"] <- "PfNamt"

    # data$text <- NA
    # textGene1 <- textGene %>% filter(genename != "RON2")
    # data$text[match(textGene$ID, data$ID)] <- textGene$genename

    ggplot(data %>% filter(significant != "LowTPM"), aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(color = "#b8ada7") +
        geom_point(data = data %>% filter(significant %in% c("up", "down")), aes(color = significant)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        scale_x_continuous(breaks = c(-8, -4, 0, 4, 8), labels = c(-8, -4, 0, 4, 8), limits = c(-8.5, 8.5)) +
        scale_color_manual(values = c(
            "up" = "#db8c6b", "down" = "#4494b0", "Rhoptry" = "purple", "DDRs" = "green"
        )) +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            plot.title = element_text(
                colour = "black", face = "bold",
                size = 14, vjust = 1, hjust = 0.5
            )
        )
    ggsave("R_6mAKDvs3D7_dotplotdata.pdf", height = 7, width = 7)
}


### 7. S KD
{
    colnames(plotdata_S_6mAKDvs3D7_TPM)[1] <- "ID"

    textGene <- data.frame(ID = c("PF3D7_0107800", "PF3D7_0505500", "PF3D7_0803400", "PF3D7_1452000"), genename = c("MRE11", "MSH6", "RAD54", "RON2"))

    data <- plotdata_S_6mAKDvs3D7_TPM
    # data$colorType <- NA
    # data$colorType[data$pvalue < 0.05 & data$log2FoldChange > 1] <- "up"
    # data$colorType[data$pvalue < 0.05 & data$log2FoldChange < -1] <- "down"
    # data$colorType[data$ID %in% Rhoptry] <- "Rhoptry"
    # data$colorType[data$ID %in% DDRs] <- "DDRs"
    # data$colorType[data$ID == "PF3D7_1350700"] <- "PfNamt"

    # data$text <- NA
    # textGene1 <- textGene %>% filter(genename != "RON2")
    # data$text[match(textGene$ID, data$ID)] <- textGene$genename

    ggplot(data %>% filter(significant != "LowTPM"), aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(color = "#b8ada7") +
        geom_point(data = data %>% filter(significant %in% c("up", "down")), aes(color = significant)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        scale_x_continuous(breaks = c(-8, -4, 0, 4, 8), labels = c(-8, -4, 0, 4, 8), limits = c(-8.5, 8.5)) +
        scale_color_manual(values = c(
            "up" = "#db8c6b", "down" = "#4494b0", "Rhoptry" = "purple", "DDRs" = "green"
        )) +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            plot.title = element_text(
                colour = "black", face = "bold",
                size = 14, vjust = 1, hjust = 0.5
            )
        )
    ggsave("plotdata_S_6mAKDvs3D7_TPM.pdf", height = 7, width = 7)
}
