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
NAMT_OEvsWT_R <- read.xlsx("NAMT-OE vs WT-R.xlsx")
colnames(NAMT_OEvsWT_R)[1] <- "ID"

interested_gene <- read.xlsx("2023 plot list for gss.xlsx")
DDRs <- c(interested_gene$HDR, interested_gene$MMEJ, interested_gene$BER, interested_gene$NER, interested_gene$MMR)
DDRs <- DDRs[!is.na(DDRs)]
Rhoptry <- interested_gene$Rhoptry
Rhoptry <- Rhoptry[!is.na(Rhoptry)]



### 4. R OE
{
    textGene <- data.frame(ID = c("PF3D7_0107800", "PF3D7_0505500", "PF3D7_0803400", "PF3D7_1452000"), genename = c("MRE11", "MSH6", "RAD54", "RON2"))

    data <- NAMT_OEvsWT_R
    data$colorType <- NA
    data$colorType[data$significant == "up"] <- "up"
    data$colorType[data$significant == "down"] <- "down"
    # data$colorType[data$ID %in% Rhoptry] <- "Rhoptry"
    data$colorType[data$ID %in% DDRs] <- "DDRs"
    # data$colorType[data$ID == "PF3D7_1350700"] <- "PfNamt"

    data$text <- NA
    textGene1 <- textGene %>% filter(genename == "RAD54")
    data$text[match(textGene1$ID, data$ID)] <- textGene1$genename

    ggplot(data %>% filter(significant != "LowTPM"), aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(color = "#b8ada7") +
        geom_point(data = data %>% filter(colorType == "up" | colorType == "down"), aes(color = colorType)) +
        geom_point(data = data %>% filter(colorType == "DDRs"), aes(color = colorType)) +
        geom_text_repel(aes(label = text), vjust = -1, size = 4) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        scale_x_continuous(breaks = c(-8, -4, 0, 4, 8), limits = c(-8.5, 8.5)) +
        scale_color_manual(values = c(
            "up" = "#db8c6b", "down" = "#4494b0", "DDRs" = "green"
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
    ggsave("NAMT_OEvsWT_R.pdf", height = 7, width = 7)
}
