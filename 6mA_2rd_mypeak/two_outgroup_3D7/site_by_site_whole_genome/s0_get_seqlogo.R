#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(seqinr)
library(ggseqlogo)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
read.data <- function() {
    data <- lapply(paste("flank_", 1:10, sep = ""), function(x) {
        result <- lapply(chroms, function(chrom) {
            result <- read.fasta(paste0(x, "/", chrom, ".peak.fa"), as.string = TRUE, seqtype = "DNA", set.attributes = F)
            result <- unlist(result)
            names(result) <- NULL
            result <- toupper(result)
            return(result)
        })
        names(result) <- chroms
        return(result)
    })
    names(data) <- paste("flank_", 1:10, sep = "")
    return(data)
}


plot_graph <- function(seq_data = NULL) {
    result <- lapply(paste("flank_", 1:10, sep = ""), function(x) {
        result <- lapply(chroms, function(chrom) {
            p <- ggplot() +
                geom_logo(seq_data[[x]][[chrom]], seq_type = "dna") +
                labs(title = paste0(x, ": ", chrom)) +
                scale_y_continuous(limits = c(0, 0.4)) +
                geom_hline(yintercept=0.05)+
                geom_hline(yintercept=0.1)+
                geom_hline(yintercept=0.2)+
                geom_hline(yintercept=0.3)+
                theme(
                    plot.title = element_text(size=20,face="bold"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    strip.background = element_blank(),
                    panel.border = element_blank(),
                    # strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
                    strip.text.x = element_text(size = 44, color = "black", face = "bold"),
                    # panel.border = element_blank(),
                    axis.line = element_line(colour = "black"),
                    axis.text.x = element_text(
                        size = 36, color = "black",
                        # angle = 45, hjust = 0.5, vjust = 0.5
                    ),
                    axis.text.y = element_text(size = 36, color = "black"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_text(size = 40, color = "black"),
                    legend.title = element_blank(),
                    legend.text = element_text(size = 30, color = "black"),
                    legend.position = "bottom",
                    plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
                    panel.spacing = unit(3, "lines")
                )
            return(p)
        })
        names(result) <- chroms
        return(result)
    })
    names(result) <- paste("flank_", 1:10, sep = "")
    return(result)
}

patchGraph <- function(graph_data = NULL) {
    final_graph <- lapply(chroms, function(chrom) {
        as.character(graph_data[["flank_1"]][[chrom]])
        cmd <- paste(paste(paste("graph_data[[\"", paste("flank_", 1:10, sep = ""), sep = ""), "\"]][[chrom]]", sep = ""), collapse = "|")
        p <- eval(str2expression(cmd))
        return(p)
    })
    names(final_graph) <- chroms
    cmd <- paste(paste(paste("final_graph[[\"", chroms, sep = ""),"\"]]",sep=""), collapse = "/")
    p <- eval(str2expression(cmd))
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
chroms <- c("Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_07_v3", "Pf3D7_08_v3", "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3", "Pf3D7_13_v3", "Pf3D7_14_v3")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern_GAWGAW_asWhole/GAWGAW/flank_nomatter_strand.remove_close_overlap")

seq_data <- read.data()
graph_data <- plot_graph(seq_data = seq_data)
patchGraph(graph_data = graph_data)
ggsave(filename="seqlogo.pdf",width=50,height=70,limitsize = FALSE)
