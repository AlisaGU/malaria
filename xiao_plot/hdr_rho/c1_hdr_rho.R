#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggalt)
library(dplyr)
# 2. functions ------------------------------------------------------------ TODO:
read_data <- function(gene = NULL) {
    wt <- read.csv(paste0("3D7-", gene, ".csv"), header = T, row.names = 1, as.is = T)
    wt <- apply(wt, 2, mean)

    kd <- read.csv(paste0("6mAKD-", gene, ".csv"), header = T, row.names = 1, as.is = T)
    kd <- apply(kd, 2, mean)

    data_for_plot <- rbind(
        data.frame(value = wt, label = "WT", position = 1:length(wt)),
        data.frame(value = kd, label = "KD", position = 1:length(kd))
    )

    data_for_plot$label <- factor(data_for_plot$label, levels = c("WT", "KD"))
    return(data_for_plot)
}

plot_as_wen <- function(data_for_plot = NULL, array_color_value = NULL,span=NULL) {
  wt<-data_for_plot %>% filter(label=="WT")
  kd<-data_for_plot %>% filter(label=="KD")
  spline_int_wt <- as.data.frame(spline(wt$position, wt$value))
  spline_int_kd <- as.data.frame(spline(kd$position, kd$value))
    p <- ggplot(data = data_for_plot, aes(x = position, y = value,group=label)) +
        # geom_point(aes(color = label))+
      geom_smooth(aes(color = label),method = "loess",span=span,se=FALSE,size=3)+
      # geom_line(data = spline_int_wt, aes(x = x, y = y),linewidth=1,color="yellow")+
      # geom_line(data = spline_int_kd, aes(x = x, y = y),linewidth=1,color="blue")+
        # geom_line(aes(color = label), linewidth = 1) +
        # geom_area(aes(fill = label), alpha = 0.1, position = "identity") +
        scale_fill_manual(values = array_color_value, guide = "none") +
        scale_x_continuous(breaks = c(0, 100, 250, 400, 500), labels = c("2kb\nupstream", "Start", "gene body", "Stop", "2kb\ndownstream")) +
        scale_color_manual(values = array_color_value) +
        labs(y = "6mA signal density") +
        theme_classic() +
        theme(
          axis.ticks.length = unit(.25, "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            # strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 44, color = "black", face = "bold"),
            # panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(
                size = 36, color = "black"
                # angle = 45, hjust = 0.5, vjust = 0.5
            ),
            axis.text.y = element_text(size = 36, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.position = c(0.85, 0.9),
            legend.text = element_text(size = 36, color = "black"),
            legend.key.width = unit(4, "cm"),
            plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
# setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/rho_hdr/wyh_ana")
setwd("I:\\projects\\malaria\\dataAndResult\\xiao_plot\\rho_hdr\\wyh_ana")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
hdr <- read_data(gene = "HDR")
rho <- read_data(gene = "Rho")
plot_as_wen(hdr, array_color_value = c("KD" = "red", "WT" = "black"),span=0.35)
ggsave("hdr.pdf",width=10,height=7)
plot_as_wen(rho, array_color_value = c("KD" = "red", "WT" = "black"),span=0.45)+
  scale_y_continuous(labels=c("1.0","1.5","2.0","2.5"),
                     breaks=c(1,1.5,2,2.5),limits=c(NA,2.6))
ggsave("rho.pdf",width=10,height=7)

