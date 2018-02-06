#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/NishaRIPseq_volcanoplot.RData")

#####load packages#####
library(calibrate)
library(ggplot2)

###setwd#####
setwd("J:/III/Waters/Group Members/Mallory/NishaRNAseq/")

#####load data for RIP seq#####
RIPseq_raw <- read.csv("VolcanoPlots/RIPseqData_organized.csv")

RIPseq_plot_data <- RIPseq_raw[,c(8,3,6,7)]
  RIPseq_plot_data$neglog10pvalue <- -log10(RIPseq_plot_data$pvalue)
  RIPseq_plot_data$neglog10padj <- -log10(RIPseq_plot_data$padj)

  
#####plot for RIP seq#####
  with(RIPseq_plot_data, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-3,6)))

  # Add colored points: red if padj<0.05, orange of log2FC>3, green if both)
  with(subset(RIPseq_plot_data, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
  with(subset(RIPseq_plot_data, abs(log2FoldChange)>3), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
  with(subset(RIPseq_plot_data, padj<.05 & abs(log2FoldChange)>3), points(log2FoldChange, -log10(padj), pch=20, col="green"))
  
  # Label points with the textxy function from the calibrate plot

  with(subset(RIPseq_plot_data, padj<.05 & abs(log2FoldChange)>3), textxy(log2FoldChange, -log10(padj), labs=GeneID, cex=.8))
  
  
####plot with ggplot for RIP seq####
  ##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
  RIPseq_plot_data$threshold = as.factor(abs(RIPseq_plot_data$log2FoldChange) > 2 & RIPseq_plot_data$padj < 0.05)
  
  ##Construct the plot object
  g <- ggplot(data=RIPseq_plot_data, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
    geom_point(alpha=0.7, size=1.75) +
    xlim(c(-2.5, 6)) + ylim(c(0, 21)) +
    xlab("log2(Fold Change)") + 
    ylab("-log10(Adjusted p-value)") +
    theme_bw() +
    theme(axis.title.y = element_text(angle = 90, size=20)) +
    theme(axis.title.x = element_text(angle = 0, size=20)) +
    theme(axis.text.y = element_text(size=16)) +
    theme(axis.text.x = element_text(size=16)) +
    theme(legend.position =  "none") +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  g
  

#####load data for RNAseq#####
RNAseq_ap2g_raw <- read.csv("VolcanoPlots/RNAseq_ap2gKO_vs_820_data.csv")
  RNAseq_ap2g_plot_data <- RNAseq_ap2g_raw[,c(1,4:8)]
  RNAseq_ap2g_plot_data <- subset(RNAseq_ap2g_raw, RNAseq_ap2g_raw$status == "OK")
  RNAseq_ap2g_plot_data <- subset(RNAseq_ap2g_plot_data, RNAseq_ap2g_plot_data$log2_foldchange != Inf &
                                    RNAseq_ap2g_plot_data$log2_foldchange != -Inf)

                              
RNAseq_alba3_raw <- read.csv("VolcanoPlots/RNAseq_alba3KO_vs_820_data.csv")
  RNAseq_alba3_plot_data <- RNAseq_alba3_raw[,c(1,4:8)]
  RNAseq_alba3_plot_data <- subset(RNAseq_alba3_raw, RNAseq_alba3_raw$status == "OK")
  RNAseq_alba3_plot_data <- subset(RNAseq_alba3_plot_data, RNAseq_alba3_plot_data$log2_foldchange != Inf &
                                    RNAseq_alba3_plot_data$log2_foldchange != -Inf)


##Highlight genes that have an absolute fold change > 2 and a q-value < 0.05
RNAseq_ap2g_plot_data$threshold = as.factor(abs(RNAseq_ap2g_plot_data$log2_foldchange) > 2 & RNAseq_ap2g_plot_data$q_value < 0.05)
  
##Highlight genes that have an absolute fold change > 2 and a q-value < 0.05
RNAseq_alba3_plot_data$threshold = as.factor(abs(RNAseq_alba3_plot_data$log2_foldchange) > 2 & RNAseq_alba3_plot_data$q_value < 0.05)
  
  
#figure out the genes that are differentially regulated in each
RNAseq_ap2g_genes <- subset(RNAseq_ap2g_plot_data, RNAseq_ap2g_plot_data$threshold == TRUE)
  RNAseq_ap2g_genes <- as.character(RNAseq_ap2g_genes$gene_id)
RNAseq_alba3_genes <- subset(RNAseq_alba3_plot_data, RNAseq_alba3_plot_data$threshold == TRUE)
  RNAseq_alba3_genes <- as.character(RNAseq_alba3_genes$gene_id)
  
RNAseq_shared_genes <- as.factor(intersect(RNAseq_ap2g_genes, RNAseq_alba3_genes))


#add colour factor for plot
RNAseq_ap2g_plot_data$colour <- NA
  for (i in 1:length(RNAseq_ap2g_plot_data$gene_id)) {
    if (RNAseq_ap2g_plot_data$gene_id[i] %in% RNAseq_shared_genes) {
      RNAseq_ap2g_plot_data$colour[i] <- "shared"
    } else {
      RNAseq_ap2g_plot_data$colour[i] <- "not shared"
    }
    if (RNAseq_ap2g_plot_data$colour[i] == "not shared") {
      RNAseq_ap2g_plot_data$colour[i] <- RNAseq_ap2g_plot_data$threshold[i]
    }
  }
RNAseq_ap2g_plot_data$colour <- as.factor(RNAseq_ap2g_plot_data$colour)

#add colour factor for plot
RNAseq_alba3_plot_data$colour <- NA
for (i in 1:length(RNAseq_alba3_plot_data$gene_id)) {
  if (RNAseq_alba3_plot_data$gene_id[i] %in% RNAseq_shared_genes) {
    RNAseq_alba3_plot_data$colour[i] <- "shared"
  } else {
    RNAseq_alba3_plot_data$colour[i] <- "not shared"
  }
  if (RNAseq_alba3_plot_data$colour[i] == "not shared") {
    RNAseq_alba3_plot_data$colour[i] <- RNAseq_alba3_plot_data$threshold[i]
  }
}
RNAseq_alba3_plot_data$colour <- as.factor(RNAseq_alba3_plot_data$colour)
  

####plot with ggplot for RNA seq####


###make function for plot###
VolcanoPlot <- function(data,gene){
    ggplot(data=data, aes(x=log2_foldchange, y=-log10(q_value), colour=colour)) +
    # geom_point(alpha=0.7, aes(size = colour)) +
    geom_point(data = subset(data, data$colour == "1"), alpha=0.9, size = 1, shape = 19) + #grey
    geom_point(data = subset(data, data$colour == "2"), alpha=0.9, size = 1.5, shape = 19) + #red - signif
    geom_point(data = subset(data, data$colour == "shared"), alpha=0.9, size = 1.5, shape = 19) + #blue shared signif
    # scale_size_manual(values = c(1,2,2)) +
    xlab(expression(paste("log"[2],"(fold change)"))) +
    ylab(expression(paste("-log"[10],"(q-value)"))) +
    scale_colour_manual(values = c("grey50","red3","blue3")) +
    scale_x_continuous(breaks = c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8), 
                       labels = c(-12,"",-10,"",-8,"",-6,"",-4,"",-2,"",0,"",2,"",4,"",6,"",8),
                       minor_breaks = c(-11,-9,-7,-5,-3,-1,1,3,5,7),
                       limits = c(-12,8)) +
    scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5,3,3.5,4),
                       labels = c(0,"",1,"",2,"",3,"",4),
                       limits = c(0,4)) +
    ggtitle(paste(gene)) +
    theme_bw() +
    theme(axis.title.y = element_text(angle = 90, size=20),
          axis.title.x = element_text(angle = 0, size=20),
          axis.text.y = element_text(size=16, colour = "black", margin = margin (0,2,0,0,"mm")),
          axis.text.x = element_text(size=16, colour = "black", margin = margin(2,0,2,0,"mm")),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.ticks.length = unit(2, "mm"),
          plot.title = element_text(size=18, hjust = 0.5),
          legend.position =  "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(5,5,5,5), "mm"))
}


RNAseq_ap2g_plot <- VolcanoPlot(data = RNAseq_ap2g_plot_data, gene = "AP2-G")
RNAseq_alba3_plot <- VolcanoPlot(data = RNAseq_alba3_plot_data, gene = "Alba3")

#####RNAseq multiplot#####
RNAseq_plot <- multiplot(RNAseq_ap2g_plot,RNAseq_alba3_plot)

png(filename = "RNAseq_volcano_plot.png", 
     pointsize = 12,
     width = 15,
     height = 20,
     units = "cm",
     res = 600,
     # quality = 600,
     bg = "white")
multiplot(RNAseq_ap2g_plot,RNAseq_alba3_plot)
dev.off()

#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/NishaRIPseq_volcanoplot.RData")
