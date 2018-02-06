#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/Brett_transcriptome_GoIs.RData")
library("ggplot2")

# Install
# install.packages("wesanderson")
# # Load
# library(wesanderson)


####set wd####
setwd("J:/III/Waters/Group Members/Mallory/BrettTranscriptomics/")


####load data####

allgene_plotdata <- read.csv("ExtractedDataForR.csv")

  split <- split(allgene_plotdata, allgene_plotdata$OrigLine)
  allgene3D7_plotdata <- split$`3D7`
  allgene7G8_plotdata <- split$`7G8`
  allgeneD10_plotdata <- split$D10
  allgeneHB3A_plotdata <- split$HB3A
  
  levels(allgene_plotdata$Line)
  unique(allgene_plotdata$Line)
  
  allgene_plotdata$Line <- factor(allgene_plotdata$Line, levels = unique(allgene_plotdata$Line))
  
  # 
  # spline_int <- data.frame()
  #   for (i in seq_along(allgene_plotdata$Line)){
  #     data <- subset(allgene_plotdata, allgene_plotdata$Line == i)
  #     cbind(spline_int, spline(data$TimePoint, data$NormGeneLevel))
  #   }
  #   
  #   
  # spline_int <- as.data.frame(spline(allgene3D7_plotdata$TimePoint, allgene3D7_plotdata$NormGeneLevel))
  
#####plots#####
  #use facet_grid() to make a grid of plots
  # + facet_grid(OrigLine ~ GeneName2)
  # The palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  
AllGene_3D7_plot <- ggplot(allgene3D7_plotdata, aes(x=TimePoint, y=NormGeneLevel)) +
    geom_point(aes(colour = factor(Line), fill =  factor(Line), shape = factor(Line)), size = 0.1) +
    geom_line(aes(colour = factor(Line)), size = 1.1) +
    # geom_line(data = spline_int, aes(x = x, y = y)) +
    # stat_smooth(aes(y=NormGeneLevel, x=TimePoint, colour = factor(SubLineNum)), formula = y ~ s(x, k = 5), method = "gam", se = FALSE) +
    scale_fill_manual(values = cbbPalette) +
    scale_shape_manual(values = c(15,19,17,23,25,8)) +
    scale_colour_manual(values = cbbPalette) +
    # xlab("Time (hpi)") + 
    # ylab("Normalized\nGene Level") + 
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.y = element_text(size=16)) +
    theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
    theme(legend.title = element_blank()) +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1,1,1,1), "cm")) +
    facet_grid(. ~ GeneName3)
  
AllGene_7G8_plot <- ggplot(allgene7G8_plotdata, aes(x=TimePoint, y=NormGeneLevel)) +
  geom_point(aes(colour = factor(Line), fill =  factor(Line), shape = factor(Line)), size = 0.1) +
  geom_line(aes(colour = factor(Line)), size = 1.1) +
  # geom_line(data = spline_int, aes(x = x, y = y)) +
  # stat_smooth(aes(y=NormGeneLevel, x=TimePoint, colour = factor(SubLineNum)), formula = y ~ s(x, k = 5), method = "gam", se = FALSE) +
  scale_fill_manual(values = cbbPalette) +
  scale_shape_manual(values = c(15,19,17,23,25,8)) +
  scale_colour_manual(values = cbbPalette) +
  # xlab("Time (hpi)") + 
  ylab("Normalized\nGene Level") + 
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  facet_grid(. ~ GeneName3)

AllGene_D10_plot <- ggplot(allgeneD10_plotdata, aes(x=TimePoint, y=NormGeneLevel)) +
  geom_point(aes(colour = factor(Line), fill =  factor(Line), shape = factor(Line)), size = 0.1) +
  geom_line(aes(colour = factor(Line)), size = 1.1) +
  # geom_line(data = spline_int, aes(x = x, y = y)) +
  # stat_smooth(aes(y=NormGeneLevel, x=TimePoint, colour = factor(SubLineNum)), formula = y ~ s(x, k = 5), method = "gam", se = FALSE) +
  scale_fill_manual(values = cbbPalette) +
  scale_shape_manual(values = c(15,19,17,23,25,8)) +
  scale_colour_manual(values = cbbPalette) +
  # xlab("Time (hpi)") + 
  # ylab("Normalized\nGene Level") + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  facet_grid(. ~ GeneName3)


AllGene_HB3A_plot <- ggplot(allgeneHB3A_plotdata, aes(x=TimePoint, y=NormGeneLevel)) +
  geom_point(aes(colour = factor(Line), fill =  factor(Line), shape = factor(Line)), size = 0.1) +
  geom_line(aes(colour = factor(Line)), size = 1.1) +
  # geom_line(data = spline_int, aes(x = x, y = y)) +
  # stat_smooth(aes(y=NormGeneLevel, x=TimePoint, colour = factor(SubLineNum)), formula = y ~ s(x, k = 5), method = "gam", se = FALSE) +
  scale_fill_manual(values = cbbPalette) +
  scale_shape_manual(values = c(15,19,17,23,25,8,43)) +
  scale_colour_manual(values = cbbPalette) +
  xlab("Time (hpi)") + 
  # ylab("Normalized\nGene Level") + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  facet_grid(. ~ GeneName3)


#21 colour palette
# tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

###for the colour palette, try using the same colours = 5,5,5,6 colours
# The palette with black:

cbbPalette3 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", 
                 "#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", 
                 "#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", 
                 "#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00")


AllGene_allplot <- ggplot(allgene_plotdata, aes(x=TimePoint, y=NormGeneLevel)) +
  geom_point(aes(colour = factor(Line), fill =  factor(Line)), size = 0.1) +
  geom_line(aes(colour = factor(Line)), size = 0.8) +
  scale_fill_manual(name = "Line", values = cbbPalette3) +
  scale_colour_manual(name = "Line", values = cbbPalette3) +
  # xlab("Time (hpi)") + 
  # ylab("Normalized Gene Level") +
  # ylab(expression(log["2"](Cy5/Cy3)))+
  labs(x = "Time (hpi)", y = expression(log["2"](Cy5/Cy3))) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title.align=0.5) +
  theme(legend.title = element_text(size=16)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  facet_grid(OrigLine ~ GeneName3)



# allplots <- multiplot(AllGene_3D7_plot, AllGene_7G8_plot, AllGene_D10_plot,AllGene_HB3A_plot, cols = 1)

ggsave("AllGene_allplot.pdf", 
       plot = AllGene_allplot, 
       units = "cm", 
       width = 25, 
       height = 25)

####extra####
##trying to figure out how to get a smooth line between points
#looks like I could use spline (https://stackoverflow.com/questions/35205795/plotting-smooth-line-through-all-data-points-maybe-polynomial-interpolation)
#but i think the problem is I'd have to calculate it separately for each line

# 
# ggplot(d) + 
#   geom_point(aes(x = hour, y = impressions, colour = cvr), size = 3) +
#   geom_line(data = spline_int, aes(x = x, y = y))


#might need to make separate DFs for all, then separate plots, then put them in a grid
  
#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/Brett_transcriptome_GoIs.RData")
