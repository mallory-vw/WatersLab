#load environment
load("J:/III/Waters/Group Members/Mallory/R/ALBA_project_heatmaps.RData")

#setwd
setwd("J:/III/Waters/Group Members/Mallory/ProteinsProject/Results_Evangelia")

######load ratio data######
library("reshape2")
library("plyr")

ratio_ALBA_BG_allres_1 <- as.data.frame(read.table("ratio4protsvsproteome1.txt", 
                                                   header=T))
ratio_ALBA_BG_allres_1 <- melt(ratio_ALBA_BG_allres_1, id.var = "Name")

AAs <- ratio_ALBA_BG_allres_1$Name[1:20]
AAs <- ordered(AAs, levels=c("E","D","R","K","H","W","F","Y","S","T","G","A","P","Q","N","C","M","V","L","I"))
levels(AAs)


ratio_ALBA_BG_allres_2 <- as.data.frame(read.table("ratio4protsvsproteome2.txt", 
                                                   header=T))
ratio_ALBA_BG_allres_2 <- melt(ratio_ALBA_BG_allres_2, id.var = "Name")

ratio_ALBA_BG_allres_3 <- as.data.frame(read.table("ratio4protsvsproteome3.txt", 
                                                   header=T))
ratio_ALBA_BG_allres_3 <- melt(ratio_ALBA_BG_allres_3, id.var = "Name")

ratio_ALBA_BG_allres_4 <- as.data.frame(read.table("ratio4protsvsproteome4.txt", 
                                                   header=T))
ratio_ALBA_BG_allres_4 <- melt(ratio_ALBA_BG_allres_4, id.var = "Name")

ratio_ALBA_BG_allres_5 <- as.data.frame(read.table("ratio4protsvsproteome5.txt", 
                                                   header=T))
ratio_ALBA_BG_allres_5 <- melt(ratio_ALBA_BG_allres_5, id.var = "Name")



#####heat maps#####
library("ggplot2")


ratio_ALBA_BG_allres_1_plot <- ggplot(ratio_ALBA_BG_allres_1, aes(Name, variable, fill = value)) +
  geom_tile() +
  ylim(levels(AAs)) +
  xlim(levels(AAs)) +
  xlab("Residue") +
  ylab("Residue") +
  scale_fill_gradient(low = "black", high = "red") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
ggsave(ratio_ALBA_BG_allres_1_plot, file = "ratio_ALBA_BG_allres_1_plot.pdf", height = 8, width = 10)

ratio_ALBA_BG_allres_2_plot <- ggplot(ratio_ALBA_BG_allres_2, aes(Name, variable, fill = value)) +
  geom_tile() +
  ylim(levels(AAs)) +
  xlim(levels(AAs)) +
  xlab("Residue") +
  ylab("Residue") +
  scale_fill_gradient(low = "black", high = "red") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
ggsave(ratio_ALBA_BG_allres_2_plot, file = "ratio_ALBA_BG_allres_2_plot.pdf", height = 8, width = 10)

ratio_ALBA_BG_allres_3_plot <- ggplot(ratio_ALBA_BG_allres_3, aes(Name, variable, fill = value)) +
  geom_tile() +
  ylim(levels(AAs)) +
  xlim(levels(AAs)) +
  xlab("Residue") +
  ylab("Residue") +
  scale_fill_gradient(low = "black", high = "red") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
ggsave(ratio_ALBA_BG_allres_3_plot, file = "ratio_ALBA_BG_allres_3_plot.pdf", height = 8, width = 10)

ratio_ALBA_BG_allres_4_plot <- ggplot(ratio_ALBA_BG_allres_4, aes(Name, variable, fill = value)) +
  geom_tile() +
  ylim(levels(AAs)) +
  xlim(levels(AAs)) +
  xlab("Residue") +
  ylab("Residue") +
  scale_fill_gradient(low = "black", high = "red") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
ggsave(ratio_ALBA_BG_allres_4_plot, file = "ratio_ALBA_BG_allres_4_plot.pdf", height = 8, width = 10)

ratio_ALBA_BG_allres_5_plot <- ggplot(ratio_ALBA_BG_allres_5, aes(Name, variable, fill = value)) +
  geom_tile() +
  ylim(levels(AAs)) +
  xlim(levels(AAs)) +
  xlab("Residue") +
  ylab("Residue") +
  scale_fill_gradient(low = "black", high = "red") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
ggsave(ratio_ALBA_BG_allres_5_plot, file = "ratio_ALBA_BG_allres_5_plot.pdf", height = 8, width = 10)



#####final save#####
save.image("J:/III/Waters/Group Members/Mallory/R/ALBA_project_heatmaps.RData")
