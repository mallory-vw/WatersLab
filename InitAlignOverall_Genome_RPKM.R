# load environment
load("J:/III/Waters/Group Members/Mallory/R/InitAlignOverall_Genome_RPKM.RData")

#load packages
library("ggplot2")
library("reshape2")

#set wd
setwd("J:/III/Waters/Group Members/Mallory/GenomeProjects/GenomesVariants/CNV")

#load indivs names
indivs <- read.table("J:/III/Waters/Group Members/Mallory/GenomeProjects/scripts/files/InitAlignOverall_indivs.txt")
indivs <- indivs[,1]

#load chr names
chr <- read.table("J:/III/Waters/Group Members/Mallory/GenomeProjects/scripts/files/InitAlignOverall_chromosomes.txt")
chr <- chr[,1]

#load chr size
chr_size <- read.table("J:/III/Waters/Group Members/Mallory/GenomeProjects/GenomesVariants/PbergheiGenome/genome.txt")
colnames(chr_size) <- c("chr","size_in_bases")
chr_size$"size_in_KB" <- chr_size$size_in_bases/1000


####load RPKM data (created in excel from bedtools multicov)####
  # created a bed file with 1000KB segments (the very last segment was not 1000KB)
  # then ran multicov using that bed file --> gave me number of alignments(reads?) overlapping each 1000kB segment
  # divided that number by the number of reads(alignments?) mapped per line
  # this gives RPKM (reads per kilobase per million mapped) --> allows for comparison across different lines while controlling for the number of reads mapped per line
RPKM <- read.csv("InitAlignOverall_multicov_RPKM.csv")
head(levels(RPKM$ContigName))
RPKM$ContigName <- factor(RPKM$ContigName, 
                          levels = RPKM$ContigName[order(RPKM$Chr, RPKM$ContigStart)])

RPKM_plot_data <- RPKM[,1:15]
RPKM_plot_data <- melt(RPKM_plot_data, id.vars=1:4)
colnames(RPKM_plot_data) <- c("Chr","ContigStart","ContigEnd","ContigName","Line","RPKM")
RPKM_plot_data$Line <- gsub('X','',RPKM_plot_data$Line)
RPKM_plot_data$LineName <- RPKM_plot_data$Line
  RPKM_plot_data$LineName <- gsub('_trimmed_aligned_sorted','',RPKM_plot_data$LineName)
  RPKM_plot_data$LineName <- as.factor(RPKM_plot_data$LineName)
head(levels(RPKM_plot_data$ContigName))


#add in range and save
#make range data
RPKM_range <- RPKM
for (i in 1:nrow(RPKM_range)) {
  RPKM_range$min[i] <- min(RPKM_range[i,5:15])
  RPKM_range$max[i] <- max(RPKM_range[i,5:15])
  RPKM_range$range[i] <- (RPKM_range$max[i])-(RPKM_range$min[i])
}

write.csv(RPKM_range, "InitAlignOverall_multicov_RPKM_range.csv")
####plots!####
#make a manhattan plot for each contig name? each chr first for each indiv?
  # not a traditional manhattan plot, just want the geom_point basically for RPKM
# hwo about a heat map?  y = count, x = contig --> how would I compare between diff indivs here?

##get plot data ready - make a list of DFs, 1 DF for each chr
RPKM_plot_data_Chr <- split(RPKM_plot_data, RPKM_plot_data$Chr)

#Create a custom colour scale based on the population factors - for all plots
myColours2 <- as.character(c("dodgerblue4","cornflowerblue", "darkgreen","olivedrab3","red4","lightcoral","orangered2","orange1","darkmagenta","mediumorchid1","grey15"))
names(myColours2) <- levels(RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_1$LineName)
colScale2 <- scale_color_manual(name = "Line", values = myColours2)

####plots coded individually#####
#PbANKA_00_v3_archived_contig_1_plot
PbANKA_00_v3_archived_contig_1_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_1, 
                                              aes(x = RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_1$ContigName, 
                                                  y = RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_1$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_1$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_00_v3_archived_contig_2_plot
PbANKA_00_v3_archived_contig_2_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_2, 
                                              aes(x = RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_2$ContigName, 
                                                  y = RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_2$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_2$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_00_v3_archived_contig_3_plot
PbANKA_00_v3_archived_contig_3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_3, 
                                              aes(x = RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_3$ContigName, 
                                                  y = RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_00_v3_archived_contig_4_plot
PbANKA_00_v3_archived_contig_4_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_4, 
                                              aes(x = RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_4$ContigName, 
                                                  y = RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_4$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_4$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_00_v3_archived_contig_5_plot
PbANKA_00_v3_archived_contig_5_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_5, 
                                              aes(x = RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_5$ContigName, 
                                                  y = RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_5$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_00_v3_archived_contig_5$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_01_v3_plot
PbANKA_01_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_01_v3, 
                                              aes(x = RPKM_plot_data_Chr$PbANKA_01_v3$ContigName, 
                                                  y = RPKM_plot_data_Chr$PbANKA_01_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_01_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_02_v3_plot
PbANKA_02_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_02_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_02_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_02_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_02_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_03_v3_plot
PbANKA_03_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_03_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_03_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_03_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_03_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_04_v3_plot
PbANKA_04_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_04_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_04_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_04_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_04_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_05_v3_plot
PbANKA_05_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_05_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_05_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_05_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_05_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_06_v3_plot
PbANKA_06_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_06_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_06_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_06_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_06_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_07_v3_plot
PbANKA_07_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_07_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_07_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_07_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_07_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_08_v3_plot
PbANKA_08_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_08_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_08_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_08_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_08_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_09_v3_plot
PbANKA_09_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_09_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_09_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_09_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_09_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_10_v3_plot
PbANKA_10_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_10_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_10_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_10_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_10_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_11_v3_plot
PbANKA_11_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_11_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_11_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_11_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_11_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_12_v3_plot
PbANKA_12_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_12_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_12_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_12_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_12_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_13_v3_plot
PbANKA_13_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_13_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_13_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_13_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_13_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_14_v3_plot
PbANKA_14_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_14_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_14_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_14_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_14_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_API_v3_plot
PbANKA_API_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_API_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_API_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_API_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_API_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_MIT_v3_plot
PbANKA_MIT_v3_plot <- ggplot(RPKM_plot_data_Chr$PbANKA_MIT_v3, 
                            aes(x = RPKM_plot_data_Chr$PbANKA_MIT_v3$ContigName, 
                                y = RPKM_plot_data_Chr$PbANKA_MIT_v3$RPKM)) +
  geom_point(aes(colour=RPKM_plot_data_Chr$PbANKA_MIT_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#####save plots#####
ggsave(PbANKA_00_v3_archived_contig_1_plot, file="PbANKA_00_v3_archived_contig_1_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_00_v3_archived_contig_2_plot, file="PbANKA_00_v3_archived_contig_2_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_00_v3_archived_contig_3_plot, file="PbANKA_00_v3_archived_contig_3_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_00_v3_archived_contig_4_plot, file="PbANKA_00_v3_archived_contig_4_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_00_v3_archived_contig_5_plot, file="PbANKA_00_v3_archived_contig_5_plot.pdf",
       height = 10, width = 14, unit = "in")




#####RPKM range > 100#####
#issues - some chr are waaaaay too big for 1 plot
#each plot can handle probably about 50-60KB and still be leigible
# but some chrs are way to big to make individual 60KB plots
# figure out a way of highlighting which chr 1KB segments to focus on
# maybe calculate the range (across all lines) for each segment
# done, now remake plots for smaller segments per chr?

#####load/make range subset data#####
#make range data
for (i in 1:nrow(RPKM)) {
  RPKM$min[i] <- min(RPKM[i,5:15])
  RPKM$max[i] <- max(RPKM[i,5:15])
  RPKM$range[i] <- (RPKM$max[i])-(RPKM$min[i])
}

#subset dataframe to range > 100
RPKM_CNV_subset <- subset(RPKM, RPKM$range > 100)

write.csv(RPKM_CNV_subset, "InitAlignOverall_multicov_RPKM_greater_100.csv")

######make plot range subset data dataframe#####
RPKM_CNV_subset_plot_data <- RPKM_CNV_subset[,1:15]
RPKM_CNV_subset_plot_data <- melt(RPKM_CNV_subset_plot_data, id.vars=1:4)
colnames(RPKM_CNV_subset_plot_data) <- c("Chr","ContigStart","ContigEnd","ContigName","Line","RPKM")
RPKM_CNV_subset_plot_data$Line <- gsub('X','',RPKM_CNV_subset_plot_data$Line)
RPKM_CNV_subset_plot_data$LineName <- RPKM_CNV_subset_plot_data$Line
RPKM_CNV_subset_plot_data$LineName <- gsub('_trimmed_aligned_sorted','',RPKM_CNV_subset_plot_data$LineName)
RPKM_CNV_subset_plot_data$LineName <- as.factor(RPKM_CNV_subset_plot_data$LineName)
head(levels(RPKM_CNV_subset_plot_data$ContigName))

####range subset plots!####
##get plot data ready - make a list of DFs, 1 DF for each chr
RPKM_CNV_subset_plot_data_Chr <- split(RPKM_CNV_subset_plot_data, RPKM_CNV_subset_plot_data$Chr)

####range subset plots coded individually#####
#PbANKA_00_v3_archived_contig_1_range_subset_plot
PbANKA_00_v3_archived_contig_1_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_1, 
                                              aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_1$ContigName, 
                                                  y = RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_1$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_1$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  colScale2

#PbANKA_00_v3_archived_contig_2_range_subset_plot
PbANKA_00_v3_archived_contig_2_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_2, 
                                              aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_2$ContigName, 
                                                  y = RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_2$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_2$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_00_v3_archived_contig_3_range_subset_plot
PbANKA_00_v3_archived_contig_3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_3, 
                                              aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_3$ContigName, 
                                                  y = RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_00_v3_archived_contig_4_range_subset_plot
PbANKA_00_v3_archived_contig_4_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_4, 
                                              aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_4$ContigName, 
                                                  y = RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_4$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_4$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_00_v3_archived_contig_5_range_subset_plot
PbANKA_00_v3_archived_contig_5_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_5, 
                                              aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_5$ContigName, 
                                                  y = RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_5$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_5$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_01_v3_range_subset_plot
PbANKA_01_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_01_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_01_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_01_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_01_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_02_v3_range_subset_plot
PbANKA_02_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_02_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_02_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_02_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_02_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_03_v3_range_subset_plot
PbANKA_03_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_03_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_03_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_03_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_03_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_04_v3_range_subset_plot
PbANKA_04_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_04_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_04_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_04_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_04_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_05_v3_range_subset_plot
PbANKA_05_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_05_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_05_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_05_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_05_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_06_v3_range_subset_plot
PbANKA_06_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_06_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_06_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_06_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_06_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_07_v3_range_subset_plot
PbANKA_07_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_07_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_07_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_07_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_07_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_08_v3_range_subset_plot
PbANKA_08_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_08_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_08_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_08_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_08_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_09_v3_range_subset_plot
PbANKA_09_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_09_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_09_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_09_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_09_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_10_v3_range_subset_plot
PbANKA_10_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_10_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_10_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_10_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_10_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_11_v3_range_subset_plot
PbANKA_11_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_11_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_11_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_11_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_11_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_12_v3_range_subset_plot
PbANKA_12_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_12_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_12_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_12_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_12_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_13_v3_range_subset_plot
PbANKA_13_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_13_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_13_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_13_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_13_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_14_v3_range_subset_plot
PbANKA_14_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_14_v3, 
                            aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_14_v3$ContigName, 
                                y = RPKM_CNV_subset_plot_data_Chr$PbANKA_14_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_14_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_API_v3_range_subset_plot
PbANKA_API_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_API_v3, 
                             aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_API_v3$ContigName, 
                                 y = RPKM_CNV_subset_plot_data_Chr$PbANKA_API_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_API_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_MIT_v3_range_subset_plot
PbANKA_MIT_v3_range_subset_plot <- ggplot(RPKM_CNV_subset_plot_data_Chr$PbANKA_MIT_v3, 
                             aes(x = RPKM_CNV_subset_plot_data_Chr$PbANKA_MIT_v3$ContigName, 
                                 y = RPKM_CNV_subset_plot_data_Chr$PbANKA_MIT_v3$RPKM)) +
  geom_point(aes(colour=RPKM_CNV_subset_plot_data_Chr$PbANKA_MIT_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2





######RPKM only GNP789,K173,820, range > 100#####

####load/make range subset line subset data#####
RPKM_GNPlines <- RPKM[,c(1:4,8:15)]

for (i in 1:nrow(RPKM_GNPlines)) {
  RPKM_GNPlines$min[i] <- min(RPKM_GNPlines[i,5:12])
  RPKM_GNPlines$max[i] <- max(RPKM_GNPlines[i,5:12])
  RPKM_GNPlines$range[i] <- (RPKM_GNPlines$max[i])-(RPKM_GNPlines$min[i])
}

####make plot range subset line subsetdata dataframe####
#subset dataframe to range > 100
RPKM_GNPlines_CNV_subset <- subset(RPKM_GNPlines, RPKM_GNPlines$range > 100)

write.csv(RPKM_GNPlines_CNV_subset, "InitAlignOverall_multicov_RPKM_GNPlines_greater_100.csv")

RPKM_GNPlines_CNV_subset_plot_data <- RPKM_GNPlines_CNV_subset[,1:12]
RPKM_GNPlines_CNV_subset_plot_data <- melt(RPKM_GNPlines_CNV_subset_plot_data, id.vars=1:4)
colnames(RPKM_GNPlines_CNV_subset_plot_data) <- c("Chr","ContigStart","ContigEnd","ContigName","Line","RPKM")
RPKM_GNPlines_CNV_subset_plot_data$Line <- gsub('X','',RPKM_GNPlines_CNV_subset_plot_data$Line)
RPKM_GNPlines_CNV_subset_plot_data$LineName <- RPKM_GNPlines_CNV_subset_plot_data$Line
RPKM_GNPlines_CNV_subset_plot_data$LineName <- gsub('_trimmed_aligned_sorted','',RPKM_GNPlines_CNV_subset_plot_data$LineName)
RPKM_GNPlines_CNV_subset_plot_data$LineName <- as.factor(RPKM_GNPlines_CNV_subset_plot_data$LineName)
head(levels(RPKM_GNPlines_CNV_subset_plot_data$ContigName))

#####range line subset plots#####
##get plot data ready - make a list of DFs, 1 DF for each chr
RPKM_GNPlines_CNV_subset_plot_data_Chr <- split(RPKM_GNPlines_CNV_subset_plot_data, RPKM_GNPlines_CNV_subset_plot_data$Chr)

####range subset plots coded individually#####
#PbANKA_00_v3_archived_contig_1_range_subset_plot
PbANKA_00_v3_archived_contig_1_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_1, 
                                                           aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_1$ContigName, 
                                                               y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_1$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_1$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  colScale2

#PbANKA_00_v3_archived_contig_2_range_subset_plot
PbANKA_00_v3_archived_contig_2_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_2, 
                                                           aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_2$ContigName, 
                                                               y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_2$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_2$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_00_v3_archived_contig_3_range_subset_plot
PbANKA_00_v3_archived_contig_3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_3, 
                                                           aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_3$ContigName, 
                                                               y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_00_v3_archived_contig_4_range_subset_plot
PbANKA_00_v3_archived_contig_4_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_4, 
                                                           aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_4$ContigName, 
                                                               y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_4$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_4$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_00_v3_archived_contig_5_range_subset_plot
PbANKA_00_v3_archived_contig_5_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_5, 
                                                           aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_5$ContigName, 
                                                               y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_5$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_00_v3_archived_contig_5$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_01_v3_range_subset_plot
PbANKA_01_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_01_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_01_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_01_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_01_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_02_v3_range_subset_plot
PbANKA_02_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_02_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_02_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_02_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_02_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_03_v3_range_subset_plot
PbANKA_03_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_03_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_03_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_03_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_03_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_04_v3_range_subset_plot
PbANKA_04_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_04_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_04_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_04_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_04_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_05_v3_range_subset_plot
PbANKA_05_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_05_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_05_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_05_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_05_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_06_v3_range_subset_plot
PbANKA_06_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_06_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_06_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_06_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_06_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_07_v3_range_subset_plot
PbANKA_07_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_07_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_07_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_07_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_07_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_08_v3_range_subset_plot
PbANKA_08_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_08_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_08_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_08_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_08_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_09_v3_range_subset_plot
PbANKA_09_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_09_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_09_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_09_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_09_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_10_v3_range_subset_plot
PbANKA_10_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_10_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_10_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_10_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_10_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_11_v3_range_subset_plot
PbANKA_11_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_11_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_11_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_11_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_11_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_12_v3_range_subset_plot
PbANKA_12_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_12_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_12_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_12_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_12_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_13_v3_range_subset_plot
PbANKA_13_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_13_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_13_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_13_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_13_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_14_v3_range_subset_plot
PbANKA_14_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_14_v3, 
                                         aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_14_v3$ContigName, 
                                             y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_14_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_14_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_API_v3_range_subset_plot
PbANKA_API_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_API_v3, 
                                          aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_API_v3$ContigName, 
                                              y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_API_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_API_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#PbANKA_MIT_v3_range_subset_plot
PbANKA_MIT_v3_range_line_subset_plot <- ggplot(RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_MIT_v3, 
                                          aes(x = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_MIT_v3$ContigName, 
                                              y = RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_MIT_v3$RPKM)) +
  geom_point(aes(colour=RPKM_GNPlines_CNV_subset_plot_data_Chr$PbANKA_MIT_v3$LineName), 
             size = 1) +
  xlab("Contig") +
  ylab("RPKM") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  colScale2

#####save range subset line subset plots#####
ggsave(PbANKA_00_v3_archived_contig_1_range_line_subset_plot, 
       file="PbANKA_00_v3_archived_contig_1_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_00_v3_archived_contig_2_range_line_subset_plot, 
       file="PbANKA_00_v3_archived_contig_2_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_00_v3_archived_contig_3_range_line_subset_plot, 
       file="PbANKA_00_v3_archived_contig_3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_00_v3_archived_contig_4_range_line_subset_plot, 
       file="PbANKA_00_v3_archived_contig_4_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_00_v3_archived_contig_5_range_line_subset_plot, 
       file="PbANKA_00_v3_archived_contig_5_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_01_v3_range_line_subset_plot, 
       file="PbANKA_01_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_02_v3_range_line_subset_plot, 
       file="PbANKA_02_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_03_v3_range_line_subset_plot, 
       file="PbANKA_03_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_04_v3_range_line_subset_plot, 
       file="PbANKA_04_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_05_v3_range_line_subset_plot, 
       file="PbANKA_05_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_06_v3_range_line_subset_plot, 
       file="PbANKA_06_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_07_v3_range_line_subset_plot, 
       file="PbANKA_07_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_08_v3_range_line_subset_plot, 
       file="PbANKA_08_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_09_v3_range_line_subset_plot, 
       file="PbANKA_09_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_10_v3_range_line_subset_plot, 
       file="PbANKA_10_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_11_v3_range_line_subset_plot, 
       file="PbANKA_11_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_12_v3_range_line_subset_plot, 
       file="PbANKA_12_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_13_v3_range_line_subset_plot, 
       file="PbANKA_13_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_14_v3_range_line_subset_plot, 
       file="PbANKA_14_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_API_v3_range_line_subset_plot, 
       file="PbANKA_API_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")

ggsave(PbANKA_MIT_v3_range_line_subset_plot, 
       file="PbANKA_MIT_v3_range_line_subset_plot.pdf",
       height = 10, width = 14, unit = "in")



####save environment####
save.image("J:/III/Waters/Group Members/Mallory/R/InitAlignOverall_Genome_RPKM.RData")
