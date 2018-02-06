# load environment
load("J:/III/Waters/Group Members/Mallory/R/InitAlignOverall_VCF_QC_plots_workspace.RData")

#load packages
library("ggplot2")
library("reshape2")
library("hexbin")

#set wd
setwd("J:/III/Waters/Group Members/Mallory/GenomeProjects/GenomesVariants/Plots")



#####LOAD DATA####
Q60_filter_qual_ao <- read.csv("InitAlignQ60_filter_QualAO.csv")
filter_qual_ao <- read.csv("InitAlign_filter_QualAO.csv")
qual_ao <- read.csv("InitAlign_QualAO.csv")

#issues with the files - some of the vars are multi allelic so it's treating AO as a factor, mega problem for the plots
#testing the filter_qual_ao by first removing the multiallelic loci and then plotting
Q60_filter_qual_ao_no_multi <- read.csv("InitAlignQ60_filter_QualAO_NoMulti.csv")
filter_qual_ao_no_multi <- read.csv("InitAlign_filter_QualAO_NoMulti.csv")
qual_ao_no_multi <- read.csv("InitAlign_QualAO_NoMulti.csv")

####plots#####
Qual_AO_Prefilter_Plot <- ggplot(qual_ao, aes(x=AO, y=QUAL)) +
  geom_bin2d(bins=100) + 
  xlab("AO") +
  ylab("QUAL") + 
  ylim(0,1600000) 

  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

Filter_Qual_AO_Plot <- ggplot(filter_qual_ao, aes(x=AO, y=QUAL)) +
    geom_bin2d() + 
    xlab("AO") +
    ylab("QUAL") + 
    ylim(0,35000) +
    xlim(0,1250) 

Q60_Filter_Qual_AO_Plot <- ggplot(Q60_filter_qual_ao, aes(x=AO, y=QUAL)) +
  geom_bin2d(bins=100) + 
  xlab("AO") +
  ylab("QUAL") + 
  ylim(0,35000) +
  xlim(0,1250)


#no multis#
Qual_AO_Prefilter_NoMulti_Plot <- ggplot(qual_ao_no_multi, aes(x=AO, y=QUAL)) +
  geom_bin2d(bins=100) + 
  xlab("AO") +
  ylab("QUAL") + 
  ylim(0,35000) +
  xlim(0,1200) 

Filter_Qual_AO_NoMulti_Plot <- ggplot(filter_qual_ao_no_multi, aes(x=AO, y=QUAL)) +
  geom_bin2d(bins=100) + 
  xlab("AO") +
  ylab("QUAL") + 
  ylim(0,35000) +
  xlim(0,100) 

Q60_Filter_Qual_AO_NoMulti_Plot <- ggplot(Q60_filter_qual_ao_no_multi, aes(x=AO, y=QUAL)) +
  geom_bin2d(bins=100) + 
  xlab("AO") +
  ylab("QUAL") + 
  ylim(0,35000) +
  xlim(0,100)



#####FINAL SAVE#####
save.image("J:/III/Waters/Group Members/Mallory/R/InitAlignOverall_VCF_QC_plots_workspace.RData")
