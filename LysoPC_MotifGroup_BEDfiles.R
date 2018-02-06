#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/LysoPC_MotifGroup_BEDfiles.RData")
library("data.table")

#####set WD#####
setwd("J:/III/Waters/Group Members/Mallory/LysoPC_project/FIRE_clusterfiles/ManualMotifConsensus/OnlyGenesWithMotifRedo16.6.2017/")


#####MotifGroup clusterfile BED file creation to create FASTA files#####
#load cluster files
MotifGroup1 <- read.table("MotifGroup1_upreg_vs_pos8_background_nodownreg_strand_Niggi_noDups.txt", sep="", header=T)
MotifGroup1b <- read.table("MotifGroup1b_upreg_vs_pos8_background_nodownreg_strand_Niggi_noDups.txt", sep="", header=T)
MotifGroup2_test1 <- read.table("MotifGroup2_test1_pos8_background_vs_downreg_strand_Niggi_noDups.txt", sep="", header=T)
MotifGroup2_test2 <- read.table("MotifGroup2_test2_upreg_vs_pos8_background_vs_downreg_strand_Niggi_noDups.txt", sep="", header=T)
MotifGroup2_test3 <- read.table("MotifGroup2_test3_upreg_motifonly_vs_pos8_background_vs_downreg_strand_Niggi_noDups.txt", sep="", header=T)
MotifGroup3 <- read.table("MotifGroup3_upreg_vs_pos8_background_nodownreg_strand_Niggi_noDups.txt", sep="", header=T)
MotifGroup4 <- read.table("MotifGroup4_upreg_vs_pos8_background_nodownreg_strand_Niggi_noDups.txt", sep="", header=T)

#load BED file
setwd("J:/III/Waters/Group Members/Mallory/LysoPC_project/BEDLists/UpTo_2kb/")
pos8_allgenes_strand_Niggi_noDups_download_BED <- read.csv("pos8_allgenes_strand_Niggi_noDups_download_BED.csv")

#want to subset pos8_allgenes_strand_Niggi_noDups_download_BED to include only the genes in the cluter file  
pos8_MotifGroup1genes_strand_Niggi_noDups_download_BED <- subset(pos8_allgenes_strand_Niggi_noDups_download_BED, 
                                                                 pos8_allgenes_strand_Niggi_noDups_download_BED$NAME %in% MotifGroup1$orf)
pos8_MotifGroup1bgenes_strand_Niggi_noDups_download_BED <- subset(pos8_allgenes_strand_Niggi_noDups_download_BED, 
                                                                  pos8_allgenes_strand_Niggi_noDups_download_BED$NAME %in% MotifGroup1b$orf)    
pos8_MotifGroup2_test1genes_strand_Niggi_noDups_download_BED <- subset(pos8_allgenes_strand_Niggi_noDups_download_BED, 
                                                                       pos8_allgenes_strand_Niggi_noDups_download_BED$NAME %in% MotifGroup2_test1$orf)
pos8_MotifGroup2_test2genes_strand_Niggi_noDups_download_BED <- subset(pos8_allgenes_strand_Niggi_noDups_download_BED, 
                                                                       pos8_allgenes_strand_Niggi_noDups_download_BED$NAME %in% MotifGroup2_test2$orf)
pos8_MotifGroup2_test3genes_strand_Niggi_noDups_download_BED <- subset(pos8_allgenes_strand_Niggi_noDups_download_BED, 
                                                                       pos8_allgenes_strand_Niggi_noDups_download_BED$NAME %in% MotifGroup2_test3$orf)
pos8_MotifGroup3genes_strand_Niggi_noDups_download_BED <- subset(pos8_allgenes_strand_Niggi_noDups_download_BED, 
                                                                 pos8_allgenes_strand_Niggi_noDups_download_BED$NAME %in% MotifGroup3$orf)
pos8_MotifGroup4genes_strand_Niggi_noDups_download_BED <- subset(pos8_allgenes_strand_Niggi_noDups_download_BED, 
                                                                 pos8_allgenes_strand_Niggi_noDups_download_BED$NAME %in% MotifGroup4$orf)

#output BED files
write.csv(pos8_MotifGroup1genes_strand_Niggi_noDups_download_BED, "pos8_MotifGroup1genes_strand_Niggi_noDups_download_BED.csv")
write.csv(pos8_MotifGroup1bgenes_strand_Niggi_noDups_download_BED, "pos8_MotifGroup1bgenes_strand_Niggi_noDups_download_BED.csv")
write.csv(pos8_MotifGroup2_test1genes_strand_Niggi_noDups_download_BED, "pos8_MotifGroup2_test1genes_strand_Niggi_noDups_download_BED.csv")
write.csv(pos8_MotifGroup2_test2genes_strand_Niggi_noDups_download_BED, "pos8_MotifGroup2_test2genes_strand_Niggi_noDups_download_BED.csv")
write.csv(pos8_MotifGroup2_test3genes_strand_Niggi_noDups_download_BED, "pos8_MotifGroup2_test3genes_strand_Niggi_noDups_download_BED.csv")
write.csv(pos8_MotifGroup3genes_strand_Niggi_noDups_download_BED, "pos8_MotifGroup3genes_strand_Niggi_noDups_download_BED.csv")
write.csv(pos8_MotifGroup4genes_strand_Niggi_noDups_download_BED, "pos8_MotifGroup4genes_strand_Niggi_noDups_download_BED.csv")

#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/LysoPC_MotifGroup_BEDfiles.RData")
