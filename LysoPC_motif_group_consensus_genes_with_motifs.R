#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/LysoPC_motif_group_consensus_genes_with_motifs.RData")

#####setwd#####
setwd("J:/III/Waters/Group Members/Mallory/LysoPC_project/GenesWithMotifs/FIRE/GeneListsFromMotifGroups")


#####load data#####
MotifGroup1_Motif1 <- read.csv("MotifGroup1_Motif1_AGTCACACTACACT.csv", header=F)
  colnames(MotifGroup1_Motif1) <- "GENEID"
MotifGroup1_Motif2 <- read.csv("MotifGroup1_Motif2_AGTCGATAAACAAGT.csv", header=F)
  colnames(MotifGroup1_Motif2) <- "GENEID"
MotifGroup1b_Motif1 <- read.csv("MotifGroup1b_Motif1_CACATCTACACTGTAGT.csv", header = F)                               
  colnames(MotifGroup1b_Motif1) <- "GENEID"
MotifGroup1b_Motif2 <- read.csv("MotifGroup1b_Motif2_AGTCACAAAGCAC.csv", header = F)
  colnames(MotifGroup1b_Motif2) <- "GENEID"
MotifGroup2_Motif1 <- read.csv("MotifGroup2_Motif1_AGTTACGTAGACGAT.csv", header = F)
  colnames(MotifGroup2_Motif1) <- "GENEID"
MotifGroup3_Motif1 <- read.csv("MotifGroup3_Motif1_TATATAGAGAG.csv", header = F)
  colnames(MotifGroup3_Motif1) <- "GENEID"
MotifGroup4_Motif1 <- read.csv("MotifGroup4_Motif1_ACCTAGTAAC.csv", header = F)
  colnames(MotifGroup4_Motif1) <- "GENEID"

  
MotifGroup_ConsensusMotifs <- read.csv("J:/III/Waters/Group Members/Mallory/LysoPC_Project/FIRE_results/FoundMotifGroup_consensusMotifs.csv", header=F)
  colnames(MotifGroup_ConsensusMotifs) <- c("Motif","RegExp")

background_noDownreg <- read.csv("J:/III/Waters/Group Members/Mallory/LysoPC_Project/GeneLists/UpTo_2kb/pos_background_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  background_noDownreg <- as.factor(background_noDownreg[,2])
downreg <- read.csv("J:/III/Waters/Group Members/Mallory/LysoPC_Project/GeneLists/downreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  downreg <- as.factor(downreg[,2])
upreg <- read.csv("J:/III/Waters/Group Members/Mallory/LysoPC_Project/GeneLists/upreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  upreg <- as.factor(upreg[,2])

#####add gene type to motif lists####
MotifGroup1_Motif1$Category <- NA
  for (i in seq_along(MotifGroup1_Motif1$GENEID)) {
    if (MotifGroup1_Motif1$GENEID[[i]] %in% upreg) {
      MotifGroup1_Motif1$Category[i] <- "upreg"
    } else if (MotifGroup1_Motif1$GENEID[i] %in% downreg) {
      MotifGroup1_Motif1$Category[i] <- "downreg"
    } else if (MotifGroup1_Motif1$GENEID[i] %in% background_noDownreg) {
      MotifGroup1_Motif1$Category[i] <- "background"
    } else {
      MotifGroup1_Motif1$Category[i] <- NA        
    }
  }
  
MotifGroup1_Motif2$Category <- NA
  for (i in seq_along(MotifGroup1_Motif2$GENEID)) {
    if (MotifGroup1_Motif2$GENEID[[i]] %in% upreg) {
      MotifGroup1_Motif2$Category[i] <- "upreg"
    } else if (MotifGroup1_Motif2$GENEID[i] %in% downreg) {
      MotifGroup1_Motif2$Category[i] <- "downreg"
    } else if (MotifGroup1_Motif2$GENEID[i] %in% background_noDownreg) {
      MotifGroup1_Motif2$Category[i] <- "background"
    } else {
      MotifGroup1_Motif2$Category[i] <- NA        
    }
  }
  
MotifGroup1b_Motif1$Category <- NA
  for (i in seq_along(MotifGroup1b_Motif1$GENEID)) {
    if (MotifGroup1b_Motif1$GENEID[[i]] %in% upreg) {
      MotifGroup1b_Motif1$Category[i] <- "upreg"
    } else if (MotifGroup1b_Motif1$GENEID[i] %in% downreg) {
      MotifGroup1b_Motif1$Category[i] <- "downreg"
    } else if (MotifGroup1b_Motif1$GENEID[i] %in% background_noDownreg) {
      MotifGroup1b_Motif1$Category[i] <- "background"
    } else {
      MotifGroup1b_Motif1$Category[i] <- NA        
    }
  }

MotifGroup1b_Motif2$Category <- NA
  for (i in seq_along(MotifGroup1b_Motif2$GENEID)) {
    if (MotifGroup1b_Motif2$GENEID[[i]] %in% upreg) {
      MotifGroup1b_Motif2$Category[i] <- "upreg"
    } else if (MotifGroup1b_Motif2$GENEID[i] %in% downreg) {
      MotifGroup1b_Motif2$Category[i] <- "downreg"
    } else if (MotifGroup1b_Motif2$GENEID[i] %in% background_noDownreg) {
      MotifGroup1b_Motif2$Category[i] <- "background"
    } else {
      MotifGroup1b_Motif2$Category[i] <- NA        
    }
  }

MotifGroup2_Motif1$Category <- NA
  for (i in seq_along(MotifGroup2_Motif1$GENEID)) {
    if (MotifGroup2_Motif1$GENEID[[i]] %in% upreg) {
      MotifGroup2_Motif1$Category[i] <- "upreg"
    } else if (MotifGroup2_Motif1$GENEID[i] %in% downreg) {
      MotifGroup2_Motif1$Category[i] <- "downreg"
    } else if (MotifGroup2_Motif1$GENEID[i] %in% background_noDownreg) {
      MotifGroup2_Motif1$Category[i] <- "background"
    } else {
      MotifGroup2_Motif1$Category[i] <- NA        
    }
  }

MotifGroup3_Motif1$Category <- NA
  for (i in seq_along(MotifGroup3_Motif1$GENEID)) {
    if (MotifGroup3_Motif1$GENEID[[i]] %in% upreg) {
      MotifGroup3_Motif1$Category[i] <- "upreg"
    } else if (MotifGroup3_Motif1$GENEID[i] %in% downreg) {
      MotifGroup3_Motif1$Category[i] <- "downreg"
    } else if (MotifGroup3_Motif1$GENEID[i] %in% background_noDownreg) {
      MotifGroup3_Motif1$Category[i] <- "background"
    } else {
      MotifGroup3_Motif1$Category[i] <- NA        
    }
  }
  
MotifGroup4_Motif1$Category <- NA
  for (i in seq_along(MotifGroup4_Motif1$GENEID)) {
    if (MotifGroup4_Motif1$GENEID[[i]] %in% upreg) {
      MotifGroup4_Motif1$Category[i] <- "upreg"
    } else if (MotifGroup4_Motif1$GENEID[i] %in% downreg) {
      MotifGroup4_Motif1$Category[i] <- "downreg"
    } else if (MotifGroup4_Motif1$GENEID[i] %in% background_noDownreg) {
      MotifGroup4_Motif1$Category[i] <- "background"
    } else {
      MotifGroup4_Motif1$Category[i] <- NA        
    }
  }

####add motif reg exp to data frame####
MotifGroup1_Motif1$Motif <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup1_Motif1"))$RegExp
MotifGroup1_Motif2$Motif <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup1_Motif2"))$RegExp
MotifGroup1b_Motif1$Motif <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup1b_Motif1"))$RegExp
MotifGroup1b_Motif2$Motif <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup1b_Motif2"))$RegExp
MotifGroup2_Motif1$Motif <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup2_Motif1"))$RegExp
MotifGroup3_Motif1$Motif <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup3_Motif1"))$RegExp
MotifGroup4_Motif1$Motif <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup4_Motif1"))$RegExp

MotifGroup1_Motif1$MotifName <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup1_Motif1"))$Motif
MotifGroup1_Motif2$MotifName <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup1_Motif2"))$Motif
MotifGroup1b_Motif1$MotifName <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup1b_Motif1"))$Motif
MotifGroup1b_Motif2$MotifName <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup1b_Motif2"))$Motif
MotifGroup2_Motif1$MotifName <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup2_Motif1"))$Motif
MotifGroup3_Motif1$MotifName <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup3_Motif1"))$Motif
MotifGroup4_Motif1$MotifName <- (subset(MotifGroup_ConsensusMotifs, MotifGroup_ConsensusMotifs$Motif == "MotifGroup4_Motif1"))$Motif

#####write out files#####
write.csv(MotifGroup1_Motif1, "MotifGroup1_Motif1_type.csv")
write.csv(MotifGroup1_Motif2, "MotifGroup1_Motif2_type.csv")
write.csv(MotifGroup1b_Motif1, "MotifGroup1b_Motif1_type.csv")
write.csv(MotifGroup1b_Motif2, "MotifGroup1b_Motif2_type.csv")
write.csv(MotifGroup2_Motif1, "MotifGroup2_Motif1_type.csv")
write.csv(MotifGroup3_Motif1, "MotifGroup3_Motif1_type.csv")
write.csv(MotifGroup4_Motif1, "MotifGroup4_Motif1_type.csv")



#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/LysoPC_motif_group_consensus_genes_with_motifs.RData")
