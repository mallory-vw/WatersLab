#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/LysoPC_motif_group_consensus_gene_list.RData")

#####setwd#####
setwd("J:/III/Waters/Group Members/Mallory/LysoPC_project/MotifClusteringConsensus/ManualConsensusGeneLists/")

#####load gene lists#####
setwd("J:/III/Waters/Group Members/Mallory/LysoPC_project/GenesWithMotifs/FIRE/GeneListsFromMotifs/")

Motif1 <- read.csv("Motif1_AGTCTTAAA.csv", header=F)
  Motif1 <- as.vector(Motif1[,1])
Motif2 <- read.csv("Motif2_ATAATAACG.csv", header=F)
  Motif2 <- as.vector(Motif2[,1])
Motif3 <- read.csv("Motif3_CACCATCTACACT.csv", header = F)
  Motif3 <- as.vector(Motif3[,1])
Motif4 <- read.csv("Motif4_CACACACTAGCACT.csv", header = F)
  Motif4 <- as.vector(Motif4[,1])
Motif5 <- read.csv("Motif5_ACCTAGTAAC.csv", header = F)
  Motif5 <- as.vector(Motif5[,1])
Motif6 <- read.csv("Motif6_ACGACAAACACTGTAT.csv", header = F)
  Motif6 <- as.vector(Motif6[,1])
Motif7 <- read.csv("Motif7_ACTAAGGTCGAT.csv", header = F)
  Motif7 <- as.vector(Motif7[,1])
Motif8 <- read.csv("Motif8_ACTCGTCTACTAACT.csv", header = F)
  Motif8 <- as.vector(Motif8[,1])
Motif9 <- read.csv("Motif9_ACTCGCTACT.csv", header = F)
  Motif9 <- as.vector(Motif9[,1])
Motif10 <- read.csv("Motif10_AGTAGTTACGACGAGAGT.csv", header = F)
  Motif10 <- as.vector(Motif10[,1])
Motif11 <- read.csv("Motif11_AGTCTACGGACT.csv", header = F)
  Motif11 <- as.vector(Motif11[,1])
Motif12 <- read.csv("Motif12_AGTCACACTCTACTCAC.csv", header = F)
  Motif12 <- as.vector(Motif12[,1])
Motif13 <- read.csv("Motif13_AGTCACACTCTACACTGT.csv", header = F)
  Motif13 <- as.vector(Motif13[,1])
Motif14 <- read.csv("Motif14_ATACACATAT.csv", header = F)
  Motif14 <- as.vector(Motif14[,1])
Motif15 <- read.csv("Motif15_ATCCGTAG.csv", header = F)
  Motif15 <- as.vector(Motif15[,1])
Motif16 <- read.csv("Motif16_ATCCGTCTACGTAACT.csv", header = F)
  Motif16 <- as.vector(Motif16[,1])
Motif17 <- read.csv("Motif17_ATTACTACGAGAAGT.csv", header = F)
  Motif17 <- as.vector(Motif17[,1])

####load gene types####  
setwd("J:/III/Waters/Group Members/Mallory/LysoPC_project/GeneLists")
  
  upreg_strand_Niggi_noDups_genelist <- read.csv("upreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
    upreg_strand_Niggi_noDups_genelist <- upreg_strand_Niggi_noDups_genelist[,2]
  downreg_strand_Niggi_noDups_genelist <- read.csv("downreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
    downreg_strand_Niggi_noDups_genelist <- downreg_strand_Niggi_noDups_genelist[,2]
  background_nodownreg_strand_Niggi_noDups_genelist <- read.csv("pos_background_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
    background_nodownreg_strand_Niggi_noDups_genelist <- background_nodownreg_strand_Niggi_noDups_genelist[,2]
  
  
  
#####add gene type to data frames#####
Motif1_Type <- as.data.frame(Motif1)
  colnames(Motif1_Type) <- c("GeneID")
  for (i in seq_along(Motif1_Type$GeneID)) {
      ifelse(Motif1_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif1_Type$GeneType[i] <- "upregulated", Motif1_Type$GeneType[i] <- NA)
      ifelse(Motif1_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif1_Type$GeneType[i] <- "downregulated", Motif1_Type$GeneType[i] <- Motif1_Type$GeneType[i])
      ifelse(Motif1_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif1_Type$GeneType[i] <- "background", Motif1_Type$GeneType[i] <- Motif1_Type$GeneType[i])
  }
  Motif1_Type$Motif <- "Motif1"
Motif2_Type <- as.data.frame(Motif2)
  colnames(Motif2_Type) <- c("GeneID")
  for (i in seq_along(Motif2_Type$GeneID)) {
    ifelse(Motif2_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif2_Type$GeneType[i] <- "upregulated", Motif2_Type$GeneType[i] <- NA)
    ifelse(Motif2_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif2_Type$GeneType[i] <- "downregulated", Motif2_Type$GeneType[i] <- Motif2_Type$GeneType[i])
    ifelse(Motif2_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif2_Type$GeneType[i] <- "background", Motif2_Type$GeneType[i] <- Motif2_Type$GeneType[i])
  }
  Motif2_Type$Motif <- "Motif2"
Motif3_Type <- as.data.frame(Motif3)
  colnames(Motif3_Type) <- c("GeneID")
  for (i in seq_along(Motif3_Type$GeneID)) {
    ifelse(Motif3_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif3_Type$GeneType[i] <- "upregulated", Motif3_Type$GeneType[i] <- NA)
    ifelse(Motif3_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif3_Type$GeneType[i] <- "downregulated", Motif3_Type$GeneType[i] <- Motif3_Type$GeneType[i])
    ifelse(Motif3_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif3_Type$GeneType[i] <- "background", Motif3_Type$GeneType[i] <- Motif3_Type$GeneType[i])
  }
  Motif3_Type$Motif <- "Motif3"
Motif4_Type <- as.data.frame(Motif4)
  colnames(Motif4_Type) <- c("GeneID")
  for (i in seq_along(Motif4_Type$GeneID)) {
    ifelse(Motif4_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif4_Type$GeneType[i] <- "upregulated", Motif4_Type$GeneType[i] <- NA)
    ifelse(Motif4_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif4_Type$GeneType[i] <- "downregulated", Motif4_Type$GeneType[i] <- Motif4_Type$GeneType[i])
    ifelse(Motif4_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif4_Type$GeneType[i] <- "background", Motif4_Type$GeneType[i] <- Motif4_Type$GeneType[i])
  }
  Motif4_Type$Motif <- "Motif4"
Motif5_Type <- as.data.frame(Motif5)
  colnames(Motif5_Type) <- c("GeneID")
  for (i in seq_along(Motif5_Type$GeneID)) {
    ifelse(Motif5_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif5_Type$GeneType[i] <- "upregulated", Motif5_Type$GeneType[i] <- NA)
    ifelse(Motif5_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif5_Type$GeneType[i] <- "downregulated", Motif5_Type$GeneType[i] <- Motif5_Type$GeneType[i])
    ifelse(Motif5_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif5_Type$GeneType[i] <- "background", Motif5_Type$GeneType[i] <- Motif5_Type$GeneType[i])
  }
  Motif5_Type$Motif <- "Motif5"
Motif6_Type <- as.data.frame(Motif6)
  colnames(Motif6_Type) <- c("GeneID")
  for (i in seq_along(Motif6_Type$GeneID)) {
    ifelse(Motif6_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif6_Type$GeneType[i] <- "upregulated", Motif6_Type$GeneType[i] <- NA)
    ifelse(Motif6_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif6_Type$GeneType[i] <- "downregulated", Motif6_Type$GeneType[i] <- Motif6_Type$GeneType[i])
    ifelse(Motif6_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif6_Type$GeneType[i] <- "background", Motif6_Type$GeneType[i] <- Motif6_Type$GeneType[i])
  }
  Motif6_Type$Motif <- "Motif6"
Motif7_Type <- as.data.frame(Motif7)
  colnames(Motif7_Type) <- c("GeneID")
  for (i in seq_along(Motif7_Type$GeneID)) {
    ifelse(Motif7_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif7_Type$GeneType[i] <- "upregulated", Motif7_Type$GeneType[i] <- NA)
    ifelse(Motif7_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif7_Type$GeneType[i] <- "downregulated", Motif7_Type$GeneType[i] <- Motif7_Type$GeneType[i])
    ifelse(Motif7_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif7_Type$GeneType[i] <- "background", Motif7_Type$GeneType[i] <- Motif7_Type$GeneType[i])
  }
  Motif7_Type$Motif <- "Motif7"
Motif8_Type <- as.data.frame(Motif8)
  colnames(Motif8_Type) <- c("GeneID")
  for (i in seq_along(Motif8_Type$GeneID)) {
    ifelse(Motif8_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif8_Type$GeneType[i] <- "upregulated", Motif8_Type$GeneType[i] <- NA)
    ifelse(Motif8_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif8_Type$GeneType[i] <- "downregulated", Motif8_Type$GeneType[i] <- Motif8_Type$GeneType[i])
    ifelse(Motif8_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif8_Type$GeneType[i] <- "background", Motif8_Type$GeneType[i] <- Motif8_Type$GeneType[i])
  }
  Motif8_Type$Motif <- "Motif8"
Motif9_Type <- as.data.frame(Motif9)
  colnames(Motif9_Type) <- c("GeneID")
  for (i in seq_along(Motif9_Type$GeneID)) {
    ifelse(Motif9_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif9_Type$GeneType[i] <- "upregulated", Motif9_Type$GeneType[i] <- NA)
    ifelse(Motif9_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif9_Type$GeneType[i] <- "downregulated", Motif9_Type$GeneType[i] <- Motif9_Type$GeneType[i])
    ifelse(Motif9_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif9_Type$GeneType[i] <- "background", Motif9_Type$GeneType[i] <- Motif9_Type$GeneType[i])
  }
  Motif9_Type$Motif <- "Motif9"
Motif10_Type <- as.data.frame(Motif10)
  colnames(Motif10_Type) <- c("GeneID")
  for (i in seq_along(Motif10_Type$GeneID)) {
    ifelse(Motif10_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif10_Type$GeneType[i] <- "upregulated", Motif10_Type$GeneType[i] <- NA)
    ifelse(Motif10_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif10_Type$GeneType[i] <- "downregulated", Motif10_Type$GeneType[i] <- Motif10_Type$GeneType[i])
    ifelse(Motif10_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif10_Type$GeneType[i] <- "background", Motif10_Type$GeneType[i] <- Motif10_Type$GeneType[i])
  }
  Motif10_Type$Motif <- "Motif10"
Motif11_Type <- as.data.frame(Motif11)
  colnames(Motif11_Type) <- c("GeneID")
  for (i in seq_along(Motif11_Type$GeneID)) {
    ifelse(Motif11_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif11_Type$GeneType[i] <- "upregulated", Motif11_Type$GeneType[i] <- NA)
    ifelse(Motif11_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif11_Type$GeneType[i] <- "downregulated", Motif11_Type$GeneType[i] <- Motif11_Type$GeneType[i])
    ifelse(Motif11_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif11_Type$GeneType[i] <- "background", Motif11_Type$GeneType[i] <- Motif11_Type$GeneType[i])
  }
  Motif11_Type$Motif <- "Motif11"
Motif12_Type <- as.data.frame(Motif12)
  colnames(Motif12_Type) <- c("GeneID")
  for (i in seq_along(Motif12_Type$GeneID)) {
    ifelse(Motif12_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif12_Type$GeneType[i] <- "upregulated", Motif12_Type$GeneType[i] <- NA)
    ifelse(Motif12_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif12_Type$GeneType[i] <- "downregulated", Motif12_Type$GeneType[i] <- Motif12_Type$GeneType[i])
    ifelse(Motif12_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif12_Type$GeneType[i] <- "background", Motif12_Type$GeneType[i] <- Motif12_Type$GeneType[i])
  }
  Motif12_Type$Motif <- "Motif12"
Motif13_Type <- as.data.frame(Motif13)
  colnames(Motif13_Type) <- c("GeneID")
  for (i in seq_along(Motif13_Type$GeneID)) {
    ifelse(Motif13_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif13_Type$GeneType[i] <- "upregulated", Motif13_Type$GeneType[i] <- NA)
    ifelse(Motif13_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif13_Type$GeneType[i] <- "downregulated", Motif13_Type$GeneType[i] <- Motif13_Type$GeneType[i])
    ifelse(Motif13_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif13_Type$GeneType[i] <- "background", Motif13_Type$GeneType[i] <- Motif13_Type$GeneType[i])
  }
  Motif13_Type$Motif <- "Motif13"
Motif14_Type <- as.data.frame(Motif14)
  colnames(Motif14_Type) <- c("GeneID")
  for (i in seq_along(Motif14_Type$GeneID)) {
    ifelse(Motif14_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif14_Type$GeneType[i] <- "upregulated", Motif14_Type$GeneType[i] <- NA)
    ifelse(Motif14_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif14_Type$GeneType[i] <- "downregulated", Motif14_Type$GeneType[i] <- Motif14_Type$GeneType[i])
    ifelse(Motif14_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif14_Type$GeneType[i] <- "background", Motif14_Type$GeneType[i] <- Motif14_Type$GeneType[i])
  }
  Motif14_Type$Motif <- "Motif14"
Motif15_Type <- as.data.frame(Motif15)
  colnames(Motif15_Type) <- c("GeneID")
  for (i in seq_along(Motif15_Type$GeneID)) {
    ifelse(Motif15_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif15_Type$GeneType[i] <- "upregulated", Motif15_Type$GeneType[i] <- NA)
    ifelse(Motif15_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif15_Type$GeneType[i] <- "downregulated", Motif15_Type$GeneType[i] <- Motif15_Type$GeneType[i])
    ifelse(Motif15_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif15_Type$GeneType[i] <- "background", Motif15_Type$GeneType[i] <- Motif15_Type$GeneType[i])
  }
  Motif15_Type$Motif <- "Motif15"
Motif16_Type <- as.data.frame(Motif16)
  colnames(Motif16_Type) <- c("GeneID")
  for (i in seq_along(Motif16_Type$GeneID)) {
    ifelse(Motif16_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif16_Type$GeneType[i] <- "upregulated", Motif16_Type$GeneType[i] <- NA)
    ifelse(Motif16_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif16_Type$GeneType[i] <- "downregulated", Motif16_Type$GeneType[i] <- Motif16_Type$GeneType[i])
    ifelse(Motif16_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif16_Type$GeneType[i] <- "background", Motif16_Type$GeneType[i] <- Motif16_Type$GeneType[i])
  }
  Motif16_Type$Motif <- "Motif16"
Motif17_Type <- as.data.frame(Motif17)
  colnames(Motif17_Type) <- c("GeneID")
  for (i in seq_along(Motif17_Type$GeneID)) {
    ifelse(Motif17_Type$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, Motif17_Type$GeneType[i] <- "upregulated", Motif17_Type$GeneType[i] <- NA)
    ifelse(Motif17_Type$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, Motif17_Type$GeneType[i] <- "downregulated", Motif17_Type$GeneType[i] <- Motif17_Type$GeneType[i])
    ifelse(Motif17_Type$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, Motif17_Type$GeneType[i] <- "background", Motif17_Type$GeneType[i] <- Motif17_Type$GeneType[i])
  }
  Motif17_Type$Motif <- "Motif17"
  
#####split into lists#####
Motif1_upreg <- as.character(subset(Motif1_Type, Motif1_Type$GeneType == "upregulated")$GeneID)
  Motif1_downreg <- as.character(subset(Motif1_Type, Motif1_Type$GeneType == "downregulated")$GeneID)
  Motif1_background <- as.character(subset(Motif1_Type, Motif1_Type$GeneType == "background")$GeneID)
Motif2_upreg <- as.character(subset(Motif2_Type, Motif2_Type$GeneType == "upregulated")$GeneID)
  Motif2_downreg <- as.character(subset(Motif2_Type, Motif2_Type$GeneType == "downregulated")$GeneID)
  Motif2_background <- as.character(subset(Motif2_Type, Motif2_Type$GeneType == "background")$GeneID)
Motif3_upreg <- as.character(subset(Motif3_Type, Motif3_Type$GeneType == "upregulated")$GeneID)
  Motif3_downreg <- as.character(subset(Motif3_Type, Motif3_Type$GeneType == "downregulated")$GeneID)
  Motif3_background <- as.character(subset(Motif3_Type, Motif3_Type$GeneType == "background")$GeneID)
Motif4_upreg <- as.character(subset(Motif4_Type, Motif4_Type$GeneType == "upregulated")$GeneID)
  Motif4_downreg <- as.character(subset(Motif4_Type, Motif4_Type$GeneType == "downregulated")$GeneID)
  Motif4_background <- as.character(subset(Motif4_Type, Motif4_Type$GeneType == "background")$GeneID)
Motif5_upreg <- as.character(subset(Motif5_Type, Motif5_Type$GeneType == "upregulated")$GeneID)
  Motif5_downreg <- as.character(subset(Motif5_Type, Motif5_Type$GeneType == "downregulated")$GeneID)
  Motif5_background <- as.character(subset(Motif5_Type, Motif5_Type$GeneType == "background")$GeneID)
Motif6_upreg <- as.character(subset(Motif6_Type, Motif6_Type$GeneType == "upregulated")$GeneID)
  Motif6_downreg <- as.character(subset(Motif6_Type, Motif6_Type$GeneType == "downregulated")$GeneID)
  Motif6_background <- as.character(subset(Motif6_Type, Motif6_Type$GeneType == "background")$GeneID)
Motif7_upreg <- as.character(subset(Motif7_Type, Motif7_Type$GeneType == "upregulated")$GeneID)
  Motif7_downreg <- as.character(subset(Motif7_Type, Motif7_Type$GeneType == "downregulated")$GeneID)
  Motif7_background <- as.character(subset(Motif7_Type, Motif7_Type$GeneType == "background")$GeneID)
Motif8_upreg <- as.character(subset(Motif8_Type, Motif8_Type$GeneType == "upregulated")$GeneID)
  Motif8_downreg <- as.character(subset(Motif8_Type, Motif8_Type$GeneType == "downregulated")$GeneID)
  Motif8_background <- as.character(subset(Motif8_Type, Motif8_Type$GeneType == "background")$GeneID)
Motif9_upreg <- as.character(subset(Motif9_Type, Motif9_Type$GeneType == "upregulated")$GeneID)
  Motif9_downreg <- as.character(subset(Motif9_Type, Motif9_Type$GeneType == "downregulated")$GeneID)
  Motif9_background <- as.character(subset(Motif9_Type, Motif9_Type$GeneType == "background")$GeneID)
Motif10_upreg <- as.character(subset(Motif10_Type, Motif10_Type$GeneType == "upregulated")$GeneID)
  Motif10_downreg <- as.character(subset(Motif10_Type, Motif10_Type$GeneType == "downregulated")$GeneID)
  Motif10_background <- as.character(subset(Motif10_Type, Motif10_Type$GeneType == "background")$GeneID)
Motif11_upreg <- as.character(subset(Motif11_Type, Motif11_Type$GeneType == "upregulated")$GeneID)
  Motif11_downreg <- as.character(subset(Motif11_Type, Motif11_Type$GeneType == "downregulated")$GeneID)
  Motif11_background <- as.character(subset(Motif11_Type, Motif11_Type$GeneType == "background")$GeneID)
Motif12_upreg <- as.character(subset(Motif12_Type, Motif12_Type$GeneType == "upregulated")$GeneID)
  Motif12_downreg <- as.character(subset(Motif12_Type, Motif12_Type$GeneType == "downregulated")$GeneID)
  Motif12_background <- as.character(subset(Motif12_Type, Motif12_Type$GeneType == "background")$GeneID)
Motif13_upreg <- as.character(subset(Motif13_Type, Motif13_Type$GeneType == "upregulated")$GeneID)
  Motif13_downreg <- as.character(subset(Motif13_Type, Motif13_Type$GeneType == "downregulated")$GeneID)
  Motif13_background <- as.character(subset(Motif13_Type, Motif13_Type$GeneType == "background")$GeneID)
Motif14_upreg <- as.character(subset(Motif14_Type, Motif14_Type$GeneType == "upregulated")$GeneID)
  Motif14_downreg <- as.character(subset(Motif14_Type, Motif14_Type$GeneType == "downregulated")$GeneID)
  Motif14_background <- as.character(subset(Motif14_Type, Motif14_Type$GeneType == "background")$GeneID)
Motif15_upreg <- as.character(subset(Motif15_Type, Motif15_Type$GeneType == "upregulated")$GeneID)
  Motif15_downreg <- as.character(subset(Motif15_Type, Motif15_Type$GeneType == "downregulated")$GeneID)
  Motif15_background <- as.character(subset(Motif15_Type, Motif15_Type$GeneType == "background")$GeneID)
Motif16_upreg <- as.character(subset(Motif16_Type, Motif16_Type$GeneType == "upregulated")$GeneID)
  Motif16_downreg <- as.character(subset(Motif16_Type, Motif16_Type$GeneType == "downregulated")$GeneID)
  Motif16_background <- as.character(subset(Motif16_Type, Motif16_Type$GeneType == "background")$GeneID)
Motif17_upreg <- as.character(subset(Motif17_Type, Motif17_Type$GeneType == "upregulated")$GeneID)
  Motif17_downreg <- as.character(subset(Motif17_Type, Motif17_Type$GeneType == "downregulated")$GeneID)
  Motif17_background <- as.character(subset(Motif17_Type, Motif17_Type$GeneType == "background")$GeneID)
  
  
  
  
  
  
#####SetGroups and VennDiagrams####
MotifGroup1 <- c("Motif3","Motif4","Motif12","Motif13","Motif14")
  MotifGroup1b <- c("Motif3","Motif4","Motif12","Motif13")
MotifGroup2 <- c("Motif8","Motif9","Motif16")
MotifGroup3 <- c("Motif10","Motif17")
MotifGroup4 <- c("Motif5","Motif11","Motif15")

library("compare")
library("VennDiagram")

# MotifGroup1 - 3,4,12,13,14
# MotifGroup2 - 8,9,16
# MotifGroup3 - 10,17
# MotifGroup4 - 11,15,5

###Motif Group 1###

MotifGroup1_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif3_Type, Motif4_Type, Motif12_Type, Motif13_Type,Motif14_Type))
  MotifGroup1_upreg <- subset(MotifGroup1_allgenes, MotifGroup1_allgenes$GeneType == "upregulated")

#compare lists between motifs#
#motif group 1 - overrepped in upreg genes --> only care about the upreg genes
  #motif3_upreg --> 142
  #motif4_upreg --> 157
  #motif12_upreg --> 148
  #motif13_upreg --> 91
  #motif14_upreg --> 259

#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup1_upreg_noDups <- MotifGroup1_upreg[!duplicated(MotifGroup1_upreg$GeneID),]
  MotifGroup1_upreg_noDups$NumMotifs <- NA
    for (i in 1:length(MotifGroup1_upreg_noDups$GeneID)){
      MotifGroup1_upreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup1_upreg,MotifGroup1_upreg$GeneID == MotifGroup1_upreg_noDups$GeneID[i]))
    }
MotifGroup1_upreg_noDups$Motif3 <- NA
for (i in 1:length(MotifGroup1_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_upreg_noDups$GeneID[i] %in% Motif3_upreg, 
         MotifGroup1_upreg_noDups$Motif3[i] <- "yes", 
         MotifGroup1_upreg_noDups$Motif3[i] <- "no")
}
MotifGroup1_upreg_noDups$Motif4 <- NA
for (i in 1:length(MotifGroup1_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_upreg_noDups$GeneID[i] %in% Motif4_upreg, 
         MotifGroup1_upreg_noDups$Motif4[i] <- "yes", 
         MotifGroup1_upreg_noDups$Motif4[i] <- "no")
}
MotifGroup1_upreg_noDups$Motif12 <- NA
for (i in 1:length(MotifGroup1_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_upreg_noDups$GeneID[i] %in% Motif12_upreg, 
         MotifGroup1_upreg_noDups$Motif12[i] <- "yes", 
         MotifGroup1_upreg_noDups$Motif12[i] <- "no")
}
MotifGroup1_upreg_noDups$Motif13 <- NA
for (i in 1:length(MotifGroup1_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_upreg_noDups$GeneID[i] %in% Motif13_upreg, 
         MotifGroup1_upreg_noDups$Motif13[i] <- "yes", 
         MotifGroup1_upreg_noDups$Motif13[i] <- "no")
}
MotifGroup1_upreg_noDups$Motif14 <- NA
for (i in 1:length(MotifGroup1_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_upreg_noDups$GeneID[i] %in% Motif14_upreg, 
         MotifGroup1_upreg_noDups$Motif14[i] <- "yes", 
         MotifGroup1_upreg_noDups$Motif14[i] <- "no")
}

##VennDiagram##
MotifGroup1_upreg_vennD_list <- list(Motif3_upreg,
                                     Motif4_upreg,
                                     Motif12_upreg,
                                     Motif13_upreg,
                                     Motif14_upreg)

venn.diagram(MotifGroup1_upreg_vennD_list, filename = "MotifGroup1_upreg_vennD.tiff", 
             category = MotifGroup1,
             fill = c("red","orange","yellow","green","blue"))
#70 genes shared between all 5

###Motif Group 1b###

MotifGroup1b_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif3_Type, Motif4_Type, Motif12_Type, Motif13_Type))
  MotifGroup1b_upreg <- subset(MotifGroup1b_allgenes, MotifGroup1b_allgenes$GeneType == "upregulated")


#compare lists between motifs#
    #motif group 1 - overrepped in upreg genes --> only care about the upreg genes
    #motif3_upreg --> 142
    #motif4_upreg --> 157
    #motif12_upreg --> 148
    #motif13_upreg --> 91

#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup1b_upreg_noDups <- MotifGroup1b_upreg[!duplicated(MotifGroup1b_upreg$GeneID),]
  MotifGroup1b_upreg_noDups$NumMotifs <- NA
    for (i in 1:length(MotifGroup1b_upreg_noDups$GeneID)){
      MotifGroup1b_upreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup1b_upreg,MotifGroup1b_upreg$GeneID == MotifGroup1b_upreg_noDups$GeneID[i]))
    }
MotifGroup1b_upreg_noDups$Motif3 <- NA
  for (i in 1:length(MotifGroup1b_upreg_noDups$GeneID)){
    ifelse(MotifGroup1b_upreg_noDups$GeneID[i] %in% Motif3_upreg, 
           MotifGroup1b_upreg_noDups$Motif3[i] <- "yes", 
           MotifGroup1b_upreg_noDups$Motif3[i] <- "no")
  }
MotifGroup1b_upreg_noDups$Motif4 <- NA
  for (i in 1:length(MotifGroup1b_upreg_noDups$GeneID)){
    ifelse(MotifGroup1b_upreg_noDups$GeneID[i] %in% Motif4_upreg, 
           MotifGroup1b_upreg_noDups$Motif4[i] <- "yes", 
           MotifGroup1b_upreg_noDups$Motif4[i] <- "no")
  }
MotifGroup1b_upreg_noDups$Motif12 <- NA
  for (i in 1:length(MotifGroup1b_upreg_noDups$GeneID)){
    ifelse(MotifGroup1b_upreg_noDups$GeneID[i] %in% Motif12_upreg, 
           MotifGroup1b_upreg_noDups$Motif12[i] <- "yes", 
           MotifGroup1b_upreg_noDups$Motif12[i] <- "no")
  }
MotifGroup1b_upreg_noDups$Motif13 <- NA
  for (i in 1:length(MotifGroup1b_upreg_noDups$GeneID)){
    ifelse(MotifGroup1b_upreg_noDups$GeneID[i] %in% Motif13_upreg, 
           MotifGroup1b_upreg_noDups$Motif13[i] <- "yes", 
           MotifGroup1b_upreg_noDups$Motif13[i] <- "no")
  }

MotifGroup1b_upreg_vennD_list <- list(Motif3_upreg,
                                      Motif4_upreg,
                                      Motif12_upreg,
                                      Motif13_upreg)

venn.diagram(MotifGroup1b_upreg_vennD_list, filename = "MotifGroup1b_upreg_vennD.tiff", 
             category = MotifGroup1b,
             fill = c("red","orange","yellow","green"))
#78 genes shared between all 4

#motif group 2 - underrep in upreg genes, maybe overrep in downreg
  #motif8_downreg --> 17
  #motif9_downreg --> 11
  #motif16_downreg --> 23

MotifGroup2_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif8_Type, Motif9_Type, Motif16_Type))
MotifGroup2_downreg <- subset(MotifGroup2_allgenes, MotifGroup2_allgenes$GeneType == "downregulated")

#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup2_downreg_noDups <- MotifGroup2_downreg[!duplicated(MotifGroup2_downreg$GeneID),]
MotifGroup2_downreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup2_downreg_noDups$GeneID)){
  MotifGroup2_downreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup2_downreg,MotifGroup2_downreg$GeneID == MotifGroup2_downreg_noDups$GeneID[i]))
}
MotifGroup2_downreg_noDups$Motif8 <- NA
for (i in 1:length(MotifGroup2_downreg_noDups$GeneID)){
  ifelse(MotifGroup2_downreg_noDups$GeneID[i] %in% Motif8_downreg, 
         MotifGroup2_downreg_noDups$Motif8[i] <- "yes", 
         MotifGroup2_downreg_noDups$Motif8[i] <- "no")
}
MotifGroup2_downreg_noDups$Motif9 <- NA
for (i in 1:length(MotifGroup2_downreg_noDups$GeneID)){
  ifelse(MotifGroup2_downreg_noDups$GeneID[i] %in% Motif9_downreg, 
         MotifGroup2_downreg_noDups$Motif9[i] <- "yes", 
         MotifGroup2_downreg_noDups$Motif9[i] <- "no")
}
MotifGroup2_downreg_noDups$Motif16 <- NA
for (i in 1:length(MotifGroup2_downreg_noDups$GeneID)){
  ifelse(MotifGroup2_downreg_noDups$GeneID[i] %in% Motif16_downreg, 
         MotifGroup2_downreg_noDups$Motif16[i] <- "yes", 
         MotifGroup2_downreg_noDups$Motif16[i] <- "no")
}

##VennDiagram##
MotifGroup2_downreg_vennD_list <- list(Motif8_downreg,
                                     Motif9_downreg,
                                     Motif16_downreg)

venn.diagram(MotifGroup2_downreg_vennD_list, filename = "MotifGroup2_downreg_vennD.tiff", 
             category = MotifGroup2,
             fill = c("red","orange","yellow"))
#7 genes shared between all 3, all of motif 8 in motif 16, 4 not in motif 8, 3 of those 4 are only motif 9

#motif group 3 --> underrep in upreg
  #if it's underrep in upreg genes, and not specifically in downreg, do i need to look at it?
  #could have negative effect, be a represor that's released?
  #motif10_upreg --> 15
  #motif17_upreg --> 17

MotifGroup3_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif10_Type, Motif17_Type))
MotifGroup3_upreg <- subset(MotifGroup3_allgenes, MotifGroup3_allgenes$GeneType == "upregulated")

#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup3_upreg_noDups <- MotifGroup3_upreg[!duplicated(MotifGroup3_upreg$GeneID),]
MotifGroup3_upreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup3_upreg_noDups$GeneID)){
  MotifGroup3_upreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup3_upreg,MotifGroup3_upreg$GeneID == MotifGroup3_upreg_noDups$GeneID[i]))
}
MotifGroup3_upreg_noDups$Motif10 <- NA
for (i in 1:length(MotifGroup3_upreg_noDups$GeneID)){
  ifelse(MotifGroup3_upreg_noDups$GeneID[i] %in% Motif10_upreg, 
         MotifGroup3_upreg_noDups$Motif10[i] <- "yes", 
         MotifGroup3_upreg_noDups$Motif10[i] <- "no")
}
MotifGroup3_upreg_noDups$Motif17 <- NA
for (i in 1:length(MotifGroup3_upreg_noDups$GeneID)){
  ifelse(MotifGroup3_upreg_noDups$GeneID[i] %in% Motif17_upreg, 
         MotifGroup3_upreg_noDups$Motif17[i] <- "yes", 
         MotifGroup3_upreg_noDups$Motif17[i] <- "no")
}

##VennDiagram##
MotifGroup3_upreg_vennD_list <- list(Motif10_upreg,
                                       Motif17_upreg)

venn.diagram(MotifGroup3_upreg_vennD_list, filename = "MotifGroup3_upreg_vennD.tiff", 
             category = MotifGroup3,
             fill = c("red","orange"))
#8 shared


#motifgroup 4 --> overrep in upreg
  #motif5_upreg --> 4 (out of just pfalc)
  #motif11_upreg --> 17
  #motif15_upreg --> 17
  
MotifGroup4_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif5_Type, Motif11_Type, Motif15_Type))
MotifGroup4_upreg <- subset(MotifGroup4_allgenes, MotifGroup4_allgenes$GeneType == "upregulated")

#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup4_upreg_noDups <- MotifGroup4_upreg[!duplicated(MotifGroup4_upreg$GeneID),]
MotifGroup4_upreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup4_upreg_noDups$GeneID)){
  MotifGroup4_upreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup4_upreg,MotifGroup4_upreg$GeneID == MotifGroup4_upreg_noDups$GeneID[i]))
}
MotifGroup4_upreg_noDups$Motif5 <- NA
for (i in 1:length(MotifGroup4_upreg_noDups$GeneID)){
  ifelse(MotifGroup4_upreg_noDups$GeneID[i] %in% Motif5_upreg, 
         MotifGroup4_upreg_noDups$Motif5[i] <- "yes", 
         MotifGroup4_upreg_noDups$Motif5[i] <- "no")
}
MotifGroup4_upreg_noDups$Motif11 <- NA
for (i in 1:length(MotifGroup4_upreg_noDups$GeneID)){
  ifelse(MotifGroup4_upreg_noDups$GeneID[i] %in% Motif11_upreg, 
         MotifGroup4_upreg_noDups$Motif11[i] <- "yes", 
         MotifGroup4_upreg_noDups$Motif11[i] <- "no")
}
MotifGroup4_upreg_noDups$Motif15 <- NA
for (i in 1:length(MotifGroup4_upreg_noDups$GeneID)){
  ifelse(MotifGroup4_upreg_noDups$GeneID[i] %in% Motif15_upreg, 
         MotifGroup4_upreg_noDups$Motif15[i] <- "yes", 
         MotifGroup4_upreg_noDups$Motif15[i] <- "no")
}

##VennDiagram##
MotifGroup4_upreg_vennD_list <- list(Motif5_upreg,
                                     Motif11_upreg,
                                     Motif15_upreg)

venn.diagram(MotifGroup4_upreg_vennD_list, filename = "MotifGroup4_upreg_vennD.tiff", 
             category = MotifGroup4,
             fill = c("red","orange","yellow"))
#complete overlap!






#####extract consensus list#####

# MotifGroup1 - 3,4,12,13,14
# MotifGroup2 - 8,9,16
# MotifGroup3 - 10,17
# MotifGroup4 - 11,15,5



MotifGroup1_genes <- Reduce(intersect, list(Motif3,Motif4,Motif12,Motif13,Motif14))
  MotifGroup1b_genes <- Reduce(intersect, list(Motif3,Motif4,Motif12,Motif13))
MotifGroup2_genes <- Reduce(intersect, list(Motif8,Motif9,Motif16))
MotifGroup3_genes <- Reduce(intersect, list(Motif10,Motif17))
MotifGroup4_genes <- Reduce(intersect, list(Motif5,Motif11,Motif15))

write.csv(MotifGroup1_genes, "MotifGroup1_genes.csv")
write.csv(MotifGroup1b_genes, "MotifGroup1b_genes.csv")
write.csv(MotifGroup2_genes, "MotifGroup2_genes.csv")
write.csv(MotifGroup3_genes, "MotifGroup3_genes.csv")
write.csv(MotifGroup4_genes, "MotifGroup4_genes.csv")

#####add in gene type#####


# MotifGroup1 - 3,4,12,13,14
# MotifGroup2 - 8,9,16
# MotifGroup3 - 10,17
# MotifGroup4 - 11,15,5



upreg_strand_Niggi_noDups_genelist <- read.csv("upreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  upreg_strand_Niggi_noDups_genelist <- upreg_strand_Niggi_noDups_genelist[,2]
downreg_strand_Niggi_noDups_genelist <- read.csv("downreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  downreg_strand_Niggi_noDups_genelist <- downreg_strand_Niggi_noDups_genelist[,2]
background_nodownreg_strand_Niggi_noDups_genelist <- read.csv("pos_background_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  background_nodownreg_strand_Niggi_noDups_genelist <- background_nodownreg_strand_Niggi_noDups_genelist[,2]

MotifGroup1_genes_types <- as.data.frame(MotifGroup1_genes)
  colnames(MotifGroup1_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup1_genes_types$GeneID)) {
  ifelse(MotifGroup1_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup1_genes_types$GeneType[i] <- "upregulated", MotifGroup1_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup1_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup1_genes_types$GeneType[i] <- "downregulated", MotifGroup1_genes_types$GeneType[i] <- MotifGroup1_genes_types$GeneType[i])
  ifelse(MotifGroup1_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup1_genes_types$GeneType[i] <- "background", MotifGroup1_genes_types$GeneType[i] <- MotifGroup1_genes_types$GeneType[i])
}

MotifGroup1b_genes_types <- as.data.frame(MotifGroup1b_genes)
colnames(MotifGroup1b_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup1b_genes_types$GeneID)) {
  ifelse(MotifGroup1b_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup1b_genes_types$GeneType[i] <- "upregulated", MotifGroup1b_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup1b_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup1b_genes_types$GeneType[i] <- "downregulated", MotifGroup1b_genes_types$GeneType[i] <- MotifGroup1b_genes_types$GeneType[i])
  ifelse(MotifGroup1b_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup1b_genes_types$GeneType[i] <- "background", MotifGroup1b_genes_types$GeneType[i] <- MotifGroup1b_genes_types$GeneType[i])
}

MotifGroup2_genes_types <- as.data.frame(MotifGroup2_genes)
  colnames(MotifGroup2_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup2_genes_types$GeneID)) {
  ifelse(MotifGroup2_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup2_genes_types$GeneType[i] <- "upregulated", MotifGroup2_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup2_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup2_genes_types$GeneType[i] <- "downregulated", MotifGroup2_genes_types$GeneType[i] <- MotifGroup2_genes_types$GeneType[i])
  ifelse(MotifGroup2_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup2_genes_types$GeneType[i] <- "background", MotifGroup2_genes_types$GeneType[i] <- MotifGroup2_genes_types$GeneType[i])
}
  
MotifGroup3_genes_types <- as.data.frame(MotifGroup3_genes)
colnames(MotifGroup3_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup3_genes_types$GeneID)) {
  ifelse(MotifGroup3_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup3_genes_types$GeneType[i] <- "upregulated", MotifGroup3_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup3_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup3_genes_types$GeneType[i] <- "downregulated", MotifGroup3_genes_types$GeneType[i] <- MotifGroup3_genes_types$GeneType[i])
  ifelse(MotifGroup3_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup3_genes_types$GeneType[i] <- "background", MotifGroup3_genes_types$GeneType[i] <- MotifGroup3_genes_types$GeneType[i])
}

MotifGroup4_genes_types <- as.data.frame(MotifGroup4_genes)
colnames(MotifGroup4_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup4_genes_types$GeneID)) {
  ifelse(MotifGroup4_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup4_genes_types$GeneType[i] <- "upregulated", MotifGroup4_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup4_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup4_genes_types$GeneType[i] <- "downregulated", MotifGroup4_genes_types$GeneType[i] <- MotifGroup4_genes_types$GeneType[i])
  ifelse(MotifGroup4_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup4_genes_types$GeneType[i] <- "background", MotifGroup4_genes_types$GeneType[i] <- MotifGroup4_genes_types$GeneType[i])
}

write.csv(MotifGroup1_genes_types, "MotifGroup1_genes_types.csv")
write.csv(MotifGroup1b_genes_types, "MotifGroup1b_genes_types.csv")
write.csv(MotifGroup2_genes_types, "MotifGroup2_genes_types.csv")
write.csv(MotifGroup3_genes_types, "MotifGroup3_genes_types.csv")
write.csv(MotifGroup4_genes_types, "MotifGroup4_genes_types.csv")



#####SetGroups and VennDiagrams 12.7.2017 MotifGroup1#####
MotifGroup1_test1 <- c("Motif4","Motif12","Motif13","Motif14")
MotifGroup1_test2 <- c("Motif3","Motif12","Motif13","Motif14")
MotifGroup1_test3 <- c("Motif3","Motif4","Motif13","Motif14")
MotifGroup1_test4 <- c("Motif3","Motif4","Motif12","Motif14")

library("compare")
library("VennDiagram")


# MotifGroup1 test 1 - 4,12,13,14
# MotifGroup1 test 2 - 3,12,13,14
# MotifGroup1 test 3 - 3,4,13,14
# MotifGroup1 test 4 - 3,4,12,14

###Motif Group 1 test1###
MotifGroup1_test1_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif4_Type, Motif12_Type, Motif13_Type, Motif14_Type))
  MotifGroup1_test1_upreg <- subset(MotifGroup1_test1_allgenes, MotifGroup1_test1_allgenes$GeneType == "upregulated")

#compare lists between motifs#
#motif group 1 - overrepped in upreg genes --> only care about the upreg genes
#Motif3_upreg --> 142
#Motif4_upreg --> 157
#Motif12_upreg --> 148
#Motif13_upreg --> 91
#Motif14_upreg --> 259

#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup1_test1_upreg_noDups <- MotifGroup1_test1_upreg[!duplicated(MotifGroup1_test1_upreg$GeneID),]
  MotifGroup1_test1_upreg_noDups$NumMotifs <- NA
    for (i in 1:length(MotifGroup1_test1_upreg_noDups$GeneID)){
      MotifGroup1_test1_upreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup1_test1_upreg,MotifGroup1_test1_upreg$GeneID == MotifGroup1_test1_upreg_noDups$GeneID[i]))
    }
MotifGroup1_test1_upreg_noDups$Motif4 <- NA
    for (i in 1:length(MotifGroup1_test1_upreg_noDups$GeneID)){
      ifelse(MotifGroup1_test1_upreg_noDups$GeneID[i] %in% Motif4_upreg, 
             MotifGroup1_test1_upreg_noDups$Motif4[i] <- "yes", 
             MotifGroup1_test1_upreg_noDups$Motif4[i] <- "no")
    }
MotifGroup1_test1_upreg_noDups$Motif12 <- NA
    for (i in 1:length(MotifGroup1_test1_upreg_noDups$GeneID)){
      ifelse(MotifGroup1_test1_upreg_noDups$GeneID[i] %in% Motif12_upreg, 
             MotifGroup1_test1_upreg_noDups$Motif12[i] <- "yes", 
             MotifGroup1_test1_upreg_noDups$Motif12[i] <- "no")
    }
MotifGroup1_test1_upreg_noDups$Motif13 <- NA
    for (i in 1:length(MotifGroup1_test1_upreg_noDups$GeneID)){
      ifelse(MotifGroup1_test1_upreg_noDups$GeneID[i] %in% Motif13_upreg, 
             MotifGroup1_test1_upreg_noDups$Motif13[i] <- "yes", 
             MotifGroup1_test1_upreg_noDups$Motif13[i] <- "no")
    }
MotifGroup1_test1_upreg_noDups$Motif14 <- NA
    for (i in 1:length(MotifGroup1_test1_upreg_noDups$GeneID)){
      ifelse(MotifGroup1_test1_upreg_noDups$GeneID[i] %in% Motif14_upreg, 
             MotifGroup1_test1_upreg_noDups$Motif14[i] <- "yes", 
             MotifGroup1_test1_upreg_noDups$Motif14[i] <- "no")
    }

##VennDiagram##
MotifGroup1_test1_upreg_vennD_list <- list(Motif4_upreg,
                                     Motif12_upreg,
                                     Motif13_upreg,
                                     Motif14_upreg)

venn.diagram(MotifGroup1_test1_upreg_vennD_list, filename = "MotifGroup1_test1_upreg_vennD.tiff", 
             category = MotifGroup1_test1,
             fill = c("orange","yellow","green","blue"))
#72 genes shared between all 4

###Motif Group 1 test2###
# MotifGroup1 test 2 - 3,12,13,14

MotifGroup1_test2_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif3_Type, Motif12_Type, Motif13_Type, Motif14_Type))
MotifGroup1_test2_upreg <- subset(MotifGroup1_test2_allgenes, MotifGroup1_test2_allgenes$GeneType == "upregulated")

#compare lists between motifs#
#motif group 1 - overrepped in upreg genes --> only care about the upreg genes
#Motif3_upreg --> 142
#Motif4_upreg --> 157
#Motif12_upreg --> 148
#Motif13_upreg --> 91
#Motif14_upreg --> 259

#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup1_test2_upreg_noDups <- MotifGroup1_test2_upreg[!duplicated(MotifGroup1_test2_upreg$GeneID),]
MotifGroup1_test2_upreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup1_test2_upreg_noDups$GeneID)){
  MotifGroup1_test2_upreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup1_test2_upreg,MotifGroup1_test2_upreg$GeneID == MotifGroup1_test2_upreg_noDups$GeneID[i]))
}
MotifGroup1_test2_upreg_noDups$Motif3 <- NA
for (i in 1:length(MotifGroup1_test2_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test2_upreg_noDups$GeneID[i] %in% Motif3_upreg, 
         MotifGroup1_test2_upreg_noDups$Motif4[i] <- "yes", 
         MotifGroup1_test2_upreg_noDups$Motif4[i] <- "no")
}
MotifGroup1_test2_upreg_noDups$Motif12 <- NA
for (i in 1:length(MotifGroup1_test2_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test2_upreg_noDups$GeneID[i] %in% Motif12_upreg, 
         MotifGroup1_test2_upreg_noDups$Motif12[i] <- "yes", 
         MotifGroup1_test2_upreg_noDups$Motif12[i] <- "no")
}
MotifGroup1_test2_upreg_noDups$Motif13 <- NA
for (i in 1:length(MotifGroup1_test2_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test2_upreg_noDups$GeneID[i] %in% Motif13_upreg, 
         MotifGroup1_test2_upreg_noDups$Motif13[i] <- "yes", 
         MotifGroup1_test2_upreg_noDups$Motif13[i] <- "no")
}
MotifGroup1_test2_upreg_noDups$Motif14 <- NA
for (i in 1:length(MotifGroup1_test2_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test2_upreg_noDups$GeneID[i] %in% Motif14_upreg, 
         MotifGroup1_test2_upreg_noDups$Motif14[i] <- "yes", 
         MotifGroup1_test2_upreg_noDups$Motif14[i] <- "no")
}

##VennDiagram##
MotifGroup1_test2_upreg_vennD_list <- list(Motif3_upreg,
                                           Motif12_upreg,
                                           Motif13_upreg,
                                           Motif14_upreg)

venn.diagram(MotifGroup1_test2_upreg_vennD_list, filename = "MotifGroup1_test2_upreg_vennD.tiff", 
             category = MotifGroup1_test2,
             fill = c("orange","yellow","green","blue"))
#72 genes shared between all 4



###Motif Group 1 test3###
# MotifGroup1 test 3 - 3,4,13,14

MotifGroup1_test3_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif3_Type, Motif4_Type, Motif13_Type, Motif14_Type))
MotifGroup1_test3_upreg <- subset(MotifGroup1_test3_allgenes, MotifGroup1_test3_allgenes$GeneType == "upregulated")

#compare lists between motifs#
#motif group 1 - overrepped in upreg genes --> only care about the upreg genes
#Motif3_upreg --> 142
#Motif4_upreg --> 157
#Motif12_upreg --> 148
#Motif13_upreg --> 91
#Motif14_upreg --> 259

#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup1_test3_upreg_noDups <- MotifGroup1_test3_upreg[!duplicated(MotifGroup1_test3_upreg$GeneID),]
MotifGroup1_test3_upreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup1_test3_upreg_noDups$GeneID)){
  MotifGroup1_test3_upreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup1_test3_upreg,MotifGroup1_test3_upreg$GeneID == MotifGroup1_test3_upreg_noDups$GeneID[i]))
}
MotifGroup1_test3_upreg_noDups$Motif3 <- NA
for (i in 1:length(MotifGroup1_test3_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test3_upreg_noDups$GeneID[i] %in% Motif3_upreg, 
         MotifGroup1_test3_upreg_noDups$Motif4[i] <- "yes", 
         MotifGroup1_test3_upreg_noDups$Motif4[i] <- "no")
}
MotifGroup1_test3_upreg_noDups$Motif4 <- NA
for (i in 1:length(MotifGroup1_test3_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test3_upreg_noDups$GeneID[i] %in% Motif4_upreg, 
         MotifGroup1_test3_upreg_noDups$Motif12[i] <- "yes", 
         MotifGroup1_test3_upreg_noDups$Motif12[i] <- "no")
}
MotifGroup1_test3_upreg_noDups$Motif13 <- NA
for (i in 1:length(MotifGroup1_test3_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test3_upreg_noDups$GeneID[i] %in% Motif13_upreg, 
         MotifGroup1_test3_upreg_noDups$Motif13[i] <- "yes", 
         MotifGroup1_test3_upreg_noDups$Motif13[i] <- "no")
}
MotifGroup1_test3_upreg_noDups$Motif14 <- NA
for (i in 1:length(MotifGroup1_test3_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test3_upreg_noDups$GeneID[i] %in% Motif14_upreg, 
         MotifGroup1_test3_upreg_noDups$Motif14[i] <- "yes", 
         MotifGroup1_test3_upreg_noDups$Motif14[i] <- "no")
}

##VennDiagram##
MotifGroup1_test3_upreg_vennD_list <- list(Motif3_upreg,
                                           Motif4_upreg,
                                           Motif13_upreg,
                                           Motif14_upreg)

venn.diagram(MotifGroup1_test3_upreg_vennD_list, filename = "MotifGroup1_test3_upreg_vennD.tiff", 
             category = MotifGroup1_test3,
             fill = c("orange","yellow","green","blue"))
#73 genes shared between all 4

###Motif Group 1 test4###
# MotifGroup1 test 4 - 3,4,12,14

MotifGroup1_test4_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif3_Type, Motif4_Type, Motif12_Type, Motif14_Type))
MotifGroup1_test4_upreg <- subset(MotifGroup1_test4_allgenes, MotifGroup1_test4_allgenes$GeneType == "upregulated")

#compare lists between motifs#
#motif group 1 - overrepped in upreg genes --> only care about the upreg genes
#Motif3_upreg --> 142
#Motif4_upreg --> 157
#Motif12_upreg --> 148
#Motif13_upreg --> 91
#Motif14_upreg --> 259

#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup1_test4_upreg_noDups <- MotifGroup1_test4_upreg[!duplicated(MotifGroup1_test4_upreg$GeneID),]
MotifGroup1_test4_upreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup1_test4_upreg_noDups$GeneID)){
  MotifGroup1_test4_upreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup1_test4_upreg,MotifGroup1_test4_upreg$GeneID == MotifGroup1_test4_upreg_noDups$GeneID[i]))
}
MotifGroup1_test4_upreg_noDups$Motif3 <- NA
for (i in 1:length(MotifGroup1_test4_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test4_upreg_noDups$GeneID[i] %in% Motif3_upreg, 
         MotifGroup1_test4_upreg_noDups$Motif4[i] <- "yes", 
         MotifGroup1_test4_upreg_noDups$Motif4[i] <- "no")
}
MotifGroup1_test4_upreg_noDups$Motif4 <- NA
for (i in 1:length(MotifGroup1_test4_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test4_upreg_noDups$GeneID[i] %in% Motif4_upreg, 
         MotifGroup1_test4_upreg_noDups$Motif12[i] <- "yes", 
         MotifGroup1_test4_upreg_noDups$Motif12[i] <- "no")
}
MotifGroup1_test4_upreg_noDups$Motif12 <- NA
for (i in 1:length(MotifGroup1_test4_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test4_upreg_noDups$GeneID[i] %in% Motif12_upreg, 
         MotifGroup1_test4_upreg_noDups$Motif13[i] <- "yes", 
         MotifGroup1_test4_upreg_noDups$Motif13[i] <- "no")
}
MotifGroup1_test4_upreg_noDups$Motif14 <- NA
for (i in 1:length(MotifGroup1_test4_upreg_noDups$GeneID)){
  ifelse(MotifGroup1_test4_upreg_noDups$GeneID[i] %in% Motif14_upreg, 
         MotifGroup1_test4_upreg_noDups$Motif14[i] <- "yes", 
         MotifGroup1_test4_upreg_noDups$Motif14[i] <- "no")
}

##VennDiagram##
MotifGroup1_test4_upreg_vennD_list <- list(Motif3_upreg,
                                           Motif4_upreg,
                                           Motif12_upreg,
                                           Motif14_upreg)

venn.diagram(MotifGroup1_test4_upreg_vennD_list, filename = "MotifGroup1_test4_upreg_vennD.tiff", 
             category = MotifGroup1_test4,
             fill = c("orange","yellow","green","blue"))
#97 genes shared between all 4



#####extract consensus list 12.7.2017 MotifGroup1#####
# MotifGroup1 test 1 - 4,12,13,14
# MotifGroup1 test 2 - 3,12,13,14
# MotifGroup1 test 3 - 3,4,13,14
# MotifGroup1 test 4 - 3,4,12,14


MotifGroup1_test1_genes <- Reduce(intersect, list(Motif4,Motif12,Motif13,Motif14))
MotifGroup1_test2_genes <- Reduce(intersect, list(Motif3,Motif12,Motif13,Motif14))
MotifGroup1_test3_genes <- Reduce(intersect, list(Motif3,Motif4,Motif13,Motif14))
MotifGroup1_test4_genes <- Reduce(intersect, list(Motif3,Motif4,Motif12,Motif14))

write.csv(MotifGroup1_test1_genes, "MotifGroup1_test1_genes.csv")
write.csv(MotifGroup1_test2_genes, "MotifGroup1_test2_genes.csv")
write.csv(MotifGroup1_test3_genes, "MotifGroup1_test3_genes.csv")
write.csv(MotifGroup1_test4_genes, "MotifGroup1_test4_genes.csv")



#####add in gene type 12.7.2017 MotifGroup1#####

# MotifGroup1 test 1 - 4,12,13,14
# MotifGroup1 test 2 - 3,12,13,14
# MotifGroup1 test 3 - 3,4,13,14
# MotifGroup1 test 4 - 3,4,12,14


MotifGroup1_test1_genes_types <- as.data.frame(MotifGroup1_test1_genes)
colnames(MotifGroup1_test1_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup1_test1_genes_types$GeneID)) {
  ifelse(MotifGroup1_test1_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup1_test1_genes_types$GeneType[i] <- "upregulated", MotifGroup1_test1_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup1_test1_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup1_test1_genes_types$GeneType[i] <- "downregulated", MotifGroup1_test1_genes_types$GeneType[i] <- MotifGroup1_test1_genes_types$GeneType[i])
  ifelse(MotifGroup1_test1_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup1_test1_genes_types$GeneType[i] <- "background", MotifGroup1_test1_genes_types$GeneType[i] <- MotifGroup1_test1_genes_types$GeneType[i])
}

MotifGroup1_test2_genes_types <- as.data.frame(MotifGroup1_test2_genes)
colnames(MotifGroup1_test2_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup1_test2_genes_types$GeneID)) {
  ifelse(MotifGroup1_test2_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup1_test2_genes_types$GeneType[i] <- "upregulated", MotifGroup1_test2_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup1_test2_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup1_test2_genes_types$GeneType[i] <- "downregulated", MotifGroup1_test2_genes_types$GeneType[i] <- MotifGroup1_test2_genes_types$GeneType[i])
  ifelse(MotifGroup1_test2_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup1_test2_genes_types$GeneType[i] <- "background", MotifGroup1_test2_genes_types$GeneType[i] <- MotifGroup1_test2_genes_types$GeneType[i])
}

MotifGroup1_test3_genes_types <- as.data.frame(MotifGroup1_test3_genes)
colnames(MotifGroup1_test3_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup1_test3_genes_types$GeneID)) {
  ifelse(MotifGroup1_test3_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup1_test3_genes_types$GeneType[i] <- "upregulated", MotifGroup1_test3_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup1_test3_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup1_test3_genes_types$GeneType[i] <- "downregulated", MotifGroup1_test3_genes_types$GeneType[i] <- MotifGroup1_test3_genes_types$GeneType[i])
  ifelse(MotifGroup1_test3_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup1_test3_genes_types$GeneType[i] <- "background", MotifGroup1_test3_genes_types$GeneType[i] <- MotifGroup1_test3_genes_types$GeneType[i])
}

MotifGroup1_test4_genes_types <- as.data.frame(MotifGroup1_test4_genes)
colnames(MotifGroup1_test4_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup1_test4_genes_types$GeneID)) {
  ifelse(MotifGroup1_test4_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup1_test4_genes_types$GeneType[i] <- "upregulated", MotifGroup1_test4_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup1_test4_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup1_test4_genes_types$GeneType[i] <- "downregulated", MotifGroup1_test4_genes_types$GeneType[i] <- MotifGroup1_test4_genes_types$GeneType[i])
  ifelse(MotifGroup1_test4_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup1_test4_genes_types$GeneType[i] <- "background", MotifGroup1_test4_genes_types$GeneType[i] <- MotifGroup1_test4_genes_types$GeneType[i])
}


write.csv(MotifGroup1_test1_genes_types, "MotifGroup1_test1_genes_types.csv")
write.csv(MotifGroup1_test2_genes_types, "MotifGroup1_test2_genes_types.csv")
write.csv(MotifGroup1_test3_genes_types, "MotifGroup1_test3_genes_types.csv")
write.csv(MotifGroup1_test4_genes_types, "MotifGroup1_test4_genes_types.csv")



#####SetGroups and VennDiagrams 13.7.2017 MotifGroup4#####
MotifGroup4_test1 <- c("Motif5","Motif15")
MotifGroup4_test2 <- c("Motif5","Motif11")
MotifGroup4_test3 <- c("Motif15","Motif11")

library("compare")
library("VennDiagram")


###Motif Group4 test1###
MotifGroup4_test1_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif5_Type, Motif15_Type))
MotifGroup4_test1_upreg <- subset(MotifGroup4_test1_allgenes, MotifGroup4_test1_allgenes$GeneType == "upregulated")

#compare lists between motifs#
#motif group 4 - overrepped in upreg genes --> only care about the upreg genes
#Motif5_upreg --> 4
#Motif15_upreg --> 17


#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup4_test1_upreg_noDups <- MotifGroup4_test1_upreg[!duplicated(MotifGroup4_test1_upreg$GeneID),]
MotifGroup4_test1_upreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup4_test1_upreg_noDups$GeneID)){
  MotifGroup4_test1_upreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup4_test1_upreg,MotifGroup4_test1_upreg$GeneID == 
                                                               MotifGroup4_test1_upreg_noDups$GeneID[i]))
}
MotifGroup4_test1_upreg_noDups$Motif5 <- NA
for (i in 1:length(MotifGroup4_test1_upreg_noDups$GeneID)){
  ifelse(MotifGroup4_test1_upreg_noDups$GeneID[i] %in% Motif5_upreg, 
         MotifGroup4_test1_upreg_noDups$Motif5[i] <- "yes", 
         MotifGroup4_test1_upreg_noDups$Motif5[i] <- "no")
}
MotifGroup4_test1_upreg_noDups$Motif15 <- NA
for (i in 1:length(MotifGroup4_test1_upreg_noDups$GeneID)){
  ifelse(MotifGroup4_test1_upreg_noDups$GeneID[i] %in% Motif15_upreg, 
         MotifGroup4_test1_upreg_noDups$Motif15[i] <- "yes", 
         MotifGroup4_test1_upreg_noDups$Motif15[i] <- "no")
}

##VennDiagram##
MotifGroup4_test1_upreg_vennD_list <- list(Motif5_upreg,
                                           Motif15_upreg)

venn.diagram(MotifGroup4_test1_upreg_vennD_list, filename = "MotifGroup4_test1_upreg_vennD.tiff", 
             category = MotifGroup4_test1,
             fill = c("blue","yellow"))
#4 genes shared between both

###Motif Group4 test2###
MotifGroup4_test2_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif5_Type, Motif11_Type))
MotifGroup4_test2_upreg <- subset(MotifGroup4_test2_allgenes, MotifGroup4_test2_allgenes$GeneType == "upregulated")

#compare lists between motifs#
#motif group 4 - overrepped in upreg genes --> only care about the upreg genes
#Motif5_upreg --> 4
#Motif11_upreg --> 17


#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup4_test2_upreg_noDups <- MotifGroup4_test2_upreg[!duplicated(MotifGroup4_test2_upreg$GeneID),]
MotifGroup4_test2_upreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup4_test2_upreg_noDups$GeneID)){
  MotifGroup4_test2_upreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup4_test2_upreg,MotifGroup4_test2_upreg$GeneID == 
                                                               MotifGroup4_test2_upreg_noDups$GeneID[i]))
}
MotifGroup4_test2_upreg_noDups$Motif5 <- NA
for (i in 1:length(MotifGroup4_test2_upreg_noDups$GeneID)){
  ifelse(MotifGroup4_test2_upreg_noDups$GeneID[i] %in% Motif5_upreg, 
         MotifGroup4_test2_upreg_noDups$Motif5[i] <- "yes", 
         MotifGroup4_test2_upreg_noDups$Motif5[i] <- "no")
}
MotifGroup4_test2_upreg_noDups$Motif11 <- NA
for (i in 1:length(MotifGroup4_test2_upreg_noDups$GeneID)){
  ifelse(MotifGroup4_test2_upreg_noDups$GeneID[i] %in% Motif11_upreg, 
         MotifGroup4_test2_upreg_noDups$Motif11[i] <- "yes", 
         MotifGroup4_test2_upreg_noDups$Motif11[i] <- "no")
}

##VennDiagram##
MotifGroup4_test2_upreg_vennD_list <- list(Motif5_upreg,
                                           Motif11_upreg)

venn.diagram(MotifGroup4_test2_upreg_vennD_list, filename = "MotifGroup4_test2_upreg_vennD.tiff", 
             category = MotifGroup4_test2,
             fill = c("blue","yellow"))
#4 genes shared between both

###Motif Group4 test3###
MotifGroup4_test3_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif11_Type, Motif15_Type))
MotifGroup4_test3_upreg <- subset(MotifGroup4_test3_allgenes, MotifGroup4_test3_allgenes$GeneType == "upregulated")

#compare lists between motifs#
#motif group 4 - overrepped in upreg genes --> only care about the upreg genes
#Motif11_upreg --> 17
#Motif15_upreg --> 17


#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup4_test3_upreg_noDups <- MotifGroup4_test3_upreg[!duplicated(MotifGroup4_test3_upreg$GeneID),]
MotifGroup4_test3_upreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup4_test3_upreg_noDups$GeneID)){
  MotifGroup4_test3_upreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup4_test3_upreg,MotifGroup4_test3_upreg$GeneID == 
                                                               MotifGroup4_test3_upreg_noDups$GeneID[i]))
}
MotifGroup4_test3_upreg_noDups$Motif11 <- NA
for (i in 1:length(MotifGroup4_test3_upreg_noDups$GeneID)){
  ifelse(MotifGroup4_test3_upreg_noDups$GeneID[i] %in% Motif11_upreg, 
         MotifGroup4_test3_upreg_noDups$Motif11[i] <- "yes", 
         MotifGroup4_test3_upreg_noDups$Motif11[i] <- "no")
}
MotifGroup4_test3_upreg_noDups$Motif15 <- NA
for (i in 1:length(MotifGroup4_test3_upreg_noDups$GeneID)){
  ifelse(MotifGroup4_test3_upreg_noDups$GeneID[i] %in% Motif15_upreg, 
         MotifGroup4_test3_upreg_noDups$Motif15[i] <- "yes", 
         MotifGroup4_test3_upreg_noDups$Motif15[i] <- "no")
}

##VennDiagram##
MotifGroup4_test3_upreg_vennD_list <- list(Motif11_upreg,
                                           Motif15_upreg)

venn.diagram(MotifGroup4_test3_upreg_vennD_list, filename = "MotifGroup4_test3_upreg_vennD.tiff", 
             category = MotifGroup4_test3,
             fill = c("blue","yellow"))
#17 genes shared between both



#####extract consensus list 12.7.2017 MotifGroup4#####
# MotifGroup4_test1 <- c("Motif5","Motif15")
# MotifGroup4_test2 <- c("Motif5","Motif11")
# MotifGroup4_test3 <- c("Motif15","Motif11")


MotifGroup4_test1_genes <- Reduce(intersect, list(Motif5,Motif15))
MotifGroup4_test2_genes <- Reduce(intersect, list(Motif5,Motif11))
MotifGroup4_test3_genes <- Reduce(intersect, list(Motif15,Motif11))

write.csv(MotifGroup4_test1_genes, "MotifGroup4_test1_genes.csv")
write.csv(MotifGroup4_test2_genes, "MotifGroup4_test2_genes.csv")
write.csv(MotifGroup4_test3_genes, "MotifGroup4_test3_genes.csv")

#####add in gene type 13.7.2017 MotifGroup4#####
# MotifGroup4_test1 <- c("Motif5","Motif15")
# MotifGroup4_test2 <- c("Motif5","Motif11")
# MotifGroup4_test3 <- c("Motif15","Motif11")

MotifGroup4_test1_genes_types <- as.data.frame(MotifGroup4_test1_genes)
colnames(MotifGroup4_test1_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup4_test1_genes_types$GeneID)) {
  ifelse(MotifGroup4_test1_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup4_test1_genes_types$GeneType[i] <- "upregulated", MotifGroup4_test1_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup4_test1_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup4_test1_genes_types$GeneType[i] <- "downregulated", MotifGroup4_test1_genes_types$GeneType[i] <- MotifGroup4_test1_genes_types$GeneType[i])
  ifelse(MotifGroup4_test1_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup4_test1_genes_types$GeneType[i] <- "background", MotifGroup4_test1_genes_types$GeneType[i] <- MotifGroup4_test1_genes_types$GeneType[i])
}

MotifGroup4_test2_genes_types <- as.data.frame(MotifGroup4_test2_genes)
colnames(MotifGroup4_test2_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup4_test2_genes_types$GeneID)) {
  ifelse(MotifGroup4_test2_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup4_test2_genes_types$GeneType[i] <- "upregulated", MotifGroup4_test2_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup4_test2_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup4_test2_genes_types$GeneType[i] <- "downregulated", MotifGroup4_test2_genes_types$GeneType[i] <- MotifGroup4_test2_genes_types$GeneType[i])
  ifelse(MotifGroup4_test2_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup4_test2_genes_types$GeneType[i] <- "background", MotifGroup4_test2_genes_types$GeneType[i] <- MotifGroup4_test2_genes_types$GeneType[i])
}

MotifGroup4_test3_genes_types <- as.data.frame(MotifGroup4_test3_genes)
colnames(MotifGroup4_test3_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup4_test3_genes_types$GeneID)) {
  ifelse(MotifGroup4_test3_genes_types$GeneID[i] %in% upreg_strand_Niggi_noDups_genelist, MotifGroup4_test3_genes_types$GeneType[i] <- "upregulated", MotifGroup4_test3_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup4_test3_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup4_test3_genes_types$GeneType[i] <- "downregulated", MotifGroup4_test3_genes_types$GeneType[i] <- MotifGroup4_test3_genes_types$GeneType[i])
  ifelse(MotifGroup4_test3_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup4_test3_genes_types$GeneType[i] <- "background", MotifGroup4_test3_genes_types$GeneType[i] <- MotifGroup4_test3_genes_types$GeneType[i])
}

write.csv(MotifGroup4_test1_genes_types, "MotifGroup4_test1_genes_types.csv")
write.csv(MotifGroup4_test2_genes_types, "MotifGroup4_test2_genes_types.csv")
write.csv(MotifGroup4_test3_genes_types, "MotifGroup4_test3_genes_types.csv")



#####SetGroups and VennDiagrams 13.7.2017 MotifGroup2#####
MotifGroup2_test1 <- c("Motif9","Motif16")
MotifGroup2_test2 <- c("Motif8","Motif16")
MotifGroup2_test3 <- c("Motif8","Motif9")

library("compare")
library("VennDiagram")


###Motif Group2 test1###
MotifGroup2_test1_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif9_Type, Motif16_Type))
MotifGroup2_test1_downreg <- subset(MotifGroup2_test1_allgenes, MotifGroup2_test1_allgenes$GeneType == "downregulated")

#compare lists between motifs#
#motif group 2 - overrepped in downreg genes --> only care about the downreg genes
#Motif9_upreg --> 10
#Motif16_upreg --> 48


#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup2_test1_downreg_noDups <- MotifGroup2_test1_downreg[!duplicated(MotifGroup2_test1_downreg$GeneID),]
MotifGroup2_test1_downreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup2_test1_downreg_noDups$GeneID)){
  MotifGroup2_test1_downreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup2_test1_downreg,MotifGroup2_test1_downreg$GeneID == 
                                                               MotifGroup2_test1_downreg_noDups$GeneID[i]))
}
MotifGroup2_test1_downreg_noDups$Motif9 <- NA
for (i in 1:length(MotifGroup2_test1_downreg_noDups$GeneID)){
  ifelse(MotifGroup2_test1_downreg_noDups$GeneID[i] %in% Motif9_downreg, 
         MotifGroup2_test1_downreg_noDups$Motif9[i] <- "yes", 
         MotifGroup2_test1_downreg_noDups$Motif9[i] <- "no")
}
MotifGroup2_test1_downreg_noDups$Motif16 <- NA
for (i in 1:length(MotifGroup2_test1_downreg_noDups$GeneID)){
  ifelse(MotifGroup2_test1_downreg_noDups$GeneID[i] %in% Motif16_downreg, 
         MotifGroup2_test1_downreg_noDups$Motif16[i] <- "yes", 
         MotifGroup2_test1_downreg_noDups$Motif16[i] <- "no")
}

##VennDiagram##
MotifGroup2_test1_downreg_vennD_list <- list(Motif9_downreg,
                                           Motif16_downreg)

venn.diagram(MotifGroup2_test1_downreg_vennD_list, filename = "MotifGroup2_test1_downreg_vennD.tiff", 
             category = MotifGroup2_test1,
             fill = c("blue","yellow"))
#8 genes shared between both

###Motif Group2test2###
MotifGroup2_test2_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif8_Type, Motif16_Type))
MotifGroup2_test2_downreg <- subset(MotifGroup2_test2_allgenes, MotifGroup2_test2_allgenes$GeneType == "downregulated")

#compare lists between motifs#
#motif group 2 - overrepped in downreg genes --> only care about the downreg genes
#Motif8_downreg --> 17
#Motif16_downreg --> 23


#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup2_test2_downreg_noDups <- MotifGroup2_test2_downreg[!duplicated(MotifGroup2_test2_downreg$GeneID),]
MotifGroup2_test2_downreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup2_test2_downreg_noDups$GeneID)){
  MotifGroup2_test2_downreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup2_test2_downreg,MotifGroup2_test2_downreg$GeneID == 
                                                               MotifGroup2_test2_downreg_noDups$GeneID[i]))
}
MotifGroup2_test2_downreg_noDups$Motif8 <- NA
for (i in 1:length(MotifGroup2_test2_downreg_noDups$GeneID)){
  ifelse(MotifGroup2_test2_downreg_noDups$GeneID[i] %in% Motif8_downreg, 
         MotifGroup2_test2_downreg_noDups$Motif8[i] <- "yes", 
         MotifGroup2_test2_downreg_noDups$Motif8[i] <- "no")
}
MotifGroup2_test2_downreg_noDups$Motif16 <- NA
for (i in 1:length(MotifGroup2_test2_downreg_noDups$GeneID)){
  ifelse(MotifGroup2_test2_downreg_noDups$GeneID[i] %in% Motif16_downreg, 
         MotifGroup2_test2_downreg_noDups$Motif16[i] <- "yes", 
         MotifGroup2_test2_downreg_noDups$Motif16[i] <- "no")
}

##VennDiagram##
MotifGroup2_test2_downreg_vennD_list <- list(Motif8_downreg,
                                           Motif16_downreg)

venn.diagram(MotifGroup2_test2_downreg_vennD_list, filename = "MotifGroup2_test2_downreg_vennD.tiff", 
             category = MotifGroup2_test2,
             fill = c("blue","yellow"))
#17 genes shared between both

###Motif Group2 test3###
MotifGroup2_test3_allgenes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(Motif8_Type, Motif9_Type))
MotifGroup2_test3_downreg <- subset(MotifGroup2_test3_allgenes, MotifGroup2_test3_allgenes$GeneType == "downregulated")

#compare lists between motifs#
#motif group 4 - overrepped in downreg genes --> only care about the downreg genes
#Motif8_downreg --> 17
#Motif9_downreg --> 11


#want a data frame with GeneID, GeneType, NumMotifs,Motif3,Motif4,Motif12,Motif13,Motif14
MotifGroup2_test3_downreg_noDups <- MotifGroup2_test3_downreg[!duplicated(MotifGroup2_test3_downreg$GeneID),]
MotifGroup2_test3_downreg_noDups$NumMotifs <- NA
for (i in 1:length(MotifGroup2_test3_downreg_noDups$GeneID)){
  MotifGroup2_test3_downreg_noDups$NumMotifs[i] <- nrow(subset(MotifGroup2_test3_downreg,MotifGroup2_test3_downreg$GeneID == 
                                                               MotifGroup2_test3_downreg_noDups$GeneID[i]))
}
MotifGroup2_test3_downreg_noDups$Motif8<- NA
for (i in 1:length(MotifGroup2_test3_downreg_noDups$GeneID)){
  ifelse(MotifGroup2_test3_downreg_noDups$GeneID[i] %in% Motif8_downreg, 
         MotifGroup2_test3_downreg_noDups$Motif8[i] <- "yes", 
         MotifGroup2_test3_downreg_noDups$Motif8[i] <- "no")
}
MotifGroup2_test3_downreg_noDups$Motif9 <- NA
for (i in 1:length(MotifGroup2_test3_downreg_noDups$GeneID)){
  ifelse(MotifGroup2_test3_downreg_noDups$GeneID[i] %in% Motif9_downreg, 
         MotifGroup2_test3_downreg_noDups$Motif9[i] <- "yes", 
         MotifGroup2_test3_downreg_noDups$Motif9[i] <- "no")
}

##VennDiagram##
MotifGroup2_test3_downreg_vennD_list <- list(Motif8_downreg,
                                           Motif9_downreg)

venn.diagram(MotifGroup2_test3_downreg_vennD_list, filename = "MotifGroup2_test3_downreg_vennD.tiff", 
             category = MotifGroup2_test3,
             fill = c("blue","yellow"))
#7 genes shared between both



#####extract consensus list 12.7.2017 MotifGroup2#####
# MotifGroup2_test1 <- c("Motif9","Motif16")
# MotifGroup2_test2 <- c("Motif8","Motif16")
# MotifGroup2_test3 <- c("Motif8","Motif9")


MotifGroup2_test1_genes <- Reduce(intersect, list(Motif9,Motif16))
MotifGroup2_test2_genes <- Reduce(intersect, list(Motif8,Motif16))
MotifGroup2_test3_genes <- Reduce(intersect, list(Motif8,Motif9))

write.csv(MotifGroup2_test1_genes, "MotifGroup2_test1_genes.csv")
write.csv(MotifGroup2_test2_genes, "MotifGroup2_test2_genes.csv")
write.csv(MotifGroup2_test3_genes, "MotifGroup2_test3_genes.csv")

#####add in gene type 13.7.2017 MotifGroup2#####
# MotifGroup2_test1 <- c("Motif9","Motif16")
# MotifGroup2_test2 <- c("Motif8","Motif16")
# MotifGroup2_test3 <- c("Motif8","Motif9")

MotifGroup2_test1_genes_types <- as.data.frame(MotifGroup2_test1_genes)
colnames(MotifGroup2_test1_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup2_test1_genes_types$GeneID)) {
  ifelse(MotifGroup2_test1_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup2_test1_genes_types$GeneType[i] <- "downregulated", MotifGroup2_test1_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup2_test1_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup2_test1_genes_types$GeneType[i] <- "downregulated", MotifGroup2_test1_genes_types$GeneType[i] <- MotifGroup2_test1_genes_types$GeneType[i])
  ifelse(MotifGroup2_test1_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup2_test1_genes_types$GeneType[i] <- "background", MotifGroup2_test1_genes_types$GeneType[i] <- MotifGroup2_test1_genes_types$GeneType[i])
}

MotifGroup2_test2_genes_types <- as.data.frame(MotifGroup2_test2_genes)
colnames(MotifGroup2_test2_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup2_test2_genes_types$GeneID)) {
  ifelse(MotifGroup2_test2_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup2_test2_genes_types$GeneType[i] <- "downregulated", MotifGroup2_test2_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup2_test2_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup2_test2_genes_types$GeneType[i] <- "downregulated", MotifGroup2_test2_genes_types$GeneType[i] <- MotifGroup2_test2_genes_types$GeneType[i])
  ifelse(MotifGroup2_test2_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup2_test2_genes_types$GeneType[i] <- "background", MotifGroup2_test2_genes_types$GeneType[i] <- MotifGroup2_test2_genes_types$GeneType[i])
}

MotifGroup2_test3_genes_types <- as.data.frame(MotifGroup2_test3_genes)
colnames(MotifGroup2_test3_genes_types) <- "GeneID"
for (i in seq_along(MotifGroup2_test3_genes_types$GeneID)) {
  ifelse(MotifGroup2_test3_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup2_test3_genes_types$GeneType[i] <- "downregulated", MotifGroup2_test3_genes_types$GeneType[i] <- NA)
  ifelse(MotifGroup2_test3_genes_types$GeneID[i] %in% downreg_strand_Niggi_noDups_genelist, MotifGroup2_test3_genes_types$GeneType[i] <- "downregulated", MotifGroup2_test3_genes_types$GeneType[i] <- MotifGroup2_test3_genes_types$GeneType[i])
  ifelse(MotifGroup2_test3_genes_types$GeneID[i] %in% background_nodownreg_strand_Niggi_noDups_genelist, MotifGroup2_test3_genes_types$GeneType[i] <- "background", MotifGroup2_test3_genes_types$GeneType[i] <- MotifGroup2_test3_genes_types$GeneType[i])
}

write.csv(MotifGroup2_test1_genes_types, "MotifGroup2_test1_genes_types.csv")
write.csv(MotifGroup2_test2_genes_types, "MotifGroup2_test2_genes_types.csv")
write.csv(MotifGroup2_test3_genes_types, "MotifGroup2_test3_genes_types.csv")












#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/LysoPC_motif_group_consensus_gene_list.RData")
