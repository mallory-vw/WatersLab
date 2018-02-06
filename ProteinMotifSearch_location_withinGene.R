#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/ProteinMotifSearch_location_withinGene.RData")


####load packages####
library(seqinr)
library(stringr)
library(stringi)
library(reshape2)

####set wd####
setwd("J:/III/Waters/Group Members/Mallory/KatieProteinMotifSearch/")

####load data####
Proteins_with_motif_TM <- read.fasta(file = "Inital_1355_gene_download.fasta",
                                     seqtype = "AA",
                                     as.string = TRUE)

gene_data <- read.csv("boolean_question_TranscriptRecordClasses_TranscriptRecordClass_Summary.csv")
  gene_data <- gene_data[,-3]
#####figure out whether last motif occurence is within 50bp of end of protein#####
MotifLocations <- str_locate_all(pattern = 'Y..[LIFMV]', Proteins_with_motif_TM)
  names(MotifLocations) <- names(Proteins_with_motif_TM)

LastMotifLocation <- lapply(MotifLocations, tail, n = 1L)
  names(LastMotifLocation) <- names(Proteins_with_motif_TM)

LastMotifLocation_df <- melt(LastMotifLocation)
  LastMotifLocation_df <- LastMotifLocation_df[,c("Var2","value","L1")]
  LastMotifLocation_df <- dcast(LastMotifLocation_df, L1 ~ Var2, value.vae = "value")
  names(LastMotifLocation_df) <- c("GeneID","MotifStart","MotifEnd")
  
LastMotifLocation_df <- merge(x = LastMotifLocation_df, y = gene_data, by.x = 'GeneID', by.y = 'source_id', all.x = TRUE, all.y=TRUE)  
  LastMotifLocation_df$MotifStart <- as.numeric(LastMotifLocation_df$MotifStart)
  LastMotifLocation_df$MotifEnd <- as.numeric(LastMotifLocation_df$MotifEnd)
  LastMotifLocation_df$Protein.Length <- as.numeric(LastMotifLocation_df$Protein.Length)
  
LastMotifLocation_df$DistanceToCTerm <- LastMotifLocation_df$Protein.Length - LastMotifLocation_df$MotifStart
  range(LastMotifLocation_df$DistanceToCTerm)

Genes_Near_Cterm <- subset(LastMotifLocation_df, LastMotifLocation_df$DistanceToCTerm <= 50)

####output list of genes with a motif near the Cterm###
write.csv(Genes_Near_Cterm, "Genes_Near_Cterm.csv")

#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/ProteinMotifSearch_location_withinGene.RData")
