#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/NishaFileMerging_Oct_workspace.RData")

library(tidyr)
library(dplyr)
library(reshape2)
library(tidyverse)

####set wd####
setwd("J:/III/Waters/Group Members/Mallory/NishaFileMerging/")

#####load data#####
Genes_with_TM_domains <- read.csv("GenesByMotifSearch_TMHMM.csv", header = TRUE)
    colnames(Genes_with_TM_domains) <- c("GeneID","TMStart","TMEnd","TMSequence","TMTopology","TranscriptID")
  TMCount<- as.data.frame(table(Genes_with_TM_domains$TranscriptID))
    colnames(TMCount) <- c("TranscriptID", "TMCount")
  Genes_with_TM_domains <- merge(Genes_with_TM_domains, TMCount, by = "TranscriptID", all.x = TRUE)[,c(2:6,1,7)]
  

Genes_with_Motifs <- read.csv("GenesByMotifSearch_SummaryV2_edit.csv", header = TRUE)
  colnames(Genes_with_Motifs) <- c("GeneID", "SourceID","GenomicLocation","ProductDescription","MotifLocations","MotifCount","Motif")
  Genes_with_Motifs$MotifLocations <- as.character(Genes_with_Motifs$MotifLocations)
  Genes_with_Motifs$Motif <- as.character(Genes_with_Motifs$Motif)
  
    # Wants:
    #   GeneID
    #   MatchLocations(split)
    #   Motif(split) (does she actually need this?) 
    #   StartMin
    #   EndMax
    #   TMSequence (does she actually need this?)
    #   
    # one excel, with the gene ID, location of the motifs and the location of transmembrane domains
    
      #start by sorting out the Motif table
      #need to split MatchLocations by ",", duplicating the rows
  
####Split Motif table where multiple motifs exist in each gene####  
#split MotifLocations into new rows, duplicating the rest of the row where necessary (each Motif has its own row)
Genes_with_Motifs_split <- Genes_with_Motifs %>% 
    mutate(MotifLocations = strsplit(MotifLocations, ",")) %>% 
    unnest(MotifLocations)
  Genes_with_Motifs_split <- Genes_with_Motifs_split[,c(1:5,7,6)]
  #remove the parentheses from MotifLocations
    Genes_with_Motifs_split$MotifLocations <- gsub("[(,), ]", "", Genes_with_Motifs_split$MotifLocations)
  #split MotifLocations into MotifStart and MotifEnd
    Genes_with_Motifs_split <- separate(Genes_with_Motifs_split, MotifLocations, sep = "-", into = c("MotifStart","MotifEnd"), remove = FALSE)
  
  
    #separately split the Motifs into a list of the 2306 motifs
    # issue here...have split is at _, which is fine for most multiple occurences, but some are so close that aren't separated by a "_"
    # need to figure out the best way to split this
  Motifs_split <- unlist(strsplit(Genes_with_Motifs$Motif, "_"))
  write.csv(Motifs_split, "Motifs_split.csv")
  

####Merge split Motif Table with TM table - TM table came split (ie each row is 1 TM domain, could be duplicate GeneIDs if a gene has >1 TM domain)####  
  
Genes_with_Motifs_TMs <- merge(Genes_with_Motifs_split, Genes_with_TM_domains, by = "GeneID")  
  Genes_with_Motifs_TMs <- Genes_with_Motifs_TMs[,c("GeneID", "ProductDescription","MotifCount","MotifStart","MotifEnd","TMCount","TMStart","TMEnd","TMSequence")]

####add a column that will list whether there is any overlap between the motifs and TMs####
Genes_with_Motifs_TMs$MotifStart <- as.integer(Genes_with_Motifs_TMs$MotifStart)
Genes_with_Motifs_TMs$MotifEnd <- as.integer(Genes_with_Motifs_TMs$MotifEnd)
Genes_with_Motifs_TMs$TMStart <- as.integer(Genes_with_Motifs_TMs$TMStart)
Genes_with_Motifs_TMs$TMEnd <- as.integer(Genes_with_Motifs_TMs$TMEnd)

Genes_with_Motifs_TMs$OverlapPresent <- NA
for (i in seq_len(nrow(Genes_with_Motifs_TMs))){
  # Genes_with_Motifs_TMs$Overlap[i] <- print(i)
  if (Genes_with_Motifs_TMs$MotifStart[i] >= Genes_with_Motifs_TMs$TMStart[i] && 
      Genes_with_Motifs_TMs$MotifEnd[i] <= Genes_with_Motifs_TMs$TMEnd[i]) {
          Genes_with_Motifs_TMs$OverlapPresent[i] <- "Overlap"
  } else {
    Genes_with_Motifs_TMs$OverlapPresent[i] <- "None"
  }
}

###output final DF#####
write.csv(Genes_with_Motifs_TMs, "Genes_with_Motifs_TMs.csv")

#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/NishaFileMerging_Oct_workspace.RData")
