#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/NishaRIPseq_up_down_ORF_FIREclusterfiles.RData")
library("data.table")

#####set WD#####
setwd("J:/III/Waters/Group Members/Mallory/NishaRNAseq/DataDownloads/")

#####load data#####
RIPseq_GoI <- read.csv("RIPseq_GoIs.csv", header=F)[,1]

PbANKA_33_ORF_genelist <- read.csv("PbANKA_33_ORF_regions_strand_genelist.csv")[,2]
PbANKA_33_upto1KB_upstream_genelist <- read.csv("PbANKA_33_upto1KB_upstream_regions_strand_genelist.csv")[,2]
PbANKA_33_upto1KB_downstream_genelist <- read.csv("PbANKA_33_upto1KB_downstream_regions_strand_genelist.csv")[,2]
PbANKA_33_1KB_upstream_genelist <- read.csv("PbANKA_33_1KB_upstream_regions_strand_genelist.csv")[,2]
PbANKA_33_1KB_downstream_genelist <- read.csv("PbANKA_33_1KB_downstream_regions_strand_genelist.csv")[,2]


####make dfs for cluster files####
PbANKA_33_ORF_clusterfile <- data.frame(orf = PbANKA_33_ORF_genelist)
      PbANKA_33_ORF_clusterfile$cluster <- NA
        for (i in seq_along(PbANKA_33_ORF_clusterfile$orf)){
          if (PbANKA_33_ORF_clusterfile$orf[i] %in% RIPseq_GoI) {
            PbANKA_33_ORF_clusterfile$cluster[i] <- 0
          } else { 
            PbANKA_33_ORF_clusterfile$cluster[i] <- 1
          }
        }
      PbANKA_33_ORF_clusterfile <- PbANKA_33_ORF_clusterfile[order(PbANKA_33_ORF_clusterfile$cluster),]

PbANKA_33_upto1KB_upstream_clusterfile <- data.frame(orf = PbANKA_33_upto1KB_upstream_genelist)
      PbANKA_33_upto1KB_upstream_clusterfile$cluster <- NA
        for (i in seq_along(PbANKA_33_upto1KB_upstream_clusterfile$orf)){
          if (PbANKA_33_upto1KB_upstream_clusterfile$orf[i] %in% RIPseq_GoI) {
            PbANKA_33_upto1KB_upstream_clusterfile$cluster[i] <- 0
          } else { 
            PbANKA_33_upto1KB_upstream_clusterfile$cluster[i] <- 1
          }
        }
      PbANKA_33_upto1KB_upstream_clusterfile <- PbANKA_33_upto1KB_upstream_clusterfile[order(PbANKA_33_upto1KB_upstream_clusterfile$cluster),]

PbANKA_33_upto1KB_downstream_clusterfile <- data.frame(orf = PbANKA_33_upto1KB_downstream_genelist)
      PbANKA_33_upto1KB_downstream_clusterfile$cluster <- NA
      for (i in seq_along(PbANKA_33_upto1KB_downstream_clusterfile$orf)){
        if (PbANKA_33_upto1KB_downstream_clusterfile$orf[i] %in% RIPseq_GoI) {
          PbANKA_33_upto1KB_downstream_clusterfile$cluster[i] <- 0
        } else { 
          PbANKA_33_upto1KB_downstream_clusterfile$cluster[i] <- 1
        }
      }
      PbANKA_33_upto1KB_downstream_clusterfile <- PbANKA_33_upto1KB_downstream_clusterfile[order(PbANKA_33_upto1KB_downstream_clusterfile$cluster),]
    
PbANKA_33_1KB_upstream_clusterfile <- data.frame(orf = PbANKA_33_1KB_upstream_genelist)
      PbANKA_33_1KB_upstream_clusterfile$cluster <- NA
      for (i in seq_along(PbANKA_33_1KB_upstream_clusterfile$orf)){
        if (PbANKA_33_1KB_upstream_clusterfile$orf[i] %in% RIPseq_GoI) {
          PbANKA_33_1KB_upstream_clusterfile$cluster[i] <- 0
        } else { 
          PbANKA_33_1KB_upstream_clusterfile$cluster[i] <- 1
        }
      }
      PbANKA_33_1KB_upstream_clusterfile <- PbANKA_33_1KB_upstream_clusterfile[order(PbANKA_33_1KB_upstream_clusterfile$cluster),]
      
PbANKA_33_1KB_downstream_clusterfile <- data.frame(orf = PbANKA_33_1KB_downstream_genelist)
      PbANKA_33_1KB_downstream_clusterfile$cluster <- NA
      for (i in seq_along(PbANKA_33_1KB_downstream_clusterfile$orf)){
        if (PbANKA_33_1KB_downstream_clusterfile$orf[i] %in% RIPseq_GoI) {
          PbANKA_33_1KB_downstream_clusterfile$cluster[i] <- 0
        } else { 
          PbANKA_33_1KB_downstream_clusterfile$cluster[i] <- 1
        }
      }
      PbANKA_33_1KB_downstream_clusterfile <- PbANKA_33_1KB_downstream_clusterfile[order(PbANKA_33_1KB_downstream_clusterfile$cluster),]
    
####export cluster files####   
      
write.csv(PbANKA_33_ORF_clusterfile,"J:/III/Waters/Group Members/Mallory/NishaRNAseq/FIREclusterfiles/PbANKA_33_ORF_clusterfile.csv")
write.csv(PbANKA_33_upto1KB_upstream_clusterfile,"J:/III/Waters/Group Members/Mallory/NishaRNAseq/FIREclusterfiles/PbANKA_33_upto1KB_upstream_clusterfile.csv")
write.csv(PbANKA_33_upto1KB_downstream_clusterfile,"J:/III/Waters/Group Members/Mallory/NishaRNAseq/FIREclusterfiles/PbANKA_33_upto1KB_downstream_clusterfile.csv")
write.csv(PbANKA_33_1KB_upstream_clusterfile,"J:/III/Waters/Group Members/Mallory/NishaRNAseq/FIREclusterfiles/PbANKA_33_1KB_upstream_clusterfile.csv")
write.csv(PbANKA_33_1KB_downstream_clusterfile,"J:/III/Waters/Group Members/Mallory/NishaRNAseq/FIREclusterfiles/PbANKA_33_1KB_downstream_clusterfile.csv")
      

      
#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/NishaRIPseq_up_down_ORF_FIREclusterfiles.RData")
