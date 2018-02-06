####load workspace#####
load("J:/III/Waters/Group Members/Mallory/NishaPullingOutGoIRobynData/GoIParsingNisha_RobynData_Workspace.RData")
#####setwd#####
setwd("J:/III/Waters/Group Members/Mallory/NishaPullingOutGoIRobynData/")

#####load data#####
NishaGeneList_IP2_IP3 <- read.csv("Alba3-Pulldownenriched_IP2_IP3.csv")
#remove last digit from Nisha's GeneIDs
for (j in seq_along(NishaGeneList_IP2_IP3$GeneID)){
  NishaGeneList_IP2_IP3$GeneID2[j] <- substr(NishaGeneList_IP2_IP3$GeneID[j],start = 0,stop = 13)
}

NishaGeneList_IP2 <- read.csv("Alba3-Pulldownenriched_IP2.csv")
#remove last digit from Nisha's GeneIDs
for (j in seq_along(NishaGeneList_IP2$GeneID)){
  NishaGeneList_IP2$GeneID2[j] <- substr(NishaGeneList_IP2$GeneID[j],start = 0,stop = 13)
}

NishaGeneList_IP3 <- read.csv("Alba3-Pulldownenriched_IP3.csv")
#remove last digit from Nisha's GeneIDs
for (j in seq_along(NishaGeneList_IP3$GeneID)){
  NishaGeneList_IP3$GeneID2[j] <- substr(NishaGeneList_IP3$GeneID[j],start = 0,stop = 13)
}

filenames <- c("collated_withCluster_PfalcOrtho_2h.csv",
                "collated_withCluster_PfalcOrtho_4h.csv",
                "collated_withCluster_PfalcOrtho_6h.csv",
                "collated_withCluster_PfalcOrtho_8h.csv",
                "collated_withCluster_PfalcOrtho_12h.csv",
                "collated_withCluster_PfalcOrtho_24h.csv",
                "collated_withCluster_PfalcOrtho_30h.csv",
                "collated_withCluster_PfalcOrtho_44h.csv")

timepoint_data_NishaSort <- lapply(filenames, read.csv, stringsAsFactors = FALSE)
names(timepoint_data_NishaSort) <- filenames

list2env(timepoint_data_NishaSort, envir=.GlobalEnv)


colnames(collated_withCluster_PfalcOrtho_2h.csv)[1] <- "OrderRobynFile"
colnames(collated_withCluster_PfalcOrtho_4h.csv)[1] <- "OrderRobynFile"
colnames(collated_withCluster_PfalcOrtho_6h.csv)[1] <- "OrderRobynFile"
colnames(collated_withCluster_PfalcOrtho_8h.csv)[1] <- "OrderRobynFile"
colnames(collated_withCluster_PfalcOrtho_12h.csv)[1] <- "OrderRobynFile"
colnames(collated_withCluster_PfalcOrtho_24h.csv)[1] <- "OrderRobynFile"
colnames(collated_withCluster_PfalcOrtho_30h.csv)[1] <- "OrderRobynFile"
colnames(collated_withCluster_PfalcOrtho_44h.csv)[1] <- "OrderRobynFile"

#####subset all data frames using NishaSort List IP2#####
collated_withCluster_PfalcOrtho_2h_NishaSubset_IP2.csv <- subset(collated_withCluster_PfalcOrtho_2h.csv, 
                                                             collated_withCluster_PfalcOrtho_2h.csv$GeneId2 %in% NishaGeneList_IP2$GeneID2)
collated_withCluster_PfalcOrtho_4h_NishaSubset_IP2.csv <- subset(collated_withCluster_PfalcOrtho_4h.csv, 
                                                             collated_withCluster_PfalcOrtho_4h.csv$GeneId2 %in% NishaGeneList_IP2$GeneID2)
collated_withCluster_PfalcOrtho_6h_NishaSubset_IP2.csv <- subset(collated_withCluster_PfalcOrtho_6h.csv, 
                                                             collated_withCluster_PfalcOrtho_6h.csv$GeneId2 %in% NishaGeneList_IP2$GeneID2)
collated_withCluster_PfalcOrtho_8h_NishaSubset_IP2.csv <- subset(collated_withCluster_PfalcOrtho_8h.csv, 
                                                             collated_withCluster_PfalcOrtho_8h.csv$GeneId2 %in% NishaGeneList_IP2$GeneID2)
collated_withCluster_PfalcOrtho_12h_NishaSubset_IP2.csv <- subset(collated_withCluster_PfalcOrtho_12h.csv, 
                                                             collated_withCluster_PfalcOrtho_12h.csv$GeneId2 %in% NishaGeneList_IP2$GeneID2)
collated_withCluster_PfalcOrtho_24h_NishaSubset_IP2.csv <- subset(collated_withCluster_PfalcOrtho_24h.csv, 
                                                             collated_withCluster_PfalcOrtho_24h.csv$GeneId2 %in% NishaGeneList_IP2$GeneID2)
collated_withCluster_PfalcOrtho_30h_NishaSubset_IP2.csv <- subset(collated_withCluster_PfalcOrtho_30h.csv, 
                                                             collated_withCluster_PfalcOrtho_30h.csv$GeneId2 %in% NishaGeneList_IP2$GeneID2)
collated_withCluster_PfalcOrtho_44h_NishaSubset_IP2.csv <- subset(collated_withCluster_PfalcOrtho_44h.csv, 
                                                             collated_withCluster_PfalcOrtho_44h.csv$GeneId2 %in% NishaGeneList_IP2$GeneID2)

######write out csvs IP2#####
write.csv(collated_withCluster_PfalcOrtho_2h_NishaSubset_IP2.csv, "collated_withCluster_PfalcOrtho_2h_NishaSubset_IP2.csv")
write.csv(collated_withCluster_PfalcOrtho_4h_NishaSubset_IP2.csv, "collated_withCluster_PfalcOrtho_4h_NishaSubset_IP2.csv")
write.csv(collated_withCluster_PfalcOrtho_6h_NishaSubset_IP2.csv, "collated_withCluster_PfalcOrtho_6h_NishaSubset_IP2.csv")
write.csv(collated_withCluster_PfalcOrtho_8h_NishaSubset_IP2.csv, "collated_withCluster_PfalcOrtho_8h_NishaSubset_IP2.csv")
write.csv(collated_withCluster_PfalcOrtho_12h_NishaSubset_IP2.csv, "collated_withCluster_PfalcOrtho_12h_NishaSubset_IP2.csv")
write.csv(collated_withCluster_PfalcOrtho_24h_NishaSubset_IP2.csv, "collated_withCluster_PfalcOrtho_24h_NishaSubset_IP2.csv")
write.csv(collated_withCluster_PfalcOrtho_30h_NishaSubset_IP2.csv, "collated_withCluster_PfalcOrtho_30h_NishaSubset_IP2.csv")
write.csv(collated_withCluster_PfalcOrtho_44h_NishaSubset_IP2.csv, "collated_withCluster_PfalcOrtho_44h_NishaSubset_IP2.csv")


#####subset all data frames using NishaSort List IP3#####
collated_withCluster_PfalcOrtho_2h_NishaSubset_IP3.csv <- subset(collated_withCluster_PfalcOrtho_2h.csv, 
                                                                 collated_withCluster_PfalcOrtho_2h.csv$GeneId2 %in% NishaGeneList_IP3$GeneID2)
collated_withCluster_PfalcOrtho_4h_NishaSubset_IP3.csv <- subset(collated_withCluster_PfalcOrtho_4h.csv, 
                                                                 collated_withCluster_PfalcOrtho_4h.csv$GeneId2 %in% NishaGeneList_IP3$GeneID2)
collated_withCluster_PfalcOrtho_6h_NishaSubset_IP3.csv <- subset(collated_withCluster_PfalcOrtho_6h.csv, 
                                                                 collated_withCluster_PfalcOrtho_6h.csv$GeneId2 %in% NishaGeneList_IP3$GeneID2)
collated_withCluster_PfalcOrtho_8h_NishaSubset_IP3.csv <- subset(collated_withCluster_PfalcOrtho_8h.csv, 
                                                                 collated_withCluster_PfalcOrtho_8h.csv$GeneId2 %in% NishaGeneList_IP3$GeneID2)
collated_withCluster_PfalcOrtho_12h_NishaSubset_IP3.csv <- subset(collated_withCluster_PfalcOrtho_12h.csv, 
                                                                  collated_withCluster_PfalcOrtho_12h.csv$GeneId2 %in% NishaGeneList_IP3$GeneID2)
collated_withCluster_PfalcOrtho_24h_NishaSubset_IP3.csv <- subset(collated_withCluster_PfalcOrtho_24h.csv, 
                                                                  collated_withCluster_PfalcOrtho_24h.csv$GeneId2 %in% NishaGeneList_IP3$GeneID2)
collated_withCluster_PfalcOrtho_30h_NishaSubset_IP3.csv <- subset(collated_withCluster_PfalcOrtho_30h.csv, 
                                                                  collated_withCluster_PfalcOrtho_30h.csv$GeneId2 %in% NishaGeneList_IP3$GeneID2)
collated_withCluster_PfalcOrtho_44h_NishaSubset_IP3.csv <- subset(collated_withCluster_PfalcOrtho_44h.csv, 
                                                                  collated_withCluster_PfalcOrtho_44h.csv$GeneId2 %in% NishaGeneList_IP3$GeneID2)

######write out csvs IP3#####
write.csv(collated_withCluster_PfalcOrtho_2h_NishaSubset_IP3.csv, "collated_withCluster_PfalcOrtho_2h_NishaSubset_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_4h_NishaSubset_IP3.csv, "collated_withCluster_PfalcOrtho_4h_NishaSubset_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_6h_NishaSubset_IP3.csv, "collated_withCluster_PfalcOrtho_6h_NishaSubset_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_8h_NishaSubset_IP3.csv, "collated_withCluster_PfalcOrtho_8h_NishaSubset_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_12h_NishaSubset_IP3.csv, "collated_withCluster_PfalcOrtho_12h_NishaSubset_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_24h_NishaSubset_IP3.csv, "collated_withCluster_PfalcOrtho_24h_NishaSubset_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_30h_NishaSubset_IP3.csv, "collated_withCluster_PfalcOrtho_30h_NishaSubset_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_44h_NishaSubset_IP3.csv, "collated_withCluster_PfalcOrtho_44h_NishaSubset_IP3.csv")

#####subset all data frames using NishaSort List IP2_IP3#####
collated_withCluster_PfalcOrtho_2h_NishaSubset_IP2_IP3.csv <- subset(collated_withCluster_PfalcOrtho_2h.csv, 
                                                                 collated_withCluster_PfalcOrtho_2h.csv$GeneId2 %in% NishaGeneList_IP2_IP3$GeneID2)
collated_withCluster_PfalcOrtho_4h_NishaSubset_IP2_IP3.csv <- subset(collated_withCluster_PfalcOrtho_4h.csv, 
                                                                 collated_withCluster_PfalcOrtho_4h.csv$GeneId2 %in% NishaGeneList_IP2_IP3$GeneID2)
collated_withCluster_PfalcOrtho_6h_NishaSubset_IP2_IP3.csv <- subset(collated_withCluster_PfalcOrtho_6h.csv, 
                                                                 collated_withCluster_PfalcOrtho_6h.csv$GeneId2 %in% NishaGeneList_IP2_IP3$GeneID2)
collated_withCluster_PfalcOrtho_8h_NishaSubset_IP2_IP3.csv <- subset(collated_withCluster_PfalcOrtho_8h.csv, 
                                                                 collated_withCluster_PfalcOrtho_8h.csv$GeneId2 %in% NishaGeneList_IP2_IP3$GeneID2)
collated_withCluster_PfalcOrtho_12h_NishaSubset_IP2_IP3.csv <- subset(collated_withCluster_PfalcOrtho_12h.csv, 
                                                                  collated_withCluster_PfalcOrtho_12h.csv$GeneId2 %in% NishaGeneList_IP2_IP3$GeneID2)
collated_withCluster_PfalcOrtho_24h_NishaSubset_IP2_IP3.csv <- subset(collated_withCluster_PfalcOrtho_24h.csv, 
                                                                  collated_withCluster_PfalcOrtho_24h.csv$GeneId2 %in% NishaGeneList_IP2_IP3$GeneID2)
collated_withCluster_PfalcOrtho_30h_NishaSubset_IP2_IP3.csv <- subset(collated_withCluster_PfalcOrtho_30h.csv, 
                                                                  collated_withCluster_PfalcOrtho_30h.csv$GeneId2 %in% NishaGeneList_IP2_IP3$GeneID2)
collated_withCluster_PfalcOrtho_44h_NishaSubset_IP2_IP3.csv <- subset(collated_withCluster_PfalcOrtho_44h.csv, 
                                                                  collated_withCluster_PfalcOrtho_44h.csv$GeneId2 %in% NishaGeneList_IP2_IP3$GeneID2)

######write out csvs IP2_IP3#####
write.csv(collated_withCluster_PfalcOrtho_2h_NishaSubset_IP2_IP3.csv, "collated_withCluster_PfalcOrtho_2h_NishaSubset_IP2_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_4h_NishaSubset_IP2_IP3.csv, "collated_withCluster_PfalcOrtho_4h_NishaSubset_IP2_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_6h_NishaSubset_IP2_IP3.csv, "collated_withCluster_PfalcOrtho_6h_NishaSubset_IP2_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_8h_NishaSubset_IP2_IP3.csv, "collated_withCluster_PfalcOrtho_8h_NishaSubset_IP2_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_12h_NishaSubset_IP2_IP3.csv, "collated_withCluster_PfalcOrtho_12h_NishaSubset_IP2_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_24h_NishaSubset_IP2_IP3.csv, "collated_withCluster_PfalcOrtho_24h_NishaSubset_IP2_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_30h_NishaSubset_IP2_IP3.csv, "collated_withCluster_PfalcOrtho_30h_NishaSubset_IP2_IP3.csv")
write.csv(collated_withCluster_PfalcOrtho_44h_NishaSubset_IP2_IP3.csv, "collated_withCluster_PfalcOrtho_44h_NishaSubset_IP2_IP3.csv")





####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/NishaPullingOutGoIRobynData/GoIParsingNisha_RobynData_Workspace.RData")
