####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/GoIParsingNisha_RobynData_2_Workspace.RData")
#####setwd#####
setwd("J:/III/Waters/Group Members/Mallory/NishaPullingOutGoIRobynData2/")

#####load data#####
NishaGeneList_RIP <- read.csv("RIP_enriched.csv")
#remove everything after a period in the GeneID (should just affect one GeneID)
NishaGeneList_RIP$GeneID2 <- gsub("\\..*",'',NishaGeneList_RIP$GeneID)
colnames(NishaGeneList_RIP) <- c("NishaGeneID","NishaGeneID2")


filenames <- c("data_2h_cluster_ortho_desc.csv",
               "data_4h_cluster_ortho_desc.csv",
               "data_6h_cluster_ortho_desc.csv",
               "data_8h_cluster_ortho_desc.csv",
               "data_12h_cluster_ortho_desc.csv",
               "data_24h_cluster_ortho_desc.csv",
               "data_30h_cluster_ortho_desc.csv",
               "data_44h_cluster_ortho_desc.csv")


timepoint_data_NishaSort <- lapply(filenames, read.csv, stringsAsFactors = FALSE)
names(timepoint_data_NishaSort) <- filenames

list2env(timepoint_data_NishaSort, envir=.GlobalEnv)


colnames(data_2h_cluster_ortho_desc.csv)[1] <- "OrderRobynFile"
colnames(data_4h_cluster_ortho_desc.csv)[1] <- "OrderRobynFile"
colnames(data_6h_cluster_ortho_desc.csv)[1] <- "OrderRobynFile"
colnames(data_8h_cluster_ortho_desc.csv)[1] <- "OrderRobynFile"
colnames(data_12h_cluster_ortho_desc.csv)[1] <- "OrderRobynFile"
colnames(data_24h_cluster_ortho_desc.csv)[1] <- "OrderRobynFile"
colnames(data_30h_cluster_ortho_desc.csv)[1] <- "OrderRobynFile"
colnames(data_44h_cluster_ortho_desc.csv)[1] <- "OrderRobynFile"


#####subset all data frames using NishaSort List RIP#####
data_2h_cluster_ortho_desc_NishaSubset_RIP.csv <- subset(data_2h_cluster_ortho_desc.csv, 
                                                                 data_2h_cluster_ortho_desc.csv$GeneID_PfalcMatch %in% NishaGeneList_RIP$NishaGeneID2)
  data_2h_cluster_ortho_desc_NishaSubset_RIP.csv <- merge(data_2h_cluster_ortho_desc_NishaSubset_RIP.csv, NishaGeneList_RIP, 
                by.x = "GeneID_PfalcMatch", 
                by.y = "NishaGeneID2", 
                all.x = T)
    data_2h_cluster_ortho_desc_NishaSubset_RIP.csv$NishaGeneID2 <- data_2h_cluster_ortho_desc_NishaSubset_RIP.csv$GeneID_PfalcMatch
    data_2h_cluster_ortho_desc_NishaSubset_RIP.csv <- data_2h_cluster_ortho_desc_NishaSubset_RIP.csv[,c(2:13,1,14:20)]

data_4h_cluster_ortho_desc_NishaSubset_RIP.csv <- subset(data_4h_cluster_ortho_desc.csv, 
                                                             data_4h_cluster_ortho_desc.csv$GeneID_PfalcMatch %in% NishaGeneList_RIP$NishaGeneID2)
    data_4h_cluster_ortho_desc_NishaSubset_RIP.csv <- merge(data_4h_cluster_ortho_desc_NishaSubset_RIP.csv, NishaGeneList_RIP, 
                                                            by.x = "GeneID_PfalcMatch", 
                                                            by.y = "NishaGeneID2", 
                                                            all.x = T)
    data_4h_cluster_ortho_desc_NishaSubset_RIP.csv$NishaGeneID2 <- data_4h_cluster_ortho_desc_NishaSubset_RIP.csv$GeneID_PfalcMatch
    data_4h_cluster_ortho_desc_NishaSubset_RIP.csv <- data_4h_cluster_ortho_desc_NishaSubset_RIP.csv[,c(2:11,1,12:18)]

data_6h_cluster_ortho_desc_NishaSubset_RIP.csv <- subset(data_6h_cluster_ortho_desc.csv, 
                                                             data_6h_cluster_ortho_desc.csv$GeneID_PfalcMatch %in% NishaGeneList_RIP$NishaGeneID2)
    data_6h_cluster_ortho_desc_NishaSubset_RIP.csv <- merge(data_6h_cluster_ortho_desc_NishaSubset_RIP.csv, NishaGeneList_RIP, 
                                                            by.x = "GeneID_PfalcMatch", 
                                                            by.y = "NishaGeneID2", 
                                                            all.x = T)
    data_6h_cluster_ortho_desc_NishaSubset_RIP.csv$NishaGeneID2 <- data_6h_cluster_ortho_desc_NishaSubset_RIP.csv$GeneID_PfalcMatch
    data_6h_cluster_ortho_desc_NishaSubset_RIP.csv <- data_6h_cluster_ortho_desc_NishaSubset_RIP.csv[,c(2:11,1,12:18)]

data_8h_cluster_ortho_desc_NishaSubset_RIP.csv <- subset(data_8h_cluster_ortho_desc.csv, 
                                                             data_8h_cluster_ortho_desc.csv$GeneID_PfalcMatch %in% NishaGeneList_RIP$NishaGeneID2)
    data_8h_cluster_ortho_desc_NishaSubset_RIP.csv <- merge(data_8h_cluster_ortho_desc_NishaSubset_RIP.csv, NishaGeneList_RIP, 
                                                            by.x = "GeneID_PfalcMatch", 
                                                            by.y = "NishaGeneID2", 
                                                            all.x = T)
    data_8h_cluster_ortho_desc_NishaSubset_RIP.csv$NishaGeneID2 <- data_8h_cluster_ortho_desc_NishaSubset_RIP.csv$GeneID_PfalcMatch
    data_8h_cluster_ortho_desc_NishaSubset_RIP.csv <- data_8h_cluster_ortho_desc_NishaSubset_RIP.csv[,c(2:11,1,12:18)]
    
data_12h_cluster_ortho_desc_NishaSubset_RIP.csv <- subset(data_12h_cluster_ortho_desc.csv, 
                                                             data_12h_cluster_ortho_desc.csv$GeneID_PfalcMatch %in% NishaGeneList_RIP$NishaGeneID2)
    data_12h_cluster_ortho_desc_NishaSubset_RIP.csv <- merge(data_12h_cluster_ortho_desc_NishaSubset_RIP.csv, NishaGeneList_RIP, 
                                                            by.x = "GeneID_PfalcMatch", 
                                                            by.y = "NishaGeneID2", 
                                                            all.x = T)
    data_12h_cluster_ortho_desc_NishaSubset_RIP.csv$NishaGeneID2 <- data_12h_cluster_ortho_desc_NishaSubset_RIP.csv$GeneID_PfalcMatch
    data_12h_cluster_ortho_desc_NishaSubset_RIP.csv <- data_12h_cluster_ortho_desc_NishaSubset_RIP.csv[,c(2:11,1,12:18)]
    
data_24h_cluster_ortho_desc_NishaSubset_RIP.csv <- subset(data_24h_cluster_ortho_desc.csv, 
                                                             data_24h_cluster_ortho_desc.csv$GeneID_PfalcMatch %in% NishaGeneList_RIP$NishaGeneID2)
    data_24h_cluster_ortho_desc_NishaSubset_RIP.csv <- merge(data_24h_cluster_ortho_desc_NishaSubset_RIP.csv, NishaGeneList_RIP, 
                                                            by.x = "GeneID_PfalcMatch", 
                                                            by.y = "NishaGeneID2", 
                                                            all.x = T)
    data_24h_cluster_ortho_desc_NishaSubset_RIP.csv$NishaGeneID2 <- data_24h_cluster_ortho_desc_NishaSubset_RIP.csv$GeneID_PfalcMatch
    data_24h_cluster_ortho_desc_NishaSubset_RIP.csv <- data_24h_cluster_ortho_desc_NishaSubset_RIP.csv[,c(2:11,1,12:18)]
    
data_30h_cluster_ortho_desc_NishaSubset_RIP.csv <- subset(data_30h_cluster_ortho_desc.csv, 
                                                             data_30h_cluster_ortho_desc.csv$GeneID_PfalcMatch %in% NishaGeneList_RIP$NishaGeneID2)
    data_30h_cluster_ortho_desc_NishaSubset_RIP.csv <- merge(data_30h_cluster_ortho_desc_NishaSubset_RIP.csv, NishaGeneList_RIP, 
                                                            by.x = "GeneID_PfalcMatch", 
                                                            by.y = "NishaGeneID2", 
                                                            all.x = T)
    data_30h_cluster_ortho_desc_NishaSubset_RIP.csv$NishaGeneID2 <- data_30h_cluster_ortho_desc_NishaSubset_RIP.csv$GeneID_PfalcMatch
    data_30h_cluster_ortho_desc_NishaSubset_RIP.csv <- data_30h_cluster_ortho_desc_NishaSubset_RIP.csv[,c(2:11,1,12:18)]
    
data_44h_cluster_ortho_desc_NishaSubset_RIP.csv <- subset(data_44h_cluster_ortho_desc.csv, 
                                                             data_44h_cluster_ortho_desc.csv$GeneID_PfalcMatch %in% NishaGeneList_RIP$NishaGeneID2)
    data_44h_cluster_ortho_desc_NishaSubset_RIP.csv <- merge(data_44h_cluster_ortho_desc_NishaSubset_RIP.csv, NishaGeneList_RIP, 
                                                            by.x = "GeneID_PfalcMatch", 
                                                            by.y = "NishaGeneID2", 
                                                            all.x = T)
    data_44h_cluster_ortho_desc_NishaSubset_RIP.csv$NishaGeneID2 <- data_44h_cluster_ortho_desc_NishaSubset_RIP.csv$GeneID_PfalcMatch
    data_44h_cluster_ortho_desc_NishaSubset_RIP.csv <- data_44h_cluster_ortho_desc_NishaSubset_RIP.csv[,c(2:11,1,12:18)]    



######write out csvs IP2#####
write.csv(data_2h_cluster_ortho_desc_NishaSubset_RIP.csv, "data_2h_cluster_ortho_desc_NishaSubset_RIP.csv")
write.csv(data_4h_cluster_ortho_desc_NishaSubset_RIP.csv, "data_4h_cluster_ortho_desc_NishaSubset_RIP.csv")
write.csv(data_6h_cluster_ortho_desc_NishaSubset_RIP.csv, "data_6h_cluster_ortho_desc_NishaSubset_RIP.csv")
write.csv(data_8h_cluster_ortho_desc_NishaSubset_RIP.csv, "data_8h_cluster_ortho_desc_NishaSubset_RIP.csv")
write.csv(data_12h_cluster_ortho_desc_NishaSubset_RIP.csv, "data_12h_cluster_ortho_desc_NishaSubset_RIP.csv")
write.csv(data_24h_cluster_ortho_desc_NishaSubset_RIP.csv, "data_24h_cluster_ortho_desc_NishaSubset_RIP.csv")
write.csv(data_30h_cluster_ortho_desc_NishaSubset_RIP.csv, "data_30h_cluster_ortho_desc_NishaSubset_RIP.csv")
write.csv(data_44h_cluster_ortho_desc_NishaSubset_RIP.csv, "data_44h_cluster_ortho_desc_NishaSubset_RIP.csv")


####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/GoIParsingNisha_RobynData_2_Workspace.RData")