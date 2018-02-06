####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/Robyn_ColumnMergingScript_workspace.RData")
#####setwd#####
setwd("J:/III/Waters/Group Members/Mallory/RobynColumnMerging/")

#####load data#####

##csvs
SOM8Clustering <- read.csv("SOM8Clustering.csv")
FcallTimepoints_2 <- read.csv("FCallTimepoints-2.csv")
FcallTimepoints_2_2 <- read.csv("FCallTimepoints-2_2.csv")


filenames <- list.files(pattern = "*h.csv")
filenames <- filenames[order(as.numeric(sub("([0-9]*).*", "\\1", filenames)))]
filenames2 <- c("data_2h.csv","data_4h.csv","data_6h.csv","data_8h.csv","data_12h.csv","data_24h.csv","data_30h.csv","data_44h.csv")

timepoint_data <- lapply(filenames, read.csv, stringsAsFactors = FALSE)
names(timepoint_data) <- filenames2


#####truncate the GeneID to remove the extra "0.1"#####
#first fix the GENEID that says this: "PBANKA_0937200.2%2CPBANKA_0937200.1"
# need to keep in mind that the API genes have longer names! 14 chars in the FCallTimepoints file
  #9 may 2017

lapply(timepoint_data, names)

#first, replicate GeneID column to GeneId2 and GeneI3
colnames(timepoint_data$data_8h.csv) <- colnames(timepoint_data$data_4h.csv)
colnames(timepoint_data$data_12h.csv) <- colnames(timepoint_data$data_4h.csv)

timepoint_data <- lapply(timepoint_data, function(x) {
  x$GeneID_OldDescMatch = x$GeneId
  x$GeneID_PfalcMatch = x$GeneId
  x
  })

##split data frames into individual
list2env(timepoint_data, envir=.GlobalEnv)


#to truncate the GeneID, but separately for the API and MITO genes
###need to search for the API rows, do those first with 15 chars, then everything else with 13
###the timepoints 2 merger file only has up to 14 chars in the API genes, not sure how that will affect the matches
####i will need to make several geneID lists
###GeneID_OldDescMatch, GeneID_PfalcMatch
###the SOM8 matching can just use the GeneID
###GeneID_OldDescMatch --> 13 chars long, 14 for API and MITO genes, MITO genes actually won't match at all and API will have duplicates --> this will amtch with FcallTimepoints
####GeneID_PfalcMatch --> 14 for regular, 15 for API and MITO --> this will match with the output of the ortholog list

#####geneID replacements#####

list <- c("API","MIT")

for (j in seq_along(data_2h.csv$GeneID_OldDescMatch)){
  ifelse(grepl("API",data_2h.csv$GeneID_OldDescMatch[j]), 
         data_2h.csv$GeneID_OldDescMatch[j] <- substr(data_2h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
         ifelse(grepl("MIT",data_2h.csv$GeneID_OldDescMatch[j]),
                data_2h.csv$GeneID_OldDescMatch[j] <- substr(data_2h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
                data_2h.csv$GeneID_OldDescMatch[j] <- substr(data_2h.csv$GeneID_OldDescMatch[j],start = 0,stop = 13)))
}
      for (j in seq_along(data_2h.csv$GeneID_PfalcMatch)){
        ifelse(grepl("API",data_2h.csv$GeneID_PfalcMatch[j]), 
               data_2h.csv$GeneID_PfalcMatch[j] <- substr(data_2h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
               ifelse(grepl("MIT",data_2h.csv$GeneID_PfalcMatch[j]),
                      data_2h.csv$GeneID_PfalcMatch[j] <- substr(data_2h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
                      data_2h.csv$GeneID_PfalcMatch[j] <- substr(data_2h.csv$GeneID_PfalcMatch[j],start = 0,stop = 14)))
      }
for (j in seq_along(data_4h.csv$GeneID_OldDescMatch)){
  ifelse(grepl("API",data_4h.csv$GeneID_OldDescMatch[j]), 
         data_4h.csv$GeneID_OldDescMatch[j] <- substr(data_4h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
         ifelse(grepl("MIT",data_4h.csv$GeneID_OldDescMatch[j]),
                data_4h.csv$GeneID_OldDescMatch[j] <- substr(data_4h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
                data_4h.csv$GeneID_OldDescMatch[j] <- substr(data_4h.csv$GeneID_OldDescMatch[j],start = 0,stop = 13)))
}
      for (j in seq_along(data_4h.csv$GeneID_PfalcMatch)){
        ifelse(grepl("API",data_4h.csv$GeneID_PfalcMatch[j]), 
               data_4h.csv$GeneID_PfalcMatch[j] <- substr(data_4h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
               ifelse(grepl("MIT",data_4h.csv$GeneID_PfalcMatch[j]),
                      data_4h.csv$GeneID_PfalcMatch[j] <- substr(data_4h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
                      data_4h.csv$GeneID_PfalcMatch[j] <- substr(data_4h.csv$GeneID_PfalcMatch[j],start = 0,stop = 14)))
      }
for (j in seq_along(data_6h.csv$GeneID_OldDescMatch)){
  ifelse(grepl("API",data_6h.csv$GeneID_OldDescMatch[j]), 
         data_6h.csv$GeneID_OldDescMatch[j] <- substr(data_6h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
         ifelse(grepl("MIT",data_6h.csv$GeneID_OldDescMatch[j]),
                data_6h.csv$GeneID_OldDescMatch[j] <- substr(data_6h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
                data_6h.csv$GeneID_OldDescMatch[j] <- substr(data_6h.csv$GeneID_OldDescMatch[j],start = 0,stop = 13)))
}
      for (j in seq_along(data_6h.csv$GeneID_PfalcMatch)){
        ifelse(grepl("API",data_6h.csv$GeneID_PfalcMatch[j]), 
               data_6h.csv$GeneID_PfalcMatch[j] <- substr(data_6h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
               ifelse(grepl("MIT",data_6h.csv$GeneID_PfalcMatch[j]),
                      data_6h.csv$GeneID_PfalcMatch[j] <- substr(data_6h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
                      data_6h.csv$GeneID_PfalcMatch[j] <- substr(data_6h.csv$GeneID_PfalcMatch[j],start = 0,stop = 14)))
      }
for (j in seq_along(data_8h.csv$GeneID_OldDescMatch)){
  ifelse(grepl("API",data_8h.csv$GeneID_OldDescMatch[j]), 
         data_8h.csv$GeneID_OldDescMatch[j] <- substr(data_8h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
         ifelse(grepl("MIT",data_8h.csv$GeneID_OldDescMatch[j]),
                data_8h.csv$GeneID_OldDescMatch[j] <- substr(data_8h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
                data_8h.csv$GeneID_OldDescMatch[j] <- substr(data_8h.csv$GeneID_OldDescMatch[j],start = 0,stop = 13)))
}
      for (j in seq_along(data_8h.csv$GeneID_PfalcMatch)){
        ifelse(grepl("API",data_8h.csv$GeneID_PfalcMatch[j]), 
               data_8h.csv$GeneID_PfalcMatch[j] <- substr(data_8h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
               ifelse(grepl("MIT",data_8h.csv$GeneID_PfalcMatch[j]),
                      data_8h.csv$GeneID_PfalcMatch[j] <- substr(data_8h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
                      data_8h.csv$GeneID_PfalcMatch[j] <- substr(data_8h.csv$GeneID_PfalcMatch[j],start = 0,stop = 14)))
      }
for (j in seq_along(data_12h.csv$GeneID_OldDescMatch)){
  ifelse(grepl("API",data_12h.csv$GeneID_OldDescMatch[j]), 
         data_12h.csv$GeneID_OldDescMatch[j] <- substr(data_12h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
         ifelse(grepl("MIT",data_12h.csv$GeneID_OldDescMatch[j]),
                data_12h.csv$GeneID_OldDescMatch[j] <- substr(data_12h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
                data_12h.csv$GeneID_OldDescMatch[j] <- substr(data_12h.csv$GeneID_OldDescMatch[j],start = 0,stop = 13)))
}
      for (j in seq_along(data_12h.csv$GeneID_PfalcMatch)){
        ifelse(grepl("API",data_12h.csv$GeneID_PfalcMatch[j]), 
               data_12h.csv$GeneID_PfalcMatch[j] <- substr(data_12h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
               ifelse(grepl("MIT",data_12h.csv$GeneID_PfalcMatch[j]),
                      data_12h.csv$GeneID_PfalcMatch[j] <- substr(data_12h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
                      data_12h.csv$GeneID_PfalcMatch[j] <- substr(data_12h.csv$GeneID_PfalcMatch[j],start = 0,stop = 14)))
      }
for (j in seq_along(data_24h.csv$GeneID_OldDescMatch)){
  ifelse(grepl("API",data_24h.csv$GeneID_OldDescMatch[j]), 
         data_24h.csv$GeneID_OldDescMatch[j] <- substr(data_24h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
         ifelse(grepl("MIT",data_24h.csv$GeneID_OldDescMatch[j]),
                data_24h.csv$GeneID_OldDescMatch[j] <- substr(data_24h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
                data_24h.csv$GeneID_OldDescMatch[j] <- substr(data_24h.csv$GeneID_OldDescMatch[j],start = 0,stop = 13)))
}
      for (j in seq_along(data_24h.csv$GeneID_PfalcMatch)){
        ifelse(grepl("API",data_24h.csv$GeneID_PfalcMatch[j]), 
               data_24h.csv$GeneID_PfalcMatch[j] <- substr(data_24h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
               ifelse(grepl("MIT",data_24h.csv$GeneID_PfalcMatch[j]),
                      data_24h.csv$GeneID_PfalcMatch[j] <- substr(data_24h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
                      data_24h.csv$GeneID_PfalcMatch[j] <- substr(data_24h.csv$GeneID_PfalcMatch[j],start = 0,stop = 14)))
      }
for (j in seq_along(data_30h.csv$GeneID_OldDescMatch)){
  ifelse(grepl("API",data_30h.csv$GeneID_OldDescMatch[j]), 
         data_30h.csv$GeneID_OldDescMatch[j] <- substr(data_30h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
         ifelse(grepl("MIT",data_30h.csv$GeneID_OldDescMatch[j]),
                data_30h.csv$GeneID_OldDescMatch[j] <- substr(data_30h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
                data_30h.csv$GeneID_OldDescMatch[j] <- substr(data_30h.csv$GeneID_OldDescMatch[j],start = 0,stop = 13)))
}
      for (j in seq_along(data_30h.csv$GeneID_PfalcMatch)){
        ifelse(grepl("API",data_30h.csv$GeneID_PfalcMatch[j]), 
               data_30h.csv$GeneID_PfalcMatch[j] <- substr(data_30h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
               ifelse(grepl("MIT",data_30h.csv$GeneID_PfalcMatch[j]),
                      data_30h.csv$GeneID_PfalcMatch[j] <- substr(data_30h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
                      data_30h.csv$GeneID_PfalcMatch[j] <- substr(data_30h.csv$GeneID_PfalcMatch[j],start = 0,stop = 14)))
      }
for (j in seq_along(data_44h.csv$GeneID_OldDescMatch)){
  ifelse(grepl("API",data_44h.csv$GeneID_OldDescMatch[j]), 
         data_44h.csv$GeneID_OldDescMatch[j] <- substr(data_44h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
         ifelse(grepl("MIT",data_44h.csv$GeneID_OldDescMatch[j]),
                data_44h.csv$GeneID_OldDescMatch[j] <- substr(data_44h.csv$GeneID_OldDescMatch[j],start = 0,stop = 14),
                data_44h.csv$GeneID_OldDescMatch[j] <- substr(data_44h.csv$GeneID_OldDescMatch[j],start = 0,stop = 13)))
}
      for (j in seq_along(data_44h.csv$GeneID_PfalcMatch)){
        ifelse(grepl("API",data_44h.csv$GeneID_PfalcMatch[j]), 
               data_44h.csv$GeneID_PfalcMatch[j] <- substr(data_44h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
               ifelse(grepl("MIT",data_44h.csv$GeneID_PfalcMatch[j]),
                      data_44h.csv$GeneID_PfalcMatch[j] <- substr(data_44h.csv$GeneID_PfalcMatch[j],start = 0,stop = 15),
                      data_44h.csv$GeneID_PfalcMatch[j] <- substr(data_44h.csv$GeneID_PfalcMatch[j],start = 0,stop = 14)))
      }

######merge in the SOM8 cluster using the GeneId#####
SOM8Clustering_merger <- SOM8Clustering[,c(2,14)]
  colnames(SOM8Clustering_merger) <- c("GeneId","SOM8Cluster")

data_2h_cluster.csv <- merge(data_2h.csv, SOM8Clustering_merger, by = "GeneId", all.x=TRUE)
  data_2h_cluster.csv <- data_2h_cluster.csv[,c(2,1,3:14)]
data_4h_cluster.csv <- merge(data_4h.csv, SOM8Clustering_merger, by = "GeneId", all.x=TRUE)
  data_4h_cluster.csv <- data_4h_cluster.csv[,c(2,1,3:12)]
data_6h_cluster.csv <- merge(data_6h.csv, SOM8Clustering_merger, by = "GeneId", all.x=TRUE)
  data_6h_cluster.csv <- data_6h_cluster.csv[,c(2,1,3:12)]
data_8h_cluster.csv <- merge(data_8h.csv, SOM8Clustering_merger, by = "GeneId", all.x=TRUE)
  data_8h_cluster.csv <- data_8h_cluster.csv[,c(2,1,3:12)]
data_12h_cluster.csv <- merge(data_12h.csv, SOM8Clustering_merger, by = "GeneId", all.x=TRUE)
  data_12h_cluster.csv <- data_12h_cluster.csv[,c(2,1,3:12)]
data_24h_cluster.csv <- merge(data_24h.csv, SOM8Clustering_merger, by = "GeneId", all.x=TRUE)
  data_24h_cluster.csv <- data_24h_cluster.csv[,c(2,1,3:12)]
data_30h_cluster.csv <- merge(data_30h.csv, SOM8Clustering_merger, by = "GeneId", all.x=TRUE)
  data_30h_cluster.csv <- data_30h_cluster.csv[,c(2,1,3:12)]
data_44h_cluster.csv <- merge(data_44h.csv, SOM8Clustering_merger, by = "GeneId", all.x=TRUE)
  data_44h_cluster.csv <- data_44h_cluster.csv[,c(2,1,3:12)]

#####figure out P falciparum orthologs#####
#take list of all 5263 genes, download the orthologs from PlasmoDB
###UPDATE 9 May 2017 - use GeneID, NOT GeneID2
# search foudn 5242 genes, 4598 Pflac genes
#load ortholog list (processed in excel) --> downloaded 9 May 2017
Pfalc_orthologGeneIds <- read.csv("GenesByOrthologs_Summary_try2.csv", stringsAsFactors = FALSE)

#merge the PfalciparumGeneID into the _merge_cluster dataframe by GeneId2 and PbergheiGeneID
#seems to automatically handle the duplicates...duplicates the rows where a Pberghei GeneID has multiple Pfalc orthologs and changes lists the different orthologs
#data_2h_cluster.csv$GeneID_PfalcMatch2 <- data_2h_cluster.csv$GeneID_PfalcMatch

data_2h_cluster_ortho.csv <- merge(data_2h_cluster.csv, Pfalc_orthologGeneIds, by.x = "GeneID_PfalcMatch", by.y = "PbergheiGeneID", all = TRUE)
  data_2h_cluster_ortho.csv <- data_2h_cluster_ortho.csv[,c(2:14,1,15:17)]
data_4h_cluster_ortho.csv <- merge(data_4h_cluster.csv, Pfalc_orthologGeneIds, by.x = "GeneID_PfalcMatch", by.y = "PbergheiGeneID", all = TRUE)
  data_4h_cluster_ortho.csv <- data_4h_cluster_ortho.csv[,c(2:12,1,13:15)]
data_6h_cluster_ortho.csv <- merge(data_6h_cluster.csv, Pfalc_orthologGeneIds, by.x = "GeneID_PfalcMatch", by.y = "PbergheiGeneID", all = TRUE)
  data_6h_cluster_ortho.csv <- data_6h_cluster_ortho.csv[,c(2:12,1,13:15)] 
data_8h_cluster_ortho.csv <- merge(data_8h_cluster.csv, Pfalc_orthologGeneIds, by.x = "GeneID_PfalcMatch", by.y = "PbergheiGeneID", all = TRUE)
  data_8h_cluster_ortho.csv <- data_8h_cluster_ortho.csv[,c(2:12,1,13:15)]
data_12h_cluster_ortho.csv <- merge(data_12h_cluster.csv, Pfalc_orthologGeneIds, by.x = "GeneID_PfalcMatch", by.y = "PbergheiGeneID", all = TRUE)
  data_12h_cluster_ortho.csv <- data_12h_cluster_ortho.csv[,c(2:12,1,13:15)]
data_24h_cluster_ortho.csv <- merge(data_24h_cluster.csv, Pfalc_orthologGeneIds, by.x = "GeneID_PfalcMatch", by.y = "PbergheiGeneID", all = TRUE)
  data_24h_cluster_ortho.csv <- data_24h_cluster_ortho.csv[,c(2:12,1,13:15)]  
data_30h_cluster_ortho.csv <- merge(data_30h_cluster.csv, Pfalc_orthologGeneIds, by.x = "GeneID_PfalcMatch", by.y = "PbergheiGeneID", all = TRUE)
  data_30h_cluster_ortho.csv <- data_30h_cluster_ortho.csv[,c(2:12,1,13:15)]
data_44h_cluster_ortho.csv <- merge(data_44h_cluster.csv, Pfalc_orthologGeneIds, by.x = "GeneID_PfalcMatch", by.y = "PbergheiGeneID", all = TRUE)
  data_44h_cluster_ortho.csv <- data_44h_cluster_ortho.csv[,c(2:12,1,13:15)]
  
######write out final merged csvs with orthos#####
write.csv(data_2h_cluster_ortho.csv, "data_2h_cluster_ortho.csv")
write.csv(data_4h_cluster_ortho.csv, "data_4h_cluster_ortho.csv")
write.csv(data_6h_cluster_ortho.csv, "data_6h_cluster_ortho.csv")
write.csv(data_8h_cluster_ortho.csv, "data_8h_cluster_ortho.csv")
write.csv(data_12h_cluster_ortho.csv, "data_12h_cluster_ortho.csv")
write.csv(data_24h_cluster_ortho.csv, "data_24h_cluster_ortho.csv")
write.csv(data_30h_cluster_ortho.csv, "data_30h_cluster_ortho.csv")
write.csv(data_44h_cluster_ortho.csv, "data_44h_cluster_ortho.csv")  
  
######merge FCallTimepoints-2 description into the timepoints data file#####
#remove extra header from FcallTimepoints_2
FcallTimepoints_2 <- FcallTimepoints_2[2:5127,]
FcallTimepoints_2_merger <- FcallTimepoints_2[,1:2]
colnames(FcallTimepoints_2_merger) <- c("GeneId","OldProductDescription")


GeneID_OldDescMatch

data_2h_cluster_ortho_desc.csv <- merge(data_2h_cluster_ortho.csv, FcallTimepoints_2_merger, by.x = "GeneID_OldDescMatch", by.y = "GeneId", all.x = TRUE)
  data_2h_cluster_ortho_desc.csv <- data_2h_cluster_ortho_desc.csv[,c(2:17,1,18)]
data_4h_cluster_ortho_desc.csv <- merge(data_4h_cluster_ortho.csv, FcallTimepoints_2_merger, by.x = "GeneID_OldDescMatch", by.y = "GeneId", all.x = TRUE)
  data_4h_cluster_ortho_desc.csv <- data_4h_cluster_ortho_desc.csv[,c(2:15,1,16)]
data_6h_cluster_ortho_desc.csv <- merge(data_6h_cluster_ortho.csv, FcallTimepoints_2_merger, by.x = "GeneID_OldDescMatch", by.y = "GeneId", all.x = TRUE)
  data_6h_cluster_ortho_desc.csv <- data_6h_cluster_ortho_desc.csv[,c(2:15,1,16)]
data_8h_cluster_ortho_desc.csv <- merge(data_8h_cluster_ortho.csv, FcallTimepoints_2_merger, by.x = "GeneID_OldDescMatch", by.y = "GeneId", all.x = TRUE)
  data_8h_cluster_ortho_desc.csv <- data_8h_cluster_ortho_desc.csv[,c(2:15,1,16)]
data_12h_cluster_ortho_desc.csv <- merge(data_12h_cluster_ortho.csv, FcallTimepoints_2_merger, by.x = "GeneID_OldDescMatch", by.y = "GeneId", all.x = TRUE)
  data_12h_cluster_ortho_desc.csv <- data_12h_cluster_ortho_desc.csv[,c(2:15,1,16)]
data_24h_cluster_ortho_desc.csv <- merge(data_24h_cluster_ortho.csv, FcallTimepoints_2_merger, by.x = "GeneID_OldDescMatch", by.y = "GeneId", all.x = TRUE)
  data_24h_cluster_ortho_desc.csv <- data_24h_cluster_ortho_desc.csv[,c(2:15,1,16)]
data_30h_cluster_ortho_desc.csv <- merge(data_30h_cluster_ortho.csv, FcallTimepoints_2_merger, by.x = "GeneID_OldDescMatch", by.y = "GeneId", all.x = TRUE)
  data_30h_cluster_ortho_desc.csv <- data_30h_cluster_ortho_desc.csv[,c(2:15,1,16)]
data_44h_cluster_ortho_desc.csv <- merge(data_44h_cluster_ortho.csv, FcallTimepoints_2_merger, by.x = "GeneID_OldDescMatch", by.y = "GeneId", all.x = TRUE)
  data_44h_cluster_ortho_desc.csv <- data_44h_cluster_ortho_desc.csv[,c(2:15,1,16)]
  
  
######write out final merged csvs with orthos amd old desc######  
write.csv(data_2h_cluster_ortho_desc.csv, "data_2h_cluster_ortho_desc.csv")
write.csv(data_4h_cluster_ortho_desc.csv, "data_4h_cluster_ortho_desc.csv")
write.csv(data_6h_cluster_ortho_desc.csv, "data_6h_cluster_ortho_desc.csv")
write.csv(data_8h_cluster_ortho_desc.csv, "data_8h_cluster_ortho_desc.csv")
write.csv(data_12h_cluster_ortho_desc.csv, "data_12h_cluster_ortho_desc.csv")
write.csv(data_24h_cluster_ortho_desc.csv, "data_24h_cluster_ortho_desc.csv")
write.csv(data_30h_cluster_ortho_desc.csv, "data_30h_cluster_ortho_desc.csv")
write.csv(data_44h_cluster_ortho_desc.csv, "data_44h_cluster_ortho_desc.csv")   

  




####save workspace####
save.image("J:/III/Waters/Group Members/Mallory/R/Robyn_ColumnMergingScript_workspace.RData")
