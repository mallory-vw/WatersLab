#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/LysoPC_genes_with_motifs.RData")

#####setwd#####
setwd("J:/III/Waters/Group Members/Mallory/LysoPC_project/GenesWithMotifs/FIRE/GeneListsFromMotifs")


#####load data#####
Motif1 <- read.csv("Motif1_AGTCTTAAA.csv", header=F)
  colnames(Motif1) <- "GENEID"
Motif2 <- read.csv("Motif2_ATAATAACG.csv", header=F)
  colnames(Motif2) <- "GENEID"
Motif3 <- read.csv("Motif3_CACCATCTACACT.csv", header = F)                               
  colnames(Motif3) <- "GENEID"
Motif4 <- read.csv("Motif4_CACACACTAGCACT.csv", header = F)
  colnames(Motif4) <- "GENEID"
Motif5 <- read.csv("Motif5_ACCTAGTAAC.csv", header = F)
  colnames(Motif5) <- "GENEID"
Motif6 <- read.csv("Motif6_ACGACAAACACTGTAT.csv", header = F)
  colnames(Motif6) <- "GENEID"
Motif7 <- read.csv("Motif7_ACTAAGGTCGAT.csv", header = F)
  colnames(Motif7) <- "GENEID"
Motif8 <- read.csv("Motif8_ACTCGTCTACTAACT.csv", header = F)
  colnames(Motif8) <- "GENEID"
Motif9 <- read.csv("Motif9_ACTCGCTACT.csv", header = F)
  colnames(Motif9) <- "GENEID"
Motif10 <- read.csv("Motif10_AGTAGTTACGACGAGAGT.csv", header = F)
  colnames(Motif10) <- "GENEID"
Motif11 <- read.csv("Motif11_AGTCTACGGACT.csv", header = F)
  colnames(Motif11) <- "GENEID"
Motif12 <- read.csv("Motif12_AGTCACACTCTACTCAC.csv", header = F)
  colnames(Motif12) <- "GENEID"
Motif13 <- read.csv("Motif13_AGTCACACTCTACACTGT.csv", header = F)
  colnames(Motif13) <- "GENEID"
Motif14 <- read.csv("Motif14_ATACACATAT.csv", header = F)
  colnames(Motif14) <- "GENEID"
Motif15 <- read.csv("Motif15_ATCCGTAG.csv", header = F)
  colnames(Motif15) <- "GENEID"  
Motif16 <- read.csv("Motif16_ATCCGTCTACGTAACT.csv", header = F)
  colnames(Motif16) <- "GENEID"  
Motif17 <- read.csv("Motif17_ATTACTACGAGAAGT.csv", header = F)
  colnames(Motif17) <- "GENEID"  
  
Motifs <- read.csv("J:/III/Waters/Group Members/Mallory/LysoPC_Project/FIRE_results/FoundMotifs.csv", header=F)
  colnames(Motifs) <- c("Motif","RegExp")

background_noDownreg <- read.csv("J:/III/Waters/Group Members/Mallory/LysoPC_Project/GeneLists/pos_background_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  background_noDownreg <- as.factor(background_noDownreg[,2])
downreg <- read.csv("J:/III/Waters/Group Members/Mallory/LysoPC_Project/GeneLists/downreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  downreg <- as.factor(downreg[,2])
upreg <- read.csv("J:/III/Waters/Group Members/Mallory/LysoPC_Project/GeneLists/upreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  upreg <- as.factor(upreg[,2])

#####add gene type to motif lists####
Motif1$Category <- NA
  for (i in seq_along(Motif1$GENEID)) {
    if (Motif1$GENEID[[i]] %in% upreg) {
      Motif1$Category[i] <- "upreg"
    } else if (Motif1$GENEID[i] %in% downreg) {
      Motif1$Category[i] <- "downreg"
      } else if (Motif1$GENEID[i] %in% background_noDownreg) {
        Motif1$Category[i] <- "background"
        } else {
          Motif1$Category[i] <- NA        
        }
  }

Motif2$Category <- NA
  for (i in seq_along(Motif2$GENEID)) {
    if (Motif2$GENEID[[i]] %in% upreg) {
      Motif2$Category[i] <- "upreg"
    } else if (Motif2$GENEID[i] %in% downreg) {
      Motif2$Category[i] <- "downreg"
    } else if (Motif2$GENEID[i] %in% background_noDownreg) {
      Motif2$Category[i] <- "background"
    } else {
      Motif2$Category[i] <- NA        
    }
  }

Motif3$Category <- NA
for (i in seq_along(Motif3$GENEID)) {
  if (Motif3$GENEID[[i]] %in% upreg) {
    Motif3$Category[i] <- "upreg"
  } else if (Motif3$GENEID[i] %in% downreg) {
    Motif3$Category[i] <- "downreg"
  } else if (Motif3$GENEID[i] %in% background_noDownreg) {
    Motif3$Category[i] <- "background"
  } else {
    Motif3$Category[i] <- NA        
  }
}

Motif4$Category <- NA
for (i in seq_along(Motif4$GENEID)) {
  if (Motif4$GENEID[[i]] %in% upreg) {
    Motif4$Category[i] <- "upreg"
  } else if (Motif4$GENEID[i] %in% downreg) {
    Motif4$Category[i] <- "downreg"
  } else if (Motif4$GENEID[i] %in% background_noDownreg) {
    Motif4$Category[i] <- "background"
  } else {
    Motif4$Category[i] <- NA        
  }
}

Motif5$Category <- NA
for (i in seq_along(Motif5$GENEID)) {
  if (Motif5$GENEID[[i]] %in% upreg) {
    Motif5$Category[i] <- "upreg"
  } else if (Motif5$GENEID[i] %in% downreg) {
    Motif5$Category[i] <- "downreg"
  } else if (Motif5$GENEID[i] %in% background_noDownreg) {
    Motif5$Category[i] <- "background"
  } else {
    Motif5$Category[i] <- NA        
  }
}

Motif6$Category <- NA
for (i in seq_along(Motif6$GENEID)) {
  if (Motif6$GENEID[[i]] %in% upreg) {
    Motif6$Category[i] <- "upreg"
  } else if (Motif6$GENEID[i] %in% downreg) {
    Motif6$Category[i] <- "downreg"
  } else if (Motif6$GENEID[i] %in% background_noDownreg) {
    Motif6$Category[i] <- "background"
  } else {
    Motif6$Category[i] <- NA        
  }
}

Motif7$Category <- NA
for (i in seq_along(Motif7$GENEID)) {
  if (Motif7$GENEID[[i]] %in% upreg) {
    Motif7$Category[i] <- "upreg"
  } else if (Motif7$GENEID[i] %in% downreg) {
    Motif7$Category[i] <- "downreg"
  } else if (Motif7$GENEID[i] %in% background_noDownreg) {
    Motif7$Category[i] <- "background"
  } else {
    Motif7$Category[i] <- NA        
  }
}

Motif8$Category <- NA
for (i in seq_along(Motif8$GENEID)) {
  if (Motif8$GENEID[[i]] %in% upreg) {
    Motif8$Category[i] <- "upreg"
  } else if (Motif8$GENEID[i] %in% downreg) {
    Motif8$Category[i] <- "downreg"
  } else if (Motif8$GENEID[i] %in% background_noDownreg) {
    Motif8$Category[i] <- "background"
  } else {
    Motif8$Category[i] <- NA        
  }
}

Motif9$Category <- NA
for (i in seq_along(Motif9$GENEID)) {
  if (Motif9$GENEID[[i]] %in% upreg) {
    Motif9$Category[i] <- "upreg"
  } else if (Motif9$GENEID[i] %in% downreg) {
    Motif9$Category[i] <- "downreg"
  } else if (Motif9$GENEID[i] %in% background_noDownreg) {
    Motif9$Category[i] <- "background"
  } else {
    Motif9$Category[i] <- NA        
  }
}

Motif10$Category <- NA
for (i in seq_along(Motif10$GENEID)) {
  if (Motif10$GENEID[[i]] %in% upreg) {
    Motif10$Category[i] <- "upreg"
  } else if (Motif10$GENEID[i] %in% downreg) {
    Motif10$Category[i] <- "downreg"
  } else if (Motif10$GENEID[i] %in% background_noDownreg) {
    Motif10$Category[i] <- "background"
  } else {
    Motif10$Category[i] <- NA        
  }
}

Motif11$Category <- NA
for (i in seq_along(Motif11$GENEID)) {
  if (Motif11$GENEID[[i]] %in% upreg) {
    Motif11$Category[i] <- "upreg"
  } else if (Motif11$GENEID[i] %in% downreg) {
    Motif11$Category[i] <- "downreg"
  } else if (Motif11$GENEID[i] %in% background_noDownreg) {
    Motif11$Category[i] <- "background"
  } else {
    Motif11$Category[i] <- NA        
  }
}

Motif12$Category <- NA
for (i in seq_along(Motif12$GENEID)) {
  if (Motif12$GENEID[[i]] %in% upreg) {
    Motif12$Category[i] <- "upreg"
  } else if (Motif12$GENEID[i] %in% downreg) {
    Motif12$Category[i] <- "downreg"
  } else if (Motif12$GENEID[i] %in% background_noDownreg) {
    Motif12$Category[i] <- "background"
  } else {
    Motif12$Category[i] <- NA        
  }
}

Motif13$Category <- NA
for (i in seq_along(Motif13$GENEID)) {
  if (Motif13$GENEID[[i]] %in% upreg) {
    Motif13$Category[i] <- "upreg"
  } else if (Motif13$GENEID[i] %in% downreg) {
    Motif13$Category[i] <- "downreg"
  } else if (Motif13$GENEID[i] %in% background_noDownreg) {
    Motif13$Category[i] <- "background"
  } else {
    Motif13$Category[i] <- NA        
  }
}

Motif14$Category <- NA
for (i in seq_along(Motif14$GENEID)) {
  if (Motif14$GENEID[[i]] %in% upreg) {
    Motif14$Category[i] <- "upreg"
  } else if (Motif14$GENEID[i] %in% downreg) {
    Motif14$Category[i] <- "downreg"
  } else if (Motif14$GENEID[i] %in% background_noDownreg) {
    Motif14$Category[i] <- "background"
  } else {
    Motif14$Category[i] <- NA        
  }
}

Motif15$Category <- NA
for (i in seq_along(Motif15$GENEID)) {
  if (Motif15$GENEID[[i]] %in% upreg) {
    Motif15$Category[i] <- "upreg"
  } else if (Motif15$GENEID[i] %in% downreg) {
    Motif15$Category[i] <- "downreg"
  } else if (Motif15$GENEID[i] %in% background_noDownreg) {
    Motif15$Category[i] <- "background"
  } else {
    Motif15$Category[i] <- NA        
  }
}

Motif16$Category <- NA
for (i in seq_along(Motif16$GENEID)) {
  if (Motif16$GENEID[[i]] %in% upreg) {
    Motif16$Category[i] <- "upreg"
  } else if (Motif16$GENEID[i] %in% downreg) {
    Motif16$Category[i] <- "downreg"
  } else if (Motif16$GENEID[i] %in% background_noDownreg) {
    Motif16$Category[i] <- "background"
  } else {
    Motif16$Category[i] <- NA        
  }
}

Motif17$Category <- NA
for (i in seq_along(Motif17$GENEID)) {
  if (Motif17$GENEID[[i]] %in% upreg) {
    Motif17$Category[i] <- "upreg"
  } else if (Motif17$GENEID[i] %in% downreg) {
    Motif17$Category[i] <- "downreg"
  } else if (Motif17$GENEID[i] %in% background_noDownreg) {
    Motif17$Category[i] <- "background"
  } else {
    Motif17$Category[i] <- NA        
  }
}

####add motif reg exp to data frame####
Motif1$Motif <- (subset(Motifs, Motifs$Motif == "Motif1"))$RegExp
Motif2$Motif <- (subset(Motifs, Motifs$Motif == "Motif2"))$RegExp
Motif3$Motif <- (subset(Motifs, Motifs$Motif == "Motif3"))$RegExp
Motif4$Motif <- (subset(Motifs, Motifs$Motif == "Motif4"))$RegExp
Motif5$Motif <- (subset(Motifs, Motifs$Motif == "Motif5"))$RegExp
Motif6$Motif <- (subset(Motifs, Motifs$Motif == "Motif6"))$RegExp
Motif7$Motif <- (subset(Motifs, Motifs$Motif == "Motif7"))$RegExp
Motif8$Motif <- (subset(Motifs, Motifs$Motif == "Motif8"))$RegExp
Motif9$Motif <- (subset(Motifs, Motifs$Motif == "Motif9"))$RegExp
Motif10$Motif <- (subset(Motifs, Motifs$Motif == "Motif10"))$RegExp
Motif11$Motif <- (subset(Motifs, Motifs$Motif == "Motif11"))$RegExp
Motif12$Motif <- (subset(Motifs, Motifs$Motif == "Motif12"))$RegExp
Motif13$Motif <- (subset(Motifs, Motifs$Motif == "Motif13"))$RegExp
Motif14$Motif <- (subset(Motifs, Motifs$Motif == "Motif14"))$RegExp
Motif15$Motif <- (subset(Motifs, Motifs$Motif == "Motif15"))$RegExp
Motif16$Motif <- (subset(Motifs, Motifs$Motif == "Motif16"))$RegExp
Motif17$Motif <- (subset(Motifs, Motifs$Motif == "Motif17"))$RegExp


Motif1$MotifName <- (subset(Motifs, Motifs$Motif == "Motif1"))$Motif
Motif2$MotifName <- (subset(Motifs, Motifs$Motif == "Motif2"))$Motif
Motif3$MotifName <- (subset(Motifs, Motifs$Motif == "Motif3"))$Motif
Motif4$MotifName <- (subset(Motifs, Motifs$Motif == "Motif4"))$Motif
Motif5$MotifName <- (subset(Motifs, Motifs$Motif == "Motif5"))$Motif
Motif6$MotifName <- (subset(Motifs, Motifs$Motif == "Motif6"))$Motif
Motif7$MotifName <- (subset(Motifs, Motifs$Motif == "Motif7"))$Motif
Motif8$MotifName <- (subset(Motifs, Motifs$Motif == "Motif8"))$Motif
Motif9$MotifName <- (subset(Motifs, Motifs$Motif == "Motif9"))$Motif
Motif10$MotifName <- (subset(Motifs, Motifs$Motif == "Motif10"))$Motif
Motif11$MotifName <- (subset(Motifs, Motifs$Motif == "Motif11"))$Motif
Motif12$MotifName <- (subset(Motifs, Motifs$Motif == "Motif12"))$Motif
Motif13$MotifName <- (subset(Motifs, Motifs$Motif == "Motif13"))$Motif
Motif14$MotifName <- (subset(Motifs, Motifs$Motif == "Motif14"))$Motif
Motif15$MotifName <- (subset(Motifs, Motifs$Motif == "Motif15"))$Motif
Motif16$MotifName <- (subset(Motifs, Motifs$Motif == "Motif16"))$Motif
Motif17$MotifName <- (subset(Motifs, Motifs$Motif == "Motif17"))$Motif

#####write out files#####
write.csv(Motif1, "Motif1_type.csv")
write.csv(Motif2, "Motif2_type.csv")
write.csv(Motif3, "Motif3_type.csv")
write.csv(Motif4, "Motif4_type.csv")
write.csv(Motif5, "Motif5_type.csv")
write.csv(Motif6, "Motif6_type.csv")
write.csv(Motif7, "Motif7_type.csv")
write.csv(Motif8, "Motif8_type.csv")
write.csv(Motif9, "Motif9_type.csv")
write.csv(Motif10, "Motif10_type.csv")
write.csv(Motif11, "Motif11_type.csv")
write.csv(Motif12, "Motif12_type.csv")
write.csv(Motif13, "Motif13_type.csv")
write.csv(Motif14, "Motif14_type.csv")
write.csv(Motif15, "Motif15_type.csv")
write.csv(Motif16, "Motif16_type.csv")
write.csv(Motif17, "Motif17_type.csv")




####combine into 1 data frame and split####
AllMotifOccurences <- rbind(Motif1,Motif2,Motif3,Motif4,Motif5,Motif6,Motif7,Motif8,Motif9,Motif10,Motif11,Motif12,Motif13,Motif14,Motif15,Motif16,Motif17)
UpregMotifOccurences <- subset(AllMotifOccurences, AllMotifOccurences$Category == "upreg")
DownregMotifOccurences <- subset(AllMotifOccurences, AllMotifOccurences$Category == "downreg")
BackgroundMotifOccurences <- subset(AllMotifOccurences, AllMotifOccurences$Category == "background")


#####output file#####
write.csv(AllMotifOccurences, "AllMotifOccurences.csv")
write.csv(UpregMotifOccurences, "UpregMotifOccurences.csv")
write.csv(DownregMotifOccurences, "DownregMotifOccurences.csv")
write.csv(BackgroundMotifOccurences, "BackgroundMotifOccurences.csv")


#####plots!#####
library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)

AllMotifOccurences$MotifName <- factor(AllMotifOccurences$MotifName, levels = c("Motif1","Motif2","Motif3",
                                                                                "Motif4","Motif5","Motif6",
                                                                                "Motif7","Motif8","Motif9",
                                                                                "Motif10","Motif11","Motif12",
                                                                                "Motif13","Motif14","Motif15","Motif16","Motif17"))
UpregMotifOccurences$MotifName <- factor(UpregMotifOccurences$MotifName, levels = c("Motif1","Motif2","Motif3",
                                                                        "Motif4","Motif5","Motif6",
                                                                        "Motif7","Motif8","Motif9",
                                                                        "Motif10","Motif11","Motif12",
                                                                        "Motif13","Motif14","Motif15","Motif16","Motif17"))
DownregMotifOccurences$MotifName <- factor(DownregMotifOccurences$MotifName, levels = c("Motif1","Motif2","Motif3",
                                                                          "Motif4","Motif5","Motif6",
                                                                          "Motif7","Motif8","Motif9",
                                                                          "Motif10","Motif11","Motif12",
                                                                          "Motif13","Motif14","Motif15","Motif16","Motif17"))
BackgroundMotifOccurences$MotifName <- factor(BackgroundMotifOccurences$MotifName, levels = c("Motif1","Motif2","Motif3",
                                                                                           "Motif4","Motif5","Motif6",
                                                                                           "Motif7","Motif8","Motif9",
                                                                                           "Motif10","Motif11","Motif12",
                                                                                           "Motif13","Motif14","Motif15","Motif16","Motif17"))


####All Motif Occurences####
###heat map
###prep data###
AllMotifGeneIDs <- unique(AllMotifOccurences$GENEID)
AllMotifMotifs <- unique(AllMotifOccurences$MotifName)

#in long format
AllMotifOccurences_countsdata <- AllMotifOccurences[,c(1,4)]

# #start by making a matrix of overlapping GeneID counts
# AllMotifOccurencesMatrix <- table(AllMotifOccurences_countsdata$GENEID, AllMotifOccurences_countsdata$MotifName)
# AllMotifOccurencesMatrix <- t(AllMotifOccurencesMatrix) %*% AllMotifOccurencesMatrix
# #that %*% is for matrix multiplication

#could also do
AllMotifOccurences_wide <- dcast(AllMotifOccurences_countsdata, GENEID ~ MotifName, value.var = "MotifName", 
                                              fun.aggregate = length)
#this gives me the wide table I'm looking for, geneID by motif with 1 or 0

#then using this crossprod, I get the counts matrix (same as AllMotifOccurencesMatrix)
AllMotifOccurencesMatrix2 <- crossprod(as.matrix(AllMotifOccurences_wide[,2:18]))

# now my data is in an easy matrix format with the overlaps between sets of genes, I also have it in long format (AllMotifOccurences)
# and wide format (AllMotifOccurences_wide)

#now i want to calculate the percentage overlap
#need to do this by percentage overlap, not number of genes
# is that even possible? the percentage of genes overlapping in Motif1 and Motif2 will be different depending how I look at it
# since they both have different numbers of genes to begin with
# see this: http://stackoverflow.com/questions/29929074/percentage-overlap-of-two-lists

# The maximum difference will be if two lists have completely different elements. 
#So we have at most n + m discrete elements, where n is size of first list and m is the size 
#of second list. One measure can be:
#   
#   2 * c / (n + m)
# 
# where c is the number of common elements. This can be calculated like this as percentage:
#   
#   200.0 * len(set(listA) & set(listB)) / (len(listA) + len(listB))

# could also do overlap coefficient?
# or jaccard index?
#intersect/(LengthA + LengthB - intersect)

#or use crossprod - will give me the different values based on the x axis?
# http://stackoverflow.com/questions/41581217/create-overlap-matrix-in-r
# this divides each value by its corresponding diagonal value (ie all of Motif 1 is divided by the Motif1 diagonal, all of 
# Motif 2 is divided by the Motif2 diagonal, etc)
AllMotifOccurencesMatrix_percent <- round(t(AllMotifOccurencesMatrix2*100/diag(AllMotifOccurencesMatrix2)), digits = 2)


###plot code###
#convert matrix to table for plot
AllMotifOccurencesCountsTable <- melt(AllMotifOccurencesMatrix2)
colnames(AllMotifOccurencesCountsTable) <- c("MotifA","MotifB","Overlap")

AllMotifOccurencesCountPlot <- ggplot(data = AllMotifOccurencesCountsTable, aes(x = MotifB, y = MotifA)) +
  geom_tile(aes(fill=Overlap))

AllMotifOccurencesPercentTable <- melt(AllMotifOccurencesMatrix_percent)
colnames(AllMotifOccurencesPercentTable) <- c("MotifA","MotifB","Overlap")
  s = subset(AllMotifOccurencesPercentTable, AllMotifOccurencesPercentTable$MotifA == AllMotifOccurencesPercentTable$MotifB)
  s$text <- diag(AllMotifOccurencesMatrix2)


AllMotifOccurencesPercentTable$PlotColours <- ifelse(AllMotifOccurencesPercentTable$MotifA == AllMotifOccurencesPercentTable$MotifB,
                                                     NA,
                                                     AllMotifOccurencesPercentTable$Overlap)
  
AllMotifOccurencesPercentPlot <- 
  ggplot(data = AllMotifOccurencesPercentTable, aes(x = MotifB, y = ordered(MotifA, levels = rev(levels(MotifA))))) +
  geom_tile(aes(fill=PlotColours)) +
  xlab("Motif") +
  ylab("Motif") +
  ggtitle("Motif Occurences: All Genes") + 
  scale_fill_gradientn(colours = c("deepskyblue4","deepskyblue1","White","red","red4"),
                       na.value = "black",
                       guide = guide_colourbar(title = "Similarity (% of genes)",
                                               title.hjust = 0.5,
                                               title.position = "top",
                                               barwidth = 10,
                                               ticks = FALSE,
                                               border = element_line(colour = "black"))) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, vjust = 0, hjust = 0),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.position = "bottom",
        legend.justification = 0.5,
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm")) + 
  geom_text(aes(label = s$text), data=s, size=6, colour="white")

#####upreg genes#####
UpregMotifGeneIDs <- unique(UpregMotifOccurences$GENEID)
UpregMotifMotifs <- unique(UpregMotifOccurences$MotifName)


#in long format
UpregMotifOccurences_countsdata <- UpregMotifOccurences[,c(1,4)]

UpregMotifOccurences_wide <- dcast(UpregMotifOccurences_countsdata, GENEID ~ MotifName, value.var = "MotifName", 
                                 fun.aggregate = length)

UpregMotifOccurencesMatrix <- crossprod(as.matrix(UpregMotifOccurences_wide[,2:17]))

UpregMotifOccurencesMatrix_percent <- round(t(UpregMotifOccurencesMatrix*100/diag(UpregMotifOccurencesMatrix)), digits = 2)



###plot code###
#convert matrix to table for plot
UpregMotifOccurencesCountsTable <- melt(UpregMotifOccurencesMatrix)
colnames(UpregMotifOccurencesCountsTable) <- c("MotifA","MotifB","Overlap")

UpregMotifOccurencesCountPlot <- ggplot(data = UpregMotifOccurencesCountsTable, aes(x = MotifB, y = MotifA)) +
  geom_tile(aes(fill=Overlap))

UpregMotifOccurencesPercentTable <- melt(UpregMotifOccurencesMatrix_percent)
colnames(UpregMotifOccurencesPercentTable) <- c("MotifA","MotifB","Overlap")
  q = subset(UpregMotifOccurencesPercentTable, UpregMotifOccurencesPercentTable$MotifA == UpregMotifOccurencesPercentTable$MotifB)
  q$text <- diag(UpregMotifOccurencesMatrix)

  
UpregMotifOccurencesPercentTable$PlotColours <- ifelse(UpregMotifOccurencesPercentTable$MotifA == UpregMotifOccurencesPercentTable$MotifB,
                                                       NA,
                                                       UpregMotifOccurencesPercentTable$Overlap)  
  
UpregMotifOccurencesPercentPlot <- 
  ggplot(data = UpregMotifOccurencesPercentTable, aes(x = MotifB, y = ordered(MotifA, levels = rev(levels(MotifA))))) +
  geom_tile(aes(fill=PlotColours)) +
  xlab("Motif") +
  ylab("Motif") +
  ggtitle("Motif Occurences: Upregulated Genes") + 
  scale_fill_gradientn(colours = c("deepskyblue4","deepskyblue1","White","red","red4"),
                       na.value = "black",
                       guide = guide_colourbar(title = "Similarity (% of genes)",
                                               title.hjust = 0.5,
                                               title.position = "top",
                                               barwidth = 10,
                                               ticks = FALSE,
                                               border = element_line(colour = "black"))) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, vjust = 0, hjust = 0),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.position = "bottom",
        legend.justification = 0.5,
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm")) + 
  geom_text(aes(label = q$text), data=q, size=6, colour="white")



#####downreg genes#####
DownregMotifGeneIDs <- unique(DownregMotifOccurences$GENEID)
DownregMotifMotifs <- unique(DownregMotifOccurences$MotifName)


#in long format
DownregMotifOccurences_countsdata <- DownregMotifOccurences[,c(1,4)]

DownregMotifOccurences_wide <- dcast(DownregMotifOccurences_countsdata, GENEID ~ MotifName, value.var = "MotifName", 
                                   fun.aggregate = length)

DownregMotifOccurencesMatrix <- crossprod(as.matrix(DownregMotifOccurences_wide[,2:14]))

DownregMotifOccurencesMatrix_percent <- round(t(DownregMotifOccurencesMatrix*100/diag(DownregMotifOccurencesMatrix)), digits = 2)



###plot code###
#convert matrix to table for plot
DownregMotifOccurencesCountsTable <- melt(DownregMotifOccurencesMatrix)
colnames(DownregMotifOccurencesCountsTable) <- c("MotifA","MotifB","Overlap")

DownregMotifOccurencesCountPlot <- ggplot(data = DownregMotifOccurencesCountsTable, aes(x = MotifB, y = MotifA)) +
  geom_tile(aes(fill=Overlap))

DownregMotifOccurencesPercentTable <- melt(DownregMotifOccurencesMatrix_percent)
colnames(DownregMotifOccurencesPercentTable) <- c("MotifA","MotifB","Overlap")
  t = subset(DownregMotifOccurencesPercentTable, DownregMotifOccurencesPercentTable$MotifA == DownregMotifOccurencesPercentTable$MotifB)
  t$text <- diag(DownregMotifOccurencesMatrix)

DownregMotifOccurencesPercentTable$PlotColours <- ifelse(DownregMotifOccurencesPercentTable$MotifA == DownregMotifOccurencesPercentTable$MotifB,
                                                         NA,
                                                         DownregMotifOccurencesPercentTable$Overlap)  
  

DownregMotifOccurencesPercentPlot <- 
  ggplot(data = DownregMotifOccurencesPercentTable, aes(x = MotifB, y = ordered(MotifA, levels = rev(levels(MotifA))))) +
  geom_tile(aes(fill=PlotColours)) +
  xlab("Motif") +
  ylab("Motif") +
  ggtitle("Motif Occurences: Downregulated Genes") + 
  scale_fill_gradientn(colours = c("deepskyblue4","deepskyblue1","White","red","red4"),
                       na.value = "black",
                       guide = guide_colourbar(title = "Similarity (% of genes)",
                                               title.hjust = 0.5,
                                               title.position = "top",
                                               barwidth = 10,
                                               ticks = FALSE,
                                               border = element_line(colour = "black"))) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, vjust = 0, hjust = 0),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.position = "bottom",
        legend.justification = 0.5,
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm")) + 
  geom_text(aes(label = t$text), data=t, size=6, colour="white")


#####background genes#####
BackgroundMotifGeneIDs <- unique(BackgroundMotifOccurences$GENEID)
BackgroundMotifMotifs <- unique(BackgroundMotifOccurences$MotifName)


#in long format
BackgroundMotifOccurences_countsdata <- BackgroundMotifOccurences[,c(1,4)]

BackgroundMotifOccurences_wide <- dcast(BackgroundMotifOccurences_countsdata, GENEID ~ MotifName, value.var = "MotifName", 
                                   fun.aggregate = length)

BackgroundMotifOccurencesMatrix <- crossprod(as.matrix(BackgroundMotifOccurences_wide[,2:18]))

BackgroundMotifOccurencesMatrix_percent <- round(t(BackgroundMotifOccurencesMatrix*100/diag(BackgroundMotifOccurencesMatrix)), digits = 2)



###plot code###
#convert matrix to table for plot
BackgroundMotifOccurencesCountsTable <- melt(BackgroundMotifOccurencesMatrix)
colnames(BackgroundMotifOccurencesCountsTable) <- c("MotifA","MotifB","Overlap")

BackgroundMotifOccurencesCountPlot <- ggplot(data = BackgroundMotifOccurencesCountsTable, aes(x = MotifB, y = MotifA)) +
  geom_tile(aes(fill=Overlap))

BackgroundMotifOccurencesPercentTable <- melt(BackgroundMotifOccurencesMatrix_percent)
colnames(BackgroundMotifOccurencesPercentTable) <- c("MotifA","MotifB","Overlap")
  y = subset(BackgroundMotifOccurencesPercentTable, BackgroundMotifOccurencesPercentTable$MotifA == BackgroundMotifOccurencesPercentTable$MotifB)
  y$text <- diag(BackgroundMotifOccurencesMatrix)

  
BackgroundMotifOccurencesPercentTable$PlotColours <- ifelse(BackgroundMotifOccurencesPercentTable$MotifA == BackgroundMotifOccurencesPercentTable$MotifB,
                                                           NA,
                                                           BackgroundMotifOccurencesPercentTable$Overlap)  
  

BackgroundMotifOccurencesPercentPlot <- 
  ggplot(data = BackgroundMotifOccurencesPercentTable, aes(x = MotifB, y = ordered(MotifA, levels = rev(levels(MotifA))))) +
  geom_tile(aes(fill=PlotColours)) +
  xlab("Motif") +
  ylab("Motif") +
  ggtitle("Motif Occurences: Background Genes") + 
  scale_fill_gradientn(colours = c("deepskyblue4","deepskyblue","White","red","red4"),
                       na.value = "black",
                       guide = guide_colourbar(title = "Similarity (% of genes)",
                                               title.hjust = 0.5,
                                               title.position = "top",
                                               barwidth = 10,
                                               ticks = FALSE,
                                               border = element_line(colour = "black"))) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, vjust = 0, hjust = 0),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.position = "bottom",
        legend.justification = 0.5,
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm")) +  
  geom_text(aes(label = y$text), data=y, size=6, colour="white")


#####save plots#####
ggsave(AllMotifOccurencesPercentPlot, file = "AllMotifOccurencesPercentPlot.pdf", height = 8.27, width = 11.69)
ggsave(UpregMotifOccurencesPercentPlot, file = "UpregMotifOccurencesPercentPlot.pdf", height = 8.27, width = 11.69)
ggsave(DownregMotifOccurencesPercentPlot, file = "DownregMotifOccurencesPercentPlot.pdf", height = 8.27, width = 11.69)
ggsave(BackgroundMotifOccurencesPercentPlot, file = "BackgroundMotifOccurencesPercentPlot.pdf", height = 8.27, width = 11.69)


#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/LysoPC_genes_with_motifs.RData")
