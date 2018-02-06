# load wd
load("J:/III/Waters/Group Members/Mallory/R/InitAlignOverall_Genome_Coverage.RData")

#set wd
setwd("J:/III/Waters/Group Members/Mallory/GenomeProjects/GenomesVariants/AlignmentQuality")

#load indivs names
indivs <- read.table("J:/III/Waters/Group Members/Mallory/GenomeProjects/scripts/files/InitAlignOverall_indivs.txt")
  indivs <- indivs[,1]
  
#load chr names
chr <- read.table("J:/III/Waters/Group Members/Mallory/GenomeProjects/scripts/files/InitAlignOverall_chromosomes.txt")
  chr <- chr[,1]
  
#load filenames in wd
filenames <- list.files()


####load data####

# was having SO many issues with this loop - problem was that indivs was a dataframe, not a factor
# once I converted it, worked perfectly
for(i in indivs)
{
  assign(paste("coverage_", i, sep=""), read.table(paste(i, "_coverage.txt", sep=""), 
                                                   col.names = c("Feature","Depth","NumBPAtDepth","SizeFeatureBP","FracFeatAtDepth")))
}

#load data into a list of data frames rather than individual ones (best practice?)
coverage_data <- lapply(filenames, read.table)
  names(coverage_data) <- paste0("coverage_",indivs)
colnames <- c("Feature","Depth","NumBPAtDepth","SizeFeatureBP","FracFeatAtDepth")
for (i in 1:length(coverage_data)){
  colnames(coverage_data[[i]]) <- colnames
}

#make a list of the coverage tables (don't really need)
coverage_tables <- as.list(names(coverage_data))


#include whole genome in the list of chr to get genomewide coverage as well as just the chr (for average coverage below)
chr_gen <- unique(coverage_233$Feature)

####calculate coverage####
average_coverage <- data.frame(matrix(NA, nrow=11, ncol=22))
colnames(average_coverage) <- chr_gen
rownames(average_coverage) <- indivs

for (j in 1:length(coverage_data)){
  for (i in chr_gen){
    subs <- subset(coverage_data[[j]], coverage_data[[j]]$Feature == i)
    average_coverage[j,i] <- subs$Depth%*%subs$NumBPAtDepth
  }
}

#reads/bp
average_coverage_bp <- data.frame(matrix(NA, nrow=11, ncol=22))
colnames(average_coverage_bp) <- chr_gen
rownames(average_coverage_bp) <- indivs

for (j in 1:length(coverage_data)){
  for (i in chr_gen){
    subs <- subset(coverage_data[[j]], coverage_data[[j]]$Feature == i)
    average_coverage_bp[j,i] <- subs$Depth%*%subs$NumBPAtDepth/mean(subs$SizeFeatureBP)
  }
}




# for different values of V1, want sum(V2*V3)/uniqueV4 (dot product is sum V2*V3)?

####output data####
write.csv(average_coverage, "InitAlignOverall_average_coverage.csv")
write.csv(average_coverage_bp, "InitAlignOverall_average_coverage_bp.csv")


#####FINAL SAVE#####
save.image("J:/III/Waters/Group Members/Mallory/R/InitAlignOverall_Genome_Coverage.RData")
