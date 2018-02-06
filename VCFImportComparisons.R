#load environment
load("J:/III/Waters/Group Members/Mallory/R/VCFImportComparisons.RData")

#install/load variant annotation package (has a readvcf tool)
library("VariantAnnotation") #load the package
library("vcfR")

#set WD
setwd("J:/III/Waters/Group Members/Mallory/GenomeProjects/GenomesVariants")

#####import line names and types#####
lines <- read.table("J:/III/Waters/Group Members/Mallory/GenomeProjects/Scripts/files/InitAlignOverall_indivs.txt")
lines <- as.list(lines[,1])

types <- read.table("J:/III/Waters/Group Members/Mallory/GenomeProjects/Scripts/files/InitAlignOverall_types.txt")
types <- as.list(types[,1])
#####Import Syn vars#####
# may need to change lines data type - need it as a list of 11
SynVCFs <- list()
for (i in lines){
  SynVCFs[[paste0(i,".Q60.sorted.filter.ann.vcf")]] <- read.vcfR(paste0("InitAlign_trimmed_aligned_sorted_solRD_GalaxyRG_GalaxyMerge_variants_",i,".Q60.sorted.filter.ann.vcf"))
}
#now do the same for all 3 types
for (j in types){
  for (i in lines){
    SynVCFs[[paste0(i,".",j,".Q60.sorted.filter.ann.vcf")]] <- read.vcfR(paste0("InitAlign_trimmed_aligned_sorted_solRD_GalaxyRG_GalaxyMerge_variants_",i,".",j,".Q60.sorted.filter.ann.vcf"))
  }
}
####Convert Syn to tidy DF and extract CHROM, POS, REF, ALT#####
#using VCFs
VCFalldata <- list()
VCFalldata <- lapply(SynVCFs, vcfR2tidy)
VCFdata <- list()
for (i in names(VCFalldata)){
  VCFdata[[i]] <- as.data.frame(VCFalldata[[i]][[1]][,c(2:3,5:6)])
  # VCFdata[[i]] <- as.data.frame(paste0("VCFalldata$",i,"$fix")[c(2:3,5:6)])
  # assign(VCFdata[[i]], paste0("VCFalldata$",i,"$fix[,1:6]"))
}
####compare Syn data frames GNPs and 820####
library("compare")
library("plyr")

#the dupsBetweenGroups function = 
#http://www.cookbook-r.com/Manipulating_data/Comparing_data_frames/

#add a column to each DF in VCFdata that lists the line
for (i in seq_along(VCFdata)){
  VCFdata[[i]]$Line <- names(VCFdata[i])
}

SynComps <- list()
SynCompsStats <- list()
####All####
SynComps$comp_820_GNPm7 <- rbind(VCFdata[["820.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                        VCFdata[["GNPm7.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_820_GNPm7$Line <- gsub(".Q60.sorted.filter.ann.vcf", "", SynComps$comp_820_GNPm7$Line)
SynComps$comp_820_GNPm7 <- cbind(SynComps$comp_820_GNPm7, InLine = dupsBetweenGroups(SynComps$comp_820_GNPm7, "Line"))
for (i in seq_along(SynComps$comp_820_GNPm7$InLine)){
  if(SynComps$comp_820_GNPm7$InLine[i] == "TRUE") {
    SynComps$comp_820_GNPm7$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_820_GNPm7$InLine[i] <- SynComps$comp_820_GNPm7$Line[i]
  }}
SynComps$comp_820_GNPm7 <- SynComps$comp_820_GNPm7[!duplicated(SynComps$comp_820_GNPm7[c("CHROM","POS")]),]
SynComps$comp_820_GNPm7 <- SynComps$comp_820_GNPm7[order(SynComps$comp_820_GNPm7$InLine),]
SynCompsStats$comp_820_GNPm7_stats <- count(SynComps$comp_820_GNPm7, 'InLine')
SynCompsStats$comp_820_GNPm7_stats$Comparison <- capture.output(cat(unique(SynComps$comp_820_GNPm7$Line), sep="_"))

SynComps$comp_820_GNPm8 <- rbind(VCFdata[["820.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                 VCFdata[["GNPm8.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_820_GNPm8$Line <- gsub(".Q60.sorted.filter.ann.vcf", "", SynComps$comp_820_GNPm8$Line)
SynComps$comp_820_GNPm8 <- cbind(SynComps$comp_820_GNPm8, InLine = dupsBetweenGroups(SynComps$comp_820_GNPm8, "Line"))
for (i in seq_along(SynComps$comp_820_GNPm8$InLine)){
  if(SynComps$comp_820_GNPm8$InLine[i] == "TRUE") {
    SynComps$comp_820_GNPm8$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_820_GNPm8$InLine[i] <- SynComps$comp_820_GNPm8$Line[i]
  }}
SynComps$comp_820_GNPm8 <- SynComps$comp_820_GNPm8[!duplicated(SynComps$comp_820_GNPm8[c("CHROM","POS")]),]
SynComps$comp_820_GNPm8 <- SynComps$comp_820_GNPm8[order(SynComps$comp_820_GNPm8$InLine),]
SynCompsStats$comp_820_GNPm8_stats <- count(SynComps$comp_820_GNPm8, 'InLine')
SynCompsStats$comp_820_GNPm8_stats$Comparison <- capture.output(cat(unique(SynComps$comp_820_GNPm8$Line), sep="_"))

SynComps$comp_820_GNPm9 <- rbind(VCFdata[["820.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                 VCFdata[["GNPm9.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_820_GNPm9$Line <- gsub(".Q60.sorted.filter.ann.vcf", "", SynComps$comp_820_GNPm9$Line)
SynComps$comp_820_GNPm9 <- cbind(SynComps$comp_820_GNPm9, InLine = dupsBetweenGroups(SynComps$comp_820_GNPm9, "Line"))
for (i in seq_along(SynComps$comp_820_GNPm9$InLine)){
  if(SynComps$comp_820_GNPm9$InLine[i] == "TRUE") {
    SynComps$comp_820_GNPm9$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_820_GNPm9$InLine[i] <- SynComps$comp_820_GNPm9$Line[i]
  }}
SynComps$comp_820_GNPm9 <- SynComps$comp_820_GNPm9[!duplicated(SynComps$comp_820_GNPm9[c("CHROM","POS")]),]
SynComps$comp_820_GNPm9 <- SynComps$comp_820_GNPm9[order(SynComps$comp_820_GNPm9$InLine),]
SynCompsStats$comp_820_GNPm9_stats <- count(SynComps$comp_820_GNPm9, 'InLine')
SynCompsStats$comp_820_GNPm9_stats$Comparison <- capture.output(cat(unique(SynComps$comp_820_GNPm9$Line), sep="_"))
####exon####
SynComps$comp_exon_820_GNPm7 <- rbind(VCFdata[["820.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                        VCFdata[["GNPm7.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_exon_820_GNPm7$Line <- gsub(".exon.Q60.sorted.filter.ann.vcf", "", SynComps$comp_exon_820_GNPm7$Line)
SynComps$comp_exon_820_GNPm7 <- cbind(SynComps$comp_exon_820_GNPm7, InLine = dupsBetweenGroups(SynComps$comp_exon_820_GNPm7, "Line"))
for (i in seq_along(SynComps$comp_exon_820_GNPm7$InLine)){
  if(SynComps$comp_exon_820_GNPm7$InLine[i] == "TRUE") {
    SynComps$comp_exon_820_GNPm7$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_exon_820_GNPm7$InLine[i] <- SynComps$comp_exon_820_GNPm7$Line[i]
  }}
SynComps$comp_exon_820_GNPm7 <- SynComps$comp_exon_820_GNPm7[!duplicated(SynComps$comp_exon_820_GNPm7[c("CHROM","POS")]),]
SynComps$comp_exon_820_GNPm7 <- SynComps$comp_exon_820_GNPm7[order(SynComps$comp_exon_820_GNPm7$InLine),]
SynCompsStats$comp_exon_820_GNPm7_stats <- count(SynComps$comp_exon_820_GNPm7, 'InLine')
SynCompsStats$comp_exon_820_GNPm7_stats$Comparison <- capture.output(cat(unique(SynComps$comp_exon_820_GNPm7$Line), sep="_"))

SynComps$comp_exon_820_GNPm8 <- rbind(VCFdata[["820.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                      VCFdata[["GNPm8.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_exon_820_GNPm8$Line <- gsub(".exon.Q60.sorted.filter.ann.vcf", "", SynComps$comp_exon_820_GNPm8$Line)
SynComps$comp_exon_820_GNPm8 <- cbind(SynComps$comp_exon_820_GNPm8, InLine = dupsBetweenGroups(SynComps$comp_exon_820_GNPm8, "Line"))
for (i in seq_along(SynComps$comp_exon_820_GNPm8$InLine)){
  if(SynComps$comp_exon_820_GNPm8$InLine[i] == "TRUE") {
    SynComps$comp_exon_820_GNPm8$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_exon_820_GNPm8$InLine[i] <- SynComps$comp_exon_820_GNPm8$Line[i]
  }}
SynComps$comp_exon_820_GNPm8 <- SynComps$comp_exon_820_GNPm8[!duplicated(SynComps$comp_exon_820_GNPm8[c("CHROM","POS")]),]
SynComps$comp_exon_820_GNPm8 <- SynComps$comp_exon_820_GNPm8[order(SynComps$comp_exon_820_GNPm8$InLine),]
SynCompsStats$comp_exon_820_GNPm8_stats <- count(SynComps$comp_exon_820_GNPm8, 'InLine')
SynCompsStats$comp_exon_820_GNPm8_stats$Comparison <- capture.output(cat(unique(SynComps$comp_exon_820_GNPm8$Line), sep="_"))

SynComps$comp_exon_820_GNPm9 <- rbind(VCFdata[["820.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                      VCFdata[["GNPm9.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_exon_820_GNPm9$Line <- gsub(".exon.Q60.sorted.filter.ann.vcf", "", SynComps$comp_exon_820_GNPm9$Line)
SynComps$comp_exon_820_GNPm9 <- cbind(SynComps$comp_exon_820_GNPm9, InLine = dupsBetweenGroups(SynComps$comp_exon_820_GNPm9, "Line"))
for (i in seq_along(SynComps$comp_exon_820_GNPm9$InLine)){
  if(SynComps$comp_exon_820_GNPm9$InLine[i] == "TRUE") {
    SynComps$comp_exon_820_GNPm9$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_exon_820_GNPm9$InLine[i] <- SynComps$comp_exon_820_GNPm9$Line[i]
  }}
SynComps$comp_exon_820_GNPm9 <- SynComps$comp_exon_820_GNPm9[!duplicated(SynComps$comp_exon_820_GNPm9[c("CHROM","POS")]),]
SynComps$comp_exon_820_GNPm9 <- SynComps$comp_exon_820_GNPm9[order(SynComps$comp_exon_820_GNPm9$InLine),]
SynCompsStats$comp_exon_820_GNPm9_stats <- count(SynComps$comp_exon_820_GNPm9, 'InLine')
SynCompsStats$comp_exon_820_GNPm9_stats$Comparison <- capture.output(cat(unique(SynComps$comp_exon_820_GNPm9$Line), sep="_"))
####CDS####
SynComps$comp_CDS_820_GNPm7 <- rbind(VCFdata[["820.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                      VCFdata[["GNPm7.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_CDS_820_GNPm7$Line <- gsub(".CDS.Q60.sorted.filter.ann.vcf", "", SynComps$comp_CDS_820_GNPm7$Line)
SynComps$comp_CDS_820_GNPm7 <- cbind(SynComps$comp_CDS_820_GNPm7, InLine = dupsBetweenGroups(SynComps$comp_CDS_820_GNPm7, "Line"))
for (i in seq_along(SynComps$comp_CDS_820_GNPm7$InLine)){
  if(SynComps$comp_CDS_820_GNPm7$InLine[i] == "TRUE") {
    SynComps$comp_CDS_820_GNPm7$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_CDS_820_GNPm7$InLine[i] <- SynComps$comp_CDS_820_GNPm7$Line[i]
  }}
SynComps$comp_CDS_820_GNPm7 <- SynComps$comp_CDS_820_GNPm7[!duplicated(SynComps$comp_CDS_820_GNPm7[c("CHROM","POS")]),]
SynComps$comp_CDS_820_GNPm7 <- SynComps$comp_CDS_820_GNPm7[order(SynComps$comp_CDS_820_GNPm7$InLine),]
SynCompsStats$comp_CDS_820_GNPm7_stats <- count(SynComps$comp_CDS_820_GNPm7, 'InLine')
SynCompsStats$comp_CDS_820_GNPm7_stats$Comparison <- capture.output(cat(unique(SynComps$comp_CDS_820_GNPm7$Line), sep="_"))

SynComps$comp_CDS_820_GNPm8 <- rbind(VCFdata[["820.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                      VCFdata[["GNPm8.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_CDS_820_GNPm8$Line <- gsub(".CDS.Q60.sorted.filter.ann.vcf", "", SynComps$comp_CDS_820_GNPm8$Line)
SynComps$comp_CDS_820_GNPm8 <- cbind(SynComps$comp_CDS_820_GNPm8, InLine = dupsBetweenGroups(SynComps$comp_CDS_820_GNPm8, "Line"))
for (i in seq_along(SynComps$comp_CDS_820_GNPm8$InLine)){
  if(SynComps$comp_CDS_820_GNPm8$InLine[i] == "TRUE") {
    SynComps$comp_CDS_820_GNPm8$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_CDS_820_GNPm8$InLine[i] <- SynComps$comp_CDS_820_GNPm8$Line[i]
  }}
SynComps$comp_CDS_820_GNPm8 <- SynComps$comp_CDS_820_GNPm8[!duplicated(SynComps$comp_CDS_820_GNPm8[c("CHROM","POS")]),]
SynComps$comp_CDS_820_GNPm8 <- SynComps$comp_CDS_820_GNPm8[order(SynComps$comp_CDS_820_GNPm8$InLine),]
SynCompsStats$comp_CDS_820_GNPm8_stats <- count(SynComps$comp_CDS_820_GNPm8, 'InLine')
SynCompsStats$comp_CDS_820_GNPm8_stats$Comparison <- capture.output(cat(unique(SynComps$comp_CDS_820_GNPm8$Line), sep="_"))

SynComps$comp_CDS_820_GNPm9 <- rbind(VCFdata[["820.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                      VCFdata[["GNPm9.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_CDS_820_GNPm9$Line <- gsub(".CDS.Q60.sorted.filter.ann.vcf", "", SynComps$comp_CDS_820_GNPm9$Line)
SynComps$comp_CDS_820_GNPm9 <- cbind(SynComps$comp_CDS_820_GNPm9, InLine = dupsBetweenGroups(SynComps$comp_CDS_820_GNPm9, "Line"))
for (i in seq_along(SynComps$comp_CDS_820_GNPm9$InLine)){
  if(SynComps$comp_CDS_820_GNPm9$InLine[i] == "TRUE") {
    SynComps$comp_CDS_820_GNPm9$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_CDS_820_GNPm9$InLine[i] <- SynComps$comp_CDS_820_GNPm9$Line[i]
  }}
SynComps$comp_CDS_820_GNPm9 <- SynComps$comp_CDS_820_GNPm9[!duplicated(SynComps$comp_CDS_820_GNPm9[c("CHROM","POS")]),]
SynComps$comp_CDS_820_GNPm9 <- SynComps$comp_CDS_820_GNPm9[order(SynComps$comp_CDS_820_GNPm9$InLine),]
SynCompsStats$comp_CDS_820_GNPm9_stats <- count(SynComps$comp_CDS_820_GNPm9, 'InLine')
SynCompsStats$comp_CDS_820_GNPm9_stats$Comparison <- capture.output(cat(unique(SynComps$comp_CDS_820_GNPm9$Line), sep="_"))
####gene####
SynComps$comp_gene_820_GNPm7 <- rbind(VCFdata[["820.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                      VCFdata[["GNPm7.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_gene_820_GNPm7$Line <- gsub(".gene.Q60.sorted.filter.ann.vcf", "", SynComps$comp_gene_820_GNPm7$Line)
SynComps$comp_gene_820_GNPm7 <- cbind(SynComps$comp_gene_820_GNPm7, InLine = dupsBetweenGroups(SynComps$comp_gene_820_GNPm7, "Line"))
for (i in seq_along(SynComps$comp_gene_820_GNPm7$InLine)){
  if(SynComps$comp_gene_820_GNPm7$InLine[i] == "TRUE") {
    SynComps$comp_gene_820_GNPm7$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_gene_820_GNPm7$InLine[i] <- SynComps$comp_gene_820_GNPm7$Line[i]
  }}
SynComps$comp_gene_820_GNPm7 <- SynComps$comp_gene_820_GNPm7[!duplicated(SynComps$comp_gene_820_GNPm7[c("CHROM","POS")]),]
SynComps$comp_gene_820_GNPm7 <- SynComps$comp_gene_820_GNPm7[order(SynComps$comp_gene_820_GNPm7$InLine),]
SynCompsStats$comp_gene_820_GNPm7_stats <- count(SynComps$comp_gene_820_GNPm7, 'InLine')
SynCompsStats$comp_gene_820_GNPm7_stats$Comparison <- capture.output(cat(unique(SynComps$comp_gene_820_GNPm7$Line), sep="_"))

SynComps$comp_gene_820_GNPm8 <- rbind(VCFdata[["820.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                      VCFdata[["GNPm8.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_gene_820_GNPm8$Line <- gsub(".gene.Q60.sorted.filter.ann.vcf", "", SynComps$comp_gene_820_GNPm8$Line)
SynComps$comp_gene_820_GNPm8 <- cbind(SynComps$comp_gene_820_GNPm8, InLine = dupsBetweenGroups(SynComps$comp_gene_820_GNPm8, "Line"))
for (i in seq_along(SynComps$comp_gene_820_GNPm8$InLine)){
  if(SynComps$comp_gene_820_GNPm8$InLine[i] == "TRUE") {
    SynComps$comp_gene_820_GNPm8$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_gene_820_GNPm8$InLine[i] <- SynComps$comp_gene_820_GNPm8$Line[i]
  }}
SynComps$comp_gene_820_GNPm8 <- SynComps$comp_gene_820_GNPm8[!duplicated(SynComps$comp_gene_820_GNPm8[c("CHROM","POS")]),]
SynComps$comp_gene_820_GNPm8 <- SynComps$comp_gene_820_GNPm8[order(SynComps$comp_gene_820_GNPm8$InLine),]
SynCompsStats$comp_gene_820_GNPm8_stats <- count(SynComps$comp_gene_820_GNPm8, 'InLine')
SynCompsStats$comp_gene_820_GNPm8_stats$Comparison <- capture.output(cat(unique(SynComps$comp_gene_820_GNPm8$Line), sep="_"))

SynComps$comp_gene_820_GNPm9 <- rbind(VCFdata[["820.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                      VCFdata[["GNPm9.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
SynComps$comp_gene_820_GNPm9$Line <- gsub(".gene.Q60.sorted.filter.ann.vcf", "", SynComps$comp_gene_820_GNPm9$Line)
SynComps$comp_gene_820_GNPm9 <- cbind(SynComps$comp_gene_820_GNPm9, InLine = dupsBetweenGroups(SynComps$comp_gene_820_GNPm9, "Line"))
for (i in seq_along(SynComps$comp_gene_820_GNPm9$InLine)){
  if(SynComps$comp_gene_820_GNPm9$InLine[i] == "TRUE") {
    SynComps$comp_gene_820_GNPm9$InLine[i] <- "BOTH"
  }  else {
    SynComps$comp_gene_820_GNPm9$InLine[i] <- SynComps$comp_gene_820_GNPm9$Line[i]
  }}
SynComps$comp_gene_820_GNPm9 <- SynComps$comp_gene_820_GNPm9[!duplicated(SynComps$comp_gene_820_GNPm9[c("CHROM","POS")]),]
SynComps$comp_gene_820_GNPm9 <- SynComps$comp_gene_820_GNPm9[order(SynComps$comp_gene_820_GNPm9$InLine),]
SynCompsStats$comp_gene_820_GNPm9_stats <- count(SynComps$comp_gene_820_GNPm9, 'InLine')
SynCompsStats$comp_gene_820_GNPm9_stats$Comparison <- capture.output(cat(unique(SynComps$comp_gene_820_GNPm9$Line), sep="_"))

####Import nonSyn vars####
nonSynVCFs <- list()
for (i in lines){
  nonSynVCFs[[paste0(i,".Q60.sorted.filter.ann.vcf")]] <- read.vcfR(paste0("InitAlign_trimmed_aligned_sorted_solRD_GalaxyRG_GalaxyMerge_nonSyn_variants_",i,".Q60.sorted.filter.ann.vcf"))
}
#now do the same for all 3 types
for (j in types){
  for (i in lines){
    nonSynVCFs[[paste0(i,".",j,".Q60.sorted.filter.ann.vcf")]] <- read.vcfR(paste0("InitAlign_trimmed_aligned_sorted_solRD_GalaxyRG_GalaxyMerge_nonSyn_variants_",i,".",j,".Q60.sorted.filter.ann.vcf"))
  }
}
####Convert nonSyn to tidy DF and extract CHROM, POS, REF, ALT#####
#using VCFs
nonSynVCFalldata <- list()
nonSynVCFalldata <- lapply(nonSynVCFs, vcfR2tidy)
nonSynVCFdata <- list()
for (i in names(nonSynVCFalldata)){
  nonSynVCFdata[[i]] <- as.data.frame(nonSynVCFalldata[[i]][[1]][,c(2:3,5:6)])
  # VCFdata[[i]] <- as.data.frame(paste0("VCFalldata$",i,"$fix")[c(2:3,5:6)])
  # assign(VCFdata[[i]], paste0("VCFalldata$",i,"$fix[,1:6]"))
}



#2 Dec
#okay, have my numbers for comparisons

#5 december
#extract list of CHROM and POS for each GNP line (not 820)


GNPm7_not820_vars <- subset(comp_820_GNPm7[,1:2], comp_820_GNPm7$InLine == "GNPm7")

GNPm8_not820_vars <- subset(comp_820_GNPm8[,1:2], comp_820_GNPm8$InLine == "GNPm8")

GNPm9_not820_vars <- subset(comp_820_GNPm9[,1:2], comp_820_GNPm9$InLine == "GNPm9")
####Compare nonSyn data frames GNPs and 820####
#add a column to each DF in VCFdata that lists the line
for (i in seq_along(nonSynVCFdata)){
  nonSynVCFdata[[i]]$Line <- names(nonSynVCFdata[i])
}

nonSynComps <- list()
nonSynCompsStats <- list()
####nonSyn All####
nonSynComps$comp_820_GNPm7 <- rbind(nonSynVCFdata[["820.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                    nonSynVCFdata[["GNPm7.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_820_GNPm7$Line <- gsub(".Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_820_GNPm7$Line)
nonSynComps$comp_820_GNPm7 <- cbind(nonSynComps$comp_820_GNPm7, InLine = dupsBetweenGroups(nonSynComps$comp_820_GNPm7, "Line"))
for (i in seq_along(nonSynComps$comp_820_GNPm7$InLine)){
  if(nonSynComps$comp_820_GNPm7$InLine[i] == "TRUE") {
    nonSynComps$comp_820_GNPm7$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_820_GNPm7$InLine[i] <- nonSynComps$comp_820_GNPm7$Line[i]
  }}
nonSynComps$comp_820_GNPm7 <- nonSynComps$comp_820_GNPm7[!duplicated(nonSynComps$comp_820_GNPm7[c("CHROM","POS")]),]
nonSynComps$comp_820_GNPm7 <- nonSynComps$comp_820_GNPm7[order(nonSynComps$comp_820_GNPm7$InLine),]
nonSynCompsStats$comp_820_GNPm7_stats <- count(nonSynComps$comp_820_GNPm7, 'InLine')
nonSynCompsStats$comp_820_GNPm7_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_820_GNPm7$Line), sep="_"))

nonSynComps$comp_820_GNPm8 <- rbind(nonSynVCFdata[["820.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                    nonSynVCFdata[["GNPm8.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_820_GNPm8$Line <- gsub(".Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_820_GNPm8$Line)
nonSynComps$comp_820_GNPm8 <- cbind(nonSynComps$comp_820_GNPm8, InLine = dupsBetweenGroups(nonSynComps$comp_820_GNPm8, "Line"))
for (i in seq_along(nonSynComps$comp_820_GNPm8$InLine)){
  if(nonSynComps$comp_820_GNPm8$InLine[i] == "TRUE") {
    nonSynComps$comp_820_GNPm8$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_820_GNPm8$InLine[i] <- nonSynComps$comp_820_GNPm8$Line[i]
  }}
nonSynComps$comp_820_GNPm8 <- nonSynComps$comp_820_GNPm8[!duplicated(nonSynComps$comp_820_GNPm8[c("CHROM","POS")]),]
nonSynComps$comp_820_GNPm8 <- nonSynComps$comp_820_GNPm8[order(nonSynComps$comp_820_GNPm8$InLine),]
nonSynCompsStats$comp_820_GNPm8_stats <- count(nonSynComps$comp_820_GNPm8, 'InLine')
nonSynCompsStats$comp_820_GNPm8_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_820_GNPm8$Line), sep="_"))

nonSynComps$comp_820_GNPm9 <- rbind(nonSynVCFdata[["820.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                    nonSynVCFdata[["GNPm9.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_820_GNPm9$Line <- gsub(".Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_820_GNPm9$Line)
nonSynComps$comp_820_GNPm9 <- cbind(nonSynComps$comp_820_GNPm9, InLine = dupsBetweenGroups(nonSynComps$comp_820_GNPm9, "Line"))
for (i in seq_along(nonSynComps$comp_820_GNPm9$InLine)){
  if(nonSynComps$comp_820_GNPm9$InLine[i] == "TRUE") {
    nonSynComps$comp_820_GNPm9$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_820_GNPm9$InLine[i] <- nonSynComps$comp_820_GNPm9$Line[i]
  }}
nonSynComps$comp_820_GNPm9 <- nonSynComps$comp_820_GNPm9[!duplicated(nonSynComps$comp_820_GNPm9[c("CHROM","POS")]),]
nonSynComps$comp_820_GNPm9 <- nonSynComps$comp_820_GNPm9[order(nonSynComps$comp_820_GNPm9$InLine),]
nonSynCompsStats$comp_820_GNPm9_stats <- count(nonSynComps$comp_820_GNPm9, 'InLine')
nonSynCompsStats$comp_820_GNPm9_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_820_GNPm9$Line), sep="_"))
####nonSyn exon####
nonSynComps$comp_exon_820_GNPm7 <- rbind(nonSynVCFdata[["820.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                         nonSynVCFdata[["GNPm7.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_exon_820_GNPm7$Line <- gsub(".exon.Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_exon_820_GNPm7$Line)
nonSynComps$comp_exon_820_GNPm7 <- cbind(nonSynComps$comp_exon_820_GNPm7, InLine = dupsBetweenGroups(nonSynComps$comp_exon_820_GNPm7, "Line"))
for (i in seq_along(nonSynComps$comp_exon_820_GNPm7$InLine)){
  if(nonSynComps$comp_exon_820_GNPm7$InLine[i] == "TRUE") {
    nonSynComps$comp_exon_820_GNPm7$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_exon_820_GNPm7$InLine[i] <- nonSynComps$comp_exon_820_GNPm7$Line[i]
  }}
nonSynComps$comp_exon_820_GNPm7 <- nonSynComps$comp_exon_820_GNPm7[!duplicated(nonSynComps$comp_exon_820_GNPm7[c("CHROM","POS")]),]
nonSynComps$comp_exon_820_GNPm7 <- nonSynComps$comp_exon_820_GNPm7[order(nonSynComps$comp_exon_820_GNPm7$InLine),]
nonSynCompsStats$comp_exon_820_GNPm7_stats <- count(nonSynComps$comp_exon_820_GNPm7, 'InLine')
nonSynCompsStats$comp_exon_820_GNPm7_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_exon_820_GNPm7$Line), sep="_"))

nonSynComps$comp_exon_820_GNPm8 <- rbind(nonSynVCFdata[["820.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                         nonSynVCFdata[["GNPm8.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_exon_820_GNPm8$Line <- gsub(".exon.Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_exon_820_GNPm8$Line)
nonSynComps$comp_exon_820_GNPm8 <- cbind(nonSynComps$comp_exon_820_GNPm8, InLine = dupsBetweenGroups(nonSynComps$comp_exon_820_GNPm8, "Line"))
for (i in seq_along(nonSynComps$comp_exon_820_GNPm8$InLine)){
  if(nonSynComps$comp_exon_820_GNPm8$InLine[i] == "TRUE") {
    nonSynComps$comp_exon_820_GNPm8$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_exon_820_GNPm8$InLine[i] <- nonSynComps$comp_exon_820_GNPm8$Line[i]
  }}
nonSynComps$comp_exon_820_GNPm8 <- nonSynComps$comp_exon_820_GNPm8[!duplicated(nonSynComps$comp_exon_820_GNPm8[c("CHROM","POS")]),]
nonSynComps$comp_exon_820_GNPm8 <- nonSynComps$comp_exon_820_GNPm8[order(nonSynComps$comp_exon_820_GNPm8$InLine),]
nonSynCompsStats$comp_exon_820_GNPm8_stats <- count(nonSynComps$comp_exon_820_GNPm8, 'InLine')
nonSynCompsStats$comp_exon_820_GNPm8_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_exon_820_GNPm8$Line), sep="_"))

nonSynComps$comp_exon_820_GNPm9 <- rbind(nonSynVCFdata[["820.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                         nonSynVCFdata[["GNPm9.exon.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_exon_820_GNPm9$Line <- gsub(".exon.Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_exon_820_GNPm9$Line)
nonSynComps$comp_exon_820_GNPm9 <- cbind(nonSynComps$comp_exon_820_GNPm9, InLine = dupsBetweenGroups(nonSynComps$comp_exon_820_GNPm9, "Line"))
for (i in seq_along(nonSynComps$comp_exon_820_GNPm9$InLine)){
  if(nonSynComps$comp_exon_820_GNPm9$InLine[i] == "TRUE") {
    nonSynComps$comp_exon_820_GNPm9$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_exon_820_GNPm9$InLine[i] <- nonSynComps$comp_exon_820_GNPm9$Line[i]
  }}
nonSynComps$comp_exon_820_GNPm9 <- nonSynComps$comp_exon_820_GNPm9[!duplicated(nonSynComps$comp_exon_820_GNPm9[c("CHROM","POS")]),]
nonSynComps$comp_exon_820_GNPm9 <- nonSynComps$comp_exon_820_GNPm9[order(nonSynComps$comp_exon_820_GNPm9$InLine),]
nonSynCompsStats$comp_exon_820_GNPm9_stats <- count(nonSynComps$comp_exon_820_GNPm9, 'InLine')
nonSynCompsStats$comp_exon_820_GNPm9_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_exon_820_GNPm9$Line), sep="_"))
####nonSyn CDS####
nonSynComps$comp_CDS_820_GNPm7 <- rbind(nonSynVCFdata[["820.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                         nonSynVCFdata[["GNPm7.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_CDS_820_GNPm7$Line <- gsub(".CDS.Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_CDS_820_GNPm7$Line)
nonSynComps$comp_CDS_820_GNPm7 <- cbind(nonSynComps$comp_CDS_820_GNPm7, InLine = dupsBetweenGroups(nonSynComps$comp_CDS_820_GNPm7, "Line"))
for (i in seq_along(nonSynComps$comp_CDS_820_GNPm7$InLine)){
  if(nonSynComps$comp_CDS_820_GNPm7$InLine[i] == "TRUE") {
    nonSynComps$comp_CDS_820_GNPm7$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_CDS_820_GNPm7$InLine[i] <- nonSynComps$comp_CDS_820_GNPm7$Line[i]
  }}
nonSynComps$comp_CDS_820_GNPm7 <- nonSynComps$comp_CDS_820_GNPm7[!duplicated(nonSynComps$comp_CDS_820_GNPm7[c("CHROM","POS")]),]
nonSynComps$comp_CDS_820_GNPm7 <- nonSynComps$comp_CDS_820_GNPm7[order(nonSynComps$comp_CDS_820_GNPm7$InLine),]
nonSynCompsStats$comp_CDS_820_GNPm7_stats <- count(nonSynComps$comp_CDS_820_GNPm7, 'InLine')
nonSynCompsStats$comp_CDS_820_GNPm7_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_CDS_820_GNPm7$Line), sep="_"))

nonSynComps$comp_CDS_820_GNPm8 <- rbind(nonSynVCFdata[["820.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                         nonSynVCFdata[["GNPm8.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_CDS_820_GNPm8$Line <- gsub(".CDS.Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_CDS_820_GNPm8$Line)
nonSynComps$comp_CDS_820_GNPm8 <- cbind(nonSynComps$comp_CDS_820_GNPm8, InLine = dupsBetweenGroups(nonSynComps$comp_CDS_820_GNPm8, "Line"))
for (i in seq_along(nonSynComps$comp_CDS_820_GNPm8$InLine)){
  if(nonSynComps$comp_CDS_820_GNPm8$InLine[i] == "TRUE") {
    nonSynComps$comp_CDS_820_GNPm8$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_CDS_820_GNPm8$InLine[i] <- nonSynComps$comp_CDS_820_GNPm8$Line[i]
  }}
nonSynComps$comp_CDS_820_GNPm8 <- nonSynComps$comp_CDS_820_GNPm8[!duplicated(nonSynComps$comp_CDS_820_GNPm8[c("CHROM","POS")]),]
nonSynComps$comp_CDS_820_GNPm8 <- nonSynComps$comp_CDS_820_GNPm8[order(nonSynComps$comp_CDS_820_GNPm8$InLine),]
nonSynCompsStats$comp_CDS_820_GNPm8_stats <- count(nonSynComps$comp_CDS_820_GNPm8, 'InLine')
nonSynCompsStats$comp_CDS_820_GNPm8_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_CDS_820_GNPm8$Line), sep="_"))

nonSynComps$comp_CDS_820_GNPm9 <- rbind(nonSynVCFdata[["820.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                         nonSynVCFdata[["GNPm9.CDS.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_CDS_820_GNPm9$Line <- gsub(".CDS.Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_CDS_820_GNPm9$Line)
nonSynComps$comp_CDS_820_GNPm9 <- cbind(nonSynComps$comp_CDS_820_GNPm9, InLine = dupsBetweenGroups(nonSynComps$comp_CDS_820_GNPm9, "Line"))
for (i in seq_along(nonSynComps$comp_CDS_820_GNPm9$InLine)){
  if(nonSynComps$comp_CDS_820_GNPm9$InLine[i] == "TRUE") {
    nonSynComps$comp_CDS_820_GNPm9$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_CDS_820_GNPm9$InLine[i] <- nonSynComps$comp_CDS_820_GNPm9$Line[i]
  }}
nonSynComps$comp_CDS_820_GNPm9 <- nonSynComps$comp_CDS_820_GNPm9[!duplicated(nonSynComps$comp_CDS_820_GNPm9[c("CHROM","POS")]),]
nonSynComps$comp_CDS_820_GNPm9 <- nonSynComps$comp_CDS_820_GNPm9[order(nonSynComps$comp_CDS_820_GNPm9$InLine),]
nonSynCompsStats$comp_CDS_820_GNPm9_stats <- count(nonSynComps$comp_CDS_820_GNPm9, 'InLine')
nonSynCompsStats$comp_CDS_820_GNPm9_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_CDS_820_GNPm9$Line), sep="_"))
####nonSyn gene####
nonSynComps$comp_gene_820_GNPm7 <- rbind(nonSynVCFdata[["820.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                         nonSynVCFdata[["GNPm7.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_gene_820_GNPm7$Line <- gsub(".gene.Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_gene_820_GNPm7$Line)
nonSynComps$comp_gene_820_GNPm7 <- cbind(nonSynComps$comp_gene_820_GNPm7, InLine = dupsBetweenGroups(nonSynComps$comp_gene_820_GNPm7, "Line"))
for (i in seq_along(nonSynComps$comp_gene_820_GNPm7$InLine)){
  if(nonSynComps$comp_gene_820_GNPm7$InLine[i] == "TRUE") {
    nonSynComps$comp_gene_820_GNPm7$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_gene_820_GNPm7$InLine[i] <- nonSynComps$comp_gene_820_GNPm7$Line[i]
  }}
nonSynComps$comp_gene_820_GNPm7 <- nonSynComps$comp_gene_820_GNPm7[!duplicated(nonSynComps$comp_gene_820_GNPm7[c("CHROM","POS")]),]
nonSynComps$comp_gene_820_GNPm7 <- nonSynComps$comp_gene_820_GNPm7[order(nonSynComps$comp_gene_820_GNPm7$InLine),]
nonSynCompsStats$comp_gene_820_GNPm7_stats <- count(nonSynComps$comp_gene_820_GNPm7, 'InLine')
nonSynCompsStats$comp_gene_820_GNPm7_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_gene_820_GNPm7$Line), sep="_"))

nonSynComps$comp_gene_820_GNPm8 <- rbind(nonSynVCFdata[["820.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                         nonSynVCFdata[["GNPm8.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_gene_820_GNPm8$Line <- gsub(".gene.Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_gene_820_GNPm8$Line)
nonSynComps$comp_gene_820_GNPm8 <- cbind(nonSynComps$comp_gene_820_GNPm8, InLine = dupsBetweenGroups(nonSynComps$comp_gene_820_GNPm8, "Line"))
for (i in seq_along(nonSynComps$comp_gene_820_GNPm8$InLine)){
  if(nonSynComps$comp_gene_820_GNPm8$InLine[i] == "TRUE") {
    nonSynComps$comp_gene_820_GNPm8$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_gene_820_GNPm8$InLine[i] <- nonSynComps$comp_gene_820_GNPm8$Line[i]
  }}
nonSynComps$comp_gene_820_GNPm8 <- nonSynComps$comp_gene_820_GNPm8[!duplicated(nonSynComps$comp_gene_820_GNPm8[c("CHROM","POS")]),]
nonSynComps$comp_gene_820_GNPm8 <- nonSynComps$comp_gene_820_GNPm8[order(nonSynComps$comp_gene_820_GNPm8$InLine),]
nonSynCompsStats$comp_gene_820_GNPm8_stats <- count(nonSynComps$comp_gene_820_GNPm8, 'InLine')
nonSynCompsStats$comp_gene_820_GNPm8_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_gene_820_GNPm8$Line), sep="_"))

nonSynComps$comp_gene_820_GNPm9 <- rbind(nonSynVCFdata[["820.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)],
                                         nonSynVCFdata[["GNPm9.gene.Q60.sorted.filter.ann.vcf"]][,c(1:2,5)])
nonSynComps$comp_gene_820_GNPm9$Line <- gsub(".gene.Q60.sorted.filter.ann.vcf", "", nonSynComps$comp_gene_820_GNPm9$Line)
nonSynComps$comp_gene_820_GNPm9 <- cbind(nonSynComps$comp_gene_820_GNPm9, InLine = dupsBetweenGroups(nonSynComps$comp_gene_820_GNPm9, "Line"))
for (i in seq_along(nonSynComps$comp_gene_820_GNPm9$InLine)){
  if(nonSynComps$comp_gene_820_GNPm9$InLine[i] == "TRUE") {
    nonSynComps$comp_gene_820_GNPm9$InLine[i] <- "BOTH"
  }  else {
    nonSynComps$comp_gene_820_GNPm9$InLine[i] <- nonSynComps$comp_gene_820_GNPm9$Line[i]
  }}
nonSynComps$comp_gene_820_GNPm9 <- nonSynComps$comp_gene_820_GNPm9[!duplicated(nonSynComps$comp_gene_820_GNPm9[c("CHROM","POS")]),]
nonSynComps$comp_gene_820_GNPm9 <- nonSynComps$comp_gene_820_GNPm9[order(nonSynComps$comp_gene_820_GNPm9$InLine),]
nonSynCompsStats$comp_gene_820_GNPm9_stats <- count(nonSynComps$comp_gene_820_GNPm9, 'InLine')
nonSynCompsStats$comp_gene_820_GNPm9_stats$Comparison <- capture.output(cat(unique(nonSynComps$comp_gene_820_GNPm9$Line), sep="_"))

#####make comp lists#####
Syn_lists <- list()
Syn_lists$Syn_GNPm7_not820 <- subset(SynComps$comp_820_GNPm7[,1:2], SynComps$comp_820_GNPm7$InLine == "GNPm7")
Syn_lists$Syn_GNPm8_not820 <- subset(SynComps$comp_820_GNPm8[,1:2], SynComps$comp_820_GNPm8$InLine == "GNPm8")
Syn_lists$Syn_GNPm9_not820 <- subset(SynComps$comp_820_GNPm9[,1:2], SynComps$comp_820_GNPm9$InLine == "GNPm9")

Syn_lists$Syn_exon_GNPm7_not820 <- subset(SynComps$comp_exon_820_GNPm7[,1:2], SynComps$comp_exon_820_GNPm7$InLine == "GNPm7")
Syn_lists$Syn_exon_GNPm8_not820 <- subset(SynComps$comp_exon_820_GNPm8[,1:2], SynComps$comp_exon_820_GNPm8$InLine == "GNPm8")
Syn_lists$Syn_exon_GNPm9_not820 <- subset(SynComps$comp_exon_820_GNPm9[,1:2], SynComps$comp_exon_820_GNPm9$InLine == "GNPm9")

Syn_lists$Syn_CDS_GNPm7_not820 <- subset(SynComps$comp_CDS_820_GNPm7[,1:2], SynComps$comp_CDS_820_GNPm7$InLine == "GNPm7")
Syn_lists$Syn_CDS_GNPm8_not820 <- subset(SynComps$comp_CDS_820_GNPm8[,1:2], SynComps$comp_CDS_820_GNPm8$InLine == "GNPm8")
Syn_lists$Syn_CDS_GNPm9_not820 <- subset(SynComps$comp_CDS_820_GNPm9[,1:2], SynComps$comp_CDS_820_GNPm9$InLine == "GNPm9")

Syn_lists$Syn_gene_GNPm7_not820 <- subset(SynComps$comp_gene_820_GNPm7[,1:2], SynComps$comp_gene_820_GNPm7$InLine == "GNPm7")
Syn_lists$Syn_gene_GNPm8_not820 <- subset(SynComps$comp_gene_820_GNPm8[,1:2], SynComps$comp_gene_820_GNPm8$InLine == "GNPm8")
Syn_lists$Syn_gene_GNPm9_not820 <- subset(SynComps$comp_gene_820_GNPm9[,1:2], SynComps$comp_gene_820_GNPm9$InLine == "GNPm9")

nonSyn_lists <- list()
nonSyn_lists$nonSyn_GNPm7_not820 <- subset(nonSynComps$comp_820_GNPm7[,1:2], nonSynComps$comp_820_GNPm7$InLine == "GNPm7")
nonSyn_lists$nonSyn_GNPm8_not820 <- subset(nonSynComps$comp_820_GNPm8[,1:2], nonSynComps$comp_820_GNPm8$InLine == "GNPm8")
nonSyn_lists$nonSyn_GNPm9_not820 <- subset(nonSynComps$comp_820_GNPm9[,1:2], nonSynComps$comp_820_GNPm9$InLine == "GNPm9")

nonSyn_lists$nonSyn_exon_GNPm7_not820 <- subset(nonSynComps$comp_exon_820_GNPm7[,1:2], nonSynComps$comp_exon_820_GNPm7$InLine == "GNPm7")
nonSyn_lists$nonSyn_exon_GNPm8_not820 <- subset(nonSynComps$comp_exon_820_GNPm8[,1:2], nonSynComps$comp_exon_820_GNPm8$InLine == "GNPm8")
nonSyn_lists$nonSyn_exon_GNPm9_not820 <- subset(nonSynComps$comp_exon_820_GNPm9[,1:2], nonSynComps$comp_exon_820_GNPm9$InLine == "GNPm9")

nonSyn_lists$nonSyn_CDS_GNPm7_not820 <- subset(nonSynComps$comp_CDS_820_GNPm7[,1:2], nonSynComps$comp_CDS_820_GNPm7$InLine == "GNPm7")
nonSyn_lists$nonSyn_CDS_GNPm8_not820 <- subset(nonSynComps$comp_CDS_820_GNPm8[,1:2], nonSynComps$comp_CDS_820_GNPm8$InLine == "GNPm8")
nonSyn_lists$nonSyn_CDS_GNPm9_not820 <- subset(nonSynComps$comp_CDS_820_GNPm9[,1:2], nonSynComps$comp_CDS_820_GNPm9$InLine == "GNPm9")

nonSyn_lists$nonSyn_gene_GNPm7_not820 <- subset(nonSynComps$comp_gene_820_GNPm7[,1:2], nonSynComps$comp_gene_820_GNPm7$InLine == "GNPm7")
nonSyn_lists$nonSyn_gene_GNPm8_not820 <- subset(nonSynComps$comp_gene_820_GNPm8[,1:2], nonSynComps$comp_gene_820_GNPm8$InLine == "GNPm8")
nonSyn_lists$nonSyn_gene_GNPm9_not820 <- subset(nonSynComps$comp_gene_820_GNPm9[,1:2], nonSynComps$comp_gene_820_GNPm9$InLine == "GNPm9")

col_names <- c("CHROM","POS")
lapply(Syn_lists, setNames, nm = col_names)
lapply(nonSyn_lists, setNames, nm = col_names)

####export comp lists#####
for (i in names(Syn_lists)){
  write.csv(Syn_lists[i], file = paste0("J:/III/Waters/Group Members/Mallory/GenomeProjects/Scripts/files/",i,".csv"))
}

for (i in names(nonSyn_lists)){
  write.table(nonSyn_lists[i], file = paste0("J:/III/Waters/Group Members/Mallory/GenomeProjects/Scripts/files/",i,".csv"))
}
#---> keep working on fixing column names here


#####final save#####
save.image("J:/III/Waters/Group Members/Mallory/R/VCFImportComparisons.RData")
