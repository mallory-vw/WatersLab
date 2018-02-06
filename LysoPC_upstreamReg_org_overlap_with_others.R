#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/LysoPC_upstreamReg_org_overlap_with_others.RData")
library("data.table")

#####set WD#####
setwd("J:/III/Waters/Group Members/Mallory/LysoPC_project")

#####load data#####
allgenes_Niggi <- read.table("all_LysoPC.txt")
  allgenes_Niggi <- allgenes_Niggi[,1]
upreg_Niggi <- read.table("UpReg_LysoPC.txt")
  upreg_Niggi <- upreg_Niggi[,1]
downreg_Niggi <- read.table("DownReg_LysoPC.txt")
  downreg_Niggi <- downreg_Niggi[,1]
background_Niggi <- read.table("Background_LysoPC.txt")
  background_Niggi <- background_Niggi[,1]
background_nodownreg_Niggi <- read.table("Background_noDownReg_LysoPC.txt")
  background_nodownreg_Niggi <- background_nodownreg_Niggi[,1]


allgenes_Niggi_noDups <- read.table("all_LysoPC_noDups.txt")
  allgenes_Niggi_noDups <- allgenes_Niggi_noDups[,1]
upreg_Niggi_noDups <- read.table("UpReg_LysoPC_noDups.txt")
  upreg_Niggi_noDups <- upreg_Niggi_noDups[,1]
downreg_Niggi_noDups <- read.table("DownReg_LysoPC_noDups.txt")
  downreg_Niggi_noDups <- downreg_Niggi_noDups[,1]
background_Niggi_noDups <- read.table("Background_LysoPC_noDups.txt")
  background_Niggi_noDups <- background_Niggi_noDups[,1]
background_nodownreg_Niggi_noDups <- read.table("Background_noDownReg_LysoPC_noDups.txt")
  background_nodownreg_Niggi_noDups <- background_nodownreg_Niggi_noDups[,1]

  
#remove split genes from all noDups lists - these are genes that are split into different transcripts (".1", ".2",etc) but the 2 transcripts 
#are in different genes lists (upreg vs downreg, etc)
# x[x != "b"]; # without elements that are "b"
# x[!x %in% B]
to_remove <- c("PF3D7_0314100","PF3D7_0613000", "PF3D7_0529400", "PF3D7_0923800", "PF3D7_0208100")
  
allgenes_Niggi_noDups_clean <- allgenes_Niggi_noDups[!allgenes_Niggi_noDups %in% to_remove]
upreg_Niggi_noDups_clean <- upreg_Niggi_noDups[!upreg_Niggi_noDups %in% to_remove]
downreg_Niggi_noDups_clean <- downreg_Niggi_noDups[!downreg_Niggi_noDups %in% to_remove]
background_Niggi_noDups_clean <- background_Niggi_noDups[!background_Niggi_noDups %in% to_remove]
background_nodownreg_Niggi_noDups_clean <- background_nodownreg_Niggi_noDups[!background_nodownreg_Niggi_noDups %in% to_remove]  
  


#####2kb load data to calculate upstream_size#####
#load in the gene info 
#the duplicates with ".1" and ".2", etc, have been removed in Excel
Pfalciparum_GFF3_allgene_2kb_strand <- read.csv("PlasmoDB_28_GeneStart_GeneStop_strand_calc.csv", header = T)

#chrom sizes
Pfalciparum_chrom_coord <- read.csv("Pfalciparum_chromosome_coord.csv")


###start with the beginning and ending of CHROM
#beginning if the first gene is on the + strand, end if the last gene is on the - strand

###label the start of each chromosome
for (i in seq_along(Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[2:5688])) {
  if (identical(Pfalciparum_GFF3_allgene_2kb_strand$SEQID[i], Pfalciparum_GFF3_allgene_2kb_strand$SEQID[i-1])) {
    Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] <- "calculate"
  } else {
    Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] <- "START"
  }
}
###label the end of each chromosome
for (i in seq_along(Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[2:5688])) {
  if (identical(Pfalciparum_GFF3_allgene_2kb_strand$SEQID[i], Pfalciparum_GFF3_allgene_2kb_strand$SEQID[i+1]) && (Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] == "calculate")) {
    Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] <- "calculate"
  } else if (Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] != "START") {
    Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] <- "END"
  }
}


###at start of CHROM is gene is on + strand, calculate upstream dist
#if gene is on - strand will go through normal calculation
for (i in seq_along(Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] == "START") && (Pfalciparum_GFF3_allgene_2kb_strand$STRAND[i] == "+")) {
    Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] <- Pfalciparum_GFF3_allgene_2kb_strand$START[i] - 1
  } 
}
###at end of CHROM if gene is on - strand, need to calculate upstream dist based on size of CHROM
#if gene is on + strand will go through normal calc
for (i in seq_along(Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] == "END") && (Pfalciparum_GFF3_allgene_2kb_strand$STRAND[i] == "-")) {
    Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] <- (Pfalciparum_chrom_coord$END[Pfalciparum_chrom_coord$CHROM == Pfalciparum_GFF3_allgene_2kb_strand$SEQID[i]]) - Pfalciparum_GFF3_allgene_2kb_strand$END[i]
  } 
}

#need to consider the other genes <2kb away from the end/start of each CHROM
# maybe just do the calculation for all of them from the CHROM ends, then replace everything else with 2001 - first just calculate the UPDIST from end of CHROM, then make UP_DIST_EDIT 2001

#for all genes on +ve strand, calculate UPDIST from beginning of CHROM -- need to remember that UPDIST needs to be one bigger to accomodate for 0-based coord
#example, to get bases 1,2, in 1 -based cood, UPSTART = 1, UPEND = 2, 0-based coord UPSTART = 0, UPEND = 2
# hence why I'm using 2001, not 2000 in UPDIST_EDIT
for (i in seq_along(Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] == "calculate") && (Pfalciparum_GFF3_allgene_2kb_strand$STRAND[i] == "+")) {
    Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] <- (Pfalciparum_GFF3_allgene_2kb_strand$START[i] - Pfalciparum_chrom_coord$START[Pfalciparum_chrom_coord$CHROM == Pfalciparum_GFF3_allgene_2kb_strand$SEQID[i]])+1
  } 
}

#for all genes on -ve strand, calculate UPDIST from end of CHROM, then for any updist > 2000, change it to 2001
for (i in seq_along(Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] == "calculate") && (Pfalciparum_GFF3_allgene_2kb_strand$STRAND[i] == "-")) {
    Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] <- ((Pfalciparum_chrom_coord$END[Pfalciparum_chrom_coord$CHROM == Pfalciparum_GFF3_allgene_2kb_strand$SEQID[i]]) - Pfalciparum_GFF3_allgene_2kb_strand$END[i])+1
  } 
}

#then move on to populating the rest of the values, ignoring the rows that have already been filled in
#change all UP_DIST values that aren't numbers to "2001" since I want 2kb upstream
for (i in seq_along(Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] == "calculate") ||
      (Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] == "START") ||
      (Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] == "END") ) {
    Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST[i] <- "2001"
  }
}
##convert values to integers
Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST <- as.integer(Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST)


#####2kb organize upstream cutoff values#####

#generate edited Up_DIST to cut off at 2KB for all UP_DIST that are already > 2000
Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST_EDIT <- ifelse(Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST > 2001, 2001, Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST)
  # remove genes with UP_DIST <= 2 (no distance between genes)
  # UPDATE, first check to see if there are any
  range(Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST_EDIT)
  # UPDATE - do not need to to this, no genes will have an UPDIST less than 2

##all genes in Niggi's list
allgenes_2kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_2kb_strand, Pfalciparum_GFF3_allgene_2kb_strand$GENEID %in% allgenes_Niggi_noDups_clean)
    #we lose 7 genes that are in allgenes_Niggi_noDups_clean (5680) but aren't in the PlasmoDB list (since allgenes_2kb_strand_Niggi_noDups_GFF3_info has 5673)
    #we lose 14 genes that are in the PlasmoDB list (5687) but not in allgenes_Niggi_noDups_clean, end up with 5673 in allgenes_2kb_strand_Niggi_noDups_GFF3_info
    
    #use allgenes_2kb_strand_Niggi_noDups_GFF3_info GENEID list moving forwards
    allgenes_2kb_strand_Niggi_noDups <- allgenes_2kb_strand_Niggi_noDups_GFF3_info$GENEID
    
##background
background_2kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_2kb_strand, Pfalciparum_GFF3_allgene_2kb_strand$GENEID %in% background_Niggi_noDups_clean)
  #we lose 7 genes that are in background_Niggi_noDups_clean (5341) but aren't in the PlasmoDB list (background_2kb_strand_Niggi_noDups_GFF3_info is 5334)

  #use background_2kb_strand_Niggi_noDups_GFF3_info GENEID list moving forwards  
  background_2kb_strand_Niggi_noDups <- background_2kb_strand_Niggi_noDups_GFF3_info$GENEID

##background no downreg  
background_nodownreg_2kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_2kb_strand, Pfalciparum_GFF3_allgene_2kb_strand$GENEID %in% background_nodownreg_Niggi_noDups_clean)
  #we lose 7 genes that are in background_nodownreg_Niggi_noDups_clean (5299) but aren't in the PlasmoDB list (background_nodownreg_strand_Niggi_noDups_GFF3_info is 5292)
    
  #use background_2kb_strand_Niggi_noDups_GFF3_info GENEID list moving forwards  
  background_nodownreg_2kb_strand_Niggi_noDups <- background_nodownreg_2kb_strand_Niggi_noDups_GFF3_info$GENEID

##upreg and downreg
  #are all pos, also none missing from either Niggi's list or the PlasmoDB list
  upreg_2kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_2kb_strand, Pfalciparum_GFF3_allgene_2kb_strand$GENEID %in% upreg_Niggi_noDups_clean)
  downreg_2kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_2kb_strand, Pfalciparum_GFF3_allgene_2kb_strand$GENEID %in% downreg_Niggi_noDups_clean)    
    
#the next step is to figure out the beginning and end of the upstream region for both the + strand genes (5' of the gene start) and
#the -ve strand genes (based on coordinate system, 3' of the gene end)
  
# +ve gene --> UPSTART = GENE_START - UP_DIST_EDIT 
#UPEND = GENE_START - 1
# -ve gene --> UPSTART = GENE_END 
#UPEND = GENE_END + (UP_DIST_EDIT - 1)
  Pfalciparum_upstream_regions_2kb_strand <- data.frame(matrix(NA, nrow = 5687))
  Pfalciparum_upstream_regions_2kb_strand$CHROM <-  Pfalciparum_GFF3_allgene_2kb_strand$SEQID
  Pfalciparum_upstream_regions_2kb_strand$GENEID <- Pfalciparum_GFF3_allgene_2kb_strand$GENEID
  Pfalciparum_upstream_regions_2kb_strand$STRAND <- Pfalciparum_GFF3_allgene_2kb_strand$STRAND
  Pfalciparum_upstream_regions_2kb_strand$GENE_START <- Pfalciparum_GFF3_allgene_2kb_strand$START
  Pfalciparum_upstream_regions_2kb_strand$GENE_END <- Pfalciparum_GFF3_allgene_2kb_strand$END
  Pfalciparum_upstream_regions_2kb_strand$UP_DIST_EDIT <- Pfalciparum_GFF3_allgene_2kb_strand$UP_DIST_EDIT
  Pfalciparum_upstream_regions_2kb_strand$UPSTART <- NA
  Pfalciparum_upstream_regions_2kb_strand$UPEND <- NA
  Pfalciparum_upstream_regions_2kb_strand <- Pfalciparum_upstream_regions_2kb_strand[,-1]
  
  
  for (i in seq_along(Pfalciparum_upstream_regions_2kb_strand$GENEID)) {
    if (Pfalciparum_upstream_regions_2kb_strand$STRAND[i] == "+") {
      Pfalciparum_upstream_regions_2kb_strand$UPSTART[i] <- Pfalciparum_upstream_regions_2kb_strand$GENE_START[i] - Pfalciparum_upstream_regions_2kb_strand$UP_DIST_EDIT[i]
      Pfalciparum_upstream_regions_2kb_strand$UPEND[i] <- Pfalciparum_upstream_regions_2kb_strand$GENE_START[i] - 1
    }
  }
  
  for (i in seq_along(Pfalciparum_upstream_regions_2kb_strand$GENEID)) {
    if (Pfalciparum_upstream_regions_2kb_strand$STRAND[i] == "-") {
      Pfalciparum_upstream_regions_2kb_strand$UPSTART[i] <- Pfalciparum_upstream_regions_2kb_strand$GENE_END[i]
      Pfalciparum_upstream_regions_2kb_strand$UPEND[i] <- Pfalciparum_upstream_regions_2kb_strand$GENE_END[i] + (Pfalciparum_upstream_regions_2kb_strand$UP_DIST_EDIT[i] - 1)
    }
  }

  
  
  
  
  
write.csv(Pfalciparum_upstream_regions_2kb_strand, "Pfalciparum_upstream_regions_2kb_strand.csv")    
 
#####2kbbsplit upstream regions into the specific gene lists#####
# allgenes_2kb_strand_Niggi_noDups_GFF3_info
# background_2kb_strand_Niggi_noDups_GFF3_info
# background_nodownreg_2kb_strand_Niggi_noDups_GFF3_info
# upreg_2kb_strand_Niggi_noDups_GFF3_info
# downreg_2kb_strand_Niggi_noDups_GFF3_info
# 
# Pfalciparum_upstream_regions_2kb_strand

##allgenes
allgenes_2kb_strand_Niggi_noDups_upstream_regions <- subset(Pfalciparum_upstream_regions_2kb_strand, Pfalciparum_upstream_regions_2kb_strand$GENEID %in% allgenes_2kb_strand_Niggi_noDups_GFF3_info$GENEID)

##background
background_2kb_strand_Niggi_noDups_upstream_regions <- subset(Pfalciparum_upstream_regions_2kb_strand, Pfalciparum_upstream_regions_2kb_strand$GENEID %in% background_2kb_strand_Niggi_noDups_GFF3_info$GENEID)

##background nodownreg
background_2kb_nodownreg_strand_Niggi_noDups_upstream_regions <- subset(Pfalciparum_upstream_regions_2kb_strand, Pfalciparum_upstream_regions_2kb_strand$GENEID %in% background_nodownreg_2kb_strand_Niggi_noDups_GFF3_info$GENEID)

##upreg and downreg
upreg_2kb_strand_Niggi_noDups_upstream_regions <- subset(Pfalciparum_upstream_regions_2kb_strand, Pfalciparum_upstream_regions_2kb_strand$GENEID %in% upreg_2kb_strand_Niggi_noDups_GFF3_info$GENEID)
downreg_2kb_strand_Niggi_noDups_upstream_regions  <- subset(Pfalciparum_upstream_regions_2kb_strand, Pfalciparum_upstream_regions_2kb_strand$GENEID %in% downreg_2kb_strand_Niggi_noDups_GFF3_info$GENEID)

#####2kb save lists of GENEids for each category#####
write.csv(upreg_2kb_strand_Niggi_noDups_upstream_regions$GENEID, "upreg_2kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
write.csv(downreg_2kb_strand_Niggi_noDups_upstream_regions$GENEID, "downreg_2kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
write.csv(background_2kb_strand_Niggi_noDups_upstream_regions$GENEID, "background_2kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
write.csv(background_2kb_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID, "background_2kb_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
write.csv(allgenes_2kb_strand_Niggi_noDups_upstream_regions$GENEID, "allgenes_2kb_strand_Niggi_noDups_upstream_regions_genelist.csv")

#####2kb create and download text for bedtools BED files#####

#same concept as previous - use the coordinates in the upstream_regions 
#dataframes to create BED files for each subset of genes
# BED columns: CHROM START END NAME

# ##lists to make
# allgenes_2kb_strand_Niggi_noDups_upstream_regions
# upreg_2kb_strand_Niggi_noDups_upstream_regions
# downreg_2kb_strand_Niggi_noDups_upstream_regions
# background_2kb_strand_Niggi_noDups_upstream_regions
# background_2kb_nodownreg_strand_Niggi_noDups_upstream_regions

# BED columns: CHROM START END NAME

upreg_2kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 339))
  upreg_2kb_strand_Niggi_noDups_download_BED$CHROM <- upreg_2kb_strand_Niggi_noDups_upstream_regions$CHROM
  #ISSUE HERE - missing CHROM in upreg_2kb_strand_Niggi_noDups_upstream_regions
  upreg_2kb_strand_Niggi_noDups_download_BED$START <- upreg_2kb_strand_Niggi_noDups_upstream_regions$UPSTART
  upreg_2kb_strand_Niggi_noDups_download_BED$END <- upreg_2kb_strand_Niggi_noDups_upstream_regions$UPEND
  upreg_2kb_strand_Niggi_noDups_download_BED$NAME <- upreg_2kb_strand_Niggi_noDups_upstream_regions$GENEID
  upreg_2kb_strand_Niggi_noDups_download_BED <- upreg_2kb_strand_Niggi_noDups_download_BED[,2:5]
    write.csv(upreg_2kb_strand_Niggi_noDups_download_BED, "upreg_2kb_strand_Niggi_noDups_download_BED.csv")

downreg_2kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 42))
  downreg_2kb_strand_Niggi_noDups_download_BED$CHROM <- downreg_2kb_strand_Niggi_noDups_upstream_regions$CHROM
  downreg_2kb_strand_Niggi_noDups_download_BED$START <- downreg_2kb_strand_Niggi_noDups_upstream_regions$UPSTART
  downreg_2kb_strand_Niggi_noDups_download_BED$END <- downreg_2kb_strand_Niggi_noDups_upstream_regions$UPEND
  downreg_2kb_strand_Niggi_noDups_download_BED$NAME <- downreg_2kb_strand_Niggi_noDups_upstream_regions$GENEID
  downreg_2kb_strand_Niggi_noDups_download_BED <- downreg_2kb_strand_Niggi_noDups_download_BED[,2:5]
    write.csv(downreg_2kb_strand_Niggi_noDups_download_BED, "downreg_2kb_strand_Niggi_noDups_download_BED.csv")

background_2kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5334))
  background_2kb_strand_Niggi_noDups_download_BED$CHROM <- background_2kb_strand_Niggi_noDups_upstream_regions$CHROM
  background_2kb_strand_Niggi_noDups_download_BED$START <- background_2kb_strand_Niggi_noDups_upstream_regions$UPSTART
  background_2kb_strand_Niggi_noDups_download_BED$END <- background_2kb_strand_Niggi_noDups_upstream_regions$UPEND
  background_2kb_strand_Niggi_noDups_download_BED$NAME <- background_2kb_strand_Niggi_noDups_upstream_regions$GENEID
  background_2kb_strand_Niggi_noDups_download_BED <- background_2kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(background_2kb_strand_Niggi_noDups_download_BED, "background_2kb_strand_Niggi_noDups_download_BED.CSV")

background_2kb_nodownreg_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5292))  
  background_2kb_nodownreg_strand_Niggi_noDups_download_BED$CHROM <- background_2kb_nodownreg_strand_Niggi_noDups_upstream_regions$CHROM
  background_2kb_nodownreg_strand_Niggi_noDups_download_BED$START <- background_2kb_nodownreg_strand_Niggi_noDups_upstream_regions$UPSTART
  background_2kb_nodownreg_strand_Niggi_noDups_download_BED$END <- background_2kb_nodownreg_strand_Niggi_noDups_upstream_regions$UPEND
  background_2kb_nodownreg_strand_Niggi_noDups_download_BED$NAME <- background_2kb_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID
  background_2kb_nodownreg_strand_Niggi_noDups_download_BED <- background_2kb_nodownreg_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(background_2kb_nodownreg_strand_Niggi_noDups_download_BED, "background_2kb_nodownreg_strand_Niggi_noDups_download_BED.csv")

allgenes_2kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5673))
  allgenes_2kb_strand_Niggi_noDups_download_BED$CHROM <- allgenes_2kb_strand_Niggi_noDups_upstream_regions$CHROM
  allgenes_2kb_strand_Niggi_noDups_download_BED$START <- allgenes_2kb_strand_Niggi_noDups_upstream_regions$UPSTART
  allgenes_2kb_strand_Niggi_noDups_download_BED$END <- allgenes_2kb_strand_Niggi_noDups_upstream_regions$UPEND
  allgenes_2kb_strand_Niggi_noDups_download_BED$NAME <- allgenes_2kb_strand_Niggi_noDups_upstream_regions$GENEID
  allgenes_2kb_strand_Niggi_noDups_download_BED <- allgenes_2kb_strand_Niggi_noDups_download_BED[,2:5]
    write.csv(allgenes_2kb_strand_Niggi_noDups_download_BED, "allgenes_2kb_strand_Niggi_noDups_download_BED.csv")





    


    
#####1kb load data to calculate upstream_size#####
#load in the gene info 
#the duplicates with ".1" and ".2", etc, have been removed in Excel
Pfalciparum_GFF3_allgene_1kb_strand <- read.csv("PlasmoDB_28_GeneStart_GeneStop_strand_calc.csv", header = T)
    

###start with the beginning and ending of CHROM
#beginning if the first gene is on the + strand, end if the last gene is on the - strand
    
###label the start of each chromosome
for (i in seq_along(Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[2:5688])) {
  if (identical(Pfalciparum_GFF3_allgene_1kb_strand$SEQID[i], Pfalciparum_GFF3_allgene_1kb_strand$SEQID[i-1])) {
    Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] <- "calculate"
  } else {
    Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] <- "START"
  }
}
###label the end of each chromosome
for (i in seq_along(Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[2:5688])) {
  if (identical(Pfalciparum_GFF3_allgene_1kb_strand$SEQID[i], Pfalciparum_GFF3_allgene_1kb_strand$SEQID[i+1]) && (Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] == "calculate")) {
    Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] <- "calculate"
  } else if (Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] != "START") {
    Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] <- "END"
  }
}
    
    
###at start of CHROM is gene is on + strand, calculate upstream dist
  #if gene is on - strand will go through normal calculation
for (i in seq_along(Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] == "START") && (Pfalciparum_GFF3_allgene_1kb_strand$STRAND[i] == "+")) {
    Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] <- Pfalciparum_GFF3_allgene_1kb_strand$START[i] - 1
  } 
}
###at end of CHROM if gene is on - strand, need to calculate upstream dist based on size of CHROM
  #if gene is on + strand will go through normal calc
for (i in seq_along(Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] == "END") && (Pfalciparum_GFF3_allgene_1kb_strand$STRAND[i] == "-")) {
    Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] <- (Pfalciparum_chrom_coord$END[Pfalciparum_chrom_coord$CHROM == Pfalciparum_GFF3_allgene_1kb_strand$SEQID[i]]) - Pfalciparum_GFF3_allgene_1kb_strand$END[i]
  } 
}
    
#need to consider the other genes <1kb away from the end/start of each CHROM
# maybe just do the calculation for all of them from the CHROM ends, then replace everything else with 2001 - first just calculate the UPDIST from end of CHROM, then make UP_DIST_EDIT 2001

#for all genes on +ve strand, calculate UPDIST from beginning of CHROM -- need to remember that UPDIST needs to be one bigger to accomodate for 0-based coord
#example, to get bases 1,2, in 1 -based cood, UPSTART = 1, UPEND = 2, 0-based coord UPSTART = 0, UPEND = 2
# hence why I'm using 1001, not 1000 in UPDIST_EDIT
for (i in seq_along(Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] == "calculate") && (Pfalciparum_GFF3_allgene_1kb_strand$STRAND[i] == "+")) {
    Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] <- (Pfalciparum_GFF3_allgene_1kb_strand$START[i] - Pfalciparum_chrom_coord$START[Pfalciparum_chrom_coord$CHROM == Pfalciparum_GFF3_allgene_1kb_strand$SEQID[i]])+1
  } 
}
    
#for all genes on -ve strand, calculate UPDIST from end of CHROM, then for any updist > 1000, change it to 1001
for (i in seq_along(Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] == "calculate") && (Pfalciparum_GFF3_allgene_1kb_strand$STRAND[i] == "-")) {
    Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] <- ((Pfalciparum_chrom_coord$END[Pfalciparum_chrom_coord$CHROM == Pfalciparum_GFF3_allgene_1kb_strand$SEQID[i]]) - Pfalciparum_GFF3_allgene_1kb_strand$END[i])+1
  } 
}
    
#then move on to populating the rest of the values, ignoring the rows that have already been filled in
#change all UP_DIST values that aren't numbers to "2001" since I want 1kb upstream
for (i in seq_along(Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] == "calculate") ||
      (Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] == "START") ||
      (Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] == "END") ) {
    Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST[i] <- "1001"
  }
}
##convert values to integers
Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST <- as.integer(Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST)
    
    
#####1kb organize upstream cutoff values#####
    
#generate edited Up_DIST to cut off at 1kb for all UP_DIST that are already > 1000
Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST_EDIT <- ifelse(Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST > 1001, 1001, Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST)
# remove genes with UP_DIST <= 2 (no distance between genes)
# UPDATE, first check to see if there are any
range(Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST_EDIT)
# UPDATE - do not need to to this, no genes will have an UPDIST less than 2
    
##all genes in Niggi's list
    allgenes_1kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_1kb_strand, Pfalciparum_GFF3_allgene_1kb_strand$GENEID %in% allgenes_Niggi_noDups_clean)
    #we lose 7 genes that are in allgenes_Niggi_noDups_clean (5680) but aren't in the PlasmoDB list (since allgenes_1kb_strand_Niggi_noDups_GFF3_info has 5673)
    #we lose 14 genes that are in the PlasmoDB list (5687) but not in allgenes_Niggi_noDups_clean, end up with 5673 in allgenes_1kb_strand_Niggi_noDups_GFF3_info
    
#use allgenes_1kb_strand_Niggi_noDups_GFF3_info GENEID list moving forwards
    allgenes_1kb_strand_Niggi_noDups <- allgenes_1kb_strand_Niggi_noDups_GFF3_info$GENEID
    
##background
    background_1kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_1kb_strand, Pfalciparum_GFF3_allgene_1kb_strand$GENEID %in% background_Niggi_noDups_clean)
    #we lose 7 genes that are in background_Niggi_noDups_clean (5341) but aren't in the PlasmoDB list (background_1kb_strand_Niggi_noDups_GFF3_info is 5334)
    
    #use background_1kb_strand_Niggi_noDups_GFF3_info GENEID list moving forwards  
    background_1kb_strand_Niggi_noDups <- background_1kb_strand_Niggi_noDups_GFF3_info$GENEID
    
##background no downreg  
    background_nodownreg_1kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_1kb_strand, Pfalciparum_GFF3_allgene_1kb_strand$GENEID %in% background_nodownreg_Niggi_noDups_clean)
    #we lose 7 genes that are in background_nodownreg_Niggi_noDups_clean (5299) but aren't in the PlasmoDB list (background_nodownreg_strand_Niggi_noDups_GFF3_info is 5292)
    
    #use background_1kb_strand_Niggi_noDups_GFF3_info GENEID list moving forwards  
    background_nodownreg_1kb_strand_Niggi_noDups <- background_nodownreg_1kb_strand_Niggi_noDups_GFF3_info$GENEID
    
##upreg and downreg
    #are all pos, also none missing from either Niggi's list or the PlasmoDB list
    upreg_1kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_1kb_strand, Pfalciparum_GFF3_allgene_1kb_strand$GENEID %in% upreg_Niggi_noDups_clean)
    downreg_1kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_1kb_strand, Pfalciparum_GFF3_allgene_1kb_strand$GENEID %in% downreg_Niggi_noDups_clean)    
    
    #the next step is to figure out the beginning and end of the upstream region for both the + strand genes (5' of the gene start) and
    #the -ve strand genes (based on coordinate system, 3' of the gene end)
    
    # +ve gene --> UPSTART = GENE_START - UP_DIST_EDIT 
    #UPEND = GENE_START - 1
    # -ve gene --> UPSTART = GENE_END 
    #UPEND = GENE_END + (UP_DIST_EDIT - 1)
    Pfalciparum_upstream_regions_1kb_strand <- data.frame(matrix(NA, nrow = 5687))
    Pfalciparum_upstream_regions_1kb_strand$CHROM <-  Pfalciparum_GFF3_allgene_1kb_strand$SEQID
    Pfalciparum_upstream_regions_1kb_strand$GENEID <- Pfalciparum_GFF3_allgene_1kb_strand$GENEID
    Pfalciparum_upstream_regions_1kb_strand$STRAND <- Pfalciparum_GFF3_allgene_1kb_strand$STRAND
    Pfalciparum_upstream_regions_1kb_strand$GENE_START <- Pfalciparum_GFF3_allgene_1kb_strand$START
    Pfalciparum_upstream_regions_1kb_strand$GENE_END <- Pfalciparum_GFF3_allgene_1kb_strand$END
    Pfalciparum_upstream_regions_1kb_strand$UP_DIST_EDIT <- Pfalciparum_GFF3_allgene_1kb_strand$UP_DIST_EDIT
    Pfalciparum_upstream_regions_1kb_strand$UPSTART <- NA
    Pfalciparum_upstream_regions_1kb_strand$UPEND <- NA
    Pfalciparum_upstream_regions_1kb_strand <- Pfalciparum_upstream_regions_1kb_strand[,-1]
    
    
    for (i in seq_along(Pfalciparum_upstream_regions_1kb_strand$GENEID)) {
      if (Pfalciparum_upstream_regions_1kb_strand$STRAND[i] == "+") {
        Pfalciparum_upstream_regions_1kb_strand$UPSTART[i] <- Pfalciparum_upstream_regions_1kb_strand$GENE_START[i] - Pfalciparum_upstream_regions_1kb_strand$UP_DIST_EDIT[i]
        Pfalciparum_upstream_regions_1kb_strand$UPEND[i] <- Pfalciparum_upstream_regions_1kb_strand$GENE_START[i] - 1
      }
    }
    
    for (i in seq_along(Pfalciparum_upstream_regions_1kb_strand$GENEID)) {
      if (Pfalciparum_upstream_regions_1kb_strand$STRAND[i] == "-") {
        Pfalciparum_upstream_regions_1kb_strand$UPSTART[i] <- Pfalciparum_upstream_regions_1kb_strand$GENE_END[i]
        Pfalciparum_upstream_regions_1kb_strand$UPEND[i] <- Pfalciparum_upstream_regions_1kb_strand$GENE_END[i] + (Pfalciparum_upstream_regions_1kb_strand$UP_DIST_EDIT[i] - 1)
      }
    }
    
    
    
    
    
    
    write.csv(Pfalciparum_upstream_regions_1kb_strand, "Pfalciparum_upstream_regions_1kb_strand.csv")    
    
#####1kb split upstream regions into the specific gene lists#####
    # allgenes_1kb_strand_Niggi_noDups_GFF3_info
    # background_1kb_strand_Niggi_noDups_GFF3_info
    # background_nodownreg_1kb_strand_Niggi_noDups_GFF3_info
    # upreg_1kb_strand_Niggi_noDups_GFF3_info
    # downreg_1kb_strand_Niggi_noDups_GFF3_info
    # 
    # Pfalciparum_upstream_regions_1kb_strand
    
    ##allgenes
    allgenes_1kb_strand_Niggi_noDups_upstream_regions <- subset(Pfalciparum_upstream_regions_1kb_strand, Pfalciparum_upstream_regions_1kb_strand$GENEID %in% allgenes_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    
    ##background
    background_1kb_strand_Niggi_noDups_upstream_regions <- subset(Pfalciparum_upstream_regions_1kb_strand, Pfalciparum_upstream_regions_1kb_strand$GENEID %in% background_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    
    ##background nodownreg
    background_1kb_nodownreg_strand_Niggi_noDups_upstream_regions <- subset(Pfalciparum_upstream_regions_1kb_strand, Pfalciparum_upstream_regions_1kb_strand$GENEID %in% background_nodownreg_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    
    ##upreg and downreg
    upreg_1kb_strand_Niggi_noDups_upstream_regions <- subset(Pfalciparum_upstream_regions_1kb_strand, Pfalciparum_upstream_regions_1kb_strand$GENEID %in% upreg_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    downreg_1kb_strand_Niggi_noDups_upstream_regions  <- subset(Pfalciparum_upstream_regions_1kb_strand, Pfalciparum_upstream_regions_1kb_strand$GENEID %in% downreg_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    
#####1kb save lists of GENEids for each category#####
    write.csv(upreg_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "upreg_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(downreg_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "downreg_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(background_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "background_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(background_1kb_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID, "background_1kb_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(allgenes_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "allgenes_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    
#####1kb create and download text for bedtools BED files#####
    
    #same concept as previous - use the coordinates in the upstream_regions 
    #dataframes to create BED files for each subset of genes
    # BED columns: CHROM START END NAME
    
    # ##lists to make
    # allgenes_1kb_strand_Niggi_noDups_upstream_regions
    # upreg_1kb_strand_Niggi_noDups_upstream_regions
    # downreg_1kb_strand_Niggi_noDups_upstream_regions
    # background_1kb_strand_Niggi_noDups_upstream_regions
    # background_1kb_nodownreg_strand_Niggi_noDups_upstream_regions
    
    # BED columns: CHROM START END NAME
    
    upreg_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 339))
    upreg_1kb_strand_Niggi_noDups_download_BED$CHROM <- upreg_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    #ISSUE HERE - missing CHROM in upreg_1kb_strand_Niggi_noDups_upstream_regions
    upreg_1kb_strand_Niggi_noDups_download_BED$START <- upreg_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    upreg_1kb_strand_Niggi_noDups_download_BED$END <- upreg_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    upreg_1kb_strand_Niggi_noDups_download_BED$NAME <- upreg_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    upreg_1kb_strand_Niggi_noDups_download_BED <- upreg_1kb_strand_Niggi_noDups_download_BED[,2:5]
    write.csv(upreg_1kb_strand_Niggi_noDups_download_BED, "upreg_1kb_strand_Niggi_noDups_download_BED.csv")
    
    downreg_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 42))
    downreg_1kb_strand_Niggi_noDups_download_BED$CHROM <- downreg_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    downreg_1kb_strand_Niggi_noDups_download_BED$START <- downreg_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    downreg_1kb_strand_Niggi_noDups_download_BED$END <- downreg_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    downreg_1kb_strand_Niggi_noDups_download_BED$NAME <- downreg_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    downreg_1kb_strand_Niggi_noDups_download_BED <- downreg_1kb_strand_Niggi_noDups_download_BED[,2:5]
    write.csv(downreg_1kb_strand_Niggi_noDups_download_BED, "downreg_1kb_strand_Niggi_noDups_download_BED.csv")
    
    background_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5334))
    background_1kb_strand_Niggi_noDups_download_BED$CHROM <- background_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    background_1kb_strand_Niggi_noDups_download_BED$START <- background_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    background_1kb_strand_Niggi_noDups_download_BED$END <- background_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    background_1kb_strand_Niggi_noDups_download_BED$NAME <- background_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    background_1kb_strand_Niggi_noDups_download_BED <- background_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(background_1kb_strand_Niggi_noDups_download_BED, "background_1kb_strand_Niggi_noDups_download_BED.CSV")
    
    background_1kb_nodownreg_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5292))  
    background_1kb_nodownreg_strand_Niggi_noDups_download_BED$CHROM <- background_1kb_nodownreg_strand_Niggi_noDups_upstream_regions$CHROM
    background_1kb_nodownreg_strand_Niggi_noDups_download_BED$START <- background_1kb_nodownreg_strand_Niggi_noDups_upstream_regions$UPSTART
    background_1kb_nodownreg_strand_Niggi_noDups_download_BED$END <- background_1kb_nodownreg_strand_Niggi_noDups_upstream_regions$UPEND
    background_1kb_nodownreg_strand_Niggi_noDups_download_BED$NAME <- background_1kb_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID
    background_1kb_nodownreg_strand_Niggi_noDups_download_BED <- background_1kb_nodownreg_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(background_1kb_nodownreg_strand_Niggi_noDups_download_BED, "background_1kb_nodownreg_strand_Niggi_noDups_download_BED.csv")
    
    allgenes_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5673))
    allgenes_1kb_strand_Niggi_noDups_download_BED$CHROM <- allgenes_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    allgenes_1kb_strand_Niggi_noDups_download_BED$START <- allgenes_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    allgenes_1kb_strand_Niggi_noDups_download_BED$END <- allgenes_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    allgenes_1kb_strand_Niggi_noDups_download_BED$NAME <- allgenes_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    allgenes_1kb_strand_Niggi_noDups_download_BED <- allgenes_1kb_strand_Niggi_noDups_download_BED[,2:5]
    write.csv(allgenes_1kb_strand_Niggi_noDups_download_BED, "allgenes_1kb_strand_Niggi_noDups_download_BED.csv")
    
    
    
    

    
    
    

    
#####2kb and 1kb Pfalc upreg and downreg#####
#load pfalc genelists
upreg_Pfalc_notrodent <- read.csv("./GeneLists/upreg_Pfalc_NOTrodent_update.csv")
    upreg_Pfalc_notrodent <- as.factor(upreg_Pfalc_notrodent[,1])

downreg_Pfalc_notrodent <- read.csv("./GeneLists/downreg_Pfalc_NOTrodent_update.csv")
    downreg_Pfalc_notrodent <- as.factor(downreg_Pfalc_notrodent[,1])    

##upreg and downreg pfalc 2kb
upreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions <- subset(Pfalciparum_upstream_regions_2kb_strand, 
                                                               Pfalciparum_upstream_regions_2kb_strand$GENEID %in% upreg_Pfalc_notrodent)
downreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions  <- subset(Pfalciparum_upstream_regions_2kb_strand, 
                                                                  Pfalciparum_upstream_regions_2kb_strand$GENEID %in% downreg_Pfalc_notrodent)
#save geneID list
write.csv(upreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions$GENEID, "upreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
write.csv(downreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions$GENEID, "downreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions_genelist.csv")

# BED columns: CHROM START END NAME

upreg_Pfalc_2kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 24))
upreg_Pfalc_2kb_strand_Niggi_noDups_download_BED$CHROM <- upreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions$CHROM
#ISSUE HERE - missing CHROM in upreg_2kb_strand_Niggi_noDups_upstream_regions
upreg_Pfalc_2kb_strand_Niggi_noDups_download_BED$START <- upreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions$UPSTART
upreg_Pfalc_2kb_strand_Niggi_noDups_download_BED$END <- upreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions$UPEND
upreg_Pfalc_2kb_strand_Niggi_noDups_download_BED$NAME <- upreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions$GENEID
upreg_Pfalc_2kb_strand_Niggi_noDups_download_BED <- upreg_Pfalc_2kb_strand_Niggi_noDups_download_BED[,2:5]
write.csv(upreg_Pfalc_2kb_strand_Niggi_noDups_download_BED, "upreg_Pfalc_2kb_strand_Niggi_noDups_download_BED.csv")


downreg_Pfalc_2kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 17))
downreg_Pfalc_2kb_strand_Niggi_noDups_download_BED$CHROM <- downreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions$CHROM
#ISSUE HERE - missing CHROM in downreg_2kb_strand_Niggi_noDups_upstream_regions
downreg_Pfalc_2kb_strand_Niggi_noDups_download_BED$START <- downreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions$UPSTART
downreg_Pfalc_2kb_strand_Niggi_noDups_download_BED$END <- downreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions$UPEND
downreg_Pfalc_2kb_strand_Niggi_noDups_download_BED$NAME <- downreg_Pfalc_2kb_strand_Niggi_noDups_upstream_regions$GENEID
downreg_Pfalc_2kb_strand_Niggi_noDups_download_BED <- downreg_Pfalc_2kb_strand_Niggi_noDups_download_BED[,2:5]
write.csv(downreg_Pfalc_2kb_strand_Niggi_noDups_download_BED, "downreg_Pfalc_2kb_strand_Niggi_noDups_download_BED.csv")


##upreg and downreg pfalc 1kb
upreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions <- subset(Pfalciparum_upstream_regions_1kb_strand, 
                                                               Pfalciparum_upstream_regions_1kb_strand$GENEID %in% upreg_Pfalc_notrodent)
downreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions  <- subset(Pfalciparum_upstream_regions_1kb_strand, 
                                                                  Pfalciparum_upstream_regions_1kb_strand$GENEID %in% downreg_Pfalc_notrodent)
#save geneID list
write.csv(upreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "upreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
write.csv(downreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "downreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")

# BED columns: CHROM START END NAME

upreg_Pfalc_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 24))
upreg_Pfalc_1kb_strand_Niggi_noDups_download_BED$CHROM <- upreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions$CHROM
#ISSUE HERE - missing CHROM in upreg_1kb_strand_Niggi_noDups_upstream_regions
upreg_Pfalc_1kb_strand_Niggi_noDups_download_BED$START <- upreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
upreg_Pfalc_1kb_strand_Niggi_noDups_download_BED$END <- upreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions$UPEND
upreg_Pfalc_1kb_strand_Niggi_noDups_download_BED$NAME <- upreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions$GENEID
upreg_Pfalc_1kb_strand_Niggi_noDups_download_BED <- upreg_Pfalc_1kb_strand_Niggi_noDups_download_BED[,2:5]
write.csv(upreg_Pfalc_1kb_strand_Niggi_noDups_download_BED, "upreg_Pfalc_1kb_strand_Niggi_noDups_download_BED.csv")


downreg_Pfalc_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 17))
downreg_Pfalc_1kb_strand_Niggi_noDups_download_BED$CHROM <- downreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions$CHROM
#ISSUE HERE - missing CHROM in downreg_1kb_strand_Niggi_noDups_upstream_regions
downreg_Pfalc_1kb_strand_Niggi_noDups_download_BED$START <- downreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
downreg_Pfalc_1kb_strand_Niggi_noDups_download_BED$END <- downreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions$UPEND
downreg_Pfalc_1kb_strand_Niggi_noDups_download_BED$NAME <- downreg_Pfalc_1kb_strand_Niggi_noDups_upstream_regions$GENEID
downreg_Pfalc_1kb_strand_Niggi_noDups_download_BED <- downreg_Pfalc_1kb_strand_Niggi_noDups_download_BED[,2:5]
write.csv(downreg_Pfalc_1kb_strand_Niggi_noDups_download_BED, "downreg_Pfalc_1kb_strand_Niggi_noDups_download_BED.csv")



    
#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/LysoPC_upstreamReg_org_overlap_with_others.RData")
  