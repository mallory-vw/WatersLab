#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/LysoPC_upstreamReg_org.RData")
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
Pfalciparum_GFF3_allgene_strand <- read.csv("PlasmoDB_28_GeneStart_GeneStop_strand_calc.csv", header = T)
#chrom sizes
Pfalciparum_chrom_coord <- read.csv("Pfalciparum_chromosome_coord.csv")


###start with the beginning and ending of CHROM
#beginning if the first gene is on the + strand, end if the last gene is on the - strand

###label the start of each chromosome
for (i in seq_along(Pfalciparum_GFF3_allgene_strand$UP_DIST[2:5688])) {
  if (identical(Pfalciparum_GFF3_allgene_strand$SEQID[i], Pfalciparum_GFF3_allgene_strand$SEQID[i-1])) {
    Pfalciparum_GFF3_allgene_strand$UP_DIST[i] <- "calculate"
  } else {
    Pfalciparum_GFF3_allgene_strand$UP_DIST[i] <- "START"
  }
}
###label the end of each chromosome
for (i in seq_along(Pfalciparum_GFF3_allgene_strand$UP_DIST[2:5688])) {
  if (identical(Pfalciparum_GFF3_allgene_strand$SEQID[i], Pfalciparum_GFF3_allgene_strand$SEQID[i+1]) && (Pfalciparum_GFF3_allgene_strand$UP_DIST[i] == "calculate")) {
    Pfalciparum_GFF3_allgene_strand$UP_DIST[i] <- "calculate"
  } else if (Pfalciparum_GFF3_allgene_strand$UP_DIST[i] != "START") {
        Pfalciparum_GFF3_allgene_strand$UP_DIST[i] <- "END"
  }
}

###at start of CHROM is gene is on + strand, calculate upstream dist
#if gene is on - strand will go through normal calculation
for (i in seq_along(Pfalciparum_GFF3_allgene_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_strand$UP_DIST[i] == "START") && (Pfalciparum_GFF3_allgene_strand$STRAND[i] == "+")) {
    Pfalciparum_GFF3_allgene_strand$UP_DIST[i] <- Pfalciparum_GFF3_allgene_strand$START[i] - 1
  } 
}
###at end of CHROM if gene is on - strand, need to calculate upstream dist based on size of CHROM
#if gene is on + strand will go through normal calc
for (i in seq_along(Pfalciparum_GFF3_allgene_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_strand$UP_DIST[i] == "END") && (Pfalciparum_GFF3_allgene_strand$STRAND[i] == "-")) {
    Pfalciparum_GFF3_allgene_strand$UP_DIST[i] <- (Pfalciparum_chrom_coord$END[Pfalciparum_chrom_coord$CHROM == Pfalciparum_GFF3_allgene_strand$SEQID[i]]) - Pfalciparum_GFF3_allgene_strand$END[i]
  } 
}
#then move on to populating the rest of the values, ignoring the rows that have already been filled in
#change all UP_DIST values that aren't numbers to "calculate"
for (i in seq_along(Pfalciparum_GFF3_allgene_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_strand$UP_DIST[i] == "calculate") ||
      (Pfalciparum_GFF3_allgene_strand$UP_DIST[i] == "START") ||
      (Pfalciparum_GFF3_allgene_strand$UP_DIST[i] == "END") ) {
    Pfalciparum_GFF3_allgene_strand$UP_DIST[i] <- "calculate"
  }
}
#change all UP_DIST values for + strand genes to "plus"
for (i in seq_along(Pfalciparum_GFF3_allgene_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_strand$UP_DIST[i] == "calculate") &&
      (Pfalciparum_GFF3_allgene_strand$STRAND[i] == "+")) {
    Pfalciparum_GFF3_allgene_strand$UP_DIST[i] <- "plus"
  }
}
#change all UP_DIST values for - strand genes to "minus"
for (i in seq_along(Pfalciparum_GFF3_allgene_strand$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_strand$UP_DIST[i] == "calculate") &&
      (Pfalciparum_GFF3_allgene_strand$STRAND[i] == "-")) {
    Pfalciparum_GFF3_allgene_strand$UP_DIST[i] <- "minus"
  }
}
#calculate upstream dist for + strand genes
#if STRAND = +, UP_DIST = START - END of line previous 
for (i in seq_along(Pfalciparum_GFF3_allgene_strand$UP_DIST)) {
  if (Pfalciparum_GFF3_allgene_strand$UP_DIST[i] == "plus") {
    Pfalciparum_GFF3_allgene_strand$UP_DIST[i] <- Pfalciparum_GFF3_allgene_strand$START[i] - 
      Pfalciparum_GFF3_allgene_strand$END[i-1]
  }
}
#calculate upstream dist for - strand genes
#if STRAND = -, UP_DIST = START of line ahead - END
for (i in seq_along(Pfalciparum_GFF3_allgene_strand$UP_DIST)) {
  if (Pfalciparum_GFF3_allgene_strand$UP_DIST[i] == "minus") {
    Pfalciparum_GFF3_allgene_strand$UP_DIST[i] <- Pfalciparum_GFF3_allgene_strand$START[i+1] - 
      Pfalciparum_GFF3_allgene_strand$END[i]
  }
}
##convert values to integers
Pfalciparum_GFF3_allgene_strand$UP_DIST <- as.integer(Pfalciparum_GFF3_allgene_strand$UP_DIST)



#####2kb organize upstream cutoff values#####

#generate edited Up_DIST to cut off at 2KB
Pfalciparum_GFF3_allgene_strand$UP_DIST_EDIT <- ifelse(Pfalciparum_GFF3_allgene_strand$UP_DIST > 2001, 2001, Pfalciparum_GFF3_allgene_strand$UP_DIST)
#remove genes with UP_DIST <= 2 (no distance between genes)
  pos_Pfalciparum_GFF3_allgene_strand <- subset(Pfalciparum_GFF3_allgene_strand, Pfalciparum_GFF3_allgene_strand$UP_DIST_EDIT >= 3)
    write.csv(pos_Pfalciparum_GFF3_allgene_strand, "pos_Pfalciparum_GFF3_allgene_strand.csv")
  neg_Pfalciparum_GFF3_allgene_strand <- subset(Pfalciparum_GFF3_allgene_strand, Pfalciparum_GFF3_allgene_strand$UP_DIST_EDIT <= 2)
    write.csv(neg_Pfalciparum_GFF3_allgene_strand, "neg_Pfalciparum_GFF3_allgene_strand.csv")

##all genes in Niggi's list
allgenes_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_strand, Pfalciparum_GFF3_allgene_strand$GENEID %in% allgenes_Niggi_noDups_clean)
  #we lose 7 genes that are in allgenes_Niggi_noDups_clean (5680) but aren't in the PlasmoDB list (since allgenes_strand_Niggi_noDups_GFF3_info has 5673)
  #we lose 14 genes that are in the PlasmoDB list (5687) but not in allgenes_Niggi_noDups_clean, end up with 5673 in allgenes_strand_Niggi_noDups_GFF3_info

  #want to remove the 53 genes with UP_DIST_EDIT <= 2
  pos_allgenes_strand_Niggi_noDups_GFF3_info <- subset(allgenes_strand_Niggi_noDups_GFF3_info, allgenes_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 3)
    pos_allgenes_strand_Niggi_noDups <- pos_allgenes_strand_Niggi_noDups_GFF3_info$GENEID  
    write.csv(pos_allgenes_strand_Niggi_noDups_GFF3_info, "pos_allgenes_strand_Niggi_noDups_GFF3_info.csv")
  neg_allgenes_strand_Niggi_noDups_GFF3_info <- subset(allgenes_strand_Niggi_noDups_GFF3_info, allgenes_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT <= 2)
    neg_allgenes_strand_Niggi_noDups <- neg_allgenes_strand_Niggi_noDups_GFF3_info$GENEID  
    write.csv(neg_allgenes_strand_Niggi_noDups_GFF3_info, "neg_allgenes_strand_Niggi_noDups_GFF3_info.csv")
    #we lose 53 genes with UP_DIST_EDIT <= 2, mostly in the mitos
  #use pos_allgenes_strand_Niggi_noDups_GFF3_info GENEID list moving forwards
  pos4_allgenes_strand_Niggi_noDups_GFF3_info <- subset(pos_allgenes_strand_Niggi_noDups_GFF3_info, pos_allgenes_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 5)
    pos4_allgenes_strand_Niggi_noDups <- pos4_allgenes_strand_Niggi_noDups_GFF3_info$GENEID
  pos5_allgenes_strand_Niggi_noDups_GFF3_info <- subset(pos_allgenes_strand_Niggi_noDups_GFF3_info, pos_allgenes_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 6)
    pos5_allgenes_strand_Niggi_noDups <- pos5_allgenes_strand_Niggi_noDups_GFF3_info$GENEID
  pos6_allgenes_strand_Niggi_noDups_GFF3_info <- subset(pos_allgenes_strand_Niggi_noDups_GFF3_info, pos_allgenes_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 7)
    pos6_allgenes_strand_Niggi_noDups <- pos6_allgenes_strand_Niggi_noDups_GFF3_info$GENEID
  pos7_allgenes_strand_Niggi_noDups_GFF3_info <- subset(pos_allgenes_strand_Niggi_noDups_GFF3_info, pos_allgenes_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 8)
    pos7_allgenes_strand_Niggi_noDups <- pos7_allgenes_strand_Niggi_noDups_GFF3_info$GENEID
  pos8_allgenes_strand_Niggi_noDups_GFF3_info <- subset(pos_allgenes_strand_Niggi_noDups_GFF3_info, pos_allgenes_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 9)
    pos8_allgenes_strand_Niggi_noDups <- pos8_allgenes_strand_Niggi_noDups_GFF3_info$GENEID
    
##background
background_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_strand, Pfalciparum_GFF3_allgene_strand$GENEID %in% background_Niggi_noDups_clean)
  #we lose 7 genes that are in background_Niggi_noDups_clean (5341) but aren't in the PlasmoDB list (background_strand_Niggi_noDups_GFF3_info is 5334)

  pos_background_strand_Niggi_noDups_GFF3_info <- subset(background_strand_Niggi_noDups_GFF3_info, background_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 3)
    pos_background_strand_Niggi_noDups <- pos_background_strand_Niggi_noDups_GFF3_info$GENEID
  neg_background_strand_Niggi_noDups_GFF3_info <- subset(background_strand_Niggi_noDups_GFF3_info, background_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT <= 2)
  #we lose 53 genes with UP_DIST_EDIT <= 2, mostly in the mitos

  #use pos_background_strand_Niggi_noDups_GFF3_info GENEID list moving forwards
  pos4_background_strand_Niggi_noDups_GFF3_info <- subset(pos_background_strand_Niggi_noDups_GFF3_info, pos_background_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 5)
    pos4_background_strand_Niggi_noDups <- pos4_background_strand_Niggi_noDups_GFF3_info$GENEID
  pos5_background_strand_Niggi_noDups_GFF3_info <- subset(pos_background_strand_Niggi_noDups_GFF3_info, pos_background_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 6)
    pos5_background_strand_Niggi_noDups <- pos5_background_strand_Niggi_noDups_GFF3_info$GENEID
  pos6_background_strand_Niggi_noDups_GFF3_info <- subset(pos_background_strand_Niggi_noDups_GFF3_info, pos_background_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 7)
    pos6_background_strand_Niggi_noDups <- pos6_background_strand_Niggi_noDups_GFF3_info$GENEID
  pos7_background_strand_Niggi_noDups_GFF3_info <- subset(pos_background_strand_Niggi_noDups_GFF3_info, pos_background_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 8)
    pos7_background_strand_Niggi_noDups <- pos7_background_strand_Niggi_noDups_GFF3_info$GENEID
  pos8_background_strand_Niggi_noDups_GFF3_info <- subset(pos_background_strand_Niggi_noDups_GFF3_info, pos_background_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 9)
    pos8_background_strand_Niggi_noDups <- pos8_background_strand_Niggi_noDups_GFF3_info$GENEID

##background no downreg  
background_nodownreg_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_strand, Pfalciparum_GFF3_allgene_strand$GENEID %in% background_nodownreg_Niggi_noDups_clean)
  #we lose 7 genes that are in background_nodownreg_Niggi_noDups_clean (5299) but aren't in the PlasmoDB list (background_nodownreg_strand_Niggi_noDups_GFF3_info is 5292)

  pos_background_nodownreg_strand_Niggi_noDups_GFF3_info <- subset(background_nodownreg_strand_Niggi_noDups_GFF3_info, background_nodownreg_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 3)
    pos_background_nodownreg_strand_Niggi_noDups <- pos_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID  
  neg_background_nodownreg_strand_Niggi_noDups_GFF3_info <- subset(background_nodownreg_strand_Niggi_noDups_GFF3_info, background_nodownreg_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT <= 2)
    #we lose 53 genes with UP_DIST_EDIT <= 2, mostly in the mitos
  
  #use pos_background_nodownreg_strand_Niggi_noDups_GFF3_info GENEID list moving forwards
  pos4_background_nodownreg_strand_Niggi_noDups_GFF3_info <- subset(pos_background_nodownreg_strand_Niggi_noDups_GFF3_info, pos_background_nodownreg_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 5)
    pos4_background_nodownreg_strand_Niggi_noDups <- pos4_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID
  pos5_background_nodownreg_strand_Niggi_noDups_GFF3_info <- subset(pos_background_nodownreg_strand_Niggi_noDups_GFF3_info, pos_background_nodownreg_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 6)
    pos5_background_nodownreg_strand_Niggi_noDups <- pos5_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID
  pos6_background_nodownreg_strand_Niggi_noDups_GFF3_info <- subset(pos_background_nodownreg_strand_Niggi_noDups_GFF3_info, pos_background_nodownreg_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 7)
    pos6_background_nodownreg_strand_Niggi_noDups <- pos6_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID
  pos7_background_nodownreg_strand_Niggi_noDups_GFF3_info <- subset(pos_background_nodownreg_strand_Niggi_noDups_GFF3_info, pos_background_nodownreg_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 8)
    pos7_background_nodownreg_strand_Niggi_noDups <- pos7_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID
  pos8_background_nodownreg_strand_Niggi_noDups_GFF3_info <- subset(pos_background_nodownreg_strand_Niggi_noDups_GFF3_info, pos_background_nodownreg_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 9)
    pos8_background_nodownreg_strand_Niggi_noDups <- pos8_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID

##upreg and downreg
  #these don't need to be split into pos and neg, are all pos, also none missing from either Niggi's list or the PlasmoDB list
  upreg_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_strand, Pfalciparum_GFF3_allgene_strand$GENEID %in% upreg_Niggi_noDups_clean)
  downreg_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_strand, Pfalciparum_GFF3_allgene_strand$GENEID %in% downreg_Niggi_noDups_clean)    
    
    
#the next step is to figure out the beginning and end of the upstream region for both the + strand genes (5' of the gene start) and
#the -ve strand genes (based on coordinate system, 3' of the gene end)
  
# +ve gene --> UPSTART = GENE_START - UP_DIST_EDIT 
              #UPEND = GENE_START - 1
# -ve gene --> UPSTART = GENE_END 
              #UPEND = GENE_END + (UP_DIST_EDIT - 1)
pos_Pfalciparum_upstream_regions_strand <- data.frame(matrix(NA, nrow = 5632))
pos_Pfalciparum_upstream_regions_strand$CHROM <-  pos_Pfalciparum_GFF3_allgene_strand$SEQID
pos_Pfalciparum_upstream_regions_strand$GENEID <- pos_Pfalciparum_GFF3_allgene_strand$GENEID
pos_Pfalciparum_upstream_regions_strand <- pos_Pfalciparum_upstream_regions_strand[,-1]
pos_Pfalciparum_upstream_regions_strand$STRAND <- pos_Pfalciparum_GFF3_allgene_strand$STRAND
pos_Pfalciparum_upstream_regions_strand$GENE_START <- pos_Pfalciparum_GFF3_allgene_strand$START
pos_Pfalciparum_upstream_regions_strand$GENE_END <- pos_Pfalciparum_GFF3_allgene_strand$END
pos_Pfalciparum_upstream_regions_strand$UP_DIST_EDIT <- pos_Pfalciparum_GFF3_allgene_strand$UP_DIST_EDIT
pos_Pfalciparum_upstream_regions_strand$UPSTART <- NA
pos_Pfalciparum_upstream_regions_strand$UPEND <- NA

for (i in seq_along(pos_Pfalciparum_upstream_regions_strand$GENEID)) {
  if (pos_Pfalciparum_upstream_regions_strand$STRAND[i] == "+") {
    pos_Pfalciparum_upstream_regions_strand$UPSTART[i] <- pos_Pfalciparum_upstream_regions_strand$GENE_START[i] - pos_Pfalciparum_upstream_regions_strand$UP_DIST_EDIT[i]
    pos_Pfalciparum_upstream_regions_strand$UPEND[i] <- pos_Pfalciparum_upstream_regions_strand$GENE_START[i] - 1
  }
}

for (i in seq_along(pos_Pfalciparum_upstream_regions_strand$GENEID)) {
  if (pos_Pfalciparum_upstream_regions_strand$STRAND[i] == "-") {
    pos_Pfalciparum_upstream_regions_strand$UPSTART[i] <- pos_Pfalciparum_upstream_regions_strand$GENE_END[i]
    pos_Pfalciparum_upstream_regions_strand$UPEND[i] <- pos_Pfalciparum_upstream_regions_strand$GENE_END[i] + (pos_Pfalciparum_upstream_regions_strand$UP_DIST_EDIT[i] - 1)
  }
}

write.csv(pos_Pfalciparum_upstream_regions_strand, "pos_Pfalciparum_upstream_regions_strand.csv")

#####2kb split upstream regions into the specific gene lists#####
##allgenes
pos_allgenes_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos_allgenes_strand_Niggi_noDups_GFF3_info$GENEID)
  pos4_allgenes_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos4_allgenes_strand_Niggi_noDups_GFF3_info$GENEID)
  pos5_allgenes_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos5_allgenes_strand_Niggi_noDups_GFF3_info$GENEID)
  pos6_allgenes_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos6_allgenes_strand_Niggi_noDups_GFF3_info$GENEID)
  pos7_allgenes_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos7_allgenes_strand_Niggi_noDups_GFF3_info$GENEID)
  pos8_allgenes_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos8_allgenes_strand_Niggi_noDups_GFF3_info$GENEID)

##background
pos_background_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos_background_strand_Niggi_noDups_GFF3_info$GENEID)
  pos4_background_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos4_background_strand_Niggi_noDups_GFF3_info$GENEID)
  pos5_background_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos5_background_strand_Niggi_noDups_GFF3_info$GENEID)
  pos6_background_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos6_background_strand_Niggi_noDups_GFF3_info$GENEID)
  pos7_background_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos7_background_strand_Niggi_noDups_GFF3_info$GENEID)
  pos8_background_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos8_background_strand_Niggi_noDups_GFF3_info$GENEID)
  
##background nodownreg
pos_background_nodownreg_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID)
  pos4_background_nodownreg_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos4_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID)
  pos5_background_nodownreg_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos5_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID)
  pos6_background_nodownreg_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos6_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID)
  pos7_background_nodownreg_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos7_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID)
  pos8_background_nodownreg_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% pos8_background_nodownreg_strand_Niggi_noDups_GFF3_info$GENEID)


##upreg and downreg
upreg_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% upreg_strand_Niggi_noDups_GFF3_info$GENEID)
downreg_strand_Niggi_noDups_upstream_regions  <- subset(pos_Pfalciparum_upstream_regions_strand, pos_Pfalciparum_upstream_regions_strand$GENEID %in% downreg_strand_Niggi_noDups_GFF3_info$GENEID)


#####2kb save lists of GENEids for each category#####
write.csv(upreg_strand_Niggi_noDups_upstream_regions$GENEID, "upreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
write.csv(downreg_strand_Niggi_noDups_upstream_regions$GENEID, "downreg_strand_Niggi_noDups_upstream_regions_genelist.csv")

write.csv(pos_background_strand_Niggi_noDups_upstream_regions$GENEID, "pos_background_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos4_background_strand_Niggi_noDups_upstream_regions$GENEID, "pos4_background_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos5_background_strand_Niggi_noDups_upstream_regions$GENEID, "pos5_background_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos6_background_strand_Niggi_noDups_upstream_regions$GENEID, "pos6_background_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos7_background_strand_Niggi_noDups_upstream_regions$GENEID, "pos7_background_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos8_background_strand_Niggi_noDups_upstream_regions$GENEID, "pos8_background_strand_Niggi_noDups_upstream_regions_genelist.csv")

write.csv(pos_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID, "pos_background_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos4_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID, "pos4_background_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos5_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID, "pos5_background_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos6_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID, "pos6_background_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos7_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID, "pos7_background_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos8_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID, "pos8_background_nodownreg_strand_Niggi_noDups_upstream_regions_genelist.csv")

write.csv(pos_allgenes_strand_Niggi_noDups_upstream_regions$GENEID, "pos_allgenes_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos4_allgenes_strand_Niggi_noDups_upstream_regions$GENEID, "pos4_allgenes_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos5_allgenes_strand_Niggi_noDups_upstream_regions$GENEID, "pos5_allgenes_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos6_allgenes_strand_Niggi_noDups_upstream_regions$GENEID, "pos6_allgenes_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos7_allgenes_strand_Niggi_noDups_upstream_regions$GENEID, "pos7_allgenes_strand_Niggi_noDups_upstream_regions_genelist.csv")
  write.csv(pos8_allgenes_strand_Niggi_noDups_upstream_regions$GENEID, "pos8_allgenes_strand_Niggi_noDups_upstream_regions_genelist.csv")
  

  

#####2kb create and download text for bedtools BED files#####
  
#same concept as previous - use the coordinates in the upstream_regions 
#dataframes to create BED files for each subset of genes
# BED columns: CHROM START END NAME

##lists to make
pos_allgenes_strand_Niggi_noDups_upstream_regions
  pos4_allgenes_strand_Niggi_noDups_upstream_regions
  pos5_allgenes_strand_Niggi_noDups_upstream_regions
  pos6_allgenes_strand_Niggi_noDups_upstream_regions
  pos7_allgenes_strand_Niggi_noDups_upstream_regions
  pos8_allgenes_strand_Niggi_noDups_upstream_regions
pos_background_strand_Niggi_noDups_upstream_regions
  pos4_background_strand_Niggi_noDups_upstream_regions
  pos5_background_strand_Niggi_noDups_upstream_regions
  pos6_background_strand_Niggi_noDups_upstream_regions
  pos7_background_strand_Niggi_noDups_upstream_regions
  pos8_background_strand_Niggi_noDups_upstream_regions
pos_background_nodownreg_strand_Niggi_noDups_upstream_regions
  pos4_background_nodownreg_strand_Niggi_noDups_upstream_regions
  pos5_background_nodownreg_strand_Niggi_noDups_upstream_regions
  pos6_background_nodownreg_strand_Niggi_noDups_upstream_regions
  pos7_background_nodownreg_strand_Niggi_noDups_upstream_regions
  pos8_background_nodownreg_strand_Niggi_noDups_upstream_regions
upreg_strand_Niggi_noDups_upstream_regions
downreg_strand_Niggi_noDups_upstream_regions



# BED columns: CHROM START END NAME

upreg_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 339))
  upreg_strand_Niggi_noDups_download_BED$CHROM <- upreg_strand_Niggi_noDups_upstream_regions$CHROM
  upreg_strand_Niggi_noDups_download_BED$START <- upreg_strand_Niggi_noDups_upstream_regions$UPSTART
  upreg_strand_Niggi_noDups_download_BED$END <- upreg_strand_Niggi_noDups_upstream_regions$UPEND
  upreg_strand_Niggi_noDups_download_BED$NAME <- upreg_strand_Niggi_noDups_upstream_regions$GENEID
  upreg_strand_Niggi_noDups_download_BED <- upreg_strand_Niggi_noDups_download_BED[,2:5]
  write.csv(upreg_strand_Niggi_noDups_download_BED, "upreg_strand_Niggi_noDups_download_BED.csv")
  
downreg_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 42))
  downreg_strand_Niggi_noDups_download_BED$CHROM <- downreg_strand_Niggi_noDups_upstream_regions$CHROM
  downreg_strand_Niggi_noDups_download_BED$START <- downreg_strand_Niggi_noDups_upstream_regions$UPSTART
  downreg_strand_Niggi_noDups_download_BED$END <- downreg_strand_Niggi_noDups_upstream_regions$UPEND
  downreg_strand_Niggi_noDups_download_BED$NAME <- downreg_strand_Niggi_noDups_upstream_regions$GENEID
  downreg_strand_Niggi_noDups_download_BED <- downreg_strand_Niggi_noDups_download_BED[,2:5]
  write.csv(downreg_strand_Niggi_noDups_download_BED, "downreg_strand_Niggi_noDups_download_BED.csv")
  
pos_background_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5281))
  pos_background_strand_Niggi_noDups_download_BED$CHROM <- pos_background_strand_Niggi_noDups_upstream_regions$CHROM
  pos_background_strand_Niggi_noDups_download_BED$START <- pos_background_strand_Niggi_noDups_upstream_regions$UPSTART
  pos_background_strand_Niggi_noDups_download_BED$END <- pos_background_strand_Niggi_noDups_upstream_regions$UPEND
  pos_background_strand_Niggi_noDups_download_BED$NAME <- pos_background_strand_Niggi_noDups_upstream_regions$GENEID
  pos_background_strand_Niggi_noDups_download_BED <- pos_background_strand_Niggi_noDups_download_BED[,2:5]  
  write.csv(pos_background_strand_Niggi_noDups_download_BED, "pos_background_strand_Niggi_noDups_download_BED.csv")
    pos4_background_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5279))
      pos4_background_strand_Niggi_noDups_download_BED$CHROM <- pos4_background_strand_Niggi_noDups_upstream_regions$CHROM
      pos4_background_strand_Niggi_noDups_download_BED$START <- pos4_background_strand_Niggi_noDups_upstream_regions$UPSTART
      pos4_background_strand_Niggi_noDups_download_BED$END <- pos4_background_strand_Niggi_noDups_upstream_regions$UPEND
      pos4_background_strand_Niggi_noDups_download_BED$NAME <- pos4_background_strand_Niggi_noDups_upstream_regions$GENEID
      pos4_background_strand_Niggi_noDups_download_BED <- pos4_background_strand_Niggi_noDups_download_BED[,2:5]  
      write.csv(pos4_background_strand_Niggi_noDups_download_BED, "pos4_background_strand_Niggi_noDups_download_BED.csv")
    pos5_background_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5278))
      pos5_background_strand_Niggi_noDups_download_BED$CHROM <- pos5_background_strand_Niggi_noDups_upstream_regions$CHROM
      pos5_background_strand_Niggi_noDups_download_BED$START <- pos5_background_strand_Niggi_noDups_upstream_regions$UPSTART
      pos5_background_strand_Niggi_noDups_download_BED$END <- pos5_background_strand_Niggi_noDups_upstream_regions$UPEND
      pos5_background_strand_Niggi_noDups_download_BED$NAME <- pos5_background_strand_Niggi_noDups_upstream_regions$GENEID
      pos5_background_strand_Niggi_noDups_download_BED <- pos5_background_strand_Niggi_noDups_download_BED[,2:5]  
      write.csv(pos5_background_strand_Niggi_noDups_download_BED, "pos5_background_strand_Niggi_noDups_download_BED.csv")
    pos6_background_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5273))
      pos6_background_strand_Niggi_noDups_download_BED$CHROM <- pos6_background_strand_Niggi_noDups_upstream_regions$CHROM
      pos6_background_strand_Niggi_noDups_download_BED$START <- pos6_background_strand_Niggi_noDups_upstream_regions$UPSTART
      pos6_background_strand_Niggi_noDups_download_BED$END <- pos6_background_strand_Niggi_noDups_upstream_regions$UPEND
      pos6_background_strand_Niggi_noDups_download_BED$NAME <- pos6_background_strand_Niggi_noDups_upstream_regions$GENEID
      pos6_background_strand_Niggi_noDups_download_BED <- pos6_background_strand_Niggi_noDups_download_BED[,2:5]  
      write.csv(pos6_background_strand_Niggi_noDups_download_BED, "pos6_background_strand_Niggi_noDups_download_BED.csv")
    pos7_background_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5270))
      pos7_background_strand_Niggi_noDups_download_BED$CHROM <- pos7_background_strand_Niggi_noDups_upstream_regions$CHROM
      pos7_background_strand_Niggi_noDups_download_BED$START <- pos7_background_strand_Niggi_noDups_upstream_regions$UPSTART
      pos7_background_strand_Niggi_noDups_download_BED$END <- pos7_background_strand_Niggi_noDups_upstream_regions$UPEND
      pos7_background_strand_Niggi_noDups_download_BED$NAME <- pos7_background_strand_Niggi_noDups_upstream_regions$GENEID
      pos7_background_strand_Niggi_noDups_download_BED <- pos7_background_strand_Niggi_noDups_download_BED[,2:5]  
      write.csv(pos7_background_strand_Niggi_noDups_download_BED, "pos7_background_strand_Niggi_noDups_download_BED.csv")
    pos8_background_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5267))
      pos8_background_strand_Niggi_noDups_download_BED$CHROM <- pos8_background_strand_Niggi_noDups_upstream_regions$CHROM
      pos8_background_strand_Niggi_noDups_download_BED$START <- pos8_background_strand_Niggi_noDups_upstream_regions$UPSTART
      pos8_background_strand_Niggi_noDups_download_BED$END <- pos8_background_strand_Niggi_noDups_upstream_regions$UPEND
      pos8_background_strand_Niggi_noDups_download_BED$NAME <- pos8_background_strand_Niggi_noDups_upstream_regions$GENEID
      pos8_background_strand_Niggi_noDups_download_BED <- pos8_background_strand_Niggi_noDups_download_BED[,2:5]  
      write.csv(pos8_background_strand_Niggi_noDups_download_BED, "pos8_background_strand_Niggi_noDups_download_BED.csv")
  
pos_background_nodownreg_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5239))  
  pos_background_nodownreg_strand_Niggi_noDups_download_BED$CHROM <- pos_background_nodownreg_strand_Niggi_noDups_upstream_regions$CHROM
  pos_background_nodownreg_strand_Niggi_noDups_download_BED$START <- pos_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPSTART
  pos_background_nodownreg_strand_Niggi_noDups_download_BED$END <- pos_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPEND
  pos_background_nodownreg_strand_Niggi_noDups_download_BED$NAME <- pos_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID
  pos_background_nodownreg_strand_Niggi_noDups_download_BED <- pos_background_nodownreg_strand_Niggi_noDups_download_BED[,2:5]  
  write.csv(pos_background_nodownreg_strand_Niggi_noDups_download_BED, "pos_background_nodownreg_strand_Niggi_noDups_download_BED.csv")
  pos4_background_nodownreg_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5237))
    pos4_background_nodownreg_strand_Niggi_noDups_download_BED$CHROM <- pos4_background_nodownreg_strand_Niggi_noDups_upstream_regions$CHROM
    pos4_background_nodownreg_strand_Niggi_noDups_download_BED$START <- pos4_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPSTART
    pos4_background_nodownreg_strand_Niggi_noDups_download_BED$END <- pos4_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPEND
    pos4_background_nodownreg_strand_Niggi_noDups_download_BED$NAME <- pos4_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID
    pos4_background_nodownreg_strand_Niggi_noDups_download_BED <- pos4_background_nodownreg_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos4_background_nodownreg_strand_Niggi_noDups_download_BED, "pos4_background_nodownreg_strand_Niggi_noDups_download_BED.csv")
  pos5_background_nodownreg_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5236))
    pos5_background_nodownreg_strand_Niggi_noDups_download_BED$CHROM <- pos5_background_nodownreg_strand_Niggi_noDups_upstream_regions$CHROM
    pos5_background_nodownreg_strand_Niggi_noDups_download_BED$START <- pos5_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPSTART
    pos5_background_nodownreg_strand_Niggi_noDups_download_BED$END <- pos5_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPEND
    pos5_background_nodownreg_strand_Niggi_noDups_download_BED$NAME <- pos5_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID
    pos5_background_nodownreg_strand_Niggi_noDups_download_BED <- pos5_background_nodownreg_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos5_background_nodownreg_strand_Niggi_noDups_download_BED, "pos5_background_nodownreg_strand_Niggi_noDups_download_BED.csv")
  pos6_background_nodownreg_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5231))
    pos6_background_nodownreg_strand_Niggi_noDups_download_BED$CHROM <- pos6_background_nodownreg_strand_Niggi_noDups_upstream_regions$CHROM
    pos6_background_nodownreg_strand_Niggi_noDups_download_BED$START <- pos6_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPSTART
    pos6_background_nodownreg_strand_Niggi_noDups_download_BED$END <- pos6_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPEND
    pos6_background_nodownreg_strand_Niggi_noDups_download_BED$NAME <- pos6_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID
    pos6_background_nodownreg_strand_Niggi_noDups_download_BED <- pos6_background_nodownreg_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos6_background_nodownreg_strand_Niggi_noDups_download_BED, "pos6_background_nodownreg_strand_Niggi_noDups_download_BED.csv")
  pos7_background_nodownreg_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5228))
    pos7_background_nodownreg_strand_Niggi_noDups_download_BED$CHROM <- pos7_background_nodownreg_strand_Niggi_noDups_upstream_regions$CHROM
    pos7_background_nodownreg_strand_Niggi_noDups_download_BED$START <- pos7_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPSTART
    pos7_background_nodownreg_strand_Niggi_noDups_download_BED$END <- pos7_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPEND
    pos7_background_nodownreg_strand_Niggi_noDups_download_BED$NAME <- pos7_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID
    pos7_background_nodownreg_strand_Niggi_noDups_download_BED <- pos7_background_nodownreg_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos7_background_nodownreg_strand_Niggi_noDups_download_BED, "pos7_background_nodownreg_strand_Niggi_noDups_download_BED.csv")
  pos8_background_nodownreg_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5225))
    pos8_background_nodownreg_strand_Niggi_noDups_download_BED$CHROM <- pos8_background_nodownreg_strand_Niggi_noDups_upstream_regions$CHROM
    pos8_background_nodownreg_strand_Niggi_noDups_download_BED$START <- pos8_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPSTART
    pos8_background_nodownreg_strand_Niggi_noDups_download_BED$END <- pos8_background_nodownreg_strand_Niggi_noDups_upstream_regions$UPEND
    pos8_background_nodownreg_strand_Niggi_noDups_download_BED$NAME <- pos8_background_nodownreg_strand_Niggi_noDups_upstream_regions$GENEID
    pos8_background_nodownreg_strand_Niggi_noDups_download_BED <- pos8_background_nodownreg_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos8_background_nodownreg_strand_Niggi_noDups_download_BED, "pos8_background_nodownreg_strand_Niggi_noDups_download_BED.csv")
  
pos_allgenes_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5620))
  pos_allgenes_strand_Niggi_noDups_download_BED$CHROM <- pos_allgenes_strand_Niggi_noDups_upstream_regions$CHROM
  pos_allgenes_strand_Niggi_noDups_download_BED$START <- pos_allgenes_strand_Niggi_noDups_upstream_regions$UPSTART
  pos_allgenes_strand_Niggi_noDups_download_BED$END <- pos_allgenes_strand_Niggi_noDups_upstream_regions$UPEND
  pos_allgenes_strand_Niggi_noDups_download_BED$NAME <- pos_allgenes_strand_Niggi_noDups_upstream_regions$GENEID
  pos_allgenes_strand_Niggi_noDups_download_BED <- pos_allgenes_strand_Niggi_noDups_download_BED[,2:5]
  write.csv(pos_allgenes_strand_Niggi_noDups_download_BED, "pos_allgenes_strand_Niggi_noDups_download_BED.csv")
  pos4_allgenes_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5618))
    pos4_allgenes_strand_Niggi_noDups_download_BED$CHROM <- pos4_allgenes_strand_Niggi_noDups_upstream_regions$CHROM
    pos4_allgenes_strand_Niggi_noDups_download_BED$START <- pos4_allgenes_strand_Niggi_noDups_upstream_regions$UPSTART
    pos4_allgenes_strand_Niggi_noDups_download_BED$END <- pos4_allgenes_strand_Niggi_noDups_upstream_regions$UPEND
    pos4_allgenes_strand_Niggi_noDups_download_BED$NAME <- pos4_allgenes_strand_Niggi_noDups_upstream_regions$GENEID
    pos4_allgenes_strand_Niggi_noDups_download_BED <- pos4_allgenes_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos4_allgenes_strand_Niggi_noDups_download_BED, "pos4_allgenes_strand_Niggi_noDups_download_BED.csv")
  pos5_allgenes_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5617))
    pos5_allgenes_strand_Niggi_noDups_download_BED$CHROM <- pos5_allgenes_strand_Niggi_noDups_upstream_regions$CHROM
    pos5_allgenes_strand_Niggi_noDups_download_BED$START <- pos5_allgenes_strand_Niggi_noDups_upstream_regions$UPSTART
    pos5_allgenes_strand_Niggi_noDups_download_BED$END <- pos5_allgenes_strand_Niggi_noDups_upstream_regions$UPEND
    pos5_allgenes_strand_Niggi_noDups_download_BED$NAME <- pos5_allgenes_strand_Niggi_noDups_upstream_regions$GENEID
    pos5_allgenes_strand_Niggi_noDups_download_BED <- pos5_allgenes_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos5_allgenes_strand_Niggi_noDups_download_BED, "pos5_allgenes_strand_Niggi_noDups_download_BED.csv")
  pos6_allgenes_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5612))
    pos6_allgenes_strand_Niggi_noDups_download_BED$CHROM <- pos6_allgenes_strand_Niggi_noDups_upstream_regions$CHROM
    pos6_allgenes_strand_Niggi_noDups_download_BED$START <- pos6_allgenes_strand_Niggi_noDups_upstream_regions$UPSTART
    pos6_allgenes_strand_Niggi_noDups_download_BED$END <- pos6_allgenes_strand_Niggi_noDups_upstream_regions$UPEND
    pos6_allgenes_strand_Niggi_noDups_download_BED$NAME <- pos6_allgenes_strand_Niggi_noDups_upstream_regions$GENEID
    pos6_allgenes_strand_Niggi_noDups_download_BED <- pos6_allgenes_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos6_allgenes_strand_Niggi_noDups_download_BED, "pos6_allgenes_strand_Niggi_noDups_download_BED.csv")
  pos7_allgenes_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5609))
    pos7_allgenes_strand_Niggi_noDups_download_BED$CHROM <- pos7_allgenes_strand_Niggi_noDups_upstream_regions$CHROM
    pos7_allgenes_strand_Niggi_noDups_download_BED$START <- pos7_allgenes_strand_Niggi_noDups_upstream_regions$UPSTART
    pos7_allgenes_strand_Niggi_noDups_download_BED$END <- pos7_allgenes_strand_Niggi_noDups_upstream_regions$UPEND
    pos7_allgenes_strand_Niggi_noDups_download_BED$NAME <- pos7_allgenes_strand_Niggi_noDups_upstream_regions$GENEID
    pos7_allgenes_strand_Niggi_noDups_download_BED <- pos7_allgenes_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos7_allgenes_strand_Niggi_noDups_download_BED, "pos7_allgenes_strand_Niggi_noDups_download_BED.csv")
  pos8_allgenes_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5606))
    pos8_allgenes_strand_Niggi_noDups_download_BED$CHROM <- pos8_allgenes_strand_Niggi_noDups_upstream_regions$CHROM
    pos8_allgenes_strand_Niggi_noDups_download_BED$START <- pos8_allgenes_strand_Niggi_noDups_upstream_regions$UPSTART
    pos8_allgenes_strand_Niggi_noDups_download_BED$END <- pos8_allgenes_strand_Niggi_noDups_upstream_regions$UPEND
    pos8_allgenes_strand_Niggi_noDups_download_BED$NAME <- pos8_allgenes_strand_Niggi_noDups_upstream_regions$GENEID
    pos8_allgenes_strand_Niggi_noDups_download_BED <- pos8_allgenes_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos8_allgenes_strand_Niggi_noDups_download_BED, "pos8_allgenes_strand_Niggi_noDups_download_BED.csv")
  
  
  


    
#####1kb load data to calculate upstream_size#####
#load in the gene info 
#the duplicates with ".1" and ".2", etc, have been removed in Excel
Pfalciparum_GFF3_allgene_strand_upto_1kb <- read.csv("PlasmoDB_28_GeneStart_GeneStop_strand_calc.csv", header = T)

###start with the beginning and ending of CHROM
    #beginning if the first gene is on the + strand, end if the last gene is on the - strand
    
###label the start of each chromosome
for (i in seq_along(Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[2:5688])) {
  if (identical(Pfalciparum_GFF3_allgene_strand_upto_1kb$SEQID[i], Pfalciparum_GFF3_allgene_strand_upto_1kb$SEQID[i-1])) {
    Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] <- "calculate"
  } else {
    Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] <- "START"
  }
}

###label the end of each chromosome
for (i in seq_along(Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[2:5688])) {
  if (identical(Pfalciparum_GFF3_allgene_strand_upto_1kb$SEQID[i], Pfalciparum_GFF3_allgene_strand_upto_1kb$SEQID[i+1]) && (Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] == "calculate")) {
    Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] <- "calculate"
  } else if (Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] != "START") {
    Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] <- "END"
  }
}
    
###at start of CHROM is gene is on + strand, calculate upstream dist
#if gene is on - strand will go through normal calculation
for (i in seq_along(Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] == "START") && (Pfalciparum_GFF3_allgene_strand_upto_1kb$STRAND[i] == "+")) {
    Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] <- Pfalciparum_GFF3_allgene_strand_upto_1kb$START[i] - 1
  } 
}
###at end of CHROM if gene is on - strand, need to calculate upstream dist based on size of CHROM
#if gene is on + strand will go through normal calc
for (i in seq_along(Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] == "END") && (Pfalciparum_GFF3_allgene_strand_upto_1kb$STRAND[i] == "-")) {
    Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] <- (Pfalciparum_chrom_coord$END[Pfalciparum_chrom_coord$CHROM == Pfalciparum_GFF3_allgene_strand_upto_1kb$SEQID[i]]) - Pfalciparum_GFF3_allgene_strand_upto_1kb$END[i]
  } 
}
#then move on to populating the rest of the values, ignoring the rows that have already been filled in
#change all UP_DIST values that aren't numbers to "calculate"
for (i in seq_along(Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] == "calculate") ||
      (Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] == "START") ||
      (Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] == "END") ) {
    Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] <- "calculate"
  }
}
#change all UP_DIST values for + strand genes to "plus"
for (i in seq_along(Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] == "calculate") &&
      (Pfalciparum_GFF3_allgene_strand_upto_1kb$STRAND[i] == "+")) {
    Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] <- "plus"
  }
}
#change all UP_DIST values for - strand genes to "minus"
for (i in seq_along(Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST)) {
  if ((Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] == "calculate") &&
      (Pfalciparum_GFF3_allgene_strand_upto_1kb$STRAND[i] == "-")) {
    Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] <- "minus"
  }
}
#calculate upstream dist for + strand genes
#if STRAND = +, UP_DIST = START - END of line previous 
for (i in seq_along(Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST)) {
  if (Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] == "plus") {
    Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] <- Pfalciparum_GFF3_allgene_strand_upto_1kb$START[i] - 
      Pfalciparum_GFF3_allgene_strand_upto_1kb$END[i-1]
  }
}
#calculate upstream dist for - strand genes
#if STRAND = -, UP_DIST = START of line ahead - END
for (i in seq_along(Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST)) {
  if (Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] == "minus") {
    Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST[i] <- Pfalciparum_GFF3_allgene_strand_upto_1kb$START[i+1] - 
      Pfalciparum_GFF3_allgene_strand_upto_1kb$END[i]
  }
}
##convert values to integers
Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST <- as.integer(Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST)

    
    
#####1kb organize upstream cutoff values#####
    
#generate edited Up_DIST to cut off at 2KB
Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST_EDIT <- ifelse(Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST > 1001, 1001, Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST)
  #remove genes with UP_DIST <= 2 (no distance between genes)
  pos_Pfalciparum_GFF3_allgene_strand_upto_1kb <- subset(Pfalciparum_GFF3_allgene_strand_upto_1kb, Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST_EDIT >= 3)
    write.csv(pos_Pfalciparum_GFF3_allgene_strand_upto_1kb, "pos_Pfalciparum_GFF3_allgene_strand_upto_1kb.csv")
  neg_Pfalciparum_GFF3_allgene_strand_upto_1kb <- subset(Pfalciparum_GFF3_allgene_strand_upto_1kb, Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST_EDIT <= 2)
    write.csv(neg_Pfalciparum_GFF3_allgene_strand_upto_1kb, "neg_Pfalciparum_GFF3_allgene_strand_upto_1kb.csv")

##all genes in Niggi's list
allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_strand_upto_1kb, Pfalciparum_GFF3_allgene_strand_upto_1kb$GENEID %in% allgenes_Niggi_noDups_clean)
  #we lose 7 genes that are in allgenes_Niggi_noDups_clean (5680) but aren't in the PlasmoDB list (since allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info has 5673)
  #we lose 14 genes that are in the PlasmoDB list (5687) but not in allgenes_Niggi_noDups_clean, end up with 5673 in allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info
    
#want to remove the 53 genes with UP_DIST_EDIT <= 2
pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info, allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 3)
pos_allgenes_upto_1kb_strand_Niggi_noDups <- pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID  
    write.csv(pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info, "pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info.csv")
neg_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info, allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT <= 2)
neg_allgenes_upto_1kb_strand_Niggi_noDups <- neg_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID  
    write.csv(neg_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info, "neg_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info.csv")
    #we lose 53 genes with UP_DIST_EDIT <= 2, mostly in the mitos
    #use pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info GENEID list moving forwards
pos4_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 5)
  pos4_allgenes_upto_1kb_strand_Niggi_noDups <- pos4_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos5_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 6)
  pos5_allgenes_upto_1kb_strand_Niggi_noDups <- pos5_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos6_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 7)
  pos6_allgenes_upto_1kb_strand_Niggi_noDups <- pos6_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos7_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 8)
  pos7_allgenes_upto_1kb_strand_Niggi_noDups <- pos7_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos8_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 9)
  pos8_allgenes_upto_1kb_strand_Niggi_noDups <- pos8_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
    
##background
background_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_strand_upto_1kb, Pfalciparum_GFF3_allgene_strand_upto_1kb$GENEID %in% background_Niggi_noDups_clean)
    #we lose 7 genes that are in background_Niggi_noDups_clean (5341) but aren't in the PlasmoDB list (background_upto_1kb_strand_Niggi_noDups_GFF3_info is 5334)
    
pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(background_upto_1kb_strand_Niggi_noDups_GFF3_info, background_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 3)
  pos_background_upto_1kb_strand_Niggi_noDups <- pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
neg_background_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(background_upto_1kb_strand_Niggi_noDups_GFF3_info, background_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT <= 2)
    #we lose 53 genes with UP_DIST_EDIT <= 2, mostly in the mitos
    
    #use pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info GENEID list moving forwards
pos4_background_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 5)
    pos4_background_upto_1kb_strand_Niggi_noDups <- pos4_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos5_background_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 6)
    pos5_background_upto_1kb_strand_Niggi_noDups <- pos5_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos6_background_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 7)
    pos6_background_upto_1kb_strand_Niggi_noDups <- pos6_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos7_background_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 8)
    pos7_background_upto_1kb_strand_Niggi_noDups <- pos7_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos8_background_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 9)
    pos8_background_upto_1kb_strand_Niggi_noDups <- pos8_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
    
##background no downreg  
  background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_strand_upto_1kb, Pfalciparum_GFF3_allgene_strand_upto_1kb$GENEID %in% background_nodownreg_Niggi_noDups_clean)
    #we lose 7 genes that are in background_nodownreg_Niggi_noDups_clean (5299) but aren't in the PlasmoDB list (background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info is 5292)
    
pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info, background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 3)
pos_background_nodownreg_upto_1kb_strand_Niggi_noDups <- pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID  
neg_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info, background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT <= 2)
    #we lose 53 genes with UP_DIST_EDIT <= 2, mostly in the mitos
    
    #use pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info GENEID list moving forwards
pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 5)
    pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups <- pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 6)
    pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups <- pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 7)
    pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups <- pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 8)
    pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups <- pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info, pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$UP_DIST_EDIT >= 9)
    pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups <- pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID
    
    ##upreg and downreg
    #these don't need to be split into pos and neg, are all pos, also none missing from either Niggi's list or the PlasmoDB list
    upreg_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_strand_upto_1kb, Pfalciparum_GFF3_allgene_strand_upto_1kb$GENEID %in% upreg_Niggi_noDups_clean)
    downreg_upto_1kb_strand_Niggi_noDups_GFF3_info <- subset(Pfalciparum_GFF3_allgene_strand_upto_1kb, Pfalciparum_GFF3_allgene_strand_upto_1kb$GENEID %in% downreg_Niggi_noDups_clean)    
    
    
    #the next step is to figure out the beginning and end of the upstream region for both the + strand genes (5' of the gene start) and
    #the -ve strand genes (based on coordinate system, 3' of the gene end)
    
# +ve gene --> UPSTART = GENE_START - UP_DIST_EDIT 
#UPEND = GENE_START - 1
# -ve gene --> UPSTART = GENE_END 
#UPEND = GENE_END + (UP_DIST_EDIT - 1)
pos_Pfalciparum_upstream_regions_upto_1kb_strand <- data.frame(matrix(NA, nrow = 5632))
pos_Pfalciparum_upstream_regions_upto_1kb_strand$CHROM <-  pos_Pfalciparum_GFF3_allgene_strand_upto_1kb$SEQID
pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID <- pos_Pfalciparum_GFF3_allgene_strand_upto_1kb$GENEID
pos_Pfalciparum_upstream_regions_upto_1kb_strand <- pos_Pfalciparum_upstream_regions_upto_1kb_strand[,-1]
pos_Pfalciparum_upstream_regions_upto_1kb_strand$STRAND <- pos_Pfalciparum_GFF3_allgene_strand_upto_1kb$STRAND
pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENE_START <- pos_Pfalciparum_GFF3_allgene_strand_upto_1kb$START
pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENE_END <- pos_Pfalciparum_GFF3_allgene_strand_upto_1kb$END
pos_Pfalciparum_upstream_regions_upto_1kb_strand$UP_DIST_EDIT <- pos_Pfalciparum_GFF3_allgene_strand_upto_1kb$UP_DIST_EDIT
pos_Pfalciparum_upstream_regions_upto_1kb_strand$UPSTART <- NA
pos_Pfalciparum_upstream_regions_upto_1kb_strand$UPEND <- NA

for (i in seq_along(pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID)) {
  if (pos_Pfalciparum_upstream_regions_upto_1kb_strand$STRAND[i] == "+") {
    pos_Pfalciparum_upstream_regions_upto_1kb_strand$UPSTART[i] <- pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENE_START[i] - pos_Pfalciparum_upstream_regions_upto_1kb_strand$UP_DIST_EDIT[i]
    pos_Pfalciparum_upstream_regions_upto_1kb_strand$UPEND[i] <- pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENE_START[i] - 1
  }
}
    
for (i in seq_along(pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID)) {
  if (pos_Pfalciparum_upstream_regions_upto_1kb_strand$STRAND[i] == "-") {
    pos_Pfalciparum_upstream_regions_upto_1kb_strand$UPSTART[i] <- pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENE_END[i]
    pos_Pfalciparum_upstream_regions_upto_1kb_strand$UPEND[i] <- pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENE_END[i] + (pos_Pfalciparum_upstream_regions_upto_1kb_strand$UP_DIST_EDIT[i] - 1)
  }
}
    
    write.csv(pos_Pfalciparum_upstream_regions_upto_1kb_strand, "pos_Pfalciparum_upstream_regions_upto_1kb_strand.csv")
    
#####1kb split upstream regions into the specific gene lists#####
    ##allgenes
    pos_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos4_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos4_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos5_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos5_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos6_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos6_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos7_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos7_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos8_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos8_allgenes_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    
    ##background
    pos_background_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos4_background_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos4_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos5_background_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos5_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos6_background_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos6_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos7_background_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos7_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos8_background_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos8_background_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    
    ##background nodownreg
    pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    
    
    ##upreg and downreg
    upreg_upto_1kb_strand_Niggi_noDups_upstream_regions <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% upreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    downreg_upto_1kb_strand_Niggi_noDups_upstream_regions  <- subset(pos_Pfalciparum_upstream_regions_upto_1kb_strand, pos_Pfalciparum_upstream_regions_upto_1kb_strand$GENEID %in% downreg_upto_1kb_strand_Niggi_noDups_GFF3_info$GENEID)
    
    
#####1kb save lists of GENEids for each category#####
    write.csv(upreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "upreg_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(downreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "downreg_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    
    write.csv(pos_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos_background_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos4_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos4_background_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos5_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos5_background_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos6_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos6_background_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos7_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos7_background_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos8_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos8_background_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    
    write.csv(pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    
    write.csv(pos_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos4_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos4_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos5_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos5_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos6_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos6_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos7_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos7_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    write.csv(pos8_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID, "pos8_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions_genelist.csv")
    
    
    
    
#####1kb create and download text for bedtools BED files#####
    
#same concept as previous - use the coordinates in the upstream_regions 
#dataframes to create BED files for each subset of genes
# BED columns: CHROM START END NAME

##lists to make
pos_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions
pos4_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions
pos5_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions
pos6_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions
pos7_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions
pos8_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions
pos_background_upto_1kb_strand_Niggi_noDups_upstream_regions
pos4_background_upto_1kb_strand_Niggi_noDups_upstream_regions
pos5_background_upto_1kb_strand_Niggi_noDups_upstream_regions
pos6_background_upto_1kb_strand_Niggi_noDups_upstream_regions
pos7_background_upto_1kb_strand_Niggi_noDups_upstream_regions
pos8_background_upto_1kb_strand_Niggi_noDups_upstream_regions
pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions
pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions
pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions
pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions
pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions
pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions
upreg_upto_1kb_strand_Niggi_noDups_upstream_regions
downreg_upto_1kb_strand_Niggi_noDups_upstream_regions
    
    
    
# BED columns: CHROM START END NAME
    
upreg_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 339))
    upreg_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- upreg_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    upreg_upto_1kb_strand_Niggi_noDups_download_BED$START <- upreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    upreg_upto_1kb_strand_Niggi_noDups_download_BED$END <- upreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    upreg_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- upreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    upreg_upto_1kb_strand_Niggi_noDups_download_BED <- upreg_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]
    write.csv(upreg_upto_1kb_strand_Niggi_noDups_download_BED, "upreg_upto_1kb_strand_Niggi_noDups_download_BED.csv")
    
downreg_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 42))
    downreg_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- downreg_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    downreg_upto_1kb_strand_Niggi_noDups_download_BED$START <- downreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    downreg_upto_1kb_strand_Niggi_noDups_download_BED$END <- downreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    downreg_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- downreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    downreg_upto_1kb_strand_Niggi_noDups_download_BED <- downreg_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]
    write.csv(downreg_upto_1kb_strand_Niggi_noDups_download_BED, "downreg_upto_1kb_strand_Niggi_noDups_download_BED.csv")
    
pos_background_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5281))
    pos_background_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos_background_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos_background_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos_background_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos_background_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos_background_upto_1kb_strand_Niggi_noDups_download_BED <- pos_background_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos_background_upto_1kb_strand_Niggi_noDups_download_BED, "pos_background_upto_1kb_strand_Niggi_noDups_download_BED.csv")
pos4_background_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5279))
    pos4_background_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos4_background_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos4_background_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos4_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos4_background_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos4_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos4_background_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos4_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos4_background_upto_1kb_strand_Niggi_noDups_download_BED <- pos4_background_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos4_background_upto_1kb_strand_Niggi_noDups_download_BED, "pos4_background_upto_1kb_strand_Niggi_noDups_download_BED.csv")
pos5_background_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5278))
    pos5_background_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos5_background_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos5_background_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos5_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos5_background_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos5_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos5_background_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos5_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos5_background_upto_1kb_strand_Niggi_noDups_download_BED <- pos5_background_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos5_background_upto_1kb_strand_Niggi_noDups_download_BED, "pos5_background_upto_1kb_strand_Niggi_noDups_download_BED.csv")
pos6_background_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5273))
    pos6_background_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos6_background_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos6_background_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos6_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos6_background_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos6_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos6_background_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos6_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos6_background_upto_1kb_strand_Niggi_noDups_download_BED <- pos6_background_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos6_background_upto_1kb_strand_Niggi_noDups_download_BED, "pos6_background_upto_1kb_strand_Niggi_noDups_download_BED.csv")
pos7_background_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5270))
    pos7_background_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos7_background_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos7_background_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos7_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos7_background_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos7_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos7_background_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos7_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos7_background_upto_1kb_strand_Niggi_noDups_download_BED <- pos7_background_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos7_background_upto_1kb_strand_Niggi_noDups_download_BED, "pos7_background_upto_1kb_strand_Niggi_noDups_download_BED.csv")
    pos8_background_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5267))
pos8_background_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos8_background_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos8_background_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos8_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos8_background_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos8_background_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos8_background_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos8_background_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos8_background_upto_1kb_strand_Niggi_noDups_download_BED <- pos8_background_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos8_background_upto_1kb_strand_Niggi_noDups_download_BED, "pos8_background_upto_1kb_strand_Niggi_noDups_download_BED.csv")
    
pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5239))  
    pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED, "pos_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED.csv")
pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5237))
    pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED, "pos4_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED.csv")
pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5236))
    pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED, "pos5_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED.csv")
pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5231))
    pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED, "pos6_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED.csv")
pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5228))
    pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED, "pos7_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED.csv")
pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5225))
    pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED <- pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED, "pos8_background_nodownreg_upto_1kb_strand_Niggi_noDups_download_BED.csv")
    
    pos_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5620))
    pos_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- pos_allgenes_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]
    write.csv(pos_allgenes_upto_1kb_strand_Niggi_noDups_download_BED, "pos_allgenes_upto_1kb_strand_Niggi_noDups_download_BED.csv")
    pos4_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5618))
    pos4_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos4_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos4_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos4_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos4_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos4_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos4_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos4_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos4_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- pos4_allgenes_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos4_allgenes_upto_1kb_strand_Niggi_noDups_download_BED, "pos4_allgenes_upto_1kb_strand_Niggi_noDups_download_BED.csv")
    pos5_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5617))
    pos5_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos5_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos5_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos5_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos5_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos5_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos5_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos5_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos5_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- pos5_allgenes_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos5_allgenes_upto_1kb_strand_Niggi_noDups_download_BED, "pos5_allgenes_upto_1kb_strand_Niggi_noDups_download_BED.csv")
    pos6_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5612))
    pos6_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos6_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos6_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos6_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos6_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos6_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos6_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos6_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos6_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- pos6_allgenes_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos6_allgenes_upto_1kb_strand_Niggi_noDups_download_BED, "pos6_allgenes_upto_1kb_strand_Niggi_noDups_download_BED.csv")
    pos7_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5609))
    pos7_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos7_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos7_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos7_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos7_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos7_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos7_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos7_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos7_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- pos7_allgenes_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos7_allgenes_upto_1kb_strand_Niggi_noDups_download_BED, "pos7_allgenes_upto_1kb_strand_Niggi_noDups_download_BED.csv")
    pos8_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- data.frame(matrix(nrow = 5606))
    pos8_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$CHROM <- pos8_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$CHROM
    pos8_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$START <- pos8_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPSTART
    pos8_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$END <- pos8_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$UPEND
    pos8_allgenes_upto_1kb_strand_Niggi_noDups_download_BED$NAME <- pos8_allgenes_upto_1kb_strand_Niggi_noDups_upstream_regions$GENEID
    pos8_allgenes_upto_1kb_strand_Niggi_noDups_download_BED <- pos8_allgenes_upto_1kb_strand_Niggi_noDups_download_BED[,2:5]  
    write.csv(pos8_allgenes_upto_1kb_strand_Niggi_noDups_download_BED, "pos8_allgenes_upto_1kb_strand_Niggi_noDups_download_BED.csv")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    

    
#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/LysoPC_upstreamReg_org.RData")
  