#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/Pberghei_upstream_region_sort.RData")
library("data.table")

#####set WD#####
setwd("J:/III/Waters/Group Members/Mallory/PbergheiUpstreamRegionsSebastian/")

#####load data#####
###to do this I've already manually parsed out the gene coordinates from the GFF3 file in excel, file has 5 columns: SEQID (CHROM), GENEID, START, END, STRAND
###also manually parsed out CHROM coord, file has 3 columns: CHROM, START, END

Pberghei_GFF3_allgene_strand_calc <- read.csv("GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv", header=T)
Pberghei_chrom_coord <- read.csv("Pberghei_chromosome_coord.csv")

#####calculate strand specific upstream dist#####
###start with the beginning and ending of CHROM
#beginning of CHROM if the first gene is on the +ve strand, end of CHROm if the last gene is on the -ve strand

###label the start of each chromosome
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST[2:5122])) {
  if (identical(Pberghei_GFF3_allgene_strand_calc$SEQID[i], Pberghei_GFF3_allgene_strand_calc$SEQID[i-1])) {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- "calculate"
  } else {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- "START"
  }
}

###label the end of each chromosome
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST[2:5120])) {
  if (identical(Pberghei_GFF3_allgene_strand_calc$SEQID[i], Pberghei_GFF3_allgene_strand_calc$SEQID[i+1]) && (Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "calculate")) {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- "calculate"
  } else {
    if (Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] != "START") {
      Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- "END"
    }
  }
}

###at start of CHROM if gene is on +ve strand, calculate upstream dist
#if gene is on -ve strand will go through normal calculation
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "START") && (Pberghei_GFF3_allgene_strand_calc$STRAND[i] == "+")) {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- Pberghei_GFF3_allgene_strand_calc$START[i] - 1
  } 
}

###at end of CHROM if gene is on -ve strand, need to calculate upstream dist based on size of CHROM
#if gene is on +ve strand will go through normal calc
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "END") && (Pberghei_GFF3_allgene_strand_calc$STRAND[i] == "-")) {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- (Pberghei_chrom_coord$END[Pberghei_chrom_coord$CHROM == Pberghei_GFF3_allgene_strand_calc$SEQID[i]]) - Pberghei_GFF3_allgene_strand_calc$END[i]
  } 
}

#then move on to populating the rest of the values, ignoring the rows that have already been filled in
#change all UP_DIST values that aren't numbers to "calculate"
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "calculate") ||
      (Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "START") ||
      (Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "END") ) {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- "calculate"
  }
}

#change all UP_DIST values for +ve strand genes to "plus"
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "calculate") &&
      (Pberghei_GFF3_allgene_strand_calc$STRAND[i] == "+")) {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- "plus"
  }
}

#change all UP_DIST values for -ve strand genes to "minus"
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "calculate") &&
      (Pberghei_GFF3_allgene_strand_calc$STRAND[i] == "-")) {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- "minus"
  }
}

#calculate upstream dist for +ve strand genes
#if STRAND = +ve, UP_DIST = START - END of line previous 
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if (Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "plus") {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- Pberghei_GFF3_allgene_strand_calc$START[i] - 
      Pberghei_GFF3_allgene_strand_calc$END[i-1]
  }
}

#calculate upstream dist for -ve strand genes
#if STRAND = -ve, UP_DIST = START of line ahead - END
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if (Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "minus") {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- Pberghei_GFF3_allgene_strand_calc$START[i+1] - 
      Pberghei_GFF3_allgene_strand_calc$END[i]
  }
}
##convert values to integers
Pberghei_GFF3_allgene_strand_calc$UP_DIST <- as.integer(Pberghei_GFF3_allgene_strand_calc$UP_DIST)

#####organize upstream cutoff values#####
#generate edited Up_DIST to cut off at 2KB
Pberghei_GFF3_allgene_strand_calc$UP_DIST_EDIT <- ifelse(Pberghei_GFF3_allgene_strand_calc$UP_DIST > 2001, 2001, Pberghei_GFF3_allgene_strand_calc$UP_DIST)

#remove genes with an UP_DIST <1 (no bases between genes)
pos_Pberghei_GFF3_allgene_strand <- subset(Pberghei_GFF3_allgene_strand_calc, Pberghei_GFF3_allgene_strand_calc$UP_DIST_EDIT >= 3)
write.csv(pos_Pberghei_GFF3_allgene_strand, "pos_Pberghei_GFF3_allgene_strand.csv")
neg_Pberghei_GFF3_allgene_strand <- subset(Pberghei_GFF3_allgene_strand_calc, Pberghei_GFF3_allgene_strand_calc$UP_DIST_EDIT <= 2)
write.csv(neg_Pberghei_GFF3_allgene_strand, "neg_Pberghei_GFF3_allgene_strand.csv")

#the next step is to figure out the beginning and end of the upstream region for both the +ve strand genes (5' of the gene start) and
#the -ve strand genes (based on coordinate system, 3' of the gene end)

# +ve gene --> UPSTART = GENE_START - UP_DIST_EDIT 
#UPEND = GENE_START - 1
# -ve gene --> UPSTART = GENE_END 
#UPEND = GENE_END + (UP_DIST_EDIT - 1)

Pberghei_upstream_regions_strand <- data.frame(matrix(NA, nrow = 5073))
Pberghei_upstream_regions_strand$CHROM <-  pos_Pberghei_GFF3_allgene_strand$SEQID
Pberghei_upstream_regions_strand$GENEID <- pos_Pberghei_GFF3_allgene_strand$GENEID
Pberghei_upstream_regions_strand <- Pberghei_upstream_regions_strand[,-1]
Pberghei_upstream_regions_strand$STRAND <- pos_Pberghei_GFF3_allgene_strand$STRAND
Pberghei_upstream_regions_strand$GENE_START <- pos_Pberghei_GFF3_allgene_strand$START
Pberghei_upstream_regions_strand$GENE_END <- pos_Pberghei_GFF3_allgene_strand$END
Pberghei_upstream_regions_strand$UP_DIST_EDIT <- pos_Pberghei_GFF3_allgene_strand$UP_DIST_EDIT
Pberghei_upstream_regions_strand$UPSTART <- NA
Pberghei_upstream_regions_strand$UPEND <- NA


for (i in seq_along(Pberghei_upstream_regions_strand$GENEID)) {
  if (Pberghei_upstream_regions_strand$STRAND[i] == "+") {
    Pberghei_upstream_regions_strand$UPSTART[i] <- Pberghei_upstream_regions_strand$GENE_START[i] - Pberghei_upstream_regions_strand$UP_DIST_EDIT[i]
    Pberghei_upstream_regions_strand$UPEND[i] <- Pberghei_upstream_regions_strand$GENE_START[i] - 1
  }
}

for (i in seq_along(Pberghei_upstream_regions_strand$GENEID)) {
  if (Pberghei_upstream_regions_strand$STRAND[i] == "-") {
    Pberghei_upstream_regions_strand$UPSTART[i] <- Pberghei_upstream_regions_strand$GENE_END[i]
    Pberghei_upstream_regions_strand$UPEND[i] <- Pberghei_upstream_regions_strand$GENE_END[i] + (Pberghei_upstream_regions_strand$UP_DIST_EDIT[i] - 1)
  }
}

write.csv(Pberghei_upstream_regions_strand$GENEID, "Pberghei_upstream_regions_strand_genelist.csv")

#####create and download text for bedtools BED files#####
# BED columns: CHROM START END NAME

Pberghei_upstream_regions_strand_download_BED <- data.frame(matrix(nrow = 5073))
Pberghei_upstream_regions_strand_download_BED$CHROM <- Pberghei_upstream_regions_strand$CHROM
Pberghei_upstream_regions_strand_download_BED$START <- Pberghei_upstream_regions_strand$UPSTART
Pberghei_upstream_regions_strand_download_BED$END <- Pberghei_upstream_regions_strand$UPEND
Pberghei_upstream_regions_strand_download_BED$NAME <- Pberghei_upstream_regions_strand$GENEID
Pberghei_upstream_regions_strand_download_BED <- Pberghei_upstream_regions_strand_download_BED[,2:5]
write.csv(Pberghei_upstream_regions_strand_download_BED, "Pberghei_upstream_regions_strand_download_BED.csv")

#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/Pberghei_upstream_region_sort.RData")