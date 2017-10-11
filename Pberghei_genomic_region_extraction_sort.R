#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/Pberghei_upstream_region_sort.RData")
library("data.table")

#####set WD#####
setwd("J:/III/Waters/Group Members/Mallory/PbergheiUpstreamRegionsSebastian/")

#####load data#####
Pberghei_GFF3_allgene_strand_calc <- read.csv("GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv", header=T)
Pberghei_GFF3_allgene_strand_calc_1000 <- read.csv("GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv", header=T)
Pberghei_chrom_coord <- read.csv("Pberghei_chromosome_coord.csv")

#####calculate strand specific upstream dist#####
###start with the beginning and ending of CHROM
#beginning if the first gene is on the + strand, end if the last gene is on the - strand

###label the start of each chromosome
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST[2:5122])) {
if (identical(Pberghei_GFF3_allgene_strand_calc$SEQID[i], Pberghei_GFF3_allgene_strand_calc$SEQID[i-1])) {
  Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- "calculate"
} else {
  Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- "START"
}
}

for (i in seq_along(Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[2:5122])) {
  if (identical(Pberghei_GFF3_allgene_strand_calc_1000$SEQID[i], Pberghei_GFF3_allgene_strand_calc_1000$SEQID[i-1])) {
    Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] <- "calculate"
  } else {
    Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] <- "START"
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

for (i in seq_along(Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[2:5120])) {
  if (identical(Pberghei_GFF3_allgene_strand_calc_1000$SEQID[i], Pberghei_GFF3_allgene_strand_calc_1000$SEQID[i+1]) && (Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] == "calculate")) {
    Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] <- "calculate"
  } else {
    if (Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] != "START") {
      Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] <- "END"
    }
  }
}

###at start of CHROM is gene is on + strand, calculate upstream dist
#if gene is on - strand will go through normal calculation
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
if ((Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "START") && (Pberghei_GFF3_allgene_strand_calc$STRAND[i] == "+")) {
  Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- Pberghei_GFF3_allgene_strand_calc$START[i] - 1
} 
}

for (i in seq_along(Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] == "START") && (Pberghei_GFF3_allgene_strand_calc_1000$STRAND[i] == "+")) {
    Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] <- Pberghei_GFF3_allgene_strand_calc_1000$START[i] - 1
  } 
}


###at end of CHROM if gene is on - strand, need to calculate upstream dist based on size of CHROM
#if gene is on + strand will go through normal calc
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "END") && (Pberghei_GFF3_allgene_strand_calc$STRAND[i] == "-")) {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- (Pberghei_chrom_coord$END[Pberghei_chrom_coord$CHROM == Pberghei_GFF3_allgene_strand_calc$SEQID[i]]) - Pberghei_GFF3_allgene_strand_calc$END[i]
  } 
}


for (i in seq_along(Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] == "END") && (Pberghei_GFF3_allgene_strand_calc_1000$STRAND[i] == "-")) {
    Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] <- (Pberghei_chrom_coord$END[Pberghei_chrom_coord$CHROM == Pberghei_GFF3_allgene_strand_calc_1000$SEQID[i]]) - Pberghei_GFF3_allgene_strand_calc_1000$END[i]
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

for (i in seq_along(Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] == "calculate") ||
      (Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] == "START") ||
      (Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] == "END") ) {
    Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] <- "calculate"
  }
}

#change all UP_DIST values for + strand genes to "plus"
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "calculate") &&
      (Pberghei_GFF3_allgene_strand_calc$STRAND[i] == "+")) {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- "plus"
}
}

for (i in seq_along(Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] == "calculate") &&
      (Pberghei_GFF3_allgene_strand_calc_1000$STRAND[i] == "+")) {
    Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] <- "plus"
  }
}

#change all UP_DIST values for - strand genes to "minus"
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "calculate") &&
      (Pberghei_GFF3_allgene_strand_calc$STRAND[i] == "-")) {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- "minus"
  }
}

for (i in seq_along(Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST)) {
  if ((Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] == "calculate") &&
      (Pberghei_GFF3_allgene_strand_calc_1000$STRAND[i] == "-")) {
    Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] <- "minus"
  }
}
#calculate upstream dist for + strand genes
#if STRAND = +, UP_DIST = START - END of line previous 
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if (Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "plus") {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- Pberghei_GFF3_allgene_strand_calc$START[i] - 
      Pberghei_GFF3_allgene_strand_calc$END[i-1]
  }
}

for (i in seq_along(Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST)) {
  if (Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] == "plus") {
    Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] <- Pberghei_GFF3_allgene_strand_calc_1000$START[i] - 
      Pberghei_GFF3_allgene_strand_calc_1000$END[i-1]
  }
}

#calculate upstream dist for - strand genes
#if STRAND = -, UP_DIST = START of line ahead - END
for (i in seq_along(Pberghei_GFF3_allgene_strand_calc$UP_DIST)) {
  if (Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] == "minus") {
    Pberghei_GFF3_allgene_strand_calc$UP_DIST[i] <- Pberghei_GFF3_allgene_strand_calc$START[i+1] - 
      Pberghei_GFF3_allgene_strand_calc$END[i]
  }
}

for (i in seq_along(Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST)) {
  if (Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] == "minus") {
    Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST[i] <- Pberghei_GFF3_allgene_strand_calc_1000$START[i+1] - 
      Pberghei_GFF3_allgene_strand_calc_1000$END[i]
  }
}

##convert values to integers
Pberghei_GFF3_allgene_strand_calc$UP_DIST <- as.integer(Pberghei_GFF3_allgene_strand_calc$UP_DIST)
Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST <- as.integer(Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST)


#####organize upstream cutoff values#####
#generate edited Up_DIST to cut off at 2KB
Pberghei_GFF3_allgene_strand_calc$UP_DIST_EDIT <- ifelse(Pberghei_GFF3_allgene_strand_calc$UP_DIST > 2001, 2001, Pberghei_GFF3_allgene_strand_calc$UP_DIST)
Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST_EDIT <- ifelse(Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST > 1001, 1001, Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST)



#remove genes with an UP_DIST <1 (no bases between genes)
pos_Pberghei_GFF3_allgene_strand <- subset(Pberghei_GFF3_allgene_strand_calc, Pberghei_GFF3_allgene_strand_calc$UP_DIST_EDIT >= 3)
  write.csv(pos_Pberghei_GFF3_allgene_strand, "pos_Pberghei_GFF3_allgene_strand.csv")
neg_Pberghei_GFF3_allgene_strand <- subset(Pberghei_GFF3_allgene_strand_calc, Pberghei_GFF3_allgene_strand_calc$UP_DIST_EDIT <= 2)
  write.csv(neg_Pberghei_GFF3_allgene_strand, "neg_Pberghei_GFF3_allgene_strand.csv")

pos_Pberghei_GFF3_allgene_strand_1000 <- subset(Pberghei_GFF3_allgene_strand_calc_1000, Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST_EDIT >= 3)
  write.csv(pos_Pberghei_GFF3_allgene_strand_1000, "pos_Pberghei_GFF3_allgene_strand_1000.csv")
neg_Pberghei_GFF3_allgene_strand_1000 <- subset(Pberghei_GFF3_allgene_strand_calc_1000, Pberghei_GFF3_allgene_strand_calc_1000$UP_DIST_EDIT <= 2)
  write.csv(neg_Pberghei_GFF3_allgene_strand_1000, "neg_Pberghei_GFF3_allgene_strand_1000.csv")  

  
  
#the next step is to figure out the beginning and end of the upstream region for both the + strand genes (5' of the gene start) and
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




Pberghei_upstream_regions_strand_1000 <- data.frame(matrix(NA, nrow = 5073))
Pberghei_upstream_regions_strand_1000$CHROM <-  pos_Pberghei_GFF3_allgene_strand_1000$SEQID
Pberghei_upstream_regions_strand_1000$GENEID <- pos_Pberghei_GFF3_allgene_strand_1000$GENEID
Pberghei_upstream_regions_strand_1000 <- Pberghei_upstream_regions_strand_1000[,-1]
Pberghei_upstream_regions_strand_1000$STRAND <- pos_Pberghei_GFF3_allgene_strand_1000$STRAND
Pberghei_upstream_regions_strand_1000$GENE_START <- pos_Pberghei_GFF3_allgene_strand_1000$START
Pberghei_upstream_regions_strand_1000$GENE_END <- pos_Pberghei_GFF3_allgene_strand_1000$END
Pberghei_upstream_regions_strand_1000$UP_DIST_EDIT <- pos_Pberghei_GFF3_allgene_strand_1000$UP_DIST_EDIT
Pberghei_upstream_regions_strand_1000$UPSTART <- NA
Pberghei_upstream_regions_strand_1000$UPEND <- NA


for (i in seq_along(Pberghei_upstream_regions_strand_1000$GENEID)) {
  if (Pberghei_upstream_regions_strand_1000$STRAND[i] == "+") {
    Pberghei_upstream_regions_strand_1000$UPSTART[i] <- Pberghei_upstream_regions_strand_1000$GENE_START[i] - Pberghei_upstream_regions_strand_1000$UP_DIST_EDIT[i]
    Pberghei_upstream_regions_strand_1000$UPEND[i] <- Pberghei_upstream_regions_strand_1000$GENE_START[i] - 1
  }
}

for (i in seq_along(Pberghei_upstream_regions_strand_1000$GENEID)) {
  if (Pberghei_upstream_regions_strand_1000$STRAND[i] == "-") {
    Pberghei_upstream_regions_strand_1000$UPSTART[i] <- Pberghei_upstream_regions_strand_1000$GENE_END[i]
    Pberghei_upstream_regions_strand_1000$UPEND[i] <- Pberghei_upstream_regions_strand_1000$GENE_END[i] + (Pberghei_upstream_regions_strand_1000$UP_DIST_EDIT[i] - 1)
  }
}

write.csv(Pberghei_upstream_regions_strand_1000$GENEID, "Pberghei_upstream_regions_strand_1000_genelist.csv")


#####create and download text for bedtools BED files#####
# BED columns: CHROM START END NAME

Pberghei_upstream_regions_strand_download_BED <- data.frame(matrix(nrow = 5073))
Pberghei_upstream_regions_strand_download_BED$CHROM <- Pberghei_upstream_regions_strand$CHROM
Pberghei_upstream_regions_strand_download_BED$START <- Pberghei_upstream_regions_strand$UPSTART
Pberghei_upstream_regions_strand_download_BED$END <- Pberghei_upstream_regions_strand$UPEND
Pberghei_upstream_regions_strand_download_BED$NAME <- Pberghei_upstream_regions_strand$GENEID
Pberghei_upstream_regions_strand_download_BED <- Pberghei_upstream_regions_strand_download_BED[,2:5]
write.csv(Pberghei_upstream_regions_strand_download_BED, "Pberghei_upstream_regions_strand_download_BED.csv")




Pberghei_upstream_regions_strand_1000_download_BED <- data.frame(matrix(nrow = 5073))
Pberghei_upstream_regions_strand_1000_download_BED$CHROM <- Pberghei_upstream_regions_strand_1000$CHROM
Pberghei_upstream_regions_strand_1000_download_BED$START <- Pberghei_upstream_regions_strand_1000$UPSTART
Pberghei_upstream_regions_strand_1000_download_BED$END <- Pberghei_upstream_regions_strand_1000$UPEND
Pberghei_upstream_regions_strand_1000_download_BED$NAME <- Pberghei_upstream_regions_strand_1000$GENEID
Pberghei_upstream_regions_strand_1000_download_BED <- Pberghei_upstream_regions_strand_1000_download_BED[,2:5]
write.csv(Pberghei_upstream_regions_strand_1000_download_BED, "Pberghei_upstream_regions_strand_1000_download_BED.csv")


#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/Pberghei_upstream_region_sort.RData")
