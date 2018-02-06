#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/NishaRIPseq_up_down_ORF_region_extract.RData")
  library("data.table")

#####set WD#####
setwd("J:/III/Waters/Group Members/Mallory/NishaRNAseq/DataDownloads/")

###load data#####
###to do this I've already manually parsed out the gene coordinates from the GFF3 file in excel, file has 5 columns: SEQID (CHROM), GENEID, START, END, STRAND
###also manually parsed out CHROM coord, file has 3 columns: CHROM, START, END

PbANKA_33_gene_coord <- read.csv("PbANKA_33_gff3_GeneCoord.csv", header=T)
PbANKA_33_chrom_coord <- read.csv("PbANKA_33_gff3_ChromCoord.csv", header=T)
RIPseq_GoI <- read.csv("RIPseq_GoIs.csv", header=F)
  RIPseq_GoI <- RIPseq_GoI[,1]
  
#####calculate strand specific up to 1KB upstream dist#####
###start with the beginning and ending of CHROM
#beginning of CHROM if the first gene is on the +ve strand, end of CHROm if the last gene is on the -ve strand

PbANKA_33_gene_upto1KB_UP_calc <- PbANKA_33_gene_coord
  PbANKA_33_gene_upto1KB_UP_calc$UP_DIST <- NA

###UP_DIST###
    ###label the start of each chromosome
    for (i in seq_along(PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[1:5245])) {
      if (identical(PbANKA_33_gene_upto1KB_UP_calc$SEQID[i], PbANKA_33_gene_upto1KB_UP_calc$SEQID[i-1])) {
        PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] <- "calculate"
      } else {
        PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] <- "START"
      }
    }
    ###label the end of each chromosome
    for (i in seq_along(PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[1:5245])) {
      if (identical(PbANKA_33_gene_upto1KB_UP_calc$SEQID[i], PbANKA_33_gene_upto1KB_UP_calc$SEQID[i+1]) && (PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] == "calculate")) {
        PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] <- "calculate"
      } else {
        if (PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] != "START") {
          PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] <- "END"
        }
      }
    }
  
  ###at start of CHROM if gene is on +ve strand, calculate upstream dist
      #if gene is on -ve strand will go through normal calculation
      for (i in seq_along(PbANKA_33_gene_upto1KB_UP_calc$UP_DIST)) {
        if ((PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] == "START") && (PbANKA_33_gene_upto1KB_UP_calc$STRAND[i] == "+")) {
          PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] <- PbANKA_33_gene_upto1KB_UP_calc$START[i] - 1
        } 
      }
  
  ###at end of CHROM if gene is on -ve strand, need to calculate upstream dist based on size of CHROM
      #if gene is on +ve strand will go through normal calc
      for (i in seq_along(PbANKA_33_gene_upto1KB_UP_calc$UP_DIST)) {
        if ((PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] == "END") && (PbANKA_33_gene_upto1KB_UP_calc$STRAND[i] == "-")) {
          PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] <- (PbANKA_33_chrom_coord$END[PbANKA_33_chrom_coord$SEQID == PbANKA_33_gene_upto1KB_UP_calc$SEQID[i]]) - PbANKA_33_gene_upto1KB_UP_calc$END[i]
        } 
      }
  
    #then move on to populating the rest of the values, ignoring the rows that have already been filled in
    #change all UP_DIST values that aren't numbers to "calculate"
    for (i in seq_along(PbANKA_33_gene_upto1KB_UP_calc$UP_DIST)) {
      if ((PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] == "calculate") ||
          (PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] == "START") ||
          (PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] == "END") ) {
        PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] <- "calculate"
      }
    }
    
    #change all UP_DIST values for +ve strand genes to "plus"
    for (i in seq_along(PbANKA_33_gene_upto1KB_UP_calc$UP_DIST)) {
      if ((PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] == "calculate") &&
          (PbANKA_33_gene_upto1KB_UP_calc$STRAND[i] == "+")) {
        PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] <- "plus"
      }
    }
    
    #change all UP_DIST values for -ve strand genes to "minus"
    for (i in seq_along(PbANKA_33_gene_upto1KB_UP_calc$UP_DIST)) {
      if ((PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] == "calculate") &&
          (PbANKA_33_gene_upto1KB_UP_calc$STRAND[i] == "-")) {
        PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] <- "minus"
      }
    }
    
    #calculate upstream dist for +ve strand genes
    #if STRAND = +ve, UP_DIST = START - END of line previous 
    for (i in seq_along(PbANKA_33_gene_upto1KB_UP_calc$UP_DIST)) {
      if (PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] == "plus") {
        PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] <- PbANKA_33_gene_upto1KB_UP_calc$START[i] - 
          PbANKA_33_gene_upto1KB_UP_calc$END[i-1]
      }
    }
    
    #calculate upstream dist for -ve strand genes
    #if STRAND = -ve, UP_DIST = START of line ahead - END
    for (i in seq_along(PbANKA_33_gene_upto1KB_UP_calc$UP_DIST)) {
      if (PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] == "minus") {
        PbANKA_33_gene_upto1KB_UP_calc$UP_DIST[i] <- PbANKA_33_gene_upto1KB_UP_calc$START[i+1] - 
          PbANKA_33_gene_upto1KB_UP_calc$END[i]
      }
    }
    ##convert values to integers
    PbANKA_33_gene_upto1KB_UP_calc$UP_DIST <- as.integer(PbANKA_33_gene_upto1KB_UP_calc$UP_DIST)
  
    
####organize up to 1KB upstream cutoff values#####
    #generate edited UP_DIST to cut off at 1KB
    PbANKA_33_gene_upto1KB_UP_calc$UP_DIST_EDIT <- ifelse(PbANKA_33_gene_upto1KB_UP_calc$UP_DIST > 1000, 1000, PbANKA_33_gene_upto1KB_UP_calc$UP_DIST)
    
    #remove genes with an UP_DIST <1 (no bases between genes)
    pos_PbANKA_33_gene_upto1KB_UP_strand <- subset(PbANKA_33_gene_upto1KB_UP_calc, PbANKA_33_gene_upto1KB_UP_calc$UP_DIST_EDIT >= 3)
      write.csv(pos_PbANKA_33_gene_upto1KB_UP_strand, "pos_PbANKA_33_gene_upto1KB_UP_strand.csv")
    
    neg_PbANKA_33_gene_upto1KB_UP_strand  <- subset(PbANKA_33_gene_upto1KB_UP_calc, PbANKA_33_gene_upto1KB_UP_calc$UP_DIST_EDIT <= 2)
      write.csv(neg_PbANKA_33_gene_upto1KB_UP_strand, "neg_PbANKA_33_gene_upto1KB_UP_strand.csv")
    
#####calculate strand specific up to 1KB downstream dist#####
      ###start with the beginning and ending of CHROM
      #beginning of CHROM if the first gene is on the -ve strand, end of CHROM if the last gene is on the +ve strand
      
      PbANKA_33_gene_upto1KB_DOWN_calc <- PbANKA_33_gene_coord
      PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST <- NA
      
    ###DOWN_DIST###
    ###label the start of each chromosome
    for (i in seq_along(PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[1:5245])) {
      if (identical(PbANKA_33_gene_upto1KB_DOWN_calc$SEQID[i], PbANKA_33_gene_upto1KB_DOWN_calc$SEQID[i-1])) {
        PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] <- "calculate"
      } else {
        PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] <- "START"
      }
    }
    ###label the end of each chromosome
    for (i in seq_along(PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[1:5245])) {
      if (identical(PbANKA_33_gene_upto1KB_DOWN_calc$SEQID[i], PbANKA_33_gene_upto1KB_DOWN_calc$SEQID[i+1]) && (PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] == "calculate")) {
        PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] <- "calculate"
      } else {
        if (PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] != "START") {
          PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] <- "END"
        }
      }
    }
    
    ###at start of CHROM if gene is on -ve strand, calculate downstream dist
    #if gene is on +ve strand will go through normal calculation
    for (i in seq_along(PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST)) {
      if ((PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] == "START") && (PbANKA_33_gene_upto1KB_DOWN_calc$STRAND[i] == "-")) {
        PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] <- PbANKA_33_gene_upto1KB_DOWN_calc$START[i] - 1
      } 
    }
    
    ###at end of CHROM if gene is on +ve strand, need to calculate downstream dist based on size of CHROM
    #if gene is on -ve strand will go through normal calc
    for (i in seq_along(PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST)) {
      if ((PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] == "END") && (PbANKA_33_gene_upto1KB_DOWN_calc$STRAND[i] == "+")) {
        PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] <- (PbANKA_33_chrom_coord$END[PbANKA_33_chrom_coord$SEQID == PbANKA_33_gene_upto1KB_DOWN_calc$SEQID[i]]) - PbANKA_33_gene_upto1KB_DOWN_calc$END[i]
      } 
    }
    
    #then move on to populating the rest of the values, ignoring the rows that have already been filled in
    #change all DOWN_DIST values that aren't numbers to "calculate"
    for (i in seq_along(PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST)) {
      if ((PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] == "calculate") ||
          (PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] == "START") ||
          (PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] == "END") ) {
        PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] <- "calculate"
      }
    }
    
    #change all DOWN_DIST values for +ve strand genes to "plus"
    for (i in seq_along(PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST)) {
      if ((PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] == "calculate") &&
          (PbANKA_33_gene_upto1KB_DOWN_calc$STRAND[i] == "+")) {
        PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] <- "plus"
      }
    }
    
    #change all DOWN_DIST values for -ve strand genes to "minus"
    for (i in seq_along(PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST)) {
      if ((PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] == "calculate") &&
          (PbANKA_33_gene_upto1KB_DOWN_calc$STRAND[i] == "-")) {
        PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] <- "minus"
      }
    }
    
    #calculate upstream dist for +ve strand genes
    #if STRAND = +ve, DOWN_DIST = START of next line - END 
    for (i in seq_along(PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST)) {
      if (PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] == "plus") {
        PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] <- PbANKA_33_gene_upto1KB_DOWN_calc$START[i+1] - PbANKA_33_gene_upto1KB_DOWN_calc$END[i]
      }
    }
    
    #calculate upstream dist for -ve strand genes
    #if STRAND = -ve, DOWN_DIST = START - END of line previous
    for (i in seq_along(PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST)) {
      if (PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] == "minus") {
        PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST[i] <- PbANKA_33_gene_upto1KB_DOWN_calc$START[i] - PbANKA_33_gene_upto1KB_DOWN_calc$END[i-1]
      }
    }
    ##convert values to integers
    PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST <- as.integer(PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST)
    
####organize up to 1KB downstream cutoff values#####
    #generate edited DOWN_DIST to cut off at 1KB
    PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST_EDIT <- ifelse(PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST > 1000, 1000, PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST)
    
    #remove genes with a DOWN_DIST <1 (no bases between genes)
    pos_PbANKA_33_gene_upto1KB_DOWN_strand <- subset(PbANKA_33_gene_upto1KB_DOWN_calc, PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST_EDIT >= 3)
      write.csv(pos_PbANKA_33_gene_upto1KB_DOWN_strand, "pos_PbANKA_33_gene_upto1KB_DOWN_strand.csv")
      
    neg_PbANKA_33_gene_strand  <- subset(PbANKA_33_gene_upto1KB_DOWN_calc, PbANKA_33_gene_upto1KB_DOWN_calc$DOWN_DIST_EDIT <= 2)
      write.csv(neg_PbANKA_33_gene_strand, "neg_PbANKA_33_gene_upto1KB_DOWN_strand.csv")
  
#####up to 1KB Upstream region coordinate organization###########      
    #the next step is to figure out the beginning and end of the upstream region for both the +ve strand genes (5' of the gene start) and
    #the -ve strand genes (based on coordinate system, 3' of the gene end)
      ##remember i have to switch to a 0-based coord system
      ##if GENE_START is 5 and UP_DIST_EDIT is 4, a 1-based system would use 1:4 as the UP_START and UP_END, a 0-based system would use 0:4
    
    # +ve gene --> UPSTART = GENE_START - (UP_DIST_EDIT + 1)
    #UPEND = GENE_START - 1
      
    # -ve gene --> UPSTART = GENE_END 
    #UPEND = GENE_END + UP_DIST_EDIT
    
    PbANKA_33_upto1KB_upstream_regions_strand <- data.frame(matrix(NA, nrow = 5198))
      PbANKA_33_upto1KB_upstream_regions_strand$CHROM <-  pos_PbANKA_33_gene_upto1KB_UP_strand$SEQID
      PbANKA_33_upto1KB_upstream_regions_strand$GENEID <- pos_PbANKA_33_gene_upto1KB_UP_strand$GENEID
      PbANKA_33_upto1KB_upstream_regions_strand <- PbANKA_33_upto1KB_upstream_regions_strand[,-1]
      PbANKA_33_upto1KB_upstream_regions_strand$STRAND <- pos_PbANKA_33_gene_upto1KB_UP_strand$STRAND
      PbANKA_33_upto1KB_upstream_regions_strand$GENE_START <- pos_PbANKA_33_gene_upto1KB_UP_strand$START
      PbANKA_33_upto1KB_upstream_regions_strand$GENE_END <- pos_PbANKA_33_gene_upto1KB_UP_strand$END
      PbANKA_33_upto1KB_upstream_regions_strand$UP_DIST_EDIT <- pos_PbANKA_33_gene_upto1KB_UP_strand$UP_DIST_EDIT
      PbANKA_33_upto1KB_upstream_regions_strand$UPSTART <- NA
      PbANKA_33_upto1KB_upstream_regions_strand$UPEND <- NA
    
    
    for (i in seq_along(PbANKA_33_upto1KB_upstream_regions_strand$GENEID)) {
      if (PbANKA_33_upto1KB_upstream_regions_strand$STRAND[i] == "+") {
        PbANKA_33_upto1KB_upstream_regions_strand$UPSTART[i] <- PbANKA_33_upto1KB_upstream_regions_strand$GENE_START[i] - (PbANKA_33_upto1KB_upstream_regions_strand$UP_DIST_EDIT[i] + 1)
        PbANKA_33_upto1KB_upstream_regions_strand$UPEND[i] <- PbANKA_33_upto1KB_upstream_regions_strand$GENE_START[i] - 1
      }
    }
    
    for (i in seq_along(PbANKA_33_upto1KB_upstream_regions_strand$GENEID)) {
      if (PbANKA_33_upto1KB_upstream_regions_strand$STRAND[i] == "-") {
        PbANKA_33_upto1KB_upstream_regions_strand$UPSTART[i] <- PbANKA_33_upto1KB_upstream_regions_strand$GENE_END[i]
        PbANKA_33_upto1KB_upstream_regions_strand$UPEND[i] <- PbANKA_33_upto1KB_upstream_regions_strand$GENE_END[i] + PbANKA_33_upto1KB_upstream_regions_strand$UP_DIST_EDIT[i]
      }
    }
    
    write.csv(PbANKA_33_upto1KB_upstream_regions_strand$GENEID, "PbANKA_33_upto1KB_upstream_regions_strand_genelist.csv")
    
#####up to 1KB downstream region coordinate organization###########      
    #the next step is to figure out the beginning and end of the downstream region for both the +ve strand genes (3' of the gene end) and
    #the -ve strand genes (based on coordinate system, 5' of the gene start)
    ##remember i have to switch to a 0-based coord system
    ##if GENE_END is 2 and DOWN_DIST_EDIT is 4, a 1-based system would use 3:6 as the DOWN_START and DOWN_END, a 0-based system would use 2:6
    
    
    # +ve gene --> DOWNSTART = GENE_END
    #DOWNEND = GENE_END + DOWN_DIST_EDIT
    
    # -ve gene --> DOWNSTART = GENE_START - (DOWN_DIST_EDIT + 1)
    #DOWNEND = GENE_START - 1

    
    PbANKA_33_upto1KB_downstream_regions_strand <- data.frame(matrix(NA, nrow = 5181))
      PbANKA_33_upto1KB_downstream_regions_strand$CHROM <-  pos_PbANKA_33_gene_upto1KB_DOWN_strand$SEQID
      PbANKA_33_upto1KB_downstream_regions_strand$GENEID <- pos_PbANKA_33_gene_upto1KB_DOWN_strand$GENEID
      PbANKA_33_upto1KB_downstream_regions_strand <- PbANKA_33_upto1KB_downstream_regions_strand[,-1]
      PbANKA_33_upto1KB_downstream_regions_strand$STRAND <- pos_PbANKA_33_gene_upto1KB_DOWN_strand$STRAND
      PbANKA_33_upto1KB_downstream_regions_strand$GENE_START <- pos_PbANKA_33_gene_upto1KB_DOWN_strand$START
      PbANKA_33_upto1KB_downstream_regions_strand$GENE_END <- pos_PbANKA_33_gene_upto1KB_DOWN_strand$END
      PbANKA_33_upto1KB_downstream_regions_strand$DOWN_DIST_EDIT <- pos_PbANKA_33_gene_upto1KB_DOWN_strand$DOWN_DIST_EDIT
      PbANKA_33_upto1KB_downstream_regions_strand$DOWNSTART <- NA
      PbANKA_33_upto1KB_downstream_regions_strand$DOWNEND <- NA
      
    
    for (i in seq_along(PbANKA_33_upto1KB_downstream_regions_strand$GENEID)) {
      if (PbANKA_33_upto1KB_downstream_regions_strand$STRAND[i] == "+") {
        PbANKA_33_upto1KB_downstream_regions_strand$DOWNSTART[i] <- PbANKA_33_upto1KB_downstream_regions_strand$GENE_END[i]
        PbANKA_33_upto1KB_downstream_regions_strand$DOWNEND[i] <- PbANKA_33_upto1KB_downstream_regions_strand$GENE_END[i] + PbANKA_33_upto1KB_downstream_regions_strand$DOWN_DIST_EDIT[i]
      }
    }
    
    for (i in seq_along(PbANKA_33_upto1KB_downstream_regions_strand$GENEID)) {
      if (PbANKA_33_upto1KB_downstream_regions_strand$STRAND[i] == "-") {
        PbANKA_33_upto1KB_downstream_regions_strand$DOWNSTART[i] <- PbANKA_33_upto1KB_downstream_regions_strand$GENE_START[i] - (PbANKA_33_upto1KB_downstream_regions_strand$DOWN_DIST_EDIT[i] + 1)
        PbANKA_33_upto1KB_downstream_regions_strand$DOWNEND[i] <- PbANKA_33_upto1KB_downstream_regions_strand$GENE_START[i] - 1
      }
    }
    
    write.csv(PbANKA_33_upto1KB_downstream_regions_strand$GENEID, "PbANKA_33_upto1KB_downstream_regions_strand_genelist.csv") 
    
    
###up to 1KB up and downstream - create and download text for bedtools BED files#####
    # BED columns: CHROM START END NAME
  
    ###upstream###
  PbANKA_33_upto1KB_upstream_regions_strand_download_BED <- data.frame(matrix(nrow = 5198))
    PbANKA_33_upto1KB_upstream_regions_strand_download_BED$CHROM <- PbANKA_33_upto1KB_upstream_regions_strand$CHROM
    PbANKA_33_upto1KB_upstream_regions_strand_download_BED$START <- PbANKA_33_upto1KB_upstream_regions_strand$UPSTART
    PbANKA_33_upto1KB_upstream_regions_strand_download_BED$END <- PbANKA_33_upto1KB_upstream_regions_strand$UPEND
    PbANKA_33_upto1KB_upstream_regions_strand_download_BED$NAME <- PbANKA_33_upto1KB_upstream_regions_strand$GENEID
    PbANKA_33_upto1KB_upstream_regions_strand_download_BED <- PbANKA_33_upto1KB_upstream_regions_strand_download_BED[,2:5]
    write.csv(PbANKA_33_upto1KB_upstream_regions_strand_download_BED, "PbANKA_33_upto1KB_upstream_regions_strand_download_BED.csv")   
    
    ###downstream###
    PbANKA_33_upto1KB_downstream_regions_strand_download_BED <- data.frame(matrix(nrow = 5181))
      PbANKA_33_upto1KB_downstream_regions_strand_download_BED$CHROM <- PbANKA_33_upto1KB_downstream_regions_strand$CHROM
      PbANKA_33_upto1KB_downstream_regions_strand_download_BED$START <- PbANKA_33_upto1KB_downstream_regions_strand$DOWNSTART
      PbANKA_33_upto1KB_downstream_regions_strand_download_BED$END <- PbANKA_33_upto1KB_downstream_regions_strand$DOWNEND
      PbANKA_33_upto1KB_downstream_regions_strand_download_BED$NAME <- PbANKA_33_upto1KB_downstream_regions_strand$GENEID
      PbANKA_33_upto1KB_downstream_regions_strand_download_BED <- PbANKA_33_upto1KB_downstream_regions_strand_download_BED[,2:5]
      write.csv(PbANKA_33_upto1KB_downstream_regions_strand_download_BED, "PbANKA_33_upto1KB_downstream_regions_strand_download_BED.csv")   
    
    
      
#####gene ORF calculations and BED file#####
      
    ##will need to subtract 1 from gene start for both +ve and -ve strand genes (since the coordinates are all based on the -ve strand)
      
      PbANKA_33_gene_ORF_calc <- PbANKA_33_gene_coord
      
      #generate edited GENE_LENGTH to cut off at 1KB
            PbANKA_33_gene_ORF_calc$GENE_LENGTH <- NA
      for (i in seq_along(PbANKA_33_gene_ORF_calc$GENE_LENGTH[1:5245])) {
        PbANKA_33_gene_ORF_calc$GENE_LENGTH[i] <- PbANKA_33_gene_upto1KB_UP_calc$END[i] - PbANKA_33_gene_upto1KB_UP_calc$START[i]
      }
      
    
      #remove genes with an GENE_LENGTH <1 
      pos_PbANKA_33_gene_ORF_strand <- subset(PbANKA_33_gene_ORF_calc, PbANKA_33_gene_ORF_calc$GENE_LENGTH >= 3)
        write.csv(pos_PbANKA_33_gene_ORF_strand, "pos_PbANKA_33_gene_ORF_strand.csv")
      
      neg_PbANKA_33_gene_ORF_strand  <- subset(PbANKA_33_gene_ORF_calc, PbANKA_33_gene_ORF_calc$GENE_LENGTH <= 2)
        write.csv(neg_PbANKA_33_gene_ORF_strand, "neg_PbANKA_33_gene_ORF_strand.csv")
      
      
      
      #the next step is to figure out the beginning and end of the ORF region for all genes
        ##will need to subtract 1 from gene start for both +ve and -ve strand genes (since the coordinates are all based on the -ve strand)
       
        #ORF_START = START - 1
        #ORF_END = END
      
      PbANKA_33_ORF_regions_strand <- data.frame(matrix(NA, nrow = 5245))
        PbANKA_33_ORF_regions_strand$CHROM <-  pos_PbANKA_33_gene_ORF_strand$SEQID
        PbANKA_33_ORF_regions_strand$GENEID <- pos_PbANKA_33_gene_ORF_strand$GENEID
        PbANKA_33_ORF_regions_strand <- PbANKA_33_ORF_regions_strand[,-1]
        PbANKA_33_ORF_regions_strand$STRAND <- pos_PbANKA_33_gene_ORF_strand$STRAND
        PbANKA_33_ORF_regions_strand$GENE_START <- pos_PbANKA_33_gene_ORF_strand$START
        PbANKA_33_ORF_regions_strand$GENE_END <- pos_PbANKA_33_gene_ORF_strand$END
        PbANKA_33_ORF_regions_strand$GENE_LENGTH <- pos_PbANKA_33_gene_ORF_strand$GENE_LENGTH
        PbANKA_33_ORF_regions_strand$ORFSTART <- NA
        PbANKA_33_ORF_regions_strand$ORFEND <- NA
        
      
      for (i in seq_along(PbANKA_33_ORF_regions_strand$GENEID)) {
        PbANKA_33_ORF_regions_strand$ORFSTART[i] <- PbANKA_33_ORF_regions_strand$GENE_START[i] - 1
        PbANKA_33_ORF_regions_strand$ORFEND[i] <- PbANKA_33_ORF_regions_strand$GENE_END[i]
      }

      write.csv(PbANKA_33_ORF_regions_strand$GENEID, "PbANKA_33_ORF_regions_strand_genelist.csv")
      
      
####make ORF BED file#####
      PbANKA_33_ORF_regions_strand_download_BED <- data.frame(matrix(nrow = 5245))
      PbANKA_33_ORF_regions_strand_download_BED$CHROM <- PbANKA_33_ORF_regions_strand$CHROM
      PbANKA_33_ORF_regions_strand_download_BED$START <- PbANKA_33_ORF_regions_strand$ORFSTART
      PbANKA_33_ORF_regions_strand_download_BED$END <- PbANKA_33_ORF_regions_strand$ORFEND
      PbANKA_33_ORF_regions_strand_download_BED$NAME <- PbANKA_33_ORF_regions_strand$GENEID
      PbANKA_33_ORF_regions_strand_download_BED <- PbANKA_33_ORF_regions_strand_download_BED[,2:5]
      write.csv(PbANKA_33_ORF_regions_strand_download_BED, "PbANKA_33_ORF_regions_strand_download_BED.csv")   

#####calculate strand specific 1kb absolute upstream dist####
  ###start with the beginning and ending of CHROM
      #beginning of CHROM if the first gene is on the +ve strand, end of CHROm if the last gene is on the -ve strand
      
      PbANKA_33_gene_1KB_UP_calc <- PbANKA_33_gene_coord
      PbANKA_33_gene_1KB_UP_calc$UP_DIST <- NA
      
      ###UP_DIST###
      
      #Calculate all upstream distances from the start of each CHROM (if on +ve strand) or end of each CHROM (if on -ve strand)
      #change any UP_DIST > 1000 to 1000
      
      
      ###if gene is on +ve strand, calculate upstream dist from the start of the CHROM - start is always at 1 so just use the GENE_START (START)
      for (i in seq_along(PbANKA_33_gene_1KB_UP_calc$UP_DIST)) {
        if (PbANKA_33_gene_1KB_UP_calc$STRAND[i] == "+") {
          PbANKA_33_gene_1KB_UP_calc$UP_DIST[i] <- PbANKA_33_gene_1KB_UP_calc$START[i] - 1
        } 
      }
      
      ###if gene is on -ve strand, need to calculate upstream dist based on size of CHROM
      for (i in seq_along(PbANKA_33_gene_1KB_UP_calc$UP_DIST)) {
        if (PbANKA_33_gene_1KB_UP_calc$STRAND[i] == "-") {
          PbANKA_33_gene_1KB_UP_calc$UP_DIST[i] <- (PbANKA_33_chrom_coord$END[PbANKA_33_chrom_coord$SEQID == PbANKA_33_gene_1KB_UP_calc$SEQID[i]]) - PbANKA_33_gene_1KB_UP_calc$END[i]
        } 
      }
      
      ##convert values to integers
      PbANKA_33_gene_1KB_UP_calc$UP_DIST <- as.integer(PbANKA_33_gene_1KB_UP_calc$UP_DIST)  
      
####organize 1KB upstream cutoff values#####
      #change any UP_DIST > 1000 to 1000
      PbANKA_33_gene_1KB_UP_calc$UP_DIST_EDIT <- ifelse(PbANKA_33_gene_1KB_UP_calc$UP_DIST > 1000, 1000, PbANKA_33_gene_1KB_UP_calc$UP_DIST)
      
      #remove genes with an UP_DIST <1 (no bases between genes)
      pos_PbANKA_33_gene_1KB_UP_strand <- subset(PbANKA_33_gene_1KB_UP_calc, PbANKA_33_gene_1KB_UP_calc$UP_DIST_EDIT >= 3)
        write.csv(pos_PbANKA_33_gene_1KB_UP_strand, "pos_PbANKA_33_gene_1KB_UP_strand.csv")
      
      neg_PbANKA_33_gene_1KB_UP_strand  <- subset(PbANKA_33_gene_1KB_UP_calc, PbANKA_33_gene_1KB_UP_calc$UP_DIST_EDIT <= 2)
        write.csv(neg_PbANKA_33_gene_1KB_UP_strand, "neg_PbANKA_33_gene_1KB_UP_strand.csv")

#####calculate strand specific 1kb absolute downstream dist####   
        ###start with the beginning and ending of CHROM
        #end of CHROM if the last gene is on the +ve strand, beginning of CHROM if the first gene is on the +ve strand
        
        PbANKA_33_gene_1KB_DOWN_calc <- PbANKA_33_gene_coord
        PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST <- NA
        
  ###DOWN_DIST###

        #Calculate all downstream distances from the start of each CHROM (if on -ve strand) or end of each CHROM (if on +ve strand)
        #change any DOWN_DIST > 1000 to 1000
        
        ###if gene is on -ve strand, calculate downstream dist from the start of the CHROM - start is always at 1 so just use the GENE_START (START)
        for (i in seq_along(PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST)) {
          if (PbANKA_33_gene_1KB_DOWN_calc$STRAND[i] == "-") {
            PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST[i] <- PbANKA_33_gene_1KB_DOWN_calc$START[i] - 1
          } 
        }
        
        ###if gene is on +ve strand, need to calculate downstream dist based on size of CHROM
        for (i in seq_along(PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST)) {
          if (PbANKA_33_gene_1KB_DOWN_calc$STRAND[i] == "+") {
            PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST[i] <- (PbANKA_33_chrom_coord$END[PbANKA_33_chrom_coord$SEQID == PbANKA_33_gene_1KB_DOWN_calc$SEQID[i]]) - PbANKA_33_gene_1KB_DOWN_calc$END[i]
          } 
        }
        
        ##convert values to integers
        PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST <- as.integer(PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST)  
      
####organize 1KB downstream cutoff values#####
        #change any DOWN_DIST > 1000 to 1000
        PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST_EDIT <- ifelse(PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST > 1000, 1000, PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST)
        
        #remove genes with an DOWN_DIST <1 (no bases between genes)
        pos_PbANKA_33_gene_1KB_DOWN_strand <- subset(PbANKA_33_gene_1KB_DOWN_calc, PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST_EDIT >= 3)
          write.csv(pos_PbANKA_33_gene_1KB_DOWN_strand, "pos_PbANKA_33_gene_1KB_DOWN_strand.csv")
        
        neg_PbANKA_33_gene_1KB_DOWN_strand  <- subset(PbANKA_33_gene_1KB_DOWN_calc, PbANKA_33_gene_1KB_DOWN_calc$DOWN_DIST_EDIT <= 2)
          write.csv(neg_PbANKA_33_gene_1KB_DOWN_strand, "neg_PbANKA_33_gene_1KB_DOWN_strand.csv")      
  
####1KB Upstream region coordinate organization###########      
          #the next step is to figure out the beginning and end of the upstream region for both the +ve strand genes (5' of the gene start) and
          #the -ve strand genes (based on coordinate system, 3' of the gene end)
          ##remember i have to switch to a 0-based coord system
          ##if GENE_START is 5 and UP_DIST_EDIT is 4, a 1-based system would use 1:4 as the UP_START and UP_END, a 0-based system would use 0:4
          
          # +ve gene --> UPSTART = GENE_START - (UP_DIST_EDIT + 1)
          #UPEND = GENE_START - 1
          
          # -ve gene --> UPSTART = GENE_END 
          #UPEND = GENE_END + UP_DIST_EDIT
          
          PbANKA_33_1KB_upstream_regions_strand <- data.frame(matrix(NA, nrow = 5243))
            PbANKA_33_1KB_upstream_regions_strand$CHROM <-  pos_PbANKA_33_gene_1KB_UP_strand$SEQID
            PbANKA_33_1KB_upstream_regions_strand$GENEID <- pos_PbANKA_33_gene_1KB_UP_strand$GENEID
            PbANKA_33_1KB_upstream_regions_strand <- PbANKA_33_1KB_upstream_regions_strand[,-1]
            PbANKA_33_1KB_upstream_regions_strand$STRAND <- pos_PbANKA_33_gene_1KB_UP_strand$STRAND
            PbANKA_33_1KB_upstream_regions_strand$GENE_START <- pos_PbANKA_33_gene_1KB_UP_strand$START
            PbANKA_33_1KB_upstream_regions_strand$GENE_END <- pos_PbANKA_33_gene_1KB_UP_strand$END
            PbANKA_33_1KB_upstream_regions_strand$UP_DIST_EDIT <- pos_PbANKA_33_gene_1KB_UP_strand$UP_DIST_EDIT
            PbANKA_33_1KB_upstream_regions_strand$UPSTART <- NA
            PbANKA_33_1KB_upstream_regions_strand$UPEND <- NA
          
          
          for (i in seq_along(PbANKA_33_1KB_upstream_regions_strand$GENEID)) {
            if (PbANKA_33_1KB_upstream_regions_strand$STRAND[i] == "+") {
              PbANKA_33_1KB_upstream_regions_strand$UPSTART[i] <- PbANKA_33_1KB_upstream_regions_strand$GENE_START[i] - (PbANKA_33_1KB_upstream_regions_strand$UP_DIST_EDIT[i] + 1)
              PbANKA_33_1KB_upstream_regions_strand$UPEND[i] <- PbANKA_33_1KB_upstream_regions_strand$GENE_START[i] - 1
            }
          }
          
          for (i in seq_along(PbANKA_33_1KB_upstream_regions_strand$GENEID)) {
            if (PbANKA_33_1KB_upstream_regions_strand$STRAND[i] == "-") {
              PbANKA_33_1KB_upstream_regions_strand$UPSTART[i] <- PbANKA_33_1KB_upstream_regions_strand$GENE_END[i]
              PbANKA_33_1KB_upstream_regions_strand$UPEND[i] <- PbANKA_33_1KB_upstream_regions_strand$GENE_END[i] + PbANKA_33_1KB_upstream_regions_strand$UP_DIST_EDIT[i]
            }
          }
          
          write.csv(PbANKA_33_1KB_upstream_regions_strand$GENEID, "PbANKA_33_1KB_upstream_regions_strand_genelist.csv")
          
####1KB downstream region coordinate organization###########      
          #the next step is to figure out the beginning and end of the downstream region for both the +ve strand genes (3' of the gene end) and
          #the -ve strand genes (based on coordinate system, 5' of the gene start)
          ##remember i have to switch to a 0-based coord system
          ##if GENE_END is 2 and DOWN_DIST_EDIT is 4, a 1-based system would use 3:6 as the DOWN_START and DOWN_END, a 0-based system would use 2:6
          
          
          # +ve gene --> DOWNSTART = GENE_END
          #DOWNEND = GENE_END + DOWN_DIST_EDIT
          
          # -ve gene --> DOWNSTART = GENE_START - (DOWN_DIST_EDIT + 1)
          #DOWNEND = GENE_START - 1
          
          
          PbANKA_33_1KB_downstream_regions_strand <- data.frame(matrix(NA, nrow = 5244))
            PbANKA_33_1KB_downstream_regions_strand$CHROM <-  pos_PbANKA_33_gene_1KB_DOWN_strand$SEQID
            PbANKA_33_1KB_downstream_regions_strand$GENEID <- pos_PbANKA_33_gene_1KB_DOWN_strand$GENEID
            PbANKA_33_1KB_downstream_regions_strand <- PbANKA_33_1KB_downstream_regions_strand[,-1]
            PbANKA_33_1KB_downstream_regions_strand$STRAND <- pos_PbANKA_33_gene_1KB_DOWN_strand$STRAND
            PbANKA_33_1KB_downstream_regions_strand$GENE_START <- pos_PbANKA_33_gene_1KB_DOWN_strand$START
            PbANKA_33_1KB_downstream_regions_strand$GENE_END <- pos_PbANKA_33_gene_1KB_DOWN_strand$END
            PbANKA_33_1KB_downstream_regions_strand$DOWN_DIST_EDIT <- pos_PbANKA_33_gene_1KB_DOWN_strand$DOWN_DIST_EDIT
            PbANKA_33_1KB_downstream_regions_strand$DOWNSTART <- NA
            PbANKA_33_1KB_downstream_regions_strand$DOWNEND <- NA
          
          
          for (i in seq_along(PbANKA_33_1KB_downstream_regions_strand$GENEID)) {
            if (PbANKA_33_1KB_downstream_regions_strand$STRAND[i] == "+") {
              PbANKA_33_1KB_downstream_regions_strand$DOWNSTART[i] <- PbANKA_33_1KB_downstream_regions_strand$GENE_END[i]
              PbANKA_33_1KB_downstream_regions_strand$DOWNEND[i] <- PbANKA_33_1KB_downstream_regions_strand$GENE_END[i] + PbANKA_33_1KB_downstream_regions_strand$DOWN_DIST_EDIT[i]
            }
          }
          
          for (i in seq_along(PbANKA_33_1KB_downstream_regions_strand$GENEID)) {
            if (PbANKA_33_1KB_downstream_regions_strand$STRAND[i] == "-") {
              PbANKA_33_1KB_downstream_regions_strand$DOWNSTART[i] <- PbANKA_33_1KB_downstream_regions_strand$GENE_START[i] - (PbANKA_33_1KB_downstream_regions_strand$DOWN_DIST_EDIT[i] + 1)
              PbANKA_33_1KB_downstream_regions_strand$DOWNEND[i] <- PbANKA_33_1KB_downstream_regions_strand$GENE_START[i] - 1
            }
          }
          
          write.csv(PbANKA_33_1KB_downstream_regions_strand$GENEID, "PbANKA_33_1KB_downstream_regions_strand_genelist.csv") 
          
          
          
          
          
          
          
### 1KB up and downstream - create and download text for bedtools BED files#####
          # BED columns: CHROM START END NAME
          
  ###upstream###
  PbANKA_33_1KB_upstream_regions_strand_download_BED <- data.frame(matrix(nrow = 5243))
    PbANKA_33_1KB_upstream_regions_strand_download_BED$CHROM <- PbANKA_33_1KB_upstream_regions_strand$CHROM
    PbANKA_33_1KB_upstream_regions_strand_download_BED$START <- PbANKA_33_1KB_upstream_regions_strand$UPSTART
    PbANKA_33_1KB_upstream_regions_strand_download_BED$END <- PbANKA_33_1KB_upstream_regions_strand$UPEND
    PbANKA_33_1KB_upstream_regions_strand_download_BED$NAME <- PbANKA_33_1KB_upstream_regions_strand$GENEID
    PbANKA_33_1KB_upstream_regions_strand_download_BED <- PbANKA_33_1KB_upstream_regions_strand_download_BED[,2:5]
    write.csv(PbANKA_33_1KB_upstream_regions_strand_download_BED, "PbANKA_33_1KB_upstream_regions_strand_download_BED.csv")   
  
  ###downstream###
  PbANKA_33_1KB_downstream_regions_strand_download_BED <- data.frame(matrix(nrow = 5244))
    PbANKA_33_1KB_downstream_regions_strand_download_BED$CHROM <- PbANKA_33_1KB_downstream_regions_strand$CHROM
    PbANKA_33_1KB_downstream_regions_strand_download_BED$START <- PbANKA_33_1KB_downstream_regions_strand$DOWNSTART
    PbANKA_33_1KB_downstream_regions_strand_download_BED$END <- PbANKA_33_1KB_downstream_regions_strand$DOWNEND
    PbANKA_33_1KB_downstream_regions_strand_download_BED$NAME <- PbANKA_33_1KB_downstream_regions_strand$GENEID
    PbANKA_33_1KB_downstream_regions_strand_download_BED <- PbANKA_33_1KB_downstream_regions_strand_download_BED[,2:5]
    write.csv(PbANKA_33_1KB_downstream_regions_strand_download_BED, "PbANKA_33_1KB_downstream_regions_strand_download_BED.csv")   
  
    
    
#####GoI BED file#####
    
PbANKA_33_Alba3GoI_upto1KB_upstream_regions_strand_download_BED <- subset(PbANKA_33_upto1KB_upstream_regions_strand_download_BED, PbANKA_33_upto1KB_upstream_regions_strand_download_BED$NAME %in% RIPseq_GoI)
  write.csv(PbANKA_33_Alba3GoI_upto1KB_upstream_regions_strand_download_BED, "J:/III/Waters/Group Members/Mallory/NishaRNAseq/BEDfiles/PbANKA_33_Alba3GoI_upto1KB_upstream_regions_strand_download_BED.csv")

PbANKA_33_Alba3GoI_upto1KB_downstream_regions_strand_download_BED <- subset(PbANKA_33_upto1KB_downstream_regions_strand_download_BED, PbANKA_33_upto1KB_downstream_regions_strand_download_BED$NAME %in% RIPseq_GoI)
  write.csv(PbANKA_33_Alba3GoI_upto1KB_downstream_regions_strand_download_BED, "J:/III/Waters/Group Members/Mallory/NishaRNAseq/BEDfiles/PbANKA_33_Alba3GoI_upto1KB_downstream_regions_strand_download_BED.csv")

PbANKA_33_Alba3GoI_ORF_regions_strand_download_BED <- subset(PbANKA_33_ORF_regions_strand_download_BED, PbANKA_33_ORF_regions_strand_download_BED$NAME %in% RIPseq_GoI)
  write.csv(PbANKA_33_Alba3GoI_ORF_regions_strand_download_BED, "J:/III/Waters/Group Members/Mallory/NishaRNAseq/BEDfiles/PbANKA_33_Alba3GoI_ORF_regions_strand_download_BED.csv")

PbANKA_33_Alba3GoI_1KB_upstream_regions_strand_download_BED <- subset(PbANKA_33_1KB_upstream_regions_strand_download_BED, PbANKA_33_1KB_upstream_regions_strand_download_BED$NAME %in% RIPseq_GoI)
  write.csv(PbANKA_33_Alba3GoI_1KB_upstream_regions_strand_download_BED, "J:/III/Waters/Group Members/Mallory/NishaRNAseq/BEDfiles/PbANKA_33_Alba3GoI_1KB_upstream_regions_strand_download_BED.csv")

PbANKA_33_Alba3GoI_1KB_downstream_regions_strand_download_BED <- subset(PbANKA_33_1KB_downstream_regions_strand_download_BED, PbANKA_33_1KB_downstream_regions_strand_download_BED$NAME %in% RIPseq_GoI)
  write.csv(PbANKA_33_Alba3GoI_1KB_downstream_regions_strand_download_BED, "J:/III/Waters/Group Members/Mallory/NishaRNAseq/BEDfiles/PbANKA_33_Alba3GoI_1KB_downstream_regions_strand_download_BED.csv")

    
#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/NishaRIPseq_up_down_ORF_region_extract.RData")
