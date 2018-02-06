#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/Pberghei_CDS_region_RPFM_calculation.RData")
require("plyr")

#####set WD#####
setwd("J:/III/Waters/Group Members/Mallory/PbergheiUpstreamRegionsSebastian/")


#####load data#####
setwd("J:/III/Waters/Group Members/Mallory/PbergheiUpstreamRegionsSebastian")
#Pberghei gff3 file
Pberghei_GFF3_CDS_strand_calc <- read.csv("GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv", header=T)
  colnames(Pberghei_GFF3_CDS_strand_calc)[3] <- "CDS_START"
  colnames(Pberghei_GFF3_CDS_strand_calc)[4] <- "CDS_END"

setwd("J:/III/Waters/Group Members/Mallory/PbergheiUpstreamRegionsSebastian/BEDfiles/")
#Pbeghei bed file
Pbeghei_CDS_regions_strand <- read.delim("Pberghei_CDS_regions_strand_download_BED.bed", header=F)
  colnames(Pbeghei_CDS_regions_strand) <- c("CHROM","BED_CDS_START","BED_CDS_END","NAME")

setwd("J:/III/Waters/Group Members/Mallory/PbergheiUpstreamRegionsSebastian/Sir2A-KO_Mapping/CDS")

#the files!
# ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes.bed
# ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes.bed
# ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes.bed
# ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes.bed
# ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes.bed
# ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes.bed
# ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes.bed
# ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes.bed
# ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes.bed
# ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes.bed

ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes <- read.delim("ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes.bed", header=F)
  colnames(ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes)[13:16] <- c("CHROM","BED_CDS_START","BED_CDS_END","NAME")
  colnames(ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes)[1] <- "READ_CHROM"
  colnames(ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes)[2] <- "OVERLAP_START"
  colnames(ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes)[3] <- "OVERLAP_END"
  colnames(ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes)[4] <- "QNAME"
  colnames(ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes)[6] <- "STRAND"
  colnames(ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes)[7] <- "READ_START"
  colnames(ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes)[8] <- "READ_END"

ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes <- read.delim("ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes.bed", header=F)
  colnames(ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes)[13:16] <- c("CHROM","BED_CDS_START","BED_CDS_END","NAME")
  colnames(ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes)[1] <- "READ_CHROM"
  colnames(ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes)[2] <- "OVERLAP_START"
  colnames(ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes)[3] <- "OVERLAP_END"
  colnames(ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes)[4] <- "QNAME"
  colnames(ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes)[6] <- "STRAND"
  colnames(ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes)[7] <- "READ_START"
  colnames(ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes)[8] <- "READ_END"  
  
ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes <- read.delim("ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes.bed", header=F)
  colnames(ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes)[13:16] <- c("CHROM","BED_CDS_START","BED_CDS_END","NAME")
  colnames(ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes)[1] <- "READ_CHROM"
  colnames(ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes)[2] <- "OVERLAP_START"
  colnames(ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes)[3] <- "OVERLAP_END"
  colnames(ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes)[4] <- "QNAME"
  colnames(ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes)[6] <- "STRAND"
  colnames(ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes)[7] <- "READ_START"
  colnames(ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes)[8] <- "READ_END"  

ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes <- read.delim("ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes.bed", header=F)
  colnames(ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes)[13:16] <- c("CHROM","BED_CDS_START","BED_CDS_END","NAME")
  colnames(ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes)[1] <- "READ_CHROM"
  colnames(ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes)[2] <- "OVERLAP_START"
  colnames(ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes)[3] <- "OVERLAP_END"
  colnames(ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes)[4] <- "QNAME"
  colnames(ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes)[6] <- "STRAND"
  colnames(ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes)[7] <- "READ_START"
  colnames(ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes)[8] <- "READ_END"  

ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes <- read.delim("ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes.bed", header=F)
  colnames(ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[13:16] <- c("CHROM","BED_CDS_START","BED_CDS_END","NAME")
  colnames(ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[1] <- "READ_CHROM"
  colnames(ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[2] <- "OVERLAP_START"
  colnames(ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[3] <- "OVERLAP_END"
  colnames(ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[4] <- "QNAME"
  colnames(ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[6] <- "STRAND"
  colnames(ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[7] <- "READ_START"
  colnames(ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[8] <- "READ_END"    
  

ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes <- read.delim("ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes.bed", header=F)
  colnames(ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes)[13:16] <- c("CHROM","BED_CDS_START","BED_CDS_END","NAME")
  colnames(ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes)[1] <- "READ_CHROM"
  colnames(ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes)[2] <- "OVERLAP_START"
  colnames(ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes)[3] <- "OVERLAP_END"
  colnames(ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes)[4] <- "QNAME"
  colnames(ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes)[6] <- "STRAND"
  colnames(ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes)[7] <- "READ_START"
  colnames(ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes)[8] <- "READ_END"
  
ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes <- read.delim("ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes.bed", header=F)
  colnames(ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes)[13:16] <- c("CHROM","BED_CDS_START","BED_CDS_END","NAME")
  colnames(ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes)[1] <- "READ_CHROM"
  colnames(ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes)[2] <- "OVERLAP_START"
  colnames(ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes)[3] <- "OVERLAP_END"
  colnames(ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes)[4] <- "QNAME"
  colnames(ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes)[6] <- "STRAND"
  colnames(ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes)[7] <- "READ_START"
  colnames(ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes)[8] <- "READ_END"  
  
ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes <- read.delim("ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes.bed", header=F)
  colnames(ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes)[13:16] <- c("CHROM","BED_CDS_START","BED_CDS_END","NAME")
  colnames(ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes)[1] <- "READ_CHROM"
  colnames(ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes)[2] <- "OVERLAP_START"
  colnames(ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes)[3] <- "OVERLAP_END"
  colnames(ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes)[4] <- "QNAME"
  colnames(ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes)[6] <- "STRAND"
  colnames(ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes)[7] <- "READ_START"
  colnames(ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes)[8] <- "READ_END"  
  
ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes <- read.delim("ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes.bed", header=F)
  colnames(ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes)[13:16] <- c("CHROM","BED_CDS_START","BED_CDS_END","NAME")
  colnames(ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes)[1] <- "READ_CHROM"
  colnames(ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes)[2] <- "OVERLAP_START"
  colnames(ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes)[3] <- "OVERLAP_END"
  colnames(ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes)[4] <- "QNAME"
  colnames(ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes)[6] <- "STRAND"
  colnames(ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes)[7] <- "READ_START"
  colnames(ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes)[8] <- "READ_END"  
  
ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes <- read.delim("ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes.bed", header=F)
  colnames(ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[13:16] <- c("CHROM","BED_CDS_START","BED_CDS_END","NAME")
  colnames(ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[1] <- "READ_CHROM"
  colnames(ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[2] <- "OVERLAP_START"
  colnames(ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[3] <- "OVERLAP_END"
  colnames(ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[4] <- "QNAME"
  colnames(ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[6] <- "STRAND"
  colnames(ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[7] <- "READ_START"
  colnames(ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes)[8] <- "READ_END"    
  
  
  
#####calculate number of reads per frag for second batch#####
  
#ATAC_507_4h_exp1_reads_strand_upstream_Pberghei_genes
ATAC_507_4h_exp1_strand_CDS_reads_per_fragment <- as.data.frame(table(ATAC_507_4h_exp1_reads_strand_CDS_Pberghei_genes$NAME))
  colnames(ATAC_507_4h_exp1_strand_CDS_reads_per_fragment) <- c("NAME", "NUMBER_READS")
  ATAC_507_4h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_507_4h_exp1_strand_CDS_reads_per_fragment, y = Pbeghei_CDS_regions_strand, by = "NAME")
    duprows <- Pbeghei_CDS_regions_strand$NAME %in%  ATAC_507_4h_exp1_strand_CDS_reads_per_fragment$NAME
    ATAC_507_4h_exp1_strand_CDS_reads_per_fragment <- rbind.fill(ATAC_507_4h_exp1_strand_CDS_reads_per_fragment, Pbeghei_CDS_regions_strand[!duprows,])
    ATAC_507_4h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_507_4h_exp1_strand_CDS_reads_per_fragment, y = Pberghei_GFF3_CDS_strand_calc[,2:5], by.x = "NAME", by.y = "GENEID", sort = F)
    ATAC_507_4h_exp1_strand_CDS_reads_per_fragment$CDS_SIZE <- ATAC_507_4h_exp1_strand_CDS_reads_per_fragment$BED_CDS_END - ATAC_507_4h_exp1_strand_CDS_reads_per_fragment$BED_CDS_START
    for (i in seq_along(ATAC_507_4h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS)){
      if (is.na(ATAC_507_4h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i])) {
        ATAC_507_4h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i] <- 0
      }
    }
    ATAC_507_4h_exp1_strand_CDS_reads_per_fragment <- ATAC_507_4h_exp1_strand_CDS_reads_per_fragment[,c(3,1,2,6,7,8,4,5,9)]
    colnames(ATAC_507_4h_exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
  write.csv(ATAC_507_4h_exp1_strand_CDS_reads_per_fragment, "ATAC_507_4h_exp1_strand_CDS_reads_per_fragment.csv")
  
#ATAC_507_10h_exp1_reads_strand_upstream_Pberghei_genes
ATAC_507_10h_exp1_strand_CDS_reads_per_fragment <- as.data.frame(table(ATAC_507_10h_exp1_reads_strand_CDS_Pberghei_genes$NAME))
  colnames(ATAC_507_10h_exp1_strand_CDS_reads_per_fragment) <- c("NAME", "NUMBER_READS")
  ATAC_507_10h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_507_10h_exp1_strand_CDS_reads_per_fragment, y = Pbeghei_CDS_regions_strand, by = "NAME")
    duprows <- Pbeghei_CDS_regions_strand$NAME %in%  ATAC_507_10h_exp1_strand_CDS_reads_per_fragment$NAME
    ATAC_507_10h_exp1_strand_CDS_reads_per_fragment <- rbind.fill(ATAC_507_10h_exp1_strand_CDS_reads_per_fragment, Pbeghei_CDS_regions_strand[!duprows,])
    ATAC_507_10h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_507_10h_exp1_strand_CDS_reads_per_fragment, y = Pberghei_GFF3_CDS_strand_calc[,2:5], by.x = "NAME", by.y = "GENEID", sort = F)
    ATAC_507_10h_exp1_strand_CDS_reads_per_fragment$CDS_SIZE <- ATAC_507_10h_exp1_strand_CDS_reads_per_fragment$BED_CDS_END - ATAC_507_10h_exp1_strand_CDS_reads_per_fragment$BED_CDS_START
    for (i in seq_along(ATAC_507_10h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS)){
      if (is.na(ATAC_507_10h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i])) {
        ATAC_507_10h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i] <- 0
      }
    }
    ATAC_507_10h_exp1_strand_CDS_reads_per_fragment <- ATAC_507_10h_exp1_strand_CDS_reads_per_fragment[,c(3,1,2,6,7,8,4,5,9)]
    colnames(ATAC_507_10h_exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
  write.csv(ATAC_507_10h_exp1_strand_CDS_reads_per_fragment, "ATAC_507_10h_exp1_strand_CDS_reads_per_fragment.csv")
  
#ATAC_507_16h_exp1_reads_strand_upstream_Pberghei_genes
ATAC_507_16h_exp1_strand_CDS_reads_per_fragment <- as.data.frame(table(ATAC_507_16h_exp1_reads_strand_CDS_Pberghei_genes$NAME))
  colnames(ATAC_507_16h_exp1_strand_CDS_reads_per_fragment) <- c("NAME", "NUMBER_READS")
  ATAC_507_16h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_507_16h_exp1_strand_CDS_reads_per_fragment, y = Pbeghei_CDS_regions_strand, by = "NAME")
    duprows <- Pbeghei_CDS_regions_strand$NAME %in%  ATAC_507_16h_exp1_strand_CDS_reads_per_fragment$NAME
    ATAC_507_16h_exp1_strand_CDS_reads_per_fragment <- rbind.fill(ATAC_507_16h_exp1_strand_CDS_reads_per_fragment, Pbeghei_CDS_regions_strand[!duprows,])
    ATAC_507_16h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_507_16h_exp1_strand_CDS_reads_per_fragment, y = Pberghei_GFF3_CDS_strand_calc[,2:5], by.x = "NAME", by.y = "GENEID", sort = F)
    ATAC_507_16h_exp1_strand_CDS_reads_per_fragment$CDS_SIZE <- ATAC_507_16h_exp1_strand_CDS_reads_per_fragment$BED_CDS_END - ATAC_507_16h_exp1_strand_CDS_reads_per_fragment$BED_CDS_START
    for (i in seq_along(ATAC_507_16h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS)){
      if (is.na(ATAC_507_16h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i])) {
        ATAC_507_16h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i] <- 0
      }
    }
    ATAC_507_16h_exp1_strand_CDS_reads_per_fragment <- ATAC_507_16h_exp1_strand_CDS_reads_per_fragment[,c(3,1,2,6,7,8,4,5,9)]
    colnames(ATAC_507_16h_exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
  write.csv(ATAC_507_16h_exp1_strand_CDS_reads_per_fragment, "ATAC_507_16h_exp1_strand_CDS_reads_per_fragment.csv")


#ATAC_507_22h_exp1_reads_strand_upstream_Pberghei_genes
ATAC_507_22h_exp1_strand_CDS_reads_per_fragment <- as.data.frame(table(ATAC_507_22h_exp1_reads_strand_CDS_Pberghei_genes$NAME))
  colnames(ATAC_507_22h_exp1_strand_CDS_reads_per_fragment) <- c("NAME", "NUMBER_READS")
  ATAC_507_22h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_507_22h_exp1_strand_CDS_reads_per_fragment, y = Pbeghei_CDS_regions_strand, by = "NAME")
    duprows <- Pbeghei_CDS_regions_strand$NAME %in%  ATAC_507_22h_exp1_strand_CDS_reads_per_fragment$NAME
    ATAC_507_22h_exp1_strand_CDS_reads_per_fragment <- rbind.fill(ATAC_507_22h_exp1_strand_CDS_reads_per_fragment, Pbeghei_CDS_regions_strand[!duprows,])
    ATAC_507_22h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_507_22h_exp1_strand_CDS_reads_per_fragment, y = Pberghei_GFF3_CDS_strand_calc[,2:5], by.x = "NAME", by.y = "GENEID", sort = F)
    ATAC_507_22h_exp1_strand_CDS_reads_per_fragment$CDS_SIZE <- ATAC_507_22h_exp1_strand_CDS_reads_per_fragment$BED_CDS_END - ATAC_507_22h_exp1_strand_CDS_reads_per_fragment$BED_CDS_START
    for (i in seq_along(ATAC_507_22h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS)){
      if (is.na(ATAC_507_22h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i])) {
        ATAC_507_22h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i] <- 0
      }
    }
    ATAC_507_22h_exp1_strand_CDS_reads_per_fragment <- ATAC_507_22h_exp1_strand_CDS_reads_per_fragment[,c(3,1,2,6,7,8,4,5,9)]
    colnames(ATAC_507_22h_exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
  write.csv(ATAC_507_22h_exp1_strand_CDS_reads_per_fragment, "ATAC_507_22h_exp1_strand_CDS_reads_per_fragment.csv")
  
  
#ATAC_507_gDNA_exp1_reads_strand_upstream_Pberghei_genes
ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment <- as.data.frame(table(ATAC_507_gDNA_exp1_reads_strand_CDS_Pberghei_genes$NAME))
  colnames(ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment) <- c("NAME", "NUMBER_READS")
  ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment, y = Pbeghei_CDS_regions_strand, by = "NAME")
    duprows <- Pbeghei_CDS_regions_strand$NAME %in%  ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment$NAME
    ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment <- rbind.fill(ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment, Pbeghei_CDS_regions_strand[!duprows,])
    ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment, y = Pberghei_GFF3_CDS_strand_calc[,2:5], by.x = "NAME", by.y = "GENEID", sort = F)
    ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment$CDS_SIZE <- ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment$BED_CDS_END - ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment$BED_CDS_START
    for (i in seq_along(ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment$NUMBER_READS)){
      if (is.na(ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i])) {
        ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i] <- 0
      }
    }
    ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment <- ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment[,c(3,1,2,6,7,8,4,5,9)]
    colnames(ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
  write.csv(ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment, "ATAC_507_gDNA_exp1_strand_CDS_reads_per_fragment.csv")
  
  

#ATAC_1022_4h_exp1_reads_strand_upstream_Pberghei_genes
ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment <- as.data.frame(table(ATAC_1022_4h_exp1_reads_strand_CDS_Pberghei_genes$NAME))
  colnames(ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment) <- c("NAME", "NUMBER_READS")
  ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment, y = Pbeghei_CDS_regions_strand, by = "NAME")
    duprows <- Pbeghei_CDS_regions_strand$NAME %in%  ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment$NAME
    ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment <- rbind.fill(ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment, Pbeghei_CDS_regions_strand[!duprows,])
    ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment, y = Pberghei_GFF3_CDS_strand_calc[,2:5], by.x = "NAME", by.y = "GENEID", sort = F)
    ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment$CDS_SIZE <- ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment$BED_CDS_END - ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment$BED_CDS_START
    for (i in seq_along(ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS)){
      if (is.na(ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i])) {
        ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i] <- 0
      }
    }
    ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment <- ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment[,c(3,1,2,6,7,8,4,5,9)]
    colnames(ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
  write.csv(ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment, "ATAC_1022_4h_exp1_strand_CDS_reads_per_fragment.csv")
  
#ATAC_1022_10h_exp1_reads_strand_upstream_Pberghei_genes
ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment <- as.data.frame(table(ATAC_1022_10h_exp1_reads_strand_CDS_Pberghei_genes$NAME))
  colnames(ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment) <- c("NAME", "NUMBER_READS")
  ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment, y = Pbeghei_CDS_regions_strand, by = "NAME")
    duprows <- Pbeghei_CDS_regions_strand$NAME %in%  ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment$NAME
    ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment <- rbind.fill(ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment, Pbeghei_CDS_regions_strand[!duprows,])
    ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment, y = Pberghei_GFF3_CDS_strand_calc[,2:5], by.x = "NAME", by.y = "GENEID", sort = F)
    ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment$CDS_SIZE <- ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment$BED_CDS_END - ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment$BED_CDS_START
    for (i in seq_along(ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS)){
      if (is.na(ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i])) {
        ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i] <- 0
      }
    }
    ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment <- ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment[,c(3,1,2,6,7,8,4,5,9)]
    colnames(ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
  write.csv(ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment, "ATAC_1022_10h_exp1_strand_CDS_reads_per_fragment.csv")
  
#ATAC_1022_16h_exp1_reads_strand_upstream_Pberghei_genes
 ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment <- as.data.frame(table(ATAC_1022_16h_exp1_reads_strand_CDS_Pberghei_genes$NAME))
  colnames(ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment) <- c("NAME", "NUMBER_READS")
  ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment, y = Pbeghei_CDS_regions_strand, by = "NAME")
    duprows <- Pbeghei_CDS_regions_strand$NAME %in%  ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment$NAME
    ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment <- rbind.fill(ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment, Pbeghei_CDS_regions_strand[!duprows,])
    ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment, y = Pberghei_GFF3_CDS_strand_calc[,2:5], by.x = "NAME", by.y = "GENEID", sort = F)
    ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment$CDS_SIZE <- ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment$BED_CDS_END - ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment$BED_CDS_START
    for (i in seq_along(ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS)){
      if (is.na(ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i])) {
        ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i] <- 0
      }
    }
    ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment <- ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment[,c(3,1,2,6,7,8,4,5,9)]
    colnames(ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
  write.csv(ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment, "ATAC_1022_16h_exp1_strand_CDS_reads_per_fragment.csv")
  
  
#ATAC_1022_22h_exp1_reads_strand_upstream_Pberghei_genes
ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment <- as.data.frame(table(ATAC_1022_22h_exp1_reads_strand_CDS_Pberghei_genes$NAME))
  colnames(ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment) <- c("NAME", "NUMBER_READS")
  ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment, y = Pbeghei_CDS_regions_strand, by = "NAME")
    duprows <- Pbeghei_CDS_regions_strand$NAME %in%  ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment$NAME
    ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment <- rbind.fill(ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment, Pbeghei_CDS_regions_strand[!duprows,])
    ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment, y = Pberghei_GFF3_CDS_strand_calc[,2:5], by.x = "NAME", by.y = "GENEID", sort = F)
    ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment$CDS_SIZE <- ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment$BED_CDS_END - ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment$BED_CDS_START
    for (i in seq_along(ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS)){
      if (is.na(ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i])) {
        ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i] <- 0
      }
    }
    ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment <- ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment[,c(3,1,2,6,7,8,4,5,9)]
    colnames(ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
  write.csv(ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment, "ATAC_1022_22h_exp1_strand_CDS_reads_per_fragment.csv")
  
  
#ATAC_1022_gDNA_exp1_reads_strand_upstream_Pberghei_genes
ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment <- as.data.frame(table(ATAC_1022_gDNA_exp1_reads_strand_CDS_Pberghei_genes$NAME))
  colnames(ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment) <- c("NAME", "NUMBER_READS")
  ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment, y = Pbeghei_CDS_regions_strand, by = "NAME")
    duprows <- Pbeghei_CDS_regions_strand$NAME %in%  ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment$NAME
    ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment <- rbind.fill(ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment, Pbeghei_CDS_regions_strand[!duprows,])
    ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment <- merge(x = ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment, y = Pberghei_GFF3_CDS_strand_calc[,2:5], by.x = "NAME", by.y = "GENEID", sort = F)
    ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment$CDS_SIZE <- ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment$BED_CDS_END - ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment$BED_CDS_START
    for (i in seq_along(ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment$NUMBER_READS)){
      if (is.na(ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i])) {
        ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment$NUMBER_READS[i] <- 0
      }
    }
    ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment <- ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment[,c(3,1,2,6,7,8,4,5,9)]
    colnames(ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
  write.csv(ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment, "ATAC_1022_gDNA_exp1_strand_CDS_reads_per_fragment.csv")
  
  
#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/Pberghei_CDS_region_RPFM_calculation.RData")
