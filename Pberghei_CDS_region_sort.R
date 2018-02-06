#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/Pberghei_CDS_region_sort.RData")
library("data.table")
library("ape")

#####set WD#####
setwd("J:/III/Waters/Group Members/Mallory/PbergheiUpstreamRegionsSebastian/")

#####load data#####
Pberghei_GFF3_allCDS_strand_calc <- read.csv("GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv")
Pberghei_chrom_coord <- read.csv("Pberghei_chromosome_coord.csv")


#####remove genes with an UP_DIST <1 (no bases in genes)#####
pos_Pberghei_GFF3_allCDS_strand <- subset(Pberghei_GFF3_allCDS_strand_calc, Pberghei_GFF3_allCDS_strand_calc$SIZE >= 3)
  write.csv(pos_Pberghei_GFF3_allCDS_strand, "pos_Pberghei_GFF3_allCDS_strand.csv")
neg_Pberghei_GFF3_allCDSstrand <- subset(Pberghei_GFF3_allCDS_strand_calc, Pberghei_GFF3_allCDS_strand_calc$SIZE <= 2)
  write.csv(neg_Pberghei_GFF3_allCDSstrand, "neg_Pberghei_GFF3_allCDS_strand.csv")


#the next step is to figure out the beginning and end of the CDS region for all genes, taking into account the BED format
#STOP = STOP
#START = START - 1
  
  #no calculations needed
  #end coordinate stays the same
  #start coordinate needs to be one base less
  #strand won't make a difference

Pberghei_CDS_regions_strand <- data.frame(matrix(NA, nrow = 13354))
Pberghei_CDS_regions_strand$CHROM <-  pos_Pberghei_GFF3_allCDS_strand$SEQID
Pberghei_CDS_regions_strand$GENEID <- pos_Pberghei_GFF3_allCDS_strand$GENEID
Pberghei_CDS_regions_strand <- Pberghei_CDS_regions_strand[,-1]
Pberghei_CDS_regions_strand$STRAND <- pos_Pberghei_GFF3_allCDS_strand$STRAND
Pberghei_CDS_regions_strand$CDS_START <- pos_Pberghei_GFF3_allCDS_strand$START
Pberghei_CDS_regions_strand$CDS_END <- pos_Pberghei_GFF3_allCDS_strand$END
Pberghei_CDS_regions_strand$CDS_BED_START <- NA
Pberghei_CDS_regions_strand$CDS_BED_END <- NA

for (i in seq_along(Pberghei_CDS_regions_strand$GENEID)) {
  Pberghei_CDS_regions_strand$CDS_BED_START[i] <- Pberghei_CDS_regions_strand$CDS_START[i] - 1
  Pberghei_CDS_regions_strand$CDS_BED_END <- Pberghei_CDS_regions_strand$CDS_END
  }

write.csv(Pberghei_CDS_regions_strand$GENEID, "Pberghei_CDS_regions_strand_genelist.csv")


####reformat it into BED file####

# BED columns: CHROM START END NAME

Pberghei_CDS_regions_strand_download_BED <- data.frame(matrix(nrow = 13354))
Pberghei_CDS_regions_strand_download_BED$CHROM <- Pberghei_CDS_regions_strand$CHROM
Pberghei_CDS_regions_strand_download_BED$START <- Pberghei_CDS_regions_strand$UPSTART
Pberghei_CDS_regions_strand_download_BED$END <- Pberghei_CDS_regions_strand$UPEND
Pberghei_CDS_regions_strand_download_BED$NAME <- Pberghei_CDS_regions_strand$GENEID
Pberghei_CDS_regions_strand_download_BED <- Pberghei_CDS_regions_strand[,2:5]
write.csv(Pberghei_CDS_regions_strand_download_BED, "Pberghei_CDS_regions_strand_download_BED.csv")




######notes#####
#GFF3 file is in 1-based coordinate system
#BED file needs to be in 0-based coordinate system


#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/Pberghei_CDS_region_sort.RData")
