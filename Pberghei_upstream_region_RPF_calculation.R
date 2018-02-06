#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/Pberghei_upstream_region_RPF_calculation.RData")
require("plyr")

#NOTE -> this calculates the number of reads, but DOES NOT include the total number of reads in the initial BAM files
#bascially it's giving absolute non-normalized counts
#Sebastian is fine with this, as he can do the normalizing himself - he normalizes to the total number of reads found in all upstream regions, etc.

#####Function#####
#this function lets you pick your GFF3 file, 
#your region BED file (1000bp upstream, 2000bp upstream, etc - CHROM, UPSTART, UPEND, NAME)
#and your reads BED file

read_bedfile_calculate_reads_per_frag <- function(filename, gff_filename, BED_filename, WD){
  setwd("J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/")
  gff <- read.csv(gff_filename, header=T)
  colnames(gff)[3] <- "GENE_START"
  colnames(gff)[4] <- "GENE_END"
  setwd("J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/GenomeBEDfiles/")
  bed <- read.delim(BED_filename, header=F)
  colnames(bed) <- c("CHROM","UPSTART","UPEND","NAME")
  setwd(WD)
  y <- read.delim(filename, header=F)
  colnames(y)[13:16] <- c("CHROM","UPSTART","UPEND","NAME")
  colnames(y)[1:4] <- c("READ_CHROM","OVERLAP_START","OVERLAP_END","QNAME")
  colnames(y)[6:8] <- c("STRAND","READ_START","READ_END")
  # return(y)
  q <- as.data.frame(table(y$NAME))
  colnames(q) <- c("NAME","NUMBER_READS")
  q <- merge(x = q, 
             y = bed, 
             by = "NAME")
  duprows <- bed$NAME %in%  q$NAME
  q <- rbind.fill(q, 
                  bed[!duprows,])
  q <- merge(x = q, 
             y = gff[,2:5], 
             by.x = "NAME", by.y = "GENEID", sort = F)
  q$UPSIZE <- q$UPEND - q$UPSTART
  for (i in seq_along(q$NUMBER_READS)){
    if (is.na(q$NUMBER_READS[i])) {
      q$NUMBER_READS[i] <- 0
    }
  }
  colnames(q)[1] <- "GENEID"
  return(q[,c(3,1,2,6,7,8,4,5,9)])
}   


##########AllReads1000bp##########
###EXP1###  
#507 
ATAC_507_4h_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_507_4h_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_4h_Exp1_strand_1000_reads_per_fragment, "ATAC_507_4h_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_507_10h_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_507_10h_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_10h_Exp1_strand_1000_reads_per_fragment, "ATAC_507_10h_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_507_16h_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_16h_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_507_16h_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_16h_Exp1_strand_1000_reads_per_fragment, "ATAC_507_16h_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_507_22h_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_22h_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_507_22h_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_22h_Exp1_strand_1000_reads_per_fragment, "ATAC_507_22h_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_507_gDNA_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gDNA_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_507_gDNA_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gDNA_Exp1_strand_1000_reads_per_fragment, "ATAC_507_gDNA_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_507_gam_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gam_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_507_gam_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gam_Exp1_strand_1000_reads_per_fragment, "ATAC_507_gam_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_507_ook_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_ook_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_507_ook_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_ook_Exp1_strand_1000_reads_per_fragment, "ATAC_507_ook_Exp1_strand_1000_reads_per_fragment.csv")

#1022 
ATAC_1022_4h_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_4h_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_1022_4h_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_4h_Exp1_strand_1000_reads_per_fragment, "ATAC_1022_4h_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_1022_10h_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_1022_10h_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_10h_Exp1_strand_1000_reads_per_fragment, "ATAC_1022_10h_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_1022_16h_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_16h_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_1022_16h_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_16h_Exp1_strand_1000_reads_per_fragment, "ATAC_1022_16h_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_1022_22h_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_22h_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_1022_22h_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_22h_Exp1_strand_1000_reads_per_fragment, "ATAC_1022_22h_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_1022_gDNA_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_gDNA_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                            gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                            BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                            WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_1022_gDNA_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_gDNA_Exp1_strand_1000_reads_per_fragment, "ATAC_1022_gDNA_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_1022_gam_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_gam_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_1022_gam_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_gam_Exp1_strand_1000_reads_per_fragment, "ATAC_1022_gam_Exp1_strand_1000_reads_per_fragment.csv")
ATAC_1022_ook_Exp1_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_ook_Exp1_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/")
colnames(ATAC_1022_ook_Exp1_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_ook_Exp1_strand_1000_reads_per_fragment, "ATAC_1022_ook_Exp1_strand_1000_reads_per_fragment.csv")    


###Exp2###  
#507 
ATAC_507_4h_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_507_4h_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_4h_Exp2_strand_1000_reads_per_fragment, "ATAC_507_4h_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_507_10h_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_507_10h_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_10h_Exp2_strand_1000_reads_per_fragment, "ATAC_507_10h_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_507_16h_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_16h_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_507_16h_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_16h_Exp2_strand_1000_reads_per_fragment, "ATAC_507_16h_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_507_22h_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_22h_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_507_22h_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_22h_Exp2_strand_1000_reads_per_fragment, "ATAC_507_22h_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_507_gDNA_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gDNA_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_507_gDNA_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gDNA_Exp2_strand_1000_reads_per_fragment, "ATAC_507_gDNA_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_507_gam_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gam_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_507_gam_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gam_Exp2_strand_1000_reads_per_fragment, "ATAC_507_gam_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_507_ook_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_ook_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_507_ook_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_ook_Exp2_strand_1000_reads_per_fragment, "ATAC_507_ook_Exp2_strand_1000_reads_per_fragment.csv")

#1022 
ATAC_1022_4h_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_4h_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_1022_4h_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_4h_Exp2_strand_1000_reads_per_fragment, "ATAC_1022_4h_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_1022_10h_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_1022_10h_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_10h_Exp2_strand_1000_reads_per_fragment, "ATAC_1022_10h_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_1022_16h_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_16h_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_1022_16h_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_16h_Exp2_strand_1000_reads_per_fragment, "ATAC_1022_16h_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_1022_22h_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_22h_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_1022_22h_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_22h_Exp2_strand_1000_reads_per_fragment, "ATAC_1022_22h_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_1022_gDNA_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_gDNA_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                            gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                            BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                            WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_1022_gDNA_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_gDNA_Exp2_strand_1000_reads_per_fragment, "ATAC_1022_gDNA_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_1022_gam_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_gam_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_1022_gam_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_gam_Exp2_strand_1000_reads_per_fragment, "ATAC_1022_gam_Exp2_strand_1000_reads_per_fragment.csv")
ATAC_1022_ook_Exp2_strand_1000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_ook_Exp2_reads_strand_1000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/")
colnames(ATAC_1022_ook_Exp2_strand_1000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_ook_Exp2_strand_1000_reads_per_fragment, "ATAC_1022_ook_Exp2_strand_1000_reads_per_fragment.csv")  

##########AllReads2000bp##########
###EXP1###  
#507 
ATAC_507_4h_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_507_4h_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_4h_Exp1_strand_2000_reads_per_fragment, "ATAC_507_4h_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_507_10h_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_507_10h_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_10h_Exp1_strand_2000_reads_per_fragment, "ATAC_507_10h_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_507_16h_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_16h_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_507_16h_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_16h_Exp1_strand_2000_reads_per_fragment, "ATAC_507_16h_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_507_22h_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_22h_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_507_22h_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_22h_Exp1_strand_2000_reads_per_fragment, "ATAC_507_22h_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_507_gDNA_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gDNA_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_507_gDNA_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gDNA_Exp1_strand_2000_reads_per_fragment, "ATAC_507_gDNA_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_507_gam_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gam_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_507_gam_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gam_Exp1_strand_2000_reads_per_fragment, "ATAC_507_gam_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_507_ook_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_ook_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_507_ook_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_ook_Exp1_strand_2000_reads_per_fragment, "ATAC_507_ook_Exp1_strand_2000_reads_per_fragment.csv")

#1022 
ATAC_1022_4h_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_4h_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_1022_4h_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_4h_Exp1_strand_2000_reads_per_fragment, "ATAC_1022_4h_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_1022_10h_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_1022_10h_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_10h_Exp1_strand_2000_reads_per_fragment, "ATAC_1022_10h_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_1022_16h_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_16h_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_1022_16h_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_16h_Exp1_strand_2000_reads_per_fragment, "ATAC_1022_16h_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_1022_22h_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_22h_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_1022_22h_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_22h_Exp1_strand_2000_reads_per_fragment, "ATAC_1022_22h_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_1022_gDNA_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_gDNA_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                            gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                            BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                            WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_1022_gDNA_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_gDNA_Exp1_strand_2000_reads_per_fragment, "ATAC_1022_gDNA_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_1022_gam_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_gam_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_1022_gam_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_gam_Exp1_strand_2000_reads_per_fragment, "ATAC_1022_gam_Exp1_strand_2000_reads_per_fragment.csv")
ATAC_1022_ook_Exp1_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_ook_Exp1_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/")
colnames(ATAC_1022_ook_Exp1_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_ook_Exp1_strand_2000_reads_per_fragment, "ATAC_1022_ook_Exp1_strand_2000_reads_per_fragment.csv")    

###Exp2###  
#507 
ATAC_507_4h_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_507_4h_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_4h_Exp2_strand_2000_reads_per_fragment, "ATAC_507_4h_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_507_10h_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_507_10h_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_10h_Exp2_strand_2000_reads_per_fragment, "ATAC_507_10h_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_507_16h_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_16h_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_507_16h_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_16h_Exp2_strand_2000_reads_per_fragment, "ATAC_507_16h_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_507_22h_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_22h_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_507_22h_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_22h_Exp2_strand_2000_reads_per_fragment, "ATAC_507_22h_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_507_gDNA_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gDNA_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_507_gDNA_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gDNA_Exp2_strand_2000_reads_per_fragment, "ATAC_507_gDNA_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_507_gam_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gam_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_507_gam_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gam_Exp2_strand_2000_reads_per_fragment, "ATAC_507_gam_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_507_ook_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_ook_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_507_ook_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_ook_Exp2_strand_2000_reads_per_fragment, "ATAC_507_ook_Exp2_strand_2000_reads_per_fragment.csv")

#1022 
ATAC_1022_4h_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_4h_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_1022_4h_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_4h_Exp2_strand_2000_reads_per_fragment, "ATAC_1022_4h_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_1022_10h_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_1022_10h_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_10h_Exp2_strand_2000_reads_per_fragment, "ATAC_1022_10h_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_1022_16h_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_16h_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_1022_16h_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_16h_Exp2_strand_2000_reads_per_fragment, "ATAC_1022_16h_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_1022_22h_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_22h_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_1022_22h_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_22h_Exp2_strand_2000_reads_per_fragment, "ATAC_1022_22h_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_1022_gDNA_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_gDNA_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                            gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                            BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                            WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_1022_gDNA_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_gDNA_Exp2_strand_2000_reads_per_fragment, "ATAC_1022_gDNA_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_1022_gam_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_gam_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_1022_gam_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_gam_Exp2_strand_2000_reads_per_fragment, "ATAC_1022_gam_Exp2_strand_2000_reads_per_fragment.csv")
ATAC_1022_ook_Exp2_strand_2000_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_ook_Exp2_reads_strand_2000_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/")
colnames(ATAC_1022_ook_Exp2_strand_2000_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_ook_Exp2_strand_2000_reads_per_fragment, "ATAC_1022_ook_Exp2_strand_2000_reads_per_fragment.csv")  

#########AllReadsCDS##########
###EXP1###  
#507 
ATAC_507_4h_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                        gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                        BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                        WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_507_4h_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_4h_Exp1_strand_CDS_reads_per_fragment, "ATAC_507_4h_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_507_10h_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_507_10h_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_10h_Exp1_strand_CDS_reads_per_fragment, "ATAC_507_10h_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_507_16h_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_16h_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_507_16h_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_16h_Exp1_strand_CDS_reads_per_fragment, "ATAC_507_16h_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_507_22h_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_22h_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_507_22h_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_22h_Exp1_strand_CDS_reads_per_fragment, "ATAC_507_22h_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_507_gDNA_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gDNA_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_507_gDNA_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gDNA_Exp1_strand_CDS_reads_per_fragment, "ATAC_507_gDNA_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_507_gam_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gam_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_507_gam_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gam_Exp1_strand_CDS_reads_per_fragment, "ATAC_507_gam_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_507_ook_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_ook_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_507_ook_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_ook_Exp1_strand_CDS_reads_per_fragment, "ATAC_507_ook_Exp1_strand_CDS_reads_per_fragment.csv")

#1022 
ATAC_1022_4h_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_4h_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_1022_4h_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_4h_Exp1_strand_CDS_reads_per_fragment, "ATAC_1022_4h_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_1022_10h_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_1022_10h_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_10h_Exp1_strand_CDS_reads_per_fragment, "ATAC_1022_10h_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_1022_16h_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_16h_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_1022_16h_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_16h_Exp1_strand_CDS_reads_per_fragment, "ATAC_1022_16h_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_1022_22h_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_22h_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_1022_22h_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_22h_Exp1_strand_CDS_reads_per_fragment, "ATAC_1022_22h_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_1022_gDNA_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_gDNA_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                           gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                           BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                           WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_1022_gDNA_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_gDNA_Exp1_strand_CDS_reads_per_fragment, "ATAC_1022_gDNA_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_1022_gam_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_gam_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_1022_gam_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_gam_Exp1_strand_CDS_reads_per_fragment, "ATAC_1022_gam_Exp1_strand_CDS_reads_per_fragment.csv")
ATAC_1022_ook_Exp1_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_ook_Exp1_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/")
colnames(ATAC_1022_ook_Exp1_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_ook_Exp1_strand_CDS_reads_per_fragment, "ATAC_1022_ook_Exp1_strand_CDS_reads_per_fragment.csv")    

###Exp2###  
#507 
ATAC_507_4h_Exp2_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                        gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                        BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                        WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/")
colnames(ATAC_507_4h_Exp2_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_4h_Exp2_strand_CDS_reads_per_fragment, "ATAC_507_4h_Exp2_strand_CDS_reads_per_fragment.csv")
ATAC_507_10h_Exp2_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/")
colnames(ATAC_507_10h_Exp2_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_10h_Exp2_strand_CDS_reads_per_fragment, "ATAC_507_10h_Exp2_strand_CDS_reads_per_fragment.csv")
ATAC_507_16h_Exp2_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_16h_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/")
colnames(ATAC_507_16h_Exp2_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_16h_Exp2_strand_CDS_reads_per_fragment, "ATAC_507_16h_Exp2_strand_CDS_reads_per_fragment.csv")
ATAC_507_22h_Exp2_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_22h_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/")
colnames(ATAC_507_22h_Exp2_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_22h_Exp2_strand_CDS_reads_per_fragment, "ATAC_507_22h_Exp2_strand_CDS_reads_per_fragment.csv")
ATAC_507_gDNA_Exp2_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gDNA_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/")
colnames(ATAC_507_gDNA_Exp2_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gDNA_Exp2_strand_CDS_reads_per_fragment, "ATAC_507_gDNA_Exp2_strand_CDS_reads_per_fragment.csv")
ATAC_507_gam_Exp2_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_gam_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/")
colnames(ATAC_507_gam_Exp2_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_gam_Exp2_strand_CDS_reads_per_fragment, "ATAC_507_gam_Exp2_strand_CDS_reads_per_fragment.csv")
ATAC_507_ook_Exp2_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_507_ook_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/")
colnames(ATAC_507_ook_Exp2_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_507_ook_Exp2_strand_CDS_reads_per_fragment, "ATAC_507_ook_Exp2_strand_CDS_reads_per_fragment.csv")

#1022 
ATAC_1022_4h_Exp2_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_4h_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                         gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                         BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                         WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/")
colnames(ATAC_1022_4h_Exp2_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_4h_Exp2_strand_CDS_reads_per_fragment, "ATAC_1022_4h_Exp2_strand_CDS_reads_per_fragment.csv")
ATAC_1022_10h_Exp2_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/")
colnames(ATAC_1022_10h_Exp2_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_10h_Exp2_strand_CDS_reads_per_fragment, "ATAC_1022_10h_Exp2_strand_CDS_reads_per_fragment.csv")
ATAC_1022_16h_Exp2_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_16h_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/")
colnames(ATAC_1022_16h_Exp2_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_16h_Exp2_strand_CDS_reads_per_fragment, "ATAC_1022_16h_Exp2_strand_CDS_reads_per_fragment.csv")
ATAC_1022_22h_Exp2_strand_CDS_reads_per_fragment <- read_bedfile_calculate_reads_per_frag("ATAC_1022_22h_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                                                          gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                                                          BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                                                          WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/")
colnames(ATAC_1022_22h_Exp2_strand_CDS_reads_per_fragment)[2] <- "GENEID"
write.csv(ATAC_1022_22h_Exp2_strand_CDS_reads_per_fragment, "ATAC_1022_22h_Exp2_strand_CDS_reads_per_fragment.csv")




write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"),"ATAC_1022_gDNA_Exp2_strand_CDS_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"),"ATAC_1022_gam_Exp2_strand_CDS_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_CDS_upstream_Pberghei_genes.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"),"ATAC_1022_ook_Exp2_strand_CDS_reads_per_fragment.csv")

###########SubNucleo1000bp###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"),"ATAC_507_4h_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"),"ATAC_507_10h_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_16h_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_22h_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"),"ATAC_507_gDNA_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",                                                                                                    
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_gam_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_ook_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_4h_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_10h_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_16h_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_22h_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_gDNA_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_gam_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_ook_Exp1_strand_1000_subNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_4h_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_10h_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_16h_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_22h_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_gDNA_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_gam_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_ook_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_4h_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_10h_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_16h_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_22h_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_gDNA_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_gam_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_1000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_ook_Exp2_strand_1000_subNucleo_reads_per_fragment.csv")
###########SubNucleo2000bp###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_4h_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_10h_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_16h_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_22h_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_gDNA_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_gam_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_ook_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_4h_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_10h_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_16h_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_22h_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_gDNA_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_gam_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_ook_Exp1_strand_2000_subNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_4h_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_10h_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_16h_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_22h_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_gDNA_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_gam_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_ook_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_4h_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_10h_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_16h_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_22h_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_gDNA_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_gam_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_2000_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_ook_Exp2_strand_2000_subNucleo_reads_per_fragment.csv")
###########SubNucleoCDS###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_4h_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_10h_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_16h_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_22h_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_gDNA_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_gam_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_ook_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_4h_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_10h_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_16h_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_22h_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_gDNA_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_gam_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_ook_Exp1_strand_CDS_subNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_4h_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_10h_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_16h_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_22h_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_gDNA_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_gam_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_ook_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_4h_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_10h_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_16h_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_22h_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_gDNA_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_gam_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_CDS_upstream_Pberghei_genes_subNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_ook_Exp2_strand_CDS_subNucleo_reads_per_fragment.csv")


###########monoNucleo1000bp###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_4h_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_10h_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_16h_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_22h_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_gDNA_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_gam_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_ook_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_4h_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_10h_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_16h_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_22h_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_gDNA_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_gam_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_ook_Exp1_strand_1000_monoNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_4h_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_10h_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_16h_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_22h_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_gDNA_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_gam_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_ook_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_4h_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_10h_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_16h_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_22h_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_gDNA_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_gam_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_1000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_ook_Exp2_strand_1000_monoNucleo_reads_per_fragment.csv")
###########monoNucleo2000bp###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_4h_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_10h_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_16h_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_22h_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_gDNA_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_gam_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_ook_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_4h_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_10h_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_16h_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_22h_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_gDNA_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_gam_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_ook_Exp1_strand_2000_monoNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_4h_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_10h_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_16h_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_22h_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                      gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                      BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                      WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_gDNA_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_gam_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_ook_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_4h_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_10h_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_16h_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_22h_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_gDNA_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_gam_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_2000_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_ook_Exp2_strand_2000_monoNucleo_reads_per_fragment.csv")
###########monoNucleoCDS###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_4h_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_10h_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_16h_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_22h_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_gDNA_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_gam_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_ook_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_4h_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_10h_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_16h_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_22h_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_gDNA_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_gam_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_ook_Exp1_strand_CDS_monoNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_4h_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_10h_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_16h_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_22h_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_gDNA_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_gam_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_ook_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_4h_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_10h_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_16h_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_22h_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_gDNA_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_gam_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_CDS_upstream_Pberghei_genes_monoNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_ook_Exp2_strand_CDS_monoNucleo_reads_per_fragment.csv")


###########interNucleo1000bp###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_4h_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_10h_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_16h_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_22h_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_gDNA_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_gam_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_ook_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_4h_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_10h_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_16h_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_22h_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_gDNA_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_gam_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_ook_Exp1_strand_1000_interNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_4h_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_10h_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_16h_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_22h_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_gDNA_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_gam_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_ook_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_4h_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_10h_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_16h_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_22h_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_gDNA_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_gam_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_1000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_ook_Exp2_strand_1000_interNucleo_reads_per_fragment.csv")
###########interNucleo2000bp###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_4h_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_10h_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_16h_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_22h_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_gDNA_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_gam_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_ook_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_4h_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_10h_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_16h_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_22h_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_gDNA_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_gam_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_ook_Exp1_strand_2000_interNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_4h_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_10h_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_16h_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_22h_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_gDNA_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_gam_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_ook_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_4h_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_10h_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_16h_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_22h_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_gDNA_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_gam_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_2000_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_ook_Exp2_strand_2000_interNucleo_reads_per_fragment.csv")

###########interNucleoCDS###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_4h_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_10h_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_16h_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_22h_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_gDNA_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_gam_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_ook_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_4h_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_10h_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_16h_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_22h_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_gDNA_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_gam_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_ook_Exp1_strand_CDS_interNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_4h_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_10h_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_16h_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_22h_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_gDNA_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_gam_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_ook_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_4h_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_10h_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_16h_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_22h_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_gDNA_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_gam_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_CDS_upstream_Pberghei_genes_interNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_ook_Exp2_strand_CDS_interNucleo_reads_per_fragment.csv")



###########diNucleo1000bp###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_4h_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_10h_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_16h_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_22h_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_gDNA_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_gam_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_507_ook_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_4h_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_10h_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_16h_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_22h_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_gDNA_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_gam_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/"), "ATAC_1022_ook_Exp1_strand_1000_diNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_4h_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_10h_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_16h_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_22h_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_gDNA_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_gam_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_507_ook_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_4h_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_10h_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_16h_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_22h_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_gDNA_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_gam_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_1000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_1000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/"), "ATAC_1022_ook_Exp2_strand_1000_diNucleo_reads_per_fragment.csv")
###########diNucleo2000bp###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_4h_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_10h_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_16h_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_22h_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_gDNA_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_gam_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_507_ook_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_4h_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_10h_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_16h_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_22h_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_gDNA_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_gam_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/"), "ATAC_1022_ook_Exp1_strand_2000_diNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_4h_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_10h_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_16h_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_22h_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_gDNA_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_gam_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_507_ook_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_4h_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_10h_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_16h_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_22h_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_gDNA_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_gam_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_2000_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gene_StartStop_strand_calc.csv",
                                                BED_filename = "Pberghei_upstream_regions_strand_2000_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/"), "ATAC_1022_ook_Exp2_strand_2000_diNucleo_reads_per_fragment.csv")
###########diNucleoCDS###########
###Exp1###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_4h_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_10h_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_16h_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_22h_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_gDNA_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                               gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                               BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                               WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_gam_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_507_ook_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_4h_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_10h_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_16h_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_22h_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_gDNA_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_gam_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp1_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/"), "ATAC_1022_ook_Exp1_strand_CDS_diNucleo_reads_per_fragment.csv")
###Exp2###
#507#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_4h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_4h_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_10h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_10h_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_16h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_16h_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_22h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_22h_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gDNA_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_gDNA_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_gam_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_gam_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_507_ook_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_507_ook_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")
#1022#
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_4h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_4h_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_10h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_10h_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_16h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_16h_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_22h_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_22h_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")  
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gDNA_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_gDNA_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_gam_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_gam_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")
write.csv(read_bedfile_calculate_reads_per_frag(filename = "ATAC_1022_ook_Exp2_reads_strand_CDS_upstream_Pberghei_genes_diNucleo.bed",
                                                gff_filename = "GeneDB_28_Pberghei_gff3_CDS_strand_calc.csv",
                                                BED_filename = "Pberghei_CDS_regions_strand_download_BED.bed",
                                                WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/"), "ATAC_1022_ook_Exp2_strand_CDS_diNucleo_reads_per_fragment.csv")
#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/Pberghei_upstream_region_RPF_calculation.RData")
