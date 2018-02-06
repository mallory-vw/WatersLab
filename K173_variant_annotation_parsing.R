#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/K173_variant_annotation_parsing.RData")

#####set wd#####
setwd("J:/III/Waters/Group Members/Mallory/GenomeProjects/GenomesVariants/AnnotationDetails")

library(data.table)

#####import data#####
K173_TrueSeq <- read.csv("K173_TrueSeq_variant_info.csv")[,1:8]
K173_TrueSeq_CDS <- read.csv("K173_TrueSeq_CDS_variant_info.csv")[,1:8]
K173_TrueSeq_gene <- read.csv("K173_TrueSeq_gene_variant_info.csv")[,1:8]
K173_TrueSeq_exon <- read.csv("K173_TrueSeq_exon_variant_info.csv")[,1:8]

K173_Kappa <- read.csv("K173_Kappa_variant_info.csv")[,1:8]
K173_Kappa_CDS <- read.csv("K173_Kappa_CDS_variant_info.csv")[,1:8]
K173_Kappa_gene <- read.csv("K173_Kappa_gene_variant_info.csv")[,1:8]
K173_Kappa_exon <- read.csv("K173_Kappa_exon_variant_info.csv")[,1:8]

K173_TrueSeq_nonSyn <- read.csv("K173_TrueSeq_nonSyn_variant_info.csv")[,1:8]
K173_TrueSeq_CDS_nonSyn <- read.csv("K173_TrueSeq_CDS_nonSyn_variant_info.csv")[,1:8]
K173_TrueSeq_gene_nonSyn <- read.csv("K173_TrueSeq_gene_nonSyn_variant_info.csv")[,1:8]
K173_TrueSeq_exon_nonSyn <- read.csv("K173_TrueSeq_exon_nonSyn_variant_info.csv")[,1:8]

K173_Kappa_nonSyn <- read.csv("K173_Kappa_nonSyn_variant_info.csv")[,1:8]
K173_Kappa_CDS_nonSyn <- read.csv("K173_Kappa_CDS_nonSyn_variant_info.csv")[,1:8]
K173_Kappa_gene_nonSyn <- read.csv("K173_Kappa_gene_nonSyn_variant_info.csv")[,1:8]
K173_Kappa_exon_nonSyn <- read.csv("K173_Kappa_exon_nonSyn_variant_info.csv")[,1:8]



K173_Kappa_gene_nonSyn_indels <- read.csv("K173_Kappa_gene_nonSyn_indels_variant_info.csv")[,1:8]
K173_Kappa_gene_nonSyn_snps <- read.csv("K173_Kappa_gene_nonSyn_snps_variant_info.csv")[,1:8]

K173_Kappa_gene_nonSyn_snps_1Alt <- subset(K173_Kappa_gene_nonSyn_snps, K173_Kappa_gene_nonSyn_snps$TYPE == "snp")
K173_Kappa_gene_nonSyn_snps_1Alt$MQM <- as.numeric(as.character(K173_Kappa_gene_nonSyn_snps_1Alt$MQM))
K173_Kappa_gene_nonSyn_snps_1Alt$MQMR <- as.numeric(as.character(K173_Kappa_gene_nonSyn_snps_1Alt$MQMR))


K173_TrueSeq_gene_nonSyn_indels <- read.csv("K173_TrueSeq_gene_nonSyn_indels_variant_info.csv")[,1:8]
K173_TrueSeq_gene_nonSyn_snps <- read.csv("K173_TrueSeq_gene_nonSyn_snps_variant_info.csv")[,1:8]

K173_TrueSeq_gene_nonSyn_snps_1Alt <- subset(K173_TrueSeq_gene_nonSyn_snps, K173_TrueSeq_gene_nonSyn_snps$TYPE == "snp")
K173_TrueSeq_gene_nonSyn_snps_1Alt$MQM <- as.numeric(as.character(K173_TrueSeq_gene_nonSyn_snps_1Alt$MQM))
K173_TrueSeq_gene_nonSyn_snps_1Alt$MQMR <- as.numeric(as.character(K173_TrueSeq_gene_nonSyn_snps_1Alt$MQMR))



#####Counts#####
K173_TrueSeq_types <- as.data.frame(table(K173_TrueSeq$TYPE))
K173_TrueSeq_CDS_types <- as.data.frame(table(K173_TrueSeq_CDS$TYPE))
K173_TrueSeq_gene_types <- as.data.frame(table(K173_TrueSeq_gene$TYPE))
K173_TrueSeq_exon_types <- as.data.frame(table(K173_TrueSeq_exon$TYPE))

K173_Kappa_types <- as.data.frame(table(K173_Kappa$TYPE))
K173_Kappa_CDS_types <- as.data.frame(table(K173_Kappa_CDS$TYPE))
K173_Kappa_gene_types <- as.data.frame(table(K173_Kappa_gene$TYPE))
K173_Kappa_exon_types <- as.data.frame(table(K173_Kappa_exon$TYPE))

K173_TrueSeq_nonSyn_types <- as.data.frame(table(K173_TrueSeq_nonSyn$TYPE))
K173_TrueSeq_CDS_nonSyn_types <- as.data.frame(table(K173_TrueSeq_CDS_nonSyn$TYPE))
K173_TrueSeq_gene_nonSyn_types <- as.data.frame(table(K173_TrueSeq_gene_nonSyn$TYPE))
K173_TrueSeq_exon_nonSyn_types <- as.data.frame(table(K173_TrueSeq_exon_nonSyn$TYPE))

K173_Kappa_nonSyn_types <- as.data.frame(table(K173_Kappa_nonSyn$TYPE))
K173_Kappa_CDS_nonSyn_types <- as.data.frame(table(K173_Kappa_CDS_nonSyn$TYPE))
K173_Kappa_gene_nonSyn_types <- as.data.frame(table(K173_Kappa_gene_nonSyn$TYPE))
K173_Kappa_exon_nonSyn_types <- as.data.frame(table(K173_Kappa_exon_nonSyn$TYPE))


####indels and snps#####
#focusing on K173 for now, making files in VCftools for indels and non indels to compare
#load the rest of the data in the same way, now have types, etc
K173_Kappa_gene_nonSyn_indels_types <- as.data.frame(table(K173_Kappa_gene_nonSyn_indels$TYPE))
K173_Kappa_gene_nonSyn_snps_types <- as.data.frame(table(K173_Kappa_gene_nonSyn_snps$TYPE))

K173_TrueSeq_gene_nonSyn_indels_types <- as.data.frame(table(K173_TrueSeq_gene_nonSyn_indels$TYPE))
K173_TrueSeq_gene_nonSyn_snps_types <- as.data.frame(table(K173_TrueSeq_gene_nonSyn_snps$TYPE))


#now I have a list of all the nonSyn gene variants in K173 Kappa, 
#as well as the vcftools parsed indels and non-indels

#for Andy and the report, focus on the variants that have only 1 ALT 
#and are the VCFtools definition of SNP

K173_Kappa_gene_nonSyn_indels_ALTs <- as.data.frame(table(K173_Kappa_gene_nonSyn_indels$NUMALT))
K173_Kappa_gene_nonSyn_snps_ALTs <- as.data.frame(table(K173_Kappa_gene_nonSyn_snps$NUMALT))

K173_TrueSeq_gene_nonSyn_indels_ALTs <- as.data.frame(table(K173_TrueSeq_gene_nonSyn_indels$NUMALT))
K173_TrueSeq_gene_nonSyn_snps_ALTs <- as.data.frame(table(K173_TrueSeq_gene_nonSyn_snps$NUMALT))

##okay, next step - parse out the ANN info for only the SNPs with 1 ALT allele
##this is pretty easy from the SNP files

K173_TrueSeq_gene_nonSyn_1ALT_SNP_list <- K173_TrueSeq_gene_nonSyn_snps_1Alt[,1:2]
K173_Kappa_gene_nonSyn_1ALT_SNP_list <- K173_Kappa_gene_nonSyn_snps_1Alt[,1:2]

TrueSeq_1ALT_SNPs_Chr <- as.data.frame(table(K173_TrueSeq_gene_nonSyn_snps_1Alt$CHROM))
Kappa_1ALT_SNPs_Chr <- as.data.frame(table(K173_Kappa_gene_nonSyn_snps_1Alt$CHROM))

TrueSeq_Indels_Chr <- as.data.frame(table(K173_TrueSeq_gene_nonSyn_indels$CHROM))
Kappa_Indels_Chr <- as.data.frame(table(K173_Kappa_gene_nonSyn_indels$CHROM))



#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/K173_variant_annotation_parsing.RData")
