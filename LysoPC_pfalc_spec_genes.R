#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/LysoPC_pfalc_spec_genes.RData")

#####set WD#####
setwd("J:/III/Waters/Group Members/Mallory/LysoPC_project/GeneLists")

# #####load gene lists#####
# #upreg
# upreg_genes <- read.csv("upreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
#   upreg_genes <- as.list(upreg_genes[,2])
# 
# upreg_Pberghei_ortho_genes <- read.csv("upreg_Pbergei_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
#   upreg_Pberghei_ortho_genes <- as.list(upreg_Pberghei_ortho_genes[,1])
# upreg_Pchabaudi_ortho_genes <- read.csv("upreg_Pchabaudi_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
#   upreg_Pchabaudi_ortho_genes <- as.list(upreg_Pchabaudi_ortho_genes[,1])
# upreg_Pyoelii17X_ortho_genes <- read.csv("upreg_Pyoelii17X_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
#   upreg_Pyoelii17X_ortho_genes <- as.list(upreg_Pyoelii17X_ortho_genes[,1])
# upreg_Pyoelii17XNL_ortho_genes <- read.csv("upreg_Pyoelii17XNL_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
#   upreg_Pyoelii17XNL_ortho_genes <- as.list(upreg_Pyoelii17XNL_ortho_genes[,1])
# upreg_PyoeliiYM_ortho_genes <- read.csv("upreg_PyoeliiYM_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
#   upreg_PyoeliiYM_ortho_genes <- as.list(upreg_PyoeliiYM_ortho_genes[,1])
# 
# #downreg
# downreg_genes <- read.csv("downreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
#   downreg_genes <- as.list(downreg_genes[,2])
#   
# downreg_Pberghei_ortho_genes <- read.csv("downreg_Pbergei_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
#   downreg_Pberghei_ortho_genes <- as.list(downreg_Pberghei_ortho_genes[,1])
# downreg_Pchabaudi_ortho_genes <- read.csv("downreg_Pchabaudi_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
#   downreg_Pchabaudi_ortho_genes <- as.list(downreg_Pchabaudi_ortho_genes[,1])
# downreg_Pyoelii17X_ortho_genes <- read.csv("downreg_Pyoelii17X_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
#   downreg_Pyoelii17X_ortho_genes <- as.list(downreg_Pyoelii17X_ortho_genes[,1])
# downreg_Pyoelii17XNL_ortho_genes <- read.csv("downreg_Pyoelii17XNL_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
#   downreg_Pyoelii17XNL_ortho_genes <- as.list(downreg_Pyoelii17XNL_ortho_genes[,1])
# downreg_PyoeliiYM_ortho_genes <- read.csv("downreg_PyoeliiYM_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
#   downreg_PyoeliiYM_ortho_genes <- as.list(downreg_PyoeliiYM_ortho_genes[,1])
# 
# 
# #####differences between lists#####
# upreg_Pfalc_specific <- list()
# upreg_Pfalc_specific <- list(upreg_Pfalc_NOTPberghei = upreg_genes[!(upreg_genes %in% upreg_Pberghei_ortho_genes )], 
#                              upreg_Pfalc_NOTPchabaudi = upreg_genes[!(upreg_genes %in% upreg_Pchabaudi_ortho_genes)],
#                              upreg_Pfalc_NOTPyoelii17X = upreg_genes[!(upreg_genes %in% upreg_Pyoelii17X_ortho_genes)],
#                              upreg_Pfalc_NOTPyoelii17XNL = upreg_genes[!(upreg_genes %in% upreg_Pyoelii17XNL_ortho_genes)],
#                              upreg_Pfalc_NOTPyoeliiYM = upreg_genes[!(upreg_genes %in% upreg_PyoeliiYM_ortho_genes)])
# 
# upreg_Pfalc_NOTrodent <- upreg_genes[!(upreg_genes %in% c(upreg_Pberghei_ortho_genes,
#                                          upreg_Pchabaudi_ortho_genes,
#                                          upreg_Pyoelii17X_ortho_genes,
#                                          upreg_Pyoelii17XNL_ortho_genes,
#                                          upreg_PyoeliiYM_ortho_genes))]
# 
# 
# downreg_Pfalc_specific <- list()
# downreg_Pfalc_specific <- list(downreg_Pfalc_NOTPberghei = downreg_genes[!(downreg_genes %in% downreg_Pberghei_ortho_genes)], 
#                              downreg_Pfalc_NOTPchabaudi = downreg_genes[!(downreg_genes %in% downreg_Pchabaudi_ortho_genes)],
#                              downreg_Pfalc_NOTPyoelii17X = downreg_genes[!(downreg_genes %in% downreg_Pyoelii17X_ortho_genes)],
#                              downreg_Pfalc_NOTPyoelii17XNL = downreg_genes[!(downreg_genes %in% downreg_Pyoelii17XNL_ortho_genes)],
#                              downreg_Pfalc_NOTPyoeliiYM = downreg_genes[!(downreg_genes %in% downreg_PyoeliiYM_ortho_genes)])
# 
# downreg_Pfalc_NOTrodent <- downreg_genes[!(downreg_genes %in% c(downreg_Pberghei_ortho_genes,
#                                                           downreg_Pchabaudi_ortho_genes,
#                                                           downreg_Pyoelii17X_ortho_genes,
#                                                           downreg_Pyoelii17XNL_ortho_genes,
#                                                           downreg_PyoeliiYM_ortho_genes))]
# 
# 
# 
# 
# 
# #####write csv of each difference list#####
# setwd("J:/III/Waters/Group Members/Mallory/LysoPC_project/GeneOrthologs")
# 
# for (i in seq_along(upreg_Pfalc_specific)) {
#   filenames = paste0(names(upreg_Pfalc_specific)[i], ".csv")
#   write.csv(upreg_Pfalc_specific[[i]], filenames)
# }
# write.csv(upreg_Pfalc_NOTrodent, "upreg_Pfalc_NOTrodent.csv")
# 
# for (i in seq_along(downreg_Pfalc_specific)) {
#   filenames = paste0(names(downreg_Pfalc_specific)[i], ".csv")
#   write.csv(downreg_Pfalc_specific[[i]], filenames)
# }
# write.csv(downreg_Pfalc_NOTrodent, "downreg_Pfalc_NOTrodent.csv")




#########REDO 14 March 2017#######
#upreg
upreg_genes <- read.csv("upreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  colnames(upreg_genes) <- c("X","GENEID")

upreg_Pberghei_ortho_genes <- read.csv("upreg_Pbergei_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
upreg_Pchabaudi_ortho_genes <- read.csv("upreg_Pchabaudi_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
upreg_Pyoelii17X_ortho_genes <- read.csv("upreg_Pyoelii17X_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
upreg_Pyoelii17XNL_ortho_genes <- read.csv("upreg_Pyoelii17XNL_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
upreg_PyoeliiYM_ortho_genes <- read.csv("upreg_PyoeliiYM_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)

upreg_all_ortho_genes <- rbind(upreg_Pberghei_ortho_genes, 
                               upreg_Pchabaudi_ortho_genes, 
                               upreg_Pyoelii17X_ortho_genes, 
                               upreg_Pyoelii17XNL_ortho_genes,
                               upreg_PyoeliiYM_ortho_genes)
  colnames(upreg_all_ortho_genes) <- "GENEID"
upreg_all_ortho_genes_unique <- as.data.frame(unique(upreg_all_ortho_genes[,1]))
  colnames(upreg_all_ortho_genes_unique) <- "GENEID"

upreg_Pfalc_NOTrodent <- as.data.frame(subset(upreg_genes, !upreg_genes$GENEID %in% upreg_all_ortho_genes_unique$GENEID))

write.csv(upreg_Pfalc_NOTrodent, "upreg_Pfalc_NOTrodent.csv")

#downreg
downreg_genes <- read.csv("downreg_strand_Niggi_noDups_upstream_regions_genelist.csv")
  colnames(downreg_genes) <- c("X","GENEID")

downreg_Pberghei_ortho_genes <- read.csv("downreg_Pbergei_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
downreg_Pchabaudi_ortho_genes <- read.csv("downreg_Pchabaudi_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
downreg_Pyoelii17X_ortho_genes <- read.csv("downreg_Pyoelii17X_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
downreg_Pyoelii17XNL_ortho_genes <- read.csv("downreg_Pyoelii17XNL_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)
downreg_PyoeliiYM_ortho_genes <- read.csv("downreg_PyoeliiYM_orthologs_Niggi_noDups_upstream_regions_genelist.csv", header=F)

downreg_all_ortho_genes <- rbind(downreg_Pberghei_ortho_genes, 
                               downreg_Pchabaudi_ortho_genes, 
                               downreg_Pyoelii17X_ortho_genes, 
                               downreg_Pyoelii17XNL_ortho_genes,
                               downreg_PyoeliiYM_ortho_genes)
  colnames(downreg_all_ortho_genes) <- "GENEID"
downreg_all_ortho_genes_unique <- as.data.frame(unique(downreg_all_ortho_genes[,1]))
  colnames(downreg_all_ortho_genes_unique) <- "GENEID"

downreg_Pfalc_NOTrodent <- as.data.frame(subset(downreg_genes, !downreg_genes$GENEID %in% downreg_all_ortho_genes_unique$GENEID))

write.csv(downreg_Pfalc_NOTrodent, "downreg_Pfalc_NOTrodent.csv")

#####REMINDER - have changed the upreg list based on MATT MARTI's data - only 24 UPREG PFALC genes now 24 May 2017#####

#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/LysoPC_pfalc_spec_genes.RData")
