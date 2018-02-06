#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/Pberghei_upstream_region_1000bp_RPF_amalgamation.RData")

#####function#####
reads_per_frag_amalgamate <- function(WD,filename_file){
  setwd(WD)
  filenames <- read.table(filename_file)[,1]
  ###tryign to figure out how to rename columns in x
  # filenames_change <- as.factor(sapply(filenames,gsub, pattern = "_reads_per_fragment.csv", replacement = ""))
  x <- read.csv(paste0(filenames[1]))[,c(2:3,5:10)]
  for (i in seq_along(filenames)){
    q <- read.csv(paste0(filenames[i]), header=T)
    # x <- cbind(q[,c(2:3,5:10)])
    # x <- cbind(q)
    colnames(q)[4] <- paste0(filenames[i])
    x <- merge(x,q[,c("GENEID",paste0(filenames[i]))], by = "GENEID")
  }
  cols <- colnames(x)
  cols <- as.list(sapply(cols,gsub,pattern = "_reads_per_fragment.csv",replacement = ""))
  colnames(x) <- cols
  return(x)
}

#####run function to amalgamate reads per fragment#####
#Exp1
write.csv(reads_per_frag_amalgamate(WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/1kbUpstream/",
                                    filename_file = "Exp1_1kb_filenames.txt"),
          "Exp1_1kb_reads_per_frag.csv")
write.csv(reads_per_frag_amalgamate(WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/2kbUpstream/",
                                    filename_file = "Exp1_2kb_filenames.txt"),
          "Exp1_2kb_reads_per_frag.csv")
write.csv(reads_per_frag_amalgamate(WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp1/CDS/",
                                    filename_file = "Exp1_CDS_filenames.txt"),
          "Exp1_CDS_reads_per_frag.csv")
#Exp2
write.csv(reads_per_frag_amalgamate(WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/1kbUpstream/",
                                    filename_file = "Exp2_1kb_filenames.txt"),
          "Exp2_1kb_reads_per_frag.csv")
write.csv(reads_per_frag_amalgamate(WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/2kbUpstream/",
                                    filename_file = "Exp2_2kb_filenames.txt"),
          "Exp2_2kb_reads_per_frag.csv")
write.csv(reads_per_frag_amalgamate(WD = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/Exp2/CDS/",
                                    filename_file = "Exp2_CDS_filenames.txt"),
          "Exp2_CDS_reads_per_frag.csv")



#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/Pberghei_upstream_region_1000bp_RPF_amalgamation.RData")
