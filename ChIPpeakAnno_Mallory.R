#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/ChIPpeakAnno.RData")


# ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# 
# biocLite("ChIPpeakAnno")
# biocLite("rtracklayer")
# biocLite("GenomicRanges")

library(rtracklayer)
library(ChIPpeakAnno)
library(dplyr)

#####Sebastian's code with newly downloaded GFF (8Dec2017)#####
#load gff (8Dec2017)#
SK_Pbgff <- import.gff("J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/PlasmoDB-28_PbergheiANKA.gff",
                    'gff',
                    version = c("3"),
                    genome = NA, 
                    colnames = NULL, 
                    which = NULL,
                    feature.type = NULL, 
                    sequenceRegionsAsSeqinfo = FALSE)

#load MACs2 peak outputs#
SK_gr_507 <- toGRanges("J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/ChipPeakAnno/FromSebastian/507_4h_MACS2_peaks.broadPeak"
                    , format=c("broadPeak"),
                    header=FALSE, 
                    comment.char="#", 
                    colNames=NULL) 

#annotate peaks using the GFF#
SK_overlaps_anno <- annotatePeakInBatch(SK_gr_507, 
                                        AnnotationData=SK_Pbgff, 
                                        output="overlapping", 
                                        maxgap=5000L)
#Convert annotations to dataframe#
SK_AnnotatedPeaks_507_4h <- as.data.frame(SK_overlaps_anno)
  #duplicate the feature column to featurenumber#
  SK_AnnotatedPeaks_507_4h$featurenumber <- as.numeric(SK_AnnotatedPeaks_507_4h$feature)
  #output as CSV#
  write.csv(SK_AnnotatedPeaks_507_4h,"J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/ChipPeakAnno/SK_AnnotatedPeaks_507_4h.csv")

#load GFF csv (8Dec2017)#
  #this is just the gff with the header removed and saved as a CSV for ease of loading into R#
SK_Pbgff_csv <- read.csv("J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/GeneDB_28_PbergheiANKA_gff.csv")
  #number the features in order#
SK_Pbgff_csv$featurenumber <- 1:nrow(SK_Pbgff_csv)
  
#Left join the annotations dataframe and GFF (only the ATTRIBUTES1, START,END, TYPE) using left_join#
SK_AnnotatedPeaks_507_4h_GeneID <- left_join(SK_AnnotatedPeaks_507_4h, SK_Pbgff_csv[,c("featurenumber","START","END","ATTRIBUTES1","TYPE")], by = "featurenumber")
      #check that the start and end coordinates match
      #this takes the START coordinate from the gff3 and the start_position from the annotated peaks
      #because we assigned the gff featurenumber based on its order in the file, there is a chance there was a mismatch between the 
      #feature listed in the Peaks and the feature we've designated with that number in the GFF
      #if the start_position from the Peaks and the START coord from the GFF (all of them in order) match, then we can be confident that 
      #the numbers we assigned to the features in the GFF file are correct based on what was assigned in the Peaks file
        identical(SK_AnnotatedPeaks_507_4h_GeneID$START, SK_AnnotatedPeaks_507_4h_GeneID$start_position)
        identical(SK_AnnotatedPeaks_507_4h_GeneID$END, SK_AnnotatedPeaks_507_4h_GeneID$end_position)
  #ATTRIBUTES1 contains the GeneIDs#
  write.csv(SK_AnnotatedPeaks_507_4h_GeneID,"J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/ChipPeakAnno/SK_AnnotatedPeaks_507_4h_GeneID.csv")

####Annotated peaks file Sebastian sent####
##compare newly generated annotated peaks to sebastian's peaks file###
  #this uses sebastian's gff file - was the annotated peaks file sebastian sent
FromSK_AnnotatedPeaks_507_4h <- read.csv("J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/ChipPeakAnno/annotatedPeakListID_MVW.csv")
  #duplicate the feature column to featurenumber#
  FromSK_AnnotatedPeaks_507_4h$featurenumber <- as.numeric(FromSK_AnnotatedPeaks_507_4h$feature)

  #Left join the annotations dataframe and GFF (only the ATTRIBUTES1, START,END, TYPE) using left_join#
  FromSK_AnnotatedPeaks_507_4h_GeneID <- left_join(FromSK_AnnotatedPeaks_507_4h, SK_Pbgff_csv[,c("featurenumber","START","END","ATTRIBUTES1","TYPE")], by = "featurenumber")
    #check that the start and end coordinates match
    #this takes the START coordinate from the gff3 and the start_position from the annotated peaks
    #because we assigned the gff featurenumber based on its order in the file, there is a chance there was a mismatch between the 
    #feature listed in the Peaks and the feature we've designated with that number in the GFF
    #if the start_position from the Peaks and the START coord from the GFF (all of them in order) match, then we can be confident that 
    #the numbers we assigned to the features in the GFF file are correct based on what was assigned in the Peaks file
      identical(FromSK_AnnotatedPeaks_507_4h_GeneID$START, FromSK_AnnotatedPeaks_507_4h_GeneID$start_position)
      identical(FromSK_AnnotatedPeaks_507_4h_GeneID$END, FromSK_AnnotatedPeaks_507_4h_GeneID$end_position)
  #ATTRIBUTES1 contains the GeneIDs#
  write.csv(FromSK_AnnotatedPeaks_507_4h_GeneID,"J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/ChipPeakAnno/FromSK_AnnotatedPeaks_507_4h_GeneID.csv")
  

  
######My code with  new GFF file (8Dec2017)#####
#load gff (8Dec2017)#
New_Pbgff <- import.gff(con="J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/PlasmoDB-28_PbergheiANKA.gff",
                        format = "gff3",
                        genome = NA, 
                        colnames = NULL, 
                        which = NULL,
                        feature.type = NULL, 
                        sequenceRegionsAsSeqinfo = FALSE)
  # #convert gff to annoGR format#
  # New_Pbgff_anno <- annoGR(New_Pbgff)
 
setwd("J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/ChipPeakAnno/")
#load MACs2 peak outputs#
macs_507_gr <- toGRanges("J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/ChipPeakAnno/FromSebastian/507_4h_MACS2_peaks.broadPeak", 
                               format=c("broadPeak"),
                               header=FALSE, 
                               comment.char="#", 
                               colNames=NULL)
#annotate peaks using the GFF#
new_overlaps_anno <- annotatePeakInBatch(macs_507_gr, 
                                       AnnotationData=New_Pbgff, 
                                       # featureType = "gene",
                                       output="overlapping", 
                                       maxgap=5000L)
  

  # View(new_overlaps_anno)
  # head(new_overlaps_anno)
  
  #issue - this doesn't give the features their proper name, but does give the chrom and coord
  #can just input the csv of the gff and do it myself
  
#Convert annotations to dataframe#
New_AnnotatedPeaks_507_4h <- as.data.frame(new_overlaps_anno)
  #duplicate the feature column to featurenumber#
  New_AnnotatedPeaks_507_4h$featurenumber <- as.numeric(New_AnnotatedPeaks_507_4h$feature)
  #output as csv#
    write.csv(New_AnnotatedPeaks_507_4h,"J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/ChipPeakAnno/New_AnnotatedPeaks_507_4h.csv")
  
#load GFF csv (8Dec2017)#
#this is just the gff with the header removed and saved as a CSV for ease of loading into R#
  New_Pbgff_csv <- read.csv("J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/GeneDB_28_PbergheiANKA_gff.csv")
  #number the features in order#
  New_Pbgff_csv$featurenumber <- 1:nrow(New_Pbgff_csv)
  #output csv with featurenumber as a csv
  write.csv(New_Pbgff_csv, "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/GeneDB_28_PbergheiANKA_gff_8Dec2017.csv")

#Left join the annotations dataframe and GFF (only the ATTRIBUTES1, START,END, TYPE) using left_join#
    #ATTRIBUTES1 contains the GeneIDs#
New_AnnotatedPeaks_507_4h_GeneID <- left_join(New_AnnotatedPeaks_507_4h, New_Pbgff_csv[,c("featurenumber","START","END","ATTRIBUTES1","TYPE")], by = "featurenumber")

  #check that the start and end coordinates match
    #this takes the START coordinate from the gff3 and the start_position from the annotated peaks
    #because we assigned the gff featurenumber based on its order in the file, there is a chance there was a mismatch between the 
    #feature listed in the Peaks and the feature we've designated with that number in the GFF
    #if the start_position from the Peaks and the START coord from the GFF (all of them in order) match, then we can be confident that 
    #the numbers we assigned to the features in the GFF file are correct based on what was assigned in the Peaks file
  identical(New_AnnotatedPeaks_507_4h_GeneID$START, New_AnnotatedPeaks_507_4h_GeneID$start_position)
  identical(New_AnnotatedPeaks_507_4h_GeneID$END, New_AnnotatedPeaks_507_4h_GeneID$end_position)
  #output csv#
  write.csv(New_AnnotatedPeaks_507_4h_GeneID,"J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/ChipPeakAnno/New_AnnotatedPeaks_507_4h_GeneID.csv")
  

  
  binOverFeature(macs_507_gr, annotationData=new_overlaps_anno,
                 radius=2000, nbins=100, FUN=c(sum, length),
                 ylab=c("score", "count"), 
                 main=c("Distribution of aggregated peak scores around TSS", 
                        "Distribution of aggregated peak numbers around TSS"))
  
  
####FILTER to only contain genes####
New_AnnotatedPeaks_507_4h_GeneID_Genes <- subset(New_AnnotatedPeaks_507_4h_GeneID,New_AnnotatedPeaks_507_4h_GeneID$TYPE == "gene")
  write.csv(New_AnnotatedPeaks_507_4h_GeneID_Genes,"J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/ChipPeakAnno/New_AnnotatedPeaks_507_4h_GeneID_Genes.csv")
  
  
  

#####NOTES#####
  # so the AnnotatedPeaks file originally sent by Sebastian contains 75987 peaks rows, 
  # when I regenerate the annotations I get 75989 peak rows using either my code or Sebastian's code
  # identical(SK_AnnotatedPeaks_507_4h,New_AnnotatedPeaks_507_4h)
  # [1] TRUE
  # > identical(SK_AnnotatedPeaks_507_4h_GeneID,New_AnnotatedPeaks_507_4h_GeneID)
  # [1] TRUE
  
  
#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/ChIPpeakAnno.RData")
