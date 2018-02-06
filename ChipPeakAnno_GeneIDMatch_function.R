ChipPeakAnno_GeneIDMatch <- function(annotations, 
                                     gff_csv = "J:/III/Waters/Group Members/Mallory/PbergheiATACseqSebastian/GeneDB_28_PbergheiANKA_gff.csv",
                                     type){
  #Function will return a data frame containing the original annotation details as well as the GeneIDs pulled from 
  #a GFF file matched by feature number
    #annotation - name of annotation object
    #gff_csv - full path (in quotes) to the GFF csv
    #type - optional - either as an individual type (ie "gene") or a list of types (ie c("gene","CDS","mRNA"))
  
  #Convert annotations to dataframe#
    anno <- as.data.frame(annotations)
      #duplicate the feature column to featurenumber#
      anno$featurenumber <- as.numeric(anno$feature)
      
  #load GFF csv#
    #this is just the gff with the header removed and saved as a CSV for ease of loading into R#
    gff <- read.csv(gff_csv)
      #number the features in order#
      gff$featurenumber <- 1:nrow(gff)
    
  #Left join the annotations dataframe and GFF (only the ATTRIBUTES1, START,END, TYPE) using left_join#
    #ATTRIBUTES1 contains the GeneIDs#
    output <- left_join(anno, gff[,c("featurenumber","START","END","ATTRIBUTES1","TYPE")], by = "featurenumber")
    
    #check that the start and end coordinates match
      #this takes the START coordinate from the gff3 and the start_position from the annotated peaks
      #because we assigned the gff featurenumber based on its order in the file, there is a chance there was a mismatch between the 
      #feature listed in the Peaks and the feature we've designated with that number in the GFF
      #if the start_position from the Peaks and the START coord from the GFF (all of them in order) match, then we can be confident that 
      #the numbers we assigned to the features in the GFF file are correct based on what was assigned in the Peaks file
    #also performs the TYPE subset if necessary
    if(missing(type)){
    if((identical(output$START, output$start_position) == TRUE) && (identical(output$END, output$end_position) == TRUE)) {
      return(output)
      } else {
      print("Error: Start and/or end coordinates of the GFF and Annotations do not match") 
      }
    } else {
    if((identical(output$START, output$start_position) == TRUE) && (identical(output$END, output$end_position) == TRUE)) {
      return(subset(output,output$TYPE %in% type))
    } else {
      print("Error: Start and/or end coordinates of the GFF and Annotations do not match") 
    }
}
}